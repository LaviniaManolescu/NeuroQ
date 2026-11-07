"""QAOA: circuit H, rz, rzz, rx, antrenat pe statevector exact (Aer) cu COBYLA."""

from __future__ import annotations

import itertools
from typing import Callable, Optional

import numpy as np
from qiskit import QuantumCircuit
from qiskit.circuit import ParameterVector
from qiskit_aer import AerSimulator
from scipy.optimize import minimize

from neuroq.core import IsingProblem, Solution
from neuroq.execution import Executor
from neuroq.solvers.base import Solver

N_STARTS = 2
TOP_SAMPLES = 10

# callback(start, iteration, cost): apelat la fiecare evaluare COBYLA, ca să se vadă live
CostCallback = Callable[[int, int, float], None]


# din notebook-ul organizatorilor (extins cu porți rz pentru câmpurile h)
def make_qaoa(
    j: np.ndarray, h: np.ndarray, p: int = 1
) -> tuple[QuantumCircuit, ParameterVector, float]:
    n = len(j)
    scale = max(1.0, float(np.abs(j).max(initial=0.0)), float(np.abs(h).max(initial=0.0)))
    theta = ParameterVector("theta", 2 * p)
    qc = QuantumCircuit(n)
    qc.h(range(n))
    for layer in range(p):
        for i in range(n):
            if h[i] != 0:
                qc.rz(2 * theta[layer] * h[i] / scale, i)
        for i, k in itertools.combinations(range(n), 2):
            if j[i, k] != 0:
                qc.rzz(2 * theta[layer] * j[i, k] / scale, i, k)
        qc.rx(2 * theta[p + layer], range(n))
    return qc, theta, scale


def basis_energies(problem: IsingProblem) -> np.ndarray:
    """Energia fiecărui indice de statevector; indicele k are qubitul i = (k >> i) & 1."""
    n = problem.n_qubits
    return np.fromiter(
        (problem.energy("".join(str((k >> i) & 1) for i in range(n))) for k in range(2**n)),
        dtype=float, count=2**n,
    )


def statevector_probabilities(circuit: QuantumCircuit, simulator: AerSimulator) -> np.ndarray:
    """Probabilitățile exacte ale stării (fără măsurători), indexate ca Statevector.

    Aer dă aceleași valori ca qiskit.quantum_info.Statevector (verificat în teste),
    dar e de ~15 ori mai rapid la 20 de qubiți; tot statevector exact, nu Estimator.
    """
    circuit = circuit.copy()
    circuit.save_statevector()
    return np.asarray(simulator.run(circuit).result().get_statevector().probabilities())


class QAOASolver(Solver):
    """Antrenează unghiurile pe statevector exact, apoi eșantionează prin Executor."""

    def __init__(
        self, executor: Executor, p: int = 1, shots: int = 2048, seed: int = 71, maxiter: int = 90,
        callback: Optional[CostCallback] = None,
    ) -> None:
        self.executor = executor
        self.callback = callback
        self.p = p
        self.shots = shots
        self.seed = seed
        self.maxiter = maxiter

    def solve(self, problem: IsingProblem) -> Solution:
        ansatz, theta, scale = make_qaoa(problem.J, problem.h, self.p)
        energies = basis_energies(problem)
        simulator = AerSimulator(method="statevector", max_parallel_threads=1)

        def bind(values: np.ndarray) -> QuantumCircuit:
            return ansatz.assign_parameters(dict(zip(theta, values)))

        # din notebook-ul organizatorilor
        def expected_energy(values: np.ndarray) -> float:
            probabilities = statevector_probabilities(bind(values), simulator)
            return float(np.dot(probabilities, energies))

        histories: list[list[float]] = []

        def recorded(values: np.ndarray) -> float:
            cost = expected_energy(values)
            histories[-1].append(cost)
            if self.callback is not None:
                self.callback(len(histories) - 1, len(histories[-1]), cost)
            return cost

        rng = np.random.default_rng(self.seed)
        starts = [rng.uniform(0, np.pi, 2 * self.p) for _ in range(N_STARTS)]
        trials = []
        for start in starts:
            histories.append([])
            trials.append(minimize(recorded, start, method="COBYLA",
                                   options={"maxiter": self.maxiter}))
        best_start = min(range(len(trials)), key=lambda s: trials[s].fun)
        trained = trials[best_start]

        measured = bind(trained.x)
        measured.measure_all()
        counts, resources = self.executor.run(measured, shots=self.shots, seed=self.seed + 1)

        scored = {bits: problem.energy(bits) for bits in counts}
        best = min(scored, key=lambda bits: (scored[bits], -counts[bits], bits))
        ranked = sorted(counts, key=lambda bits: (-counts[bits], bits))[:TOP_SAMPLES]
        return Solution(
            bitstring=best,
            energy=scored[best],
            samples=counts,
            resources=resources,
            metadata={
                "p": self.p,
                "angles": trained.x.tolist(),
                "expected_energy": float(trained.fun),
                "nfev": int(sum(t.nfev for t in trials)),
                "scale": scale,
                "seed": self.seed,
                "cost_history": histories,
                "best_start": best_start,
                "top_samples": [
                    {"bitstring": b, "count": int(counts[b]), "energy": scored[b]} for b in ranked
                ],
            },
        )
