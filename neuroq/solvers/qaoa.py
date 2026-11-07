"""QAOA: circuit H, rz, rzz, rx, antrenat pe Statevector cu COBYLA."""

from __future__ import annotations

import itertools

import numpy as np
from qiskit import QuantumCircuit
from qiskit.circuit import ParameterVector
from qiskit.quantum_info import Statevector
from scipy.optimize import minimize

from neuroq.core import IsingProblem, Solution
from neuroq.execution import Executor
from neuroq.solvers.base import Solver

N_STARTS = 2


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
    return np.array([
        problem.energy("".join(str((k >> i) & 1) for i in range(n))) for k in range(2**n)
    ])


class QAOASolver(Solver):
    """Antrenează unghiurile pe Statevector, apoi eșantionează prin Executor."""

    def __init__(
        self, executor: Executor, p: int = 1, shots: int = 2048, seed: int = 71, maxiter: int = 90
    ) -> None:
        self.executor = executor
        self.p = p
        self.shots = shots
        self.seed = seed
        self.maxiter = maxiter

    def solve(self, problem: IsingProblem) -> Solution:
        ansatz, theta, scale = make_qaoa(problem.J, problem.h, self.p)
        energies = basis_energies(problem)

        def bind(values: np.ndarray) -> QuantumCircuit:
            return ansatz.assign_parameters(dict(zip(theta, values)))

        # din notebook-ul organizatorilor
        def expected_energy(values: np.ndarray) -> float:
            probabilities = Statevector.from_instruction(bind(values)).probabilities()
            return float(np.dot(probabilities, energies))

        rng = np.random.default_rng(self.seed)
        starts = [rng.uniform(0, np.pi, 2 * self.p) for _ in range(N_STARTS)]
        trials = [
            minimize(expected_energy, start, method="COBYLA", options={"maxiter": self.maxiter})
            for start in starts
        ]
        trained = min(trials, key=lambda result: result.fun)

        measured = bind(trained.x)
        measured.measure_all()
        counts, resources = self.executor.run(measured, shots=self.shots, seed=self.seed + 1)

        scored = {bits: problem.energy(bits) for bits in counts}
        best = min(scored, key=lambda bits: (scored[bits], -counts[bits], bits))
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
            },
        )
