"""QAOA (H, rz, rzz, rx) și warm-start QAOA (ry, rz, rzz, mixer ry-rz-ry), antrenate pe
statevector exact (Aer) cu COBYLA."""

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
from neuroq.solvers.warm_start import DEFAULT_EPSILON, regularize, relax

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


def make_warm_start_qaoa(
    j: np.ndarray, h: np.ndarray, c: np.ndarray, p: int = 1
) -> tuple[QuantumCircuit, ParameterVector, float]:
    """Warm-start QAOA (Egger et al., 2021): pornește din probabilitățile c ale relaxării.

    Starea inițială e ry(θᵢ)|0⟩ cu θᵢ = 2·arcsin(√cᵢ), deci P(bitul i = 1) = cᵢ. Stratul de
    cost e același ca la make_qaoa; mixerul ry(−θᵢ) rz(−2β) ry(θᵢ) are starea inițială ca
    stare proprie.
    """
    n = len(j)
    scale = max(1.0, float(np.abs(j).max(initial=0.0)), float(np.abs(h).max(initial=0.0)))
    angles = 2 * np.arcsin(np.sqrt(np.asarray(c, dtype=float)))
    theta = ParameterVector("theta", 2 * p)
    qc = QuantumCircuit(n)
    for i in range(n):
        qc.ry(float(angles[i]), i)
    for layer in range(p):
        for i in range(n):
            if h[i] != 0:
                qc.rz(2 * theta[layer] * h[i] / scale, i)
        for i, k in itertools.combinations(range(n), 2):
            if j[i, k] != 0:
                qc.rzz(2 * theta[layer] * j[i, k] / scale, i, k)
        for i in range(n):
            qc.ry(-float(angles[i]), i)
            qc.rz(-2 * theta[p + layer], i)
            qc.ry(float(angles[i]), i)
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


def interp_angles(angles: np.ndarray, p_old: int) -> np.ndarray:
    """INTERP (Zhou et al., 2018): unghiurile antrenate la p_old straturi → start pentru p_old + 1.

    Schema e cea din make_qaoa: [γ₁..γp, β₁..βp]. γ și β se interpolează liniar fiecare,
    x_nou[i] = (i/p)·x[i−1] + ((p−i)/p)·x[i], i = 0..p, cu x[−1] = x[p] = 0.
    """
    angles = np.asarray(angles, dtype=float)
    if p_old < 1 or angles.shape != (2 * p_old,):
        raise ValueError(f"pentru p = {p_old} trebuie {2 * p_old} unghiuri, nu {angles.shape}")
    halves = []
    for part in (angles[:p_old], angles[p_old:]):
        padded = np.concatenate(([0.0], part, [0.0]))  # padded[i + 1] = part[i]
        halves.append([(i / p_old) * padded[i] + ((p_old - i) / p_old) * padded[i + 1]
                       for i in range(p_old + 1)])
    return np.concatenate(halves)


class QAOASolver(Solver):
    """Antrenează unghiurile pe statevector exact, apoi eșantionează prin Executor.

    `initial_angles` (2·p valori) înlocuiește primul start aleator (de exemplu cu rezultatul
    INTERP de la p−1); al doilea start rămâne cel aleator, cu același seed.
    """

    def __init__(
        self, executor: Executor, p: int = 1, shots: int = 2048, seed: int = 71, maxiter: int = 90,
        callback: Optional[CostCallback] = None, warm_start: bool = False,
        epsilon: float = DEFAULT_EPSILON, relax_seed: Optional[int] = None,
        initial_angles: Optional[np.ndarray] = None,
    ) -> None:
        if initial_angles is not None and len(initial_angles) != 2 * p:
            raise ValueError(f"initial_angles trebuie să aibă {2 * p} valori pentru p = {p}")
        self.initial_angles = initial_angles
        self.executor = executor
        self.callback = callback
        self.warm_start = warm_start
        self.epsilon = epsilon
        self.relax_seed = seed if relax_seed is None else relax_seed
        self.p = p
        self.shots = shots
        self.seed = seed
        self.maxiter = maxiter

    def solve(self, problem: IsingProblem) -> Solution:
        warm_c = relaxed = None
        if self.warm_start:
            relaxed = relax(problem, self.relax_seed)
            warm_c = regularize(relaxed, self.epsilon)
            ansatz, theta, scale = make_warm_start_qaoa(problem.J, problem.h, warm_c, self.p)
        else:
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
        if self.initial_angles is not None:
            starts[0] = np.asarray(self.initial_angles, dtype=float)
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
                # c = înclinările regularizate (cele din circuit); c_relaxed = relaxarea brută
                "warm_start": (None if warm_c is None else {
                    "c": warm_c.tolist(), "c_relaxed": relaxed.tolist(), "epsilon": self.epsilon}),
                "best_start": best_start,
                "top_samples": [
                    {"bitstring": b, "count": int(counts[b]), "energy": scored[b]} for b in ranked
                ],
            },
        )


def train_layers(
    problem: IsingProblem, executor: Executor, p_max: int, **options
) -> list[Solution]:
    """Antrenează p = 1, 2, ..., p_max în ordine; fiecare pornește INTERP din cel dinainte.

    `options` merg la QAOASolver (seed, maxiter, warm_start, callback...). Întoarce câte o
    Solution pentru fiecare p.
    """
    solutions: list[Solution] = []
    for p in range(1, p_max + 1):
        initial = interp_angles(solutions[-1].metadata["angles"], p - 1) if solutions else None
        solutions.append(QAOASolver(executor, p=p, initial_angles=initial, **options).solve(problem))
    return solutions


def trained_circuit(
    problem: IsingProblem, angles: list[float], p: int = 1,
    warm_c: Optional[list[float]] = None,
) -> QuantumCircuit:
    """Circuitul QAOA (sau WS-QAOA, cu `warm_c`) cu unghiurile antrenate și măsurători."""
    if warm_c is None:
        ansatz, theta, _ = make_qaoa(problem.J, problem.h, p)
    else:
        ansatz, theta, _ = make_warm_start_qaoa(problem.J, problem.h, np.asarray(warm_c), p)
    circuit = ansatz.assign_parameters(dict(zip(theta, angles)))
    circuit.measure_all()
    return circuit
