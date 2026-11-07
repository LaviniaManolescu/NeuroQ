"""Relaxarea continuă pentru warm-start QAOA (Egger et al., 2021). Fără Qiskit.

Fiecare bit devine o probabilitate c ∈ [0, 1] (c = 1 înseamnă bitul 1), cu z = 1 − 2c.
Cu h = 0, punctul c = 0,5 peste tot e staționar (oglinda A↔B), deci pornim din mai multe
puncte aleatorii cu seed fix: soluția relaxată alege o orientare, ca fixarea unui qubit.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import minimize

from neuroq.core import IsingProblem

N_STARTS = 10
DEFAULT_EPSILON = 0.25


def relaxed_energy(problem: IsingProblem, c: np.ndarray) -> float:
    """Energia pentru biți continui: offset + Σ h·z + Σ J·z·z, cu z = 1 − 2c.

    Pe șiruri 0/1 coincide cu IsingProblem.energy (verificat în teste); există separat
    doar pentru că IsingProblem.energy primește șiruri de biți, nu valori continue.
    """
    z = 1.0 - 2.0 * np.asarray(c, dtype=float)
    return float(problem.offset + problem.h @ z + z @ problem.J @ z)


def _gradient(problem: IsingProblem, c: np.ndarray) -> np.ndarray:
    z = 1.0 - 2.0 * c
    return -2.0 * (problem.h + (problem.J + problem.J.T) @ z)


def relax(problem: IsingProblem, seed: int, n_starts: int = N_STARTS) -> np.ndarray:
    """Minimul relaxării din `n_starts` porniri aleatorii cu seed fix (L-BFGS-B pe [0, 1]ⁿ)."""
    rng = np.random.default_rng(seed)
    bounds = [(0.0, 1.0)] * problem.n_qubits
    best_c, best_value = None, np.inf
    for _ in range(n_starts):
        result = minimize(
            lambda c: relaxed_energy(problem, c), rng.uniform(0.0, 1.0, problem.n_qubits),
            jac=lambda c: _gradient(problem, c), method="L-BFGS-B", bounds=bounds,
        )
        if result.fun < best_value - 1e-12:
            best_c, best_value = np.clip(result.x, 0.0, 1.0), float(result.fun)
    return best_c


def regularize(c: np.ndarray, epsilon: float = DEFAULT_EPSILON) -> np.ndarray:
    """Taie valorile în [ε, 1 − ε], ca circuitul să nu rămână blocat în șirul de start."""
    if not 0.0 <= epsilon <= 0.5:
        raise ValueError("epsilon trebuie să fie între 0 și 0.5")
    return np.clip(np.asarray(c, dtype=float), epsilon, 1.0 - epsilon)
