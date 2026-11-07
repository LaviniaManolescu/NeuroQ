"""Solver clasic puternic: arbore de acoperire cu paritate, căutare locală, recoacere simulată."""

from __future__ import annotations

import numpy as np

from neuroq.core import IsingProblem, Solution
from neuroq.solvers import _classical_algos as algos
from neuroq.solvers.base import Solver

STRATEGIES = ("greedy_tree", "greedy_ls", "anneal", "best", "exact")
MAX_EXACT_QUBITS = 24  # exact_vectorized enumeră 2^(n-1) stări în numpy


def evidence_from_ising(problem: IsingProblem) -> tuple[np.ndarray, np.ndarray]:
    """Matricile (a, d) cerute de algoritmi: a = max(-2J, 0), d = max(2J, 0).

    Contează doar diferența a - d = -2J și existența muchiei; energia raportată
    se calculează oricum cu IsingProblem.energy.
    """
    return np.maximum(-2 * problem.J, 0.0), np.maximum(2 * problem.J, 0.0)


class ClassicalSolver(Solver):
    """Rezolvă un IsingProblem fără câmpuri h cu algoritmii clasici, cu seed fix."""

    def __init__(
        self, strategy: str = "best", seed: int = 0, restarts: int = 20, steps: int = 4000
    ) -> None:
        if strategy not in STRATEGIES:
            raise ValueError(f"strategy trebuie să fie una din {STRATEGIES}, nu {strategy!r}")
        self.strategy = strategy
        self.seed = seed
        self.restarts = restarts
        self.steps = steps

    def solve(self, problem: IsingProblem) -> Solution:
        if np.any(problem.h != 0):
            raise ValueError(
                "ClassicalSolver lucrează doar pe probleme fără câmpuri h; "
                "aplică-l înainte de fix_variable"
            )
        if self.strategy == "exact" and problem.n_qubits > MAX_EXACT_QUBITS:
            raise ValueError(
                f"{problem.n_qubits} qubiți depășesc limita de {MAX_EXACT_QUBITS} pentru enumerare"
            )
        a, d = evidence_from_ising(problem)
        x, info = self._run(a, d)
        bitstring = "".join(str(int(b)) for b in x)
        return Solution(
            bitstring=bitstring,
            energy=problem.energy(bitstring),
            metadata={"strategy": self.strategy, "seed": self.seed, **info},
        )

    def _run(self, a: np.ndarray, d: np.ndarray) -> tuple[np.ndarray, dict]:
        if self.strategy == "exact":
            x, _ = algos.exact_vectorized(a, d)
            return x, {"exact_backend": "exact_vectorized"}
        if self.strategy == "best":
            x, _, info = algos.solve_classical(
                a, d, restarts=self.restarts, steps=self.steps, seed=self.seed
            )
            return x, {
                "winner": info["winner"],
                "components": int(info["components"]),
                # costurile interne ale algoritmilor pe (a, d) adaptat, nu energii IsingProblem
                "adapted_scores": {k: float(v) for k, v in info["scores"].items()},
            }
        if self.strategy == "anneal":
            x = algos.local_search(algos.anneal(a, d, steps=self.steps, seed=self.seed), a, d)
            return x, {}
        x, components = algos.greedy_tree(a, d)
        if self.strategy == "greedy_ls":
            x = algos.local_search(x, a, d)
        return x, {"components": int(components)}
