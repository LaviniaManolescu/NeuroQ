"""Referința exactă: enumerează toate șirurile."""

from __future__ import annotations

import itertools

import numpy as np

from neuroq.core import IsingProblem, Solution
from neuroq.solvers.base import Solver


class ExactSolver(Solver):
    """Enumeră toate cele 2^n șiruri; doar pentru probleme mici."""

    def __init__(self, max_qubits: int = 16) -> None:
        self.max_qubits = max_qubits

    def solve(self, problem: IsingProblem) -> Solution:
        n = problem.n_qubits
        if n > self.max_qubits:
            raise ValueError(f"{n} qubiți depășesc limita de {self.max_qubits} pentru enumerare")
        bitstrings = ["".join(bits) for bits in itertools.product("01", repeat=n)]
        energies = np.array([problem.energy(b) for b in bitstrings])
        best = float(energies.min())
        optimal = [b for b, e in zip(bitstrings, energies) if np.isclose(e, best)]
        return Solution(
            bitstring=optimal[0],
            energy=best,
            metadata={"optimal_bitstrings": optimal, "n_states": len(bitstrings)},
        )
