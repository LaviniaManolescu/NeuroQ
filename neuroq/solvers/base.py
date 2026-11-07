"""Clasa de bază pentru solverele pe Ising."""

from __future__ import annotations

from abc import ABC, abstractmethod

from neuroq.core import IsingProblem, Solution


class Solver(ABC):
    """Caută șirul de energie minimă al unui IsingProblem."""

    @abstractmethod
    def solve(self, problem: IsingProblem) -> Solution:
        """Întoarce soluția, cu șirul în ordinea proiectului."""
