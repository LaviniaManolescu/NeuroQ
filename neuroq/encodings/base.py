"""Clasa de bază pentru codificări."""

from __future__ import annotations

from abc import ABC, abstractmethod

from neuroq.core import IsingProblem, ReadMatrix


class Encoding(ABC):
    """Transformă o matrice de citiri într-o problemă Ising."""

    @abstractmethod
    def encode(self, matrix: ReadMatrix) -> IsingProblem:
        """Construiește IsingProblem, cu decoder-ul legat de această matrice."""
