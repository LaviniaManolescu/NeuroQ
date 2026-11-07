"""Codificarea pe poziții (metoda organizatorilor): un qubit pe SNP."""

from __future__ import annotations

import itertools
from typing import Sequence

import numpy as np

from neuroq.core import MISSING, Haplotypes, IsingProblem, ReadMatrix
from neuroq.encodings.base import Encoding

SUSPECT = -1


# din notebook-ul organizatorilor
def build_evidence(reads: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Contoarele A („la fel”) și D („opus”) pe perechi de poziții, triunghiul de sus."""
    n = reads.shape[1]
    a = np.zeros((n, n))
    d = np.zeros((n, n))
    for row in reads:
        observed = np.flatnonzero(row >= 0)
        for i, j in itertools.combinations(observed, 2):
            if row[i] == row[j]:
                a[i, j] += 1
            else:
                d[i, j] += 1
    return a, d


# din notebook-ul organizatorilor
def ising_from_evidence(a: np.ndarray, d: np.ndarray) -> tuple[np.ndarray, float]:
    """J = (D − A) / 2 și offsetul Σ (A + D) / 2."""
    return (d - a) / 2, float((a + d).sum() / 2)


# din notebook-ul organizatorilor
def pairwise_cost(x: Sequence[int], a: np.ndarray, d: np.ndarray) -> float:
    """Numărul de voturi încălcate; folosit doar ca verificare independentă a energiei."""
    total = 0.0
    for i, j in itertools.combinations(range(len(x)), 2):
        total += a[i, j] if x[i] != x[j] else d[i, j]
    return float(total)


def assign_reads(values: np.ndarray, copy_a: np.ndarray, copy_b: np.ndarray) -> np.ndarray:
    """Fiecare citire merge la copia cu mai puține nepotriviri; egalitate = suspectă."""
    assignment = np.empty(values.shape[0], dtype=int)
    for r, row in enumerate(values):
        seen = row != MISSING
        miss_a = int((row[seen] != copy_a[seen]).sum())
        miss_b = int((row[seen] != copy_b[seen]).sum())
        assignment[r] = 0 if miss_a < miss_b else 1 if miss_b < miss_a else SUSPECT
    return assignment


class PositionEncoding(Encoding):
    """Un qubit pe poziție; șirul măsurat este copia A, copia B este complementul."""

    def encode(self, matrix: ReadMatrix) -> IsingProblem:
        a, d = build_evidence(matrix.values)
        J, offset = ising_from_evidence(a, d)
        values = matrix.values.copy()

        def decode(bitstring: str) -> Haplotypes:
            copy_a = np.array([int(c) for c in bitstring])
            copy_b = 1 - copy_a
            return Haplotypes(copy_a, copy_b, assign_reads(values, copy_a, copy_b))

        return IsingProblem(J=J, offset=offset, variable_kind="position", decoder=decode)
