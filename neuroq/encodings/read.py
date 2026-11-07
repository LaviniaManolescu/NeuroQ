"""Codificarea pe citiri (metoda noastră): un qubit pe citire."""

from __future__ import annotations

import itertools

import numpy as np

from neuroq.core import MISSING, Haplotypes, IsingProblem, ReadMatrix
from neuroq.encodings.base import Encoding

UNPHASED = -1  # coloană fără majoritate în niciuna dintre echipe
SUSPECT = -1


def build_read_evidence(m: ReadMatrix) -> tuple[np.ndarray, np.ndarray]:
    """Pentru fiecare pereche de citiri i < j: potriviri a și contradicții d pe SNP-urile comune."""
    values = m.values
    n = values.shape[0]
    a = np.zeros((n, n))
    d = np.zeros((n, n))
    for i, j in itertools.combinations(range(n), 2):
        both = (values[i] != MISSING) & (values[j] != MISSING)
        same = values[i, both] == values[j, both]
        a[i, j] = np.count_nonzero(same)
        d[i, j] = np.count_nonzero(~same)
    return a, d


def team_majority(values: np.ndarray) -> np.ndarray:
    """Majoritatea pe fiecare coloană a unei echipe de citiri; egalitate sau lipsă = -1."""
    ones = (values == 1).sum(axis=0)
    zeros = (values == 0).sum(axis=0)
    return np.where(ones > zeros, 1, np.where(zeros > ones, 0, UNPHASED))


def _fill_from_other(own: np.ndarray, other: np.ndarray) -> np.ndarray:
    """Unde echipa proprie e nedecisă, ia complementul celeilalte, dacă aceea e decisă."""
    filled = own.copy()
    gap = (own == UNPHASED) & (other != UNPHASED)
    filled[gap] = 1 - other[gap]
    return filled


def _mismatches(row: np.ndarray, copy: np.ndarray) -> int:
    known = (row != MISSING) & (copy != UNPHASED)
    return int(np.count_nonzero(row[known] != copy[known]))


class ReadEncoding(Encoding):
    """Un qubit pe citire; bitul spune echipa (0 = copia A, 1 = copia B)."""

    def encode(self, matrix: ReadMatrix) -> IsingProblem:
        a, d = build_read_evidence(matrix)
        J = (d - a) / 2
        offset = float((a + d).sum() / 2)
        values = matrix.values.copy()

        def decode(bitstring: str) -> Haplotypes:
            teams = np.array([int(c) for c in bitstring])
            major_a = team_majority(values[teams == 0])
            major_b = team_majority(values[teams == 1])
            copy_a = _fill_from_other(major_a, major_b)
            copy_b = _fill_from_other(major_b, major_a)
            copies = (copy_a, copy_b)
            assignment = np.array([
                SUSPECT
                if _mismatches(row, copies[1 - team]) < _mismatches(row, copies[team])
                else team
                for row, team in zip(values, teams)
            ])
            return Haplotypes(copy_a, copy_b, assignment)

        return IsingProblem(J=J, offset=offset, variable_kind="read", decoder=decode)
