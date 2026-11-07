"""Baseline clasic: lipirea lacomă a citirilor, direct pe ReadMatrix."""

from __future__ import annotations

from collections import deque

import numpy as np

from neuroq.core import MISSING, Haplotypes, ReadMatrix

UNKNOWN = -1
COPY_A, COPY_B = 0, 1


class GreedyPhaser:
    """Pune citirile pe rând în copia cu care se potrivesc mai bine; depinde de ordine.

    Nu moștenește Solver, pentru că lucrează pe ReadMatrix, nu pe IsingProblem.
    """

    def phase(self, m: ReadMatrix, order_seed: int) -> Haplotypes:
        values = m.values
        n_reads, n_snps = values.shape
        copies = np.full((2, n_snps), UNKNOWN)
        assignment = np.full(n_reads, UNKNOWN)
        if n_reads == 0:
            return Haplotypes(copies[COPY_A], copies[COPY_B], assignment)

        order = np.random.default_rng(order_seed).permutation(n_reads)
        queue = deque(int(r) for r in order)
        deferred = 0
        first = True
        while queue:
            r = queue.popleft()
            row = values[r]
            overlap = (row != MISSING) & (copies[COPY_A] != UNKNOWN)
            if not first and not overlap.any() and deferred < len(queue):
                queue.append(r)  # fără coloane comune: se amână la coadă
                deferred += 1
                continue
            deferred = 0
            first = False
            team = self._closest_copy(row, copies, overlap)
            assignment[r] = team
            self._extend(row, copies, team)
        return Haplotypes(copies[COPY_A], copies[COPY_B], assignment)

    @staticmethod
    def _closest_copy(row: np.ndarray, copies: np.ndarray, overlap: np.ndarray) -> int:
        """Copia cu mai puține nepotriviri; egalitate (inclusiv fără suprapunere) → A."""
        miss_a = np.count_nonzero(row[overlap] != copies[COPY_A, overlap])
        miss_b = np.count_nonzero(row[overlap] != copies[COPY_B, overlap])
        return COPY_B if miss_b < miss_a else COPY_A

    @staticmethod
    def _extend(row: np.ndarray, copies: np.ndarray, team: int) -> None:
        """Completează coloanele încă necunoscute cu citirea, iar cealaltă copie cu complementul."""
        new = (row != MISSING) & (copies[team] == UNKNOWN)
        copies[team, new] = row[new]
        copies[1 - team, new] = 1 - row[new]
