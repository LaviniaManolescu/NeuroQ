"""Sursa sintetică: exemplul fix din notebook sau un generator simplu de citiri."""

from __future__ import annotations

import numpy as np

from neuroq.core import MISSING, SourceData

# din notebook-ul organizatorilor: perechea ascunsă și cele 9 citiri (READS)
NOTEBOOK_PAIR = np.array([[0, 1, 1, 0, 1], [1, 0, 0, 1, 0]])
NOTEBOOK_READS = np.array([
    [0, 1, 1, -1, -1], [1, 0, 0, -1, -1],
    [-1, 1, 1, 0, -1], [-1, 0, 0, 1, -1],
    [-1, -1, 1, 0, 1], [-1, -1, 0, 1, 0],
    [0, -1, -1, 0, 1], [1, -1, -1, 1, 0],
    [0, 0, -1, -1, -1],
])

SNP_SPACING = 1000  # distanța, în litere, dintre două SNP-uri sintetice


class SyntheticSource:
    """Generează o pereche ascunsă de haplotipuri și citiri din ea, cu seed fix."""

    def __init__(
        self,
        n_snps: int = 5,
        n_reads: int = 9,
        snps_per_read: int = 3,
        error_rate: float = 0.0,
        seed: int = 0,
    ) -> None:
        if not 2 <= snps_per_read <= n_snps:
            raise ValueError("snps_per_read trebuie să fie între 2 și n_snps")
        if not 0.0 <= error_rate <= 1.0:
            raise ValueError("error_rate trebuie să fie între 0 și 1")
        self.n_snps = n_snps
        self.n_reads = n_reads
        self.snps_per_read = snps_per_read
        self.error_rate = error_rate
        self.seed = seed
        self._notebook = False

    @classmethod
    def notebook_example(cls) -> SyntheticSource:
        """Exemplul fix din notebook: 5 SNP-uri, 9 citiri, A=01101, B=10010."""
        source = cls(n_snps=5, n_reads=9)
        source._notebook = True
        return source

    def load(self) -> SourceData:
        if self._notebook:
            values = NOTEBOOK_READS.copy()
            pair = NOTEBOOK_PAIR.copy()
            technology = "notebook"
            read_starts = np.array([self._first_covered(row) for row in values])
        else:
            values, pair, read_starts = self._generate()
            technology = "synthetic"
        n_reads, n_snps = values.shape
        return SourceData(
            snp_positions=np.arange(1, n_snps + 1) * SNP_SPACING,
            ref_alleles=["A"] * n_snps,
            alt_alleles=["G"] * n_snps,
            read_ids=[f"read{i + 1}" for i in range(n_reads)],
            read_starts=read_starts,
            values=values,
            technology=technology,
            true_haplotypes=(pair[0], pair[1]),
        )

    def _generate(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        rng = np.random.default_rng(self.seed)
        copy_a = rng.integers(0, 2, self.n_snps)
        pair = np.array([copy_a, 1 - copy_a])
        values = np.full((self.n_reads, self.n_snps), MISSING)
        starts = rng.integers(0, self.n_snps - self.snps_per_read + 1, self.n_reads)
        copies = rng.integers(0, 2, self.n_reads)
        for r, (start, copy) in enumerate(zip(starts, copies)):
            span = slice(start, start + self.snps_per_read)
            segment = pair[copy, span].copy()
            flips = rng.random(self.snps_per_read) < self.error_rate
            segment[flips] = 1 - segment[flips]
            values[r, span] = segment
        read_starts = (starts + 1) * SNP_SPACING
        return values, pair, read_starts

    @staticmethod
    def _first_covered(row: np.ndarray) -> int:
        return int((np.flatnonzero(row != MISSING)[0] + 1) * SNP_SPACING)
