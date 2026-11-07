"""Fereastra de SNP-uri și downsampling-ul citirilor."""

from __future__ import annotations

import numpy as np

from neuroq.core import MISSING, ReadMatrix, SelectionParams, SourceData

MIN_SNPS_PER_READ = 2


class Downsampler:
    """Păstrează aleator, cu seed fix, un număr țintă de citiri."""

    def __init__(self, target_reads: int, seed: int) -> None:
        if target_reads < 1:
            raise ValueError("target_reads trebuie să fie cel puțin 1")
        self.target_reads = target_reads
        self.seed = seed

    def apply(self, matrix: ReadMatrix) -> ReadMatrix:
        n_reads = matrix.values.shape[0]
        if n_reads <= self.target_reads:
            return matrix
        rng = np.random.default_rng(self.seed)
        keep = np.sort(rng.choice(n_reads, size=self.target_reads, replace=False))
        return ReadMatrix(
            values=matrix.values[keep],
            read_ids=[matrix.read_ids[i] for i in keep],
            snp_ids=list(matrix.snp_ids),
            technology=matrix.technology,
            special_labels=dict(matrix.special_labels),
            true_haplotypes=matrix.true_haplotypes,
        )


class WindowSelector:
    """Taie o fereastră de SNP-uri și păstrează citirile care acoperă cel puțin două."""

    def select(self, data: SourceData, params: SelectionParams) -> ReadMatrix:
        n_total = data.values.shape[1]
        start, stop = params.first_snp, params.first_snp + params.n_snps
        if params.n_snps < 2 or start < 0 or stop > n_total:
            raise ValueError(f"fereastra [{start}, {stop}) nu încape în {n_total} SNP-uri")

        window = data.values[:, start:stop]
        keep = np.flatnonzero((window != MISSING).sum(axis=1) >= MIN_SNPS_PER_READ)
        truth = None
        if data.true_haplotypes is not None:
            truth = (data.true_haplotypes[0][start:stop], data.true_haplotypes[1][start:stop])

        matrix = ReadMatrix(
            values=window[keep],
            read_ids=[data.read_ids[i] for i in keep],
            snp_ids=[str(p) for p in data.snp_positions[start:stop]],
            technology=data.technology,
            special_labels={
                name: col - start
                for name, col in data.special_labels.items()
                if start <= col < stop
            },
            true_haplotypes=truth,
        )
        if params.target_reads is not None:
            matrix = Downsampler(params.target_reads, params.seed).apply(matrix)
        return matrix
