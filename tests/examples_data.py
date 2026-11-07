"""Exemple mici cu răspuns cunoscut, folosite de mai multe teste."""

import numpy as np

from neuroq.core import MISSING, ReadMatrix
from neuroq.sources.synthetic import NOTEBOOK_READS

X = "100111010"
Y = "011000101"
READ_LENGTH = 3


def matrix_from(values: np.ndarray, technology: str = "test") -> ReadMatrix:
    n_reads, n_snps = values.shape
    return ReadMatrix(values=values, read_ids=[f"r{i + 1}" for i in range(n_reads)],
                      snp_ids=[str(i + 1) for i in range(n_snps)], technology=technology)


def notebook_matrix() -> ReadMatrix:
    return matrix_from(NOTEBOOK_READS.copy(), technology="notebook")


def xy_matrix() -> ReadMatrix:
    """14 citiri fără erori, câte 3 poziții consecutive, alternând X și Y (R1=X[1-3], R2=Y[1-3], ...)."""
    n = len(X)
    rows = []
    for start in range(n - READ_LENGTH + 1):
        for copy in (X, Y):
            row = np.full(n, MISSING)
            row[start:start + READ_LENGTH] = [int(c) for c in copy[start:start + READ_LENGTH]]
            rows.append(row)
    return matrix_from(np.array(rows))
