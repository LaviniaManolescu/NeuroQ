"""Teste pentru modulul sources."""

import numpy as np

from neuroq.sources import SyntheticSource


def test_notebook_example_matches_notebook():
    data = SyntheticSource.notebook_example().load()
    assert data.values.shape == (9, 5)
    assert "".join(map(str, data.true_haplotypes[0])) == "01101"
    assert "".join(map(str, data.true_haplotypes[1])) == "10010"


def test_generator_is_reproducible_and_clean_reads_match_truth():
    a = SyntheticSource(n_snps=8, n_reads=12, snps_per_read=3, seed=4).load()
    b = SyntheticSource(n_snps=8, n_reads=12, snps_per_read=3, seed=4).load()
    assert np.array_equal(a.values, b.values)
    assert ((a.values != -1).sum(axis=1) == 3).all()
    copy_a, copy_b = a.true_haplotypes
    assert np.array_equal(copy_b, 1 - copy_a)
    for row in a.values:
        seen = row != -1
        assert (row[seen] == copy_a[seen]).all() or (row[seen] == copy_b[seen]).all()
