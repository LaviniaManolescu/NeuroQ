"""Teste pentru modulul selection."""

import numpy as np

from neuroq.core import SelectionParams
from neuroq.selection import Downsampler, WindowSelector
from neuroq.sources import SyntheticSource


def test_full_window_keeps_all_notebook_reads():
    data = SyntheticSource.notebook_example().load()
    matrix = WindowSelector().select(data, SelectionParams(first_snp=0, n_snps=5))
    assert matrix.values.shape == (9, 5)
    assert np.array_equal(matrix.values, data.values)


def test_window_drops_reads_with_fewer_than_two_snps():
    data = SyntheticSource.notebook_example().load()
    matrix = WindowSelector().select(data, SelectionParams(first_snp=0, n_snps=2))
    # doar citirile 1, 2 și 9 acoperă SNP-urile 1 și 2
    assert matrix.read_ids == ["read1", "read2", "read9"]
    assert "".join(map(str, matrix.true_haplotypes[0])) == "01"


def test_downsampling_is_reproducible():
    data = SyntheticSource(n_snps=6, n_reads=30, seed=1).load()
    params = SelectionParams(first_snp=0, n_snps=6, target_reads=10, seed=5)
    first = WindowSelector().select(data, params)
    second = WindowSelector().select(data, params)
    assert len(first.read_ids) == 10
    assert first.read_ids == second.read_ids
    other = Downsampler(10, seed=6).apply(WindowSelector().select(
        data, SelectionParams(first_snp=0, n_snps=6)))
    assert other.read_ids != first.read_ids
