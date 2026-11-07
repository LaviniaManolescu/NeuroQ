"""Teste pentru modulul sources."""

import numpy as np
import pytest

from neuroq.sources import SyntheticSource
from neuroq.sources.synthetic import ILLUMINA, NANOPORE, TechnologyProfile


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



def _measured_error(source) -> float:
    """Erori estimate: pe fiecare citire, nepotrivirile cu copia cea mai apropiată."""
    copy_a, copy_b = source.true_haplotypes
    seen = source.values != -1
    errors = sum(min((row[s] != copy_a[s]).sum(), (row[s] != copy_b[s]).sum())
                 for row, s in zip(source.values, seen))
    return errors / seen.sum()


def test_illumina_profile_statistics():
    source = SyntheticSource.from_profile(ILLUMINA, n_snps=2000, seed=11).load()
    snps_per_read = (source.values != -1).sum(axis=1)
    assert 2.0 <= snps_per_read.mean() <= 3.0
    assert 0.0003 <= _measured_error(source) <= 0.003


def test_nanopore_profile_statistics():
    source = SyntheticSource.from_profile(NANOPORE, n_snps=300, seed=11).load()
    assert 0.04 <= _measured_error(source) <= 0.06
    lengths = NANOPORE.read_lengths(np.random.default_rng(0), 20_000)
    assert 10_000 <= np.median(lengths) <= 50_000
    assert lengths.min() >= 1_000 and lengths.max() <= 100_000
    assert (source.values != -1).sum(axis=1).mean() > 10  # citiri lungi, multe SNP-uri


@pytest.mark.parametrize("profile", [ILLUMINA, NANOPORE], ids=["illumina", "nanopore"])
def test_profile_is_reproducible(profile):
    first = SyntheticSource.from_profile(profile, n_snps=40, seed=3).load()
    second = SyntheticSource.from_profile(profile, n_snps=40, seed=3).load()
    assert np.array_equal(first.values, second.values)
    assert np.array_equal(first.read_starts, second.read_starts)


@pytest.mark.parametrize("profile", [ILLUMINA, NANOPORE], ids=["illumina", "nanopore"])
def test_truth_is_consistent_with_error_free_reads(profile):
    clean = TechnologyProfile(**{**profile.__dict__, "error_rate": 0.0})
    source = SyntheticSource.from_profile(clean, n_snps=40, seed=5).load()
    copy_a, copy_b = source.true_haplotypes
    assert np.array_equal(copy_b, 1 - copy_a)
    for row in source.values:
        seen = row != -1
        assert (row[seen] == copy_a[seen]).all() or (row[seen] == copy_b[seen]).all()


def test_truth_depends_only_on_n_snps_and_seed():
    illumina = SyntheticSource.from_profile(ILLUMINA, n_snps=30, seed=9).load()
    nanopore = SyntheticSource.from_profile(NANOPORE, n_snps=30, seed=9).load()
    for x, y in zip(illumina.true_haplotypes, nanopore.true_haplotypes):
        assert np.array_equal(x, y)
    assert np.array_equal(illumina.snp_positions, nanopore.snp_positions)
    assert illumina.values.shape != nanopore.values.shape  # citirile diferă
