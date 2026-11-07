"""Teste pentru modulul evaluation."""

import numpy as np
import pytest

from neuroq.evaluation import accuracy, switch_errors, top_samples

TRUTH = np.array([0, 1, 1, 0, 1])


def test_mirror_counts_as_correct():
    assert accuracy(1 - TRUTH, TRUTH) == 1.0
    assert switch_errors(1 - TRUTH, TRUTH) == 0


def test_one_switch_in_the_middle():
    flipped = np.array([0, 1, 0, 1, 0])  # primele două ca A, restul ca B
    assert switch_errors(flipped, TRUTH) == 1
    assert accuracy(flipped, TRUTH) == pytest.approx(0.6)


def test_top_samples_sorted_by_count():
    top = top_samples({"00": 3, "11": 9, "01": 5}, k=2)
    assert top == [{"bitstring": "11", "count": 9}, {"bitstring": "01", "count": 5}]


def test_unphased_snp_counts_as_error_in_both_orientations():
    copy = np.array([0, 1, -1, 0, 1])
    assert accuracy(copy, TRUTH) == pytest.approx(0.8)
    assert accuracy(np.array([1, 0, -1, 1, 0]), TRUTH) == pytest.approx(0.8)
    assert switch_errors(copy, TRUTH) == 0
