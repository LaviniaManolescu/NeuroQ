"""Teste pentru modulul encodings."""

import itertools

import numpy as np
import pytest

from neuroq.core import ReadMatrix
from neuroq.encodings import PositionEncoding
from neuroq.encodings.position import build_evidence, pairwise_cost
from neuroq.sources.synthetic import NOTEBOOK_READS


def notebook_matrix() -> ReadMatrix:
    n_reads, n_snps = NOTEBOOK_READS.shape
    return ReadMatrix(values=NOTEBOOK_READS, read_ids=[f"r{i}" for i in range(n_reads)],
                      snp_ids=[str(i) for i in range(n_snps)], technology="notebook")


def test_tiny_evidence_from_notebook():
    ta, td = build_evidence(np.array([[0, 1, 1], [1, 0, 0], [0, 0, -1]]))
    assert np.array_equal(ta, [[0, 1, 0], [0, 0, 2], [0, 0, 0]])
    assert np.array_equal(td, [[0, 2, 2], [0, 0, 0], [0, 0, 0]])


def test_energy_equals_pairwise_cost_on_all_32_strings():
    problem = PositionEncoding().encode(notebook_matrix())
    a, d = build_evidence(NOTEBOOK_READS)
    for x in itertools.product([0, 1], repeat=5):
        bits = "".join(map(str, x))
        assert problem.energy(bits) == pytest.approx(pairwise_cost(x, a, d))


def test_notebook_offset_and_cost():
    problem = PositionEncoding().encode(notebook_matrix())
    assert problem.offset == pytest.approx(12.5)
    assert problem.energy("01101") == pytest.approx(1.0)


def test_decoder_gives_copies_and_flags_contradicting_read():
    haps = PositionEncoding().encode(notebook_matrix()).decoder("01101")
    assert "".join(map(str, haps.copy_a)) == "01101"
    assert "".join(map(str, haps.copy_b)) == "10010"
    assert haps.read_assignment.tolist() == [0, 1, 0, 1, 0, 1, 0, 1, -1]
