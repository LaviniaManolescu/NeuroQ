"""Teste pentru modulul encodings."""

import itertools

import numpy as np
import pytest

from examples_data import X, Y, notebook_matrix, xy_matrix
from neuroq.encodings import PositionEncoding, ReadEncoding
from neuroq.encodings.position import build_evidence, pairwise_cost
from neuroq.encodings.read import build_read_evidence
from neuroq.solvers import ExactSolver
from neuroq.sources.synthetic import NOTEBOOK_READS


def text(copy) -> str:
    return "".join(str(int(v)) for v in copy)


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


def test_read_evidence_mirrors_position_evidence_on_transpose():
    m = notebook_matrix()
    a_pos, d_pos = build_evidence(m.values.T)
    a_read, d_read = build_read_evidence(m)
    assert np.array_equal(a_pos, a_read) and np.array_equal(d_pos, d_read)


def test_read_encoding_on_notebook_reconstructs_both_copies():
    problem = ReadEncoding().encode(notebook_matrix())
    assert problem.variable_kind == "read" and problem.n_qubits == 9
    haps = problem.decoder(ExactSolver().solve(problem).bitstring)
    assert {text(haps.copy_a), text(haps.copy_b)} == {"01101", "10010"}


@pytest.mark.parametrize("encoding", [PositionEncoding(), ReadEncoding()])
def test_both_encodings_recover_x_and_y_with_zero_cost(encoding):
    problem = encoding.encode(xy_matrix())
    solution = ExactSolver().solve(problem)
    assert solution.energy == pytest.approx(0.0)
    haps = problem.decoder(solution.bitstring)
    assert {text(haps.copy_a), text(haps.copy_b)} == {X, Y}
