"""Teste pentru modulul solvers."""

import numpy as np
import pytest

from examples_data import X, Y, notebook_matrix, xy_matrix
from neuroq.core import IsingProblem, SelectionParams
from neuroq.encodings import PositionEncoding
from neuroq.execution import IdealSimulator
from neuroq.selection import WindowSelector
from neuroq.solvers import ExactSolver, GreedyPhaser, QAOASolver
from neuroq.solvers.qaoa import make_qaoa
from neuroq.sources import SyntheticSource


def notebook_problem() -> IsingProblem:
    data = SyntheticSource.notebook_example().load()
    matrix = WindowSelector().select(data, SelectionParams(first_snp=0, n_snps=5))
    return PositionEncoding().encode(matrix)


def test_exact_finds_both_mirror_optima_with_cost_one():
    solution = ExactSolver().solve(notebook_problem())
    assert solution.bitstring == "01101"
    assert solution.energy == pytest.approx(1.0)
    assert set(solution.metadata["optimal_bitstrings"]) == {"01101", "10010"}


def test_exact_refuses_too_many_qubits():
    with pytest.raises(ValueError):
        ExactSolver(max_qubits=4).solve(notebook_problem())


def test_qaoa_best_observed_cost_is_one():
    solution = QAOASolver(IdealSimulator()).solve(notebook_problem())
    assert solution.energy == pytest.approx(1.0)
    assert solution.bitstring in {"01101", "10010"}
    assert sum(solution.samples.values()) == 2048


def test_circuit_matches_notebook_shape_and_adds_rz_for_h():
    problem = notebook_problem()
    qc, theta, _ = make_qaoa(problem.J, problem.h, p=1)
    ops = qc.count_ops()
    assert ops.get("rzz", 0) == np.count_nonzero(problem.J)
    assert ops.get("rx", 0) == 5 and ops.get("h", 0) == 5 and "rz" not in ops
    qc_h, _, _ = make_qaoa(problem.J, np.array([0.5, 0, 0, 0, 0]), p=1)
    assert qc_h.count_ops().get("rz", 0) == 1


def _text(copy) -> str:
    return "".join(str(int(v)) for v in copy)


def test_greedy_returns_valid_haplotypes():
    m = notebook_matrix()
    haps = GreedyPhaser().phase(m, order_seed=0)
    assert haps.copy_a.shape == haps.copy_b.shape == (5,)
    assert len(haps.read_assignment) == 9
    assert set(haps.read_assignment.tolist()) <= {0, 1}
    known = haps.copy_a != -1
    assert np.array_equal(haps.copy_b[known], 1 - haps.copy_a[known])


def test_greedy_is_reproducible_for_the_same_seed():
    first = GreedyPhaser().phase(notebook_matrix(), order_seed=3)
    second = GreedyPhaser().phase(notebook_matrix(), order_seed=3)
    assert np.array_equal(first.copy_a, second.copy_a)
    assert np.array_equal(first.read_assignment, second.read_assignment)


@pytest.mark.parametrize("seed", range(10))
def test_greedy_recovers_x_and_y_on_clean_reads(seed):
    haps = GreedyPhaser().phase(xy_matrix(), order_seed=seed)
    assert {_text(haps.copy_a), _text(haps.copy_b)} == {X, Y}
