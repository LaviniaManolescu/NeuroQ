"""Teste pentru modulul solvers."""

import numpy as np
import pytest

from neuroq.core import IsingProblem, SelectionParams
from neuroq.encodings import PositionEncoding
from neuroq.execution import IdealSimulator
from neuroq.selection import WindowSelector
from neuroq.solvers import ExactSolver, QAOASolver
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
