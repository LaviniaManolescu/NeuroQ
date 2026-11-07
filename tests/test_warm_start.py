"""Teste pentru relaxarea continuă folosită de warm-start QAOA."""

import itertools

import numpy as np
import pytest

from examples_data import notebook_matrix
from neuroq.encodings import PositionEncoding, ReadEncoding
from neuroq.solvers import ExactSolver
from neuroq.solvers.warm_start import regularize, relax, relaxed_energy


def notebook_problem():
    return PositionEncoding().encode(notebook_matrix())


def test_relaxed_energy_equals_ising_energy_on_bitstrings():
    for problem in (notebook_problem(), ReadEncoding().encode(notebook_matrix())):
        for bits in itertools.product([0, 1], repeat=problem.n_qubits):
            text = "".join(map(str, bits))
            assert relaxed_energy(problem, np.array(bits)) == pytest.approx(problem.energy(text))


def test_relax_is_reproducible_and_breaks_the_mirror():
    problem = notebook_problem()
    first, second = relax(problem, seed=71), relax(problem, seed=71)
    assert np.array_equal(first, second)
    assert ((first >= 0) & (first <= 1)).all()
    assert not np.allclose(first, 0.5)  # alege o orientare, nu rămâne în punctul simetric
    rounded = "".join(str(int(v > 0.5)) for v in first)
    assert problem.energy(rounded) == pytest.approx(ExactSolver().solve(problem).energy)


def test_regularize_clips_to_epsilon():
    assert regularize(np.array([0.0, 0.4, 1.0]), 0.25).tolist() == [0.25, 0.4, 0.75]
    with pytest.raises(ValueError):
        regularize(np.array([0.5]), 0.6)


def test_warm_start_circuit_uses_ry_mixer_and_no_hadamard():
    from neuroq.solvers.qaoa import make_warm_start_qaoa

    problem = notebook_problem()
    c = regularize(relax(problem, 71))
    qc, _, _ = make_warm_start_qaoa(problem.J, problem.h, c, p=1)
    ops = qc.count_ops()
    assert "h" not in ops and "rx" not in ops
    assert ops["rzz"] == np.count_nonzero(problem.J) and ops["ry"] == 3 * problem.n_qubits


def test_initial_state_matches_relaxed_probabilities():
    from qiskit.quantum_info import Statevector

    from neuroq.solvers.qaoa import make_warm_start_qaoa

    problem = notebook_problem()
    c = regularize(relax(problem, 71))
    qc, theta, _ = make_warm_start_qaoa(problem.J, problem.h, c, p=1)
    probs = Statevector.from_instruction(qc.assign_parameters([0.0, 0.0])).probabilities()
    expected = np.array([np.prod([c[i] if (k >> i) & 1 else 1 - c[i] for i in range(len(c))])
                         for k in range(2 ** len(c))])
    assert np.allclose(probs, expected)


def test_ws_qaoa_finds_cost_one_on_notebook():
    from neuroq.execution import IdealSimulator
    from neuroq.solvers import QAOASolver

    solution = QAOASolver(IdealSimulator(), warm_start=True).solve(notebook_problem())
    assert solution.energy == pytest.approx(1.0)
    assert solution.metadata["warm_start"]["epsilon"] == 0.25
    assert len(solution.metadata["warm_start"]["c"]) == 5
