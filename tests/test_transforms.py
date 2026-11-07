"""Teste pentru modulul transforms: fixarea unui qubit."""

import itertools

import pytest

from examples_data import notebook_matrix
from neuroq.encodings import PositionEncoding, ReadEncoding
from neuroq.execution import IdealSimulator
from neuroq.solvers import ExactSolver, QAOASolver
from neuroq.solvers.qaoa import make_qaoa
from neuroq.transforms import fix_variable


def notebook_problem():
    return PositionEncoding().encode(notebook_matrix())


@pytest.mark.parametrize("index,value", [(0, 0), (0, 1), (2, 0), (4, 1)])
def test_reduced_energy_equals_original_with_bit_reinserted(index, value):
    original = notebook_problem()
    reduced = fix_variable(original, index, value)
    assert reduced.n_qubits == 4
    for bits in itertools.product("01", repeat=4):
        s = "".join(bits)
        full = s[:index] + str(value) + s[index:]
        assert reduced.energy(s) == pytest.approx(original.energy(full))


def test_fixed_problem_has_a_single_optimum():
    reduced = fix_variable(notebook_problem())
    solution = ExactSolver().solve(reduced)
    assert solution.metadata["optimal_bitstrings"] == ["1101"]
    assert solution.energy == pytest.approx(1.0)
    assert reduced.fixed_variables == {0: 0}
    copy_a = "".join(map(str, reduced.decoder(solution.bitstring).copy_a))
    assert copy_a == "01101"


def test_chained_fixes_record_original_indices():
    twice = fix_variable(fix_variable(notebook_problem(), 0, 0), 1, 1)
    assert twice.fixed_variables == {0: 0, 2: 1}
    assert twice.energy("101") == pytest.approx(notebook_problem().energy("01101"))


def test_fix_removes_only_the_mirror_on_read_encoding():
    # Citirea 9 contrazice datele și poate sta în oricare echipă cu același cost:
    # 4 optime fără fixare (2 × oglindă), 2 după fixare. Fixarea scoate doar oglinda.
    original = ReadEncoding().encode(notebook_matrix())
    reduced = fix_variable(original)
    full = ExactSolver().solve(original)
    fixed = ExactSolver().solve(reduced)
    assert fixed.energy == pytest.approx(full.energy)
    assert len(fixed.metadata["optimal_bitstrings"]) * 2 == len(full.metadata["optimal_bitstrings"])
    for bits in fixed.metadata["optimal_bitstrings"]:
        assert "".join(map(str, reduced.decoder(bits).copy_a)) == "01101"


def test_qaoa_runs_on_fixed_problem_with_rz_gates():
    reduced = fix_variable(notebook_problem())
    qc, _, _ = make_qaoa(reduced.J, reduced.h, p=1)
    assert qc.count_ops().get("rz", 0) == sum(reduced.h != 0) > 0
    solution = QAOASolver(IdealSimulator(), shots=2048, seed=71).solve(reduced)
    assert sum(solution.samples.values()) == 2048
    assert all(len(b) == 4 for b in solution.samples)
    assert solution.energy == pytest.approx(1.0)  # stabil pe seed-urile 71–75
