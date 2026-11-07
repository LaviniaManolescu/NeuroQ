"""Teste pentru contractele din neuroq.core."""

import json
from pathlib import Path

import numpy as np
import pytest

from neuroq.core import IsingProblem, ReadMatrix, ResultRow

SAMPLE_PATH = Path(__file__).resolve().parent.parent / "examples" / "results_sample.json"

# Exemplul din notebook-ul organizatorilor: 5 SNP-uri, h = 0, c0 = 12.5.
# Indicii din notebook sunt de la 1; în cod sunt de la 0.
NOTEBOOK_J = {
    (1, 2): 0.5, (1, 3): 1.0, (1, 4): -1.0, (1, 5): 1.0,
    (2, 3): -2.0, (2, 4): 1.0, (2, 5): 0.0,
    (3, 4): 2.0, (3, 5): -1.0,
    (4, 5): 2.0,
}
NOTEBOOK_OFFSET = 12.5


def notebook_problem(offset: float) -> IsingProblem:
    J = np.zeros((5, 5))
    for (i, j), value in NOTEBOOK_J.items():
        J[i - 1, j - 1] = value
    return IsingProblem(J=J, offset=offset)


def test_energy_without_offset_on_notebook_example():
    assert notebook_problem(0.0).energy("01101") == pytest.approx(-11.5)


def test_cost_with_offset_on_notebook_example():
    assert notebook_problem(NOTEBOOK_OFFSET).energy("01101") == pytest.approx(1.0)


def test_01101_is_a_minimum_on_notebook_example():
    problem = notebook_problem(NOTEBOOK_OFFSET)
    energies = {format(k, "05b"): problem.energy(format(k, "05b")) for k in range(32)}
    assert min(energies.values()) == pytest.approx(energies["01101"])


def test_mirror_symmetry():
    problem = notebook_problem(NOTEBOOK_OFFSET)
    assert problem.energy("10010") == pytest.approx(problem.energy("01101"))


def test_leftmost_bit_is_qubit_zero():
    # Un singur câmp h pe qubitul 0: doar primul caracter din stânga contează.
    problem = IsingProblem(J=np.zeros((3, 3)), h=np.array([1.0, 0.0, 0.0]))
    assert problem.energy("100") == pytest.approx(-1.0)
    assert problem.energy("001") == pytest.approx(1.0)


def test_energy_rejects_bad_bitstrings():
    problem = notebook_problem(0.0)
    with pytest.raises(ValueError):
        problem.energy("0110")
    with pytest.raises(ValueError):
        problem.energy("01a01")


def test_ising_rejects_lower_triangle():
    J = np.zeros((2, 2))
    J[1, 0] = 1.0
    with pytest.raises(ValueError):
        IsingProblem(J=J)


def test_read_matrix_rejects_bad_values():
    with pytest.raises(ValueError):
        ReadMatrix(values=np.array([[0, 2]]), read_ids=["r1"], snp_ids=["s1", "s2"],
                   technology="illumina")


def test_sample_results_file_matches_result_row():
    data = json.loads(SAMPLE_PATH.read_text(encoding="utf-8"))
    assert data["is_sample"] is True
    assert data["schema_version"] == "0.1"
    assert len(data["rows"]) == 3

    for raw in data["rows"]:
        row = ResultRow.from_dict(raw)
        assert row.to_dict() == raw

        assert all(isinstance(v, (int, float)) and not isinstance(v, bool)
                   for v in row.metrics.values())

        haps = row.details["haplotypes"]
        assert len(haps["copy_a"]) == len(haps["copy_b"]) == row.instance["n_snps"]
        assert len(haps["read_assignment"]) == row.instance["n_reads"]
        assert set(haps["read_assignment"]) <= {-1, 0, 1}

        counts = [s["count"] for s in row.details["top_samples"]]
        assert counts == sorted(counts, reverse=True)
        for sample in row.details["top_samples"]:
            assert len(sample["bitstring"]) == row.resources.n_qubits
