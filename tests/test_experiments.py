"""Teste pentru modulul experiments: drumul complet și comparația pe ferestre mici."""

import json

import pytest

from neuroq.experiments.notebook_example import run_notebook_example, write_results
from neuroq.encodings import ReadEncoding
from neuroq.experiments.resource_comparison import (
    illumina_case, nanopore_case, notebook_case, resources_row, run_comparison, write_outputs,
)
from schema_check import validate_results_file


def test_notebook_pipeline_finds_01101_and_writes_valid_results(tmp_path):
    rows = run_notebook_example()
    path = tmp_path / "results.json"
    write_results(rows, path)
    loaded = validate_results_file(path, is_sample=False)

    exact, qaoa = loaded
    assert exact.method == "exact" and qaoa.method == "qaoa"
    assert exact.details["haplotypes"]["copy_a"] == "01101"
    assert exact.metrics["pair_cost"] == pytest.approx(1.0)
    assert exact.metrics["accuracy"] == 1.0
    assert qaoa.metrics["pair_cost"] == pytest.approx(1.0)
    assert qaoa.details["haplotypes"]["copy_a"] in {"01101", "10010"}
    assert 0 < qaoa.metrics["optimum_frequency"] <= 1


def test_comparison_on_notebook_writes_valid_json_and_js(tmp_path):
    rows = run_comparison([notebook_case()])
    json_path, js_path = write_outputs(rows, tmp_path)
    loaded = validate_results_file(json_path, is_sample=False)

    assert len(loaded) == 9  # greedy + 2 codificări × (exact, classical, qaoa, qaoa fixat)
    greedy = [r for r in loaded if r.method == "greedy"]
    assert len(greedy) == 1 and greedy[0].encoding is None and greedy[0].resources is None
    exact_pos = next(r for r in loaded if r.method == "exact" and r.encoding == "position")
    assert exact_pos.details["haplotypes"]["copy_a"] == "01101"
    assert exact_pos.metrics["pair_cost"] == pytest.approx(1.0)
    fixed = [r for r in loaded if r.transforms == ["fix_variable"]]
    assert {r.resources.n_qubits for r in fixed} == {4, 8}

    js = js_path.read_text(encoding="utf-8")
    assert js.startswith("window.RESULTS = ") and js.rstrip().endswith(";")
    assert json.loads(js[len("window.RESULTS = "):].rstrip().rstrip(";")) == json.loads(
        json_path.read_text(encoding="utf-8"))


def test_resources_row_counts_rzz_and_passes_schema(tmp_path):
    case = notebook_case()
    row = resources_row(case, ReadEncoding(), {"skipped_reason": "test"})
    problem = ReadEncoding().encode(case.matrix)
    assert row.method == "resources" and row.encoding == "read"
    assert row.metrics["n_variables"] == 9
    assert row.metrics["rzz_per_layer"] == row.metrics["n_couplings"] == int((problem.J != 0).sum())
    assert row.resources.two_qubit_gates_before == row.metrics["rzz_per_layer"]
    json_path, _ = write_outputs([row], tmp_path)
    assert validate_results_file(json_path, is_sample=False)[0].method == "resources"


def test_synthetic_cases_split_full_and_resource_encodings():
    illumina, nanopore = illumina_case(), nanopore_case()
    assert illumina.full == ("position",) and set(illumina.resources_only) == {"read"}
    assert illumina.resources_only["read"]["previous_measurement"]["qaoa_seconds"] == 1246
    assert nanopore.full == ("read",) and set(nanopore.resources_only) == {"position"}
    assert illumina.matrix.values.shape == (24, 10)
    assert nanopore.matrix.values.shape == (9, 18)
