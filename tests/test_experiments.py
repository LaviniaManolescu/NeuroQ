"""Teste pentru modulul experiments: drumul complet pe exemplul din notebook."""

import pytest

from neuroq.experiments.notebook_example import run_notebook_example, write_results
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
