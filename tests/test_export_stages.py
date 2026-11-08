"""Teste pentru exportul etapelor către interfață (results/stages.js)."""

import json
import shutil
import subprocess

import pytest

from neuroq.experiments.export_stages import JS_PREFIX, export_stages

INSTANCES = {"notebook", "illumina", "nanopore"}
FIELDS = {"metadata", "evidence", "a", "d", "ising", "cost_history", "best_start",
          "top_samples", "all_samples", "result", "haplotypes", "read_assignment", "metrics", "circuit_png"}


def _reject_constant(name):
    raise ValueError(f"constantă nepermisă în JSON: {name}")


def load_stages(path):
    text = path.read_text(encoding="utf-8")
    assert text.startswith(JS_PREFIX) and text.rstrip().endswith(";")
    return json.loads(text[len(JS_PREFIX):].rstrip().rstrip(";"), parse_constant=_reject_constant)


@pytest.fixture(scope="module")
def exported(tmp_path_factory):
    out = tmp_path_factory.mktemp("stages")
    return out, export_stages(out)


def test_stages_js_has_all_instances_and_fields(exported):
    _, path = exported
    stages = load_stages(path)
    assert set(stages) == INSTANCES
    for name, entry in stages.items():
        assert set(entry) == FIELDS, name
        assert {"instance", "encoding", "n_qubits", "n_reads", "n_snps",
                "qaoa_seconds"} <= set(entry["metadata"])
        assert set(entry["metrics"]) == {"accuracy", "switch_errors", "mec_cost"}
        ev = entry["evidence"]
        assert len(ev["values"]) == entry["metadata"]["n_reads"] == len(entry["read_assignment"])
        assert all(v in (0, 1, None) for row in ev["values"] for v in row)
        for read in entry["read_assignment"]:
            assert read["team"] in (-1, 0, 1)
            assert len(read["snps"]) == len(read["values"]) == len(read["errors"])
    assert stages["notebook"]["metadata"]["encoding"] == "position"
    assert stages["nanopore"]["metadata"]["encoding"] == "read"


def test_notebook_entry_matches_the_lab(exported):
    _, path = exported
    nb = load_stages(path)["notebook"]
    assert nb["result"]["bitstring"] in {"01101", "10010"}
    assert nb["result"]["energy"] == pytest.approx(1.0)
    assert set(nb["haplotypes"]["decoded"]) == {"01101", "10010"}
    assert nb["ising"]["offset"] == pytest.approx(12.5) and len(nb["ising"]["edges"]) == 9


def test_all_samples_is_the_full_distribution(exported):
    import inspect

    from neuroq.solvers import QAOASolver

    shots = inspect.signature(QAOASolver.__init__).parameters["shots"].default
    _, path = exported
    for name, entry in load_stages(path).items():
        samples = entry["all_samples"]
        assert sum(s["count"] for s in samples) == shots, name
        values = [int(s["bitstring"], 2) for s in samples]
        assert values == sorted(set(values)), f"{name}: ordine binară"
        assert all(s["count"] >= 1 and len(s["bitstring"]) == entry["metadata"]["n_qubits"]
                   for s in samples)
        by_count = sorted(samples, key=lambda s: (-s["count"], s["bitstring"]))
        assert entry["top_samples"] == by_count[:len(entry["top_samples"])], name
        chosen = next(s for s in samples if s["bitstring"] == entry["result"]["bitstring"])
        assert chosen["energy"] == pytest.approx(entry["result"]["energy"])


def test_circuit_images_exist(exported):
    out, path = exported
    for entry in load_stages(path).values():
        assert entry["circuit_png"] and (out / entry["circuit_png"]).is_file()


def test_export_is_reproducible(exported, tmp_path):
    _, first = exported
    second = export_stages(tmp_path)

    def without_timing(path):
        stages = load_stages(path)
        for entry in stages.values():
            entry["metadata"].pop("qaoa_seconds")
        return stages

    assert without_timing(first) == without_timing(second)


def test_node_can_parse_stages_js(exported):
    node = shutil.which("node")
    if node is None:
        pytest.skip("node nu este instalat")
    _, path = exported
    subprocess.run([node, "--check", str(path)], check=True)


def test_inputs_round_trip_with_instance_metadata(tmp_path):
    import numpy as np

    from neuroq.experiments.export_stages import export_inputs
    from neuroq.experiments.stages import load_case
    from neuroq.sources import load_source_json

    for path in export_inputs(tmp_path):
        name = path.stem
        raw = json.loads(path.read_text(encoding="utf-8"))
        assert raw["metadata"]["instance"] == name
        loaded = load_source_json(path)
        assert np.array_equal(loaded.values, load_case(name).matrix.values)
