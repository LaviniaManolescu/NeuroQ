"""Teste pentru exportul etapelor către interfață (results/stages.js)."""

import json
import shutil
import subprocess

import pytest

from neuroq.experiments.export_stages import JS_PREFIX, export_stages

INSTANCES = {"notebook", "illumina", "nanopore"}
FIELDS = {"metadata", "evidence", "a", "d", "ising", "cost_history", "best_start",
          "top_samples", "result", "haplotypes", "read_assignment", "metrics", "circuit_png", "variants"}
VARIANT_KEYS = {"qaoa", "qaoa_fixed", "ws_qaoa"}
VARIANT_FIELDS = {"n_qubits", "n_rzz", "ising", "cost_history", "best_start", "circuit_png",
                  "all_samples", "top_samples", "result", "haplotypes", "read_counts", "metrics",
                  "qaoa_seconds", "confidence"}
EXTRA_FIELDS = {"qaoa": set(), "qaoa_fixed": {"fixed_variables"},
                "ws_qaoa": {"c", "c_regularized", "epsilon"}}


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
            for variant in entry["variants"].values():
                variant.pop("qaoa_seconds")
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


# ---- variantele de QAOA (variants) ---------------------------------------------------------------

def _mirror(bits):
    return bits.translate(str.maketrans("01", "10"))


def _energy(ising, bits):
    """Energia unui șir pe o problemă exportată: offset + Σ h·z + Σ J·z·z (doar pentru verificare)."""
    z = [1 - 2 * int(b) for b in bits]
    total = ising["offset"] + sum(h * zi for h, zi in zip(ising["h"], z))
    return total + sum(e["J"] * z[e["i"]] * z[e["j"]] for e in ising["edges"])


def test_every_instance_has_the_three_variants_with_their_fields(exported):
    _, path = exported
    for name, entry in load_stages(path).items():
        assert set(entry["variants"]) == VARIANT_KEYS, name
        for key, variant in entry["variants"].items():
            assert set(variant) == VARIANT_FIELDS | EXTRA_FIELDS[key], (name, key)


def test_qaoa_variant_is_the_instance_level_result(exported):
    _, path = exported
    for name, entry in load_stages(path).items():
        variant = entry["variants"]["qaoa"]
        for field in ("cost_history", "best_start", "top_samples", "haplotypes", "metrics",
                      "circuit_png"):
            assert variant[field] == entry[field], (name, field)
        assert variant["result"]["full_bitstring"] == entry["result"]["bitstring"]
        assert variant["result"]["energy"] == pytest.approx(entry["result"]["energy"])
        assert variant["n_qubits"] == entry["metadata"]["n_qubits"]


def test_variant_sizes_strings_and_samples(exported):
    import inspect

    from neuroq.solvers import QAOASolver

    shots = inspect.signature(QAOASolver.__init__).parameters["shots"].default
    _, path = exported
    for name, entry in load_stages(path).items():
        n = entry["metadata"]["n_qubits"]
        for key, variant in entry["variants"].items():
            expected = n - 1 if key == "qaoa_fixed" else n
            assert variant["n_qubits"] == expected, (name, key)
            samples = variant["all_samples"]
            assert sum(s["count"] for s in samples) == shots, (name, key)
            values = [int(s["bitstring"], 2) for s in samples]
            assert values == sorted(set(values)), (name, key, "ordine binară")
            assert all(len(s["bitstring"]) == expected for s in samples), (name, key)
            by_count = sorted(samples, key=lambda s: (-s["count"], s["bitstring"]))
            assert variant["top_samples"] == by_count[:len(variant["top_samples"])], (name, key)
            chosen = next(s for s in samples if s["bitstring"] == variant["result"]["bitstring"])
            assert chosen["energy"] == pytest.approx(variant["result"]["energy"])
            assert len(variant["result"]["full_bitstring"]) == n


def test_confidence_adds_the_mirror_except_when_a_qubit_is_fixed(exported):
    _, path = exported
    for name, entry in load_stages(path).items():
        for key, variant in entry["variants"].items():
            counts = {s["bitstring"]: s["count"] for s in variant["all_samples"]}
            chosen = variant["result"]["bitstring"]
            hits = counts[chosen] + (0 if key == "qaoa_fixed" else counts.get(_mirror(chosen), 0))
            assert variant["confidence"] == pytest.approx(hits / sum(counts.values())), (name, key)
            assert 0 < variant["confidence"] <= 1


def test_n_rzz_ising_and_read_counts(exported):
    _, path = exported
    for name, entry in load_stages(path).items():
        for key, variant in entry["variants"].items():
            ising = variant["ising"]
            assert variant["n_rzz"] == len(ising["edges"]) > 0, (name, key)
            assert len(ising["h"]) == variant["n_qubits"]
            assert any(ising["h"]) == (key == "qaoa_fixed"), (name, key)  # h apare doar după fixare
            counts = variant["read_counts"]
            assert counts["a"] + counts["b"] + counts["suspect"] == entry["metadata"]["n_reads"]


def test_fixed_variant_reinserts_the_fixed_bit_and_keeps_the_energy(exported):
    _, path = exported
    for name, entry in load_stages(path).items():
        variant = entry["variants"]["qaoa_fixed"]
        fixed = {int(i): b for i, b in variant["fixed_variables"].items()}
        full, short = variant["result"]["full_bitstring"], variant["result"]["bitstring"]
        assert fixed and all(full[i] == str(b) for i, b in fixed.items()), name
        assert "".join(ch for i, ch in enumerate(full) if i not in fixed) == short
        original = {**entry["ising"], "h": [0.0] * entry["metadata"]["n_qubits"]}
        assert _energy(original, full) == pytest.approx(variant["result"]["energy"]), name


def test_warm_start_variant_exports_both_tilts(exported):
    _, path = exported
    for name, entry in load_stages(path).items():
        variant = entry["variants"]["ws_qaoa"]
        c, regularized, eps = variant["c"], variant["c_regularized"], variant["epsilon"]
        assert len(c) == len(regularized) == variant["n_qubits"]
        assert all(0 <= v <= 1 for v in c), name
        assert all(eps <= v <= 1 - eps for v in regularized), name
        assert regularized == pytest.approx([min(max(v, eps), 1 - eps) for v in c]), name


def test_variant_circuit_images_exist(exported):
    out, path = exported
    for name, entry in load_stages(path).items():
        paths = {variant["circuit_png"] for variant in entry["variants"].values()}
        assert len(paths) == 3 and None not in paths, name  # câte o imagine separată pe variantă
        assert all((out / p).is_file() for p in paths), name


def test_notebook_variants_find_cost_one_with_the_same_split_as_01101(exported):
    _, path = exported
    variants = load_stages(path)["notebook"]["variants"]
    for key, variant in variants.items():
        assert variant["result"]["energy"] == pytest.approx(1.0), key
        assert variant["result"]["full_bitstring"] in {"01101", "10010"}, key
        assert set(variant["haplotypes"]["decoded"]) == {"01101", "10010"}, key
    assert variants["qaoa_fixed"]["fixed_variables"] == {"0": 0}
    assert variants["qaoa_fixed"]["result"]["full_bitstring"] == "01101"


def test_a_variant_that_cannot_be_computed_is_left_out_with_its_reason(tmp_path, monkeypatch, capsys):
    from neuroq.experiments import export_stages as module
    from neuroq.experiments.stages import load_case, run_stages

    real = module.run_stages

    def fail_on_fix(case, encoding, **options):
        if options.get("fix"):
            raise ValueError("nu se poate fixa pe această instanță")
        return real(case, encoding, **options)

    monkeypatch.setattr(module, "run_stages", fail_on_fix)
    base = real(load_case("notebook"), "position")
    variants, skipped = module.build_variants("notebook", "position", base, tmp_path)
    assert set(variants) == {"qaoa", "ws_qaoa"}
    assert skipped == [("qaoa_fixed", "nu se poate fixa pe această instanță")]
    assert "variantă omisă pe notebook: qaoa_fixed" in capsys.readouterr().err
