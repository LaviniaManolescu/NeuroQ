"""Teste pentru demonstrație: etapele (doar calcul) și scriptul care salvează figurile."""

import matplotlib

matplotlib.use("Agg")

import pytest  # noqa: E402

from neuroq.experiments.demo import main  # noqa: E402
from neuroq.experiments.stages import load_case, resolve_encoding, run_stages  # noqa: E402

EXPECTED = ["01_citiri", "02_a_d", "03_J", "04_graf", "05_circuit", "06_cost_qaoa",
            "07_top10", "08_graf_solutie", "09_haplotipuri"]


def test_run_stages_on_notebook_finds_01101():
    calls = []
    result = run_stages(load_case("notebook"), "position", callback=lambda *c: calls.append(c))
    assert result.problem.energy(result.full_bitstring) == pytest.approx(1.0)
    assert {result.full_bitstring} <= {"01101", "10010"}
    assert result.metrics["accuracy"] == 1.0 and result.other_encoding is None
    assert len(calls) == sum(len(h) for h in result.cost_history)
    assert len(result.top_samples) == 10


def test_run_stages_with_fix_reinserts_the_fixed_bit():
    result = run_stages(load_case("notebook"), "read", fix=True)
    assert result.solved.n_qubits == 8 and len(result.full_bitstring) == 9
    assert result.full_bitstring[0] == "0"


def test_large_instances_only_run_the_small_encoding():
    assert resolve_encoding("illumina", None) == "position"
    assert resolve_encoding("nanopore", None) == "read"
    with pytest.raises(ValueError, match="doar resursele"):
        resolve_encoding("illumina", "citiri")


def test_demo_creates_expected_figures(tmp_path):
    saved = main(["--instance", "notebook", "--out", str(tmp_path)])
    names = sorted(p.stem for p in saved)
    assert names == sorted(f"{e}_pozitii" for e in EXPECTED)
    assert all(p.parent == tmp_path / "demo_notebook" and p.stat().st_size > 0 for p in saved)


def test_run_stages_with_warm_start_draws_ry_circuit():
    result = run_stages(load_case("notebook"), "position", warm=True)
    assert result.warm and not result.fixed
    assert "ry" in result.circuit.count_ops() and "h" not in result.circuit.count_ops()
    assert result.problem.energy(result.full_bitstring) == pytest.approx(1.0)


def test_instance_all_runs_the_three_instances_with_a_summary(tmp_path, capsys):
    from neuroq.experiments.stages import INSTANCES

    saved = main(["--instance", "all", "--out", str(tmp_path)])
    assert len(saved) == 3 * len(EXPECTED)
    for name in INSTANCES:
        folder = tmp_path / f"demo_{name}"
        assert len(list(folder.glob("*.png"))) == len(EXPECTED), name
        assert any("_top10_" in p.name for p in folder.iterdir()), name
    out = capsys.readouterr().out
    assert "Rezumat: 3 instanțe" in out
    for name in INSTANCES:
        assert name in out.split("Rezumat")[1]


def test_instance_all_ignores_encoding_with_a_note(tmp_path, capsys):
    saved = main(["--instance", "all", "--encoding", "citiri", "--out", str(tmp_path)])
    assert len(saved) == 3 * len(EXPECTED)
    assert "--encoding se ignoră" in capsys.readouterr().out
