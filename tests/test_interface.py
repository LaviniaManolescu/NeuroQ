"""Teste pentru interfață: tot ce referă GenoMEtric.html există, iar datele se citesc fără fetch."""

import json
import re
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
INTERFACE = ROOT / "neuroq" / "interface"
PAGE = INTERFACE / "GenoMEtric.html"
RESULTS = ROOT / "results"
REFERENCE = re.compile(r'<(?:script|link)[^>]*?(?:src|href)="([^"]+)"')


def page_references() -> list[str]:
    return REFERENCE.findall(PAGE.read_text(encoding="utf-8"))


def own_scripts() -> list[Path]:
    return [p for p in INTERFACE.rglob("*.js")]


def test_every_referenced_file_exists():
    refs = page_references()
    assert refs, "GenoMEtric.html nu referă niciun fișier"
    missing = [r for r in refs if not (INTERFACE / r).resolve().is_file()]
    assert not missing, f"lipsesc: {missing}"


def test_page_loads_results_and_stages_from_results_folder():
    resolved = {(INTERFACE / r).resolve() for r in page_references()}
    assert (RESULTS / "results.js").resolve() in resolved
    assert (RESULTS / "stages.js").resolve() in resolved


def test_every_circuit_image_in_stages_exists():
    text = (RESULTS / "stages.js").read_text(encoding="utf-8")
    stages = json.loads(text[len("window.STAGES = "):].rstrip().rstrip(";"))
    for name, entry in stages.items():
        if entry["circuit_png"]:
            assert (RESULTS / entry["circuit_png"]).is_file(), name


def test_committed_stages_js_has_all_samples():
    text = (RESULTS / "stages.js").read_text(encoding="utf-8")
    stages = json.loads(text[len("window.STAGES = "):].rstrip().rstrip(";"))
    for name, entry in stages.items():
        assert sum(s["count"] for s in entry["all_samples"]) == 2048, name


def test_scripts_work_from_file_urls():
    # fără fetch și fără resurse externe: pagina trebuie să meargă deschisă din fișier
    for script in own_scripts():
        source = script.read_text(encoding="utf-8")
        assert "fetch(" not in source, script.name
        assert not re.search(r"https?://(?!www\.w3\.org/2000/svg)", source), script.name


def test_scripts_add_only_the_neuroq_global():
    for script in own_scripts():
        source = script.read_text(encoding="utf-8")
        assert source.lstrip().startswith(("//", "(function")), script.name
        assert "window.NeuroQ" in source, script.name


def test_node_check_scripts():
    node = shutil.which("node")
    if node is None:
        pytest.skip("node nu este instalat")
    for script in own_scripts():
        subprocess.run([node, "--check", str(script)], check=True)


# ---- Histograma etapei 5 pe cazuri limită, rulată în browserul headless -------------------------
# Datele sunt construite aici, nu din results/; JS-ul real rulează în Edge/Chrome (node lipsește).

import html as html_lib
import random

from neuroq.interface.screenshots import find_browser

HARNESS = ROOT / "tests" / "histogram_harness.html"
READS = 4


def _mirror(bits: str) -> str:
    return bits.translate(str.maketrans("01", "10"))


def _stage(samples: list[tuple[str, int]], chosen: str) -> dict:
    """Etapă minimă, cu exact câmpurile citite de views.top (nu vine din datele reale)."""
    n_bits = len(chosen)
    reads = [{"read_id": f"r{i}", "team": 0, "snps": [i % n_bits], "values": [i % 2],
              "errors": [False]} for i in range(READS)]
    return {
        "metadata": {"instance": "test", "encoding": "position", "n_reads": READS,
                     "n_snps": n_bits, "n_qubits": n_bits},
        "evidence": {"snp_ids": [str(1000 + i) for i in range(n_bits)]},
        "ising": {"edges": [], "offset": 0.0},
        "result": {"bitstring": chosen, "energy": 0.0},
        "read_assignment": reads,
        "all_samples": [{"bitstring": b, "count": c, "energy": float(c % 7)} for b, c in samples],
    }


def _edge_cases() -> dict[str, tuple[list[tuple[str, int]], str]]:
    rng = random.Random(20261108)

    def strings(n_bits: int, how_many: int, first_bit: str = "") -> list[str]:
        pool = [format(k, f"0{n_bits}b") for k in range(2 ** n_bits)]
        return rng.sample([b for b in pool if b.startswith(first_bit)], how_many)

    many = strings(11, 320)                                  # peste 300 de șiruri distincte, toate rare
    many_samples = [(b, rng.randint(3, 9)) for b in many]
    frequent = strings(10, 90, first_bit="0")                # 90 de șiruri, fiecare ≥ pragul de rar
    one_side = strings(6, 8, first_bit="0")                  # oglinzile au toate 0 măsurători
    return {
        "a_doua_siruri_oglinzi": ([("0110", 1000), ("1001", 1048)], "1001"),
        "a_doua_siruri_fara_oglinzi": ([("0110", 700), ("0001", 300)], "0110"),
        "b_peste_300_de_siruri": (many_samples, many[0]),
        "b_multe_siruri_frecvente": ([(b, 30) for b in frequent], frequent[0]),
        "c_toate_oglinzile_nemasurate": ([(b, 40 + i) for i, b in enumerate(one_side)], one_side[0]),
        "d_un_singur_qubit": ([("0", 1100), ("1", 948)], "0"),
        "d_un_singur_qubit_un_singur_sir": ([("0", 2048)], "0"),
    }


EDGE_CASES = _edge_cases()


@pytest.fixture(scope="module")
def histogram_runs(tmp_path_factory):
    browser = find_browser()
    if browser is None:
        pytest.skip("nu am găsit Edge sau Chrome pentru a rula JS-ul interfeței")
    folder = tmp_path_factory.mktemp("harness")
    cases = [{"name": name, "stage": _stage(samples, chosen)}
             for name, (samples, chosen) in EDGE_CASES.items()]
    page = folder / "harness.html"
    page.write_text(HARNESS.read_text(encoding="utf-8")
                    .replace("__ROOT__", INTERFACE.as_uri()).replace("__CASES__", json.dumps(cases)),
                    encoding="utf-8")
    done = subprocess.run(
        [browser, "--headless=new", "--disable-gpu", "--no-first-run",
         f"--user-data-dir={folder / 'profile'}", "--virtual-time-budget=5000", "--dump-dom",
         page.as_uri()], capture_output=True, timeout=180, check=True)
    dom = done.stdout.decode("utf-8")
    match = re.search(r'<pre id="out">(.*?)</pre>', dom, re.S)
    assert match, dom[:400]
    return {run["name"]: run for run in json.loads(html_lib.unescape(match.group(1)))}


def _expected_bars(samples, chosen, limits) -> list[str]:
    """Regula, scrisă independent de JS: frecvente (sau toate, dacă sunt puține), șirul ales și
    oglinda fiecărei bare afișate."""
    total = sum(c for _, c in samples)
    if len(samples) <= limits["maxBars"]:
        base = {b for b, _ in samples}
    else:
        base = {b for b, c in samples if c >= limits["rareShare"] * total} | {chosen}
    wanted = base | {_mirror(b) for b in base} | {chosen, _mirror(chosen)}
    return sorted(wanted)


@pytest.mark.parametrize("name", list(EDGE_CASES))
def test_histogram_edge_case(name, histogram_runs):
    samples, chosen = EDGE_CASES[name]
    run = histogram_runs[name]
    assert run["error"] is None, run["error"]

    shots = sum(c for _, c in samples)
    bars = run["bars"]
    bits = [b for b, _ in bars]
    measured = {b: c for b, c in bars if c > 0}

    # suma barelor plus „+N rare” este numărul de măsurători
    assert run["total"] == shots
    assert sum(c for _, c in bars) + run["hiddenShots"] == shots
    assert run["hiddenStrings"] == len(samples) - len(measured)
    assert all(measured[b] == dict(samples)[b] for b in measured)

    # regula: exact șirurile așteptate, în ordine binară, fără dubluri
    assert bits == _expected_bars(samples, chosen, run["limits"])
    assert len(set(bits)) == len(bits)

    # oglinda oricărei bare afișate există, iar șirul ales și oglinda lui sunt mereu pe axă
    assert all(_mirror(b) in bits for b in bits)
    assert chosen in bits and _mirror(chosen) in bits

    # cele trei numere ale perechii selectate dau numărul de citiri
    assert sum(run["cards"]) == READS


def test_edge_cases_cover_the_requested_situations(histogram_runs):
    limits = histogram_runs["b_peste_300_de_siruri"]["limits"]
    assert len(EDGE_CASES["a_doua_siruri_oglinzi"][0]) == 2
    assert len(EDGE_CASES["b_peste_300_de_siruri"][0]) > 300
    assert len(EDGE_CASES["b_peste_300_de_siruri"][0]) > limits["maxBars"]
    assert histogram_runs["b_peste_300_de_siruri"]["hiddenStrings"] > 300 - 4  # aproape toate rămân „rare”
    empty = [b for b, c in histogram_runs["c_toate_oglinzile_nemasurate"]["bars"] if c == 0]
    assert len(empty) == len(EDGE_CASES["c_toate_oglinzile_nemasurate"][0])
    assert {len(b) for b, _ in histogram_runs["d_un_singur_qubit"]["bars"]} == {1}
    dense = histogram_runs["b_multe_siruri_frecvente"]
    assert len(dense["bars"]) > limits["maxBars"] and dense["hiddenStrings"] == 0
