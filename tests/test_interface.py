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
