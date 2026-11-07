"""Capturi ale etapelor interfeței cu un browser headless (Edge sau Chrome), pentru verificare.

Rulare: python neuroq/interface/screenshots.py [folder]
Fiecare captură deschide GenoMEtric.html cu ?instance=...&step=...&frame=end, deci animațiile
apar în cadrul final. Nu e parte din pipeline; doar ajută la verificarea vizuală.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAGE = HERE / "GenoMEtric.html"
DEFAULT_OUT = HERE.parent.parent / "results" / "figures" / "interface"
INSTANCES = ("notebook", "illumina", "nanopore")
STEPS = range(1, 8)
CANDIDATES = (
    r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
    r"C:\Program Files\Google\Chrome\Application\chrome.exe",
)
WINDOW = "1600,1400"


def find_browser() -> str | None:
    for name in ("msedge", "chrome", "chromium"):
        if found := shutil.which(name):
            return found
    return next((c for c in CANDIDATES if Path(c).exists()), None)


def capture(browser: str, instance: str, step: int, out: Path, profile: Path) -> Path:
    path = out / f"{instance}_{step}.png"
    url = PAGE.as_uri() + f"?instance={instance}&step={step}&frame=end"
    subprocess.run([
        browser, "--headless=new", "--disable-gpu", "--hide-scrollbars", "--no-first-run",
        f"--user-data-dir={profile}", f"--window-size={WINDOW}", "--virtual-time-budget=4000",
        f"--screenshot={path}", url,
    ], check=True, capture_output=True, timeout=60)
    return path


def main() -> None:
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_OUT
    browser = find_browser()
    if browser is None:
        sys.exit("nu am găsit Edge sau Chrome pentru capturi")
    out.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as profile:
        for instance in INSTANCES:
            for step in STEPS:
                print(capture(browser, instance, step, out, Path(profile)))


if __name__ == "__main__":
    main()
