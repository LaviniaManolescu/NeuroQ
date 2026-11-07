"""Notebook-ul de demonstrație rulează de sus în jos, fără erori, în mai puțin de un minut."""

import os
import sys
import time
from pathlib import Path

import nbformat
from nbclient import NotebookClient

NOTEBOOK = Path(__file__).resolve().parent.parent / "notebooks" / "demo.ipynb"
MAX_SECONDS = 60


def test_demo_notebook_executes(monkeypatch):
    # kernelul python3 pornește `python` din PATH: punem întâi interpretorul care rulează testele
    monkeypatch.setenv("PATH", str(Path(sys.executable).parent) + os.pathsep + os.environ["PATH"])
    nb = nbformat.read(NOTEBOOK, as_version=4)
    start = time.perf_counter()
    NotebookClient(nb, timeout=120, kernel_name="python3",
                   resources={"metadata": {"path": str(NOTEBOOK.parent)}}).execute()
    assert time.perf_counter() - start < MAX_SECONDS
    images = sum("image/png" in out.get("data", {})
                 for cell in nb.cells if cell.cell_type == "code" for out in cell.outputs)
    assert images >= 8
