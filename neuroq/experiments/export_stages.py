"""Exportă etapele demonstrației pentru interfață: results/stages.js (window.STAGES = {...}).

Refolosește run_stages; aici doar se convertesc rezultatele în JSON simplu (fără numpy,
fără NaN) și se desenează circuitul. Fișierul se scrie atomic (.tmp + os.replace).
Rulare: python -m neuroq.experiments.export_stages
"""

from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path
from typing import Any, Optional

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from neuroq.core import MISSING, ReadMatrix, SourceData  # noqa: E402
from neuroq.experiments.stages import StagesResult, load_case, run_stages  # noqa: E402
from neuroq.sources import save_reads_data_json  # noqa: E402
from neuroq.viz import plot_circuit, save_figure  # noqa: E402

# Instanțele exportate și codificarea care încape pe fiecare (fără fixare).
EXPORTS = (("notebook", "position"), ("illumina", "position"), ("nanopore", "read"))
OUTPUT_DIR = Path("results")
CIRCUIT_DIR = Path("figures") / "stages"  # relativ la OUTPUT_DIR
JS_PREFIX = "window.STAGES = "


def _bits(copy) -> str:
    return "".join("-" if v == MISSING else str(int(v)) for v in copy)


def _matrix(values: np.ndarray) -> list[list[float]]:
    return [[float(v) for v in row] for row in np.asarray(values)]


def _evidence(result: StagesResult) -> dict[str, Any]:
    m = result.matrix
    positions = [int(s) for s in m.snp_ids]
    spans = []
    for row in m.values:
        seen = np.flatnonzero(row != MISSING)
        spans.append([positions[seen[0]], positions[seen[-1]]] if len(seen) else None)
    return {
        "values": [[None if v == MISSING else int(v) for v in row] for row in m.values],
        "read_ids": list(m.read_ids),
        "snp_ids": list(m.snp_ids),
        "snp_positions_bp": positions,
        # intervalul acoperit pe SNP-uri; capetele reale ale citirilor nu sunt păstrate
        "read_spans_bp": spans,
    }


def _read_assignment(result: StagesResult) -> list[dict[str, Any]]:
    haps = result.haplotypes
    copies = (np.asarray(haps.copy_a), np.asarray(haps.copy_b))
    rows = []
    for read_id, row, team in zip(result.matrix.read_ids, result.matrix.values,
                                  haps.read_assignment):
        seen = np.flatnonzero(row != MISSING)
        team = int(team)
        if team in (0, 1):
            reference = copies[team]
        else:  # suspectă: comparăm cu copia cea mai apropiată
            misses = [int(np.sum(row[seen] != c[seen])) for c in copies]
            reference = copies[int(np.argmin(misses))]
        rows.append({
            "read_id": read_id,
            "team": team,
            "snps": [int(i) for i in seen],
            "values": [int(row[i]) for i in seen],
            "errors": [bool(reference[i] != MISSING and row[i] != reference[i]) for i in seen],
        })
    return rows


def _circuit_png(result: StagesResult, out_dir: Path) -> Optional[str]:
    relative = CIRCUIT_DIR / f"{result.instance}_circuit.png"
    fig = plot_circuit(result.circuit)
    if fig is None:
        return None
    save_figure(fig, out_dir / relative)
    plt.close(fig)
    return relative.as_posix()


def stage_entry(result: StagesResult, out_dir: Path) -> dict[str, Any]:
    """Tot ce afișează interfața pentru o instanță, ca JSON simplu."""
    m, problem, haps = result.matrix, result.problem, result.haplotypes
    J = problem.J
    truth = m.true_haplotypes
    return {
        "metadata": {
            "instance": result.instance,
            "encoding": result.encoding,
            "n_qubits": int(result.solved.n_qubits),
            "n_reads": len(m.read_ids),
            "n_snps": len(m.snp_ids),
            "qaoa_seconds": round(float(result.qaoa_seconds), 3),
        },
        "evidence": _evidence(result),
        "a": _matrix(result.a),
        "d": _matrix(result.d),
        "ising": {
            "edges": [{"i": int(i), "j": int(j), "J": float(J[i, j])}
                      for i, j in zip(*np.nonzero(J))],
            "offset": float(problem.offset),
        },
        "cost_history": [[float(c) for c in h] for h in result.cost_history],
        "best_start": int(result.best_start),
        "top_samples": [{"bitstring": s["bitstring"], "count": int(s["count"]),
                         "energy": float(s["energy"])} for s in result.top_samples],
        # distribuția completă din QAOASolver (samples), fără nicio rulare nouă; ordine binară
        "all_samples": [{"bitstring": bits, "count": int(count),
                         "energy": float(result.solved.energy(bits))}
                        for bits, count in sorted(result.qaoa.samples.items(),
                                                  key=lambda item: int(item[0], 2))],
        "result": {"bitstring": result.full_bitstring,
                   "energy": float(problem.energy(result.full_bitstring))},
        "haplotypes": {
            "decoded": [_bits(haps.copy_a), _bits(haps.copy_b)],
            "truth": [_bits(truth[0]), _bits(truth[1])] if truth is not None else None,
        },
        "read_assignment": _read_assignment(result),
        "metrics": {key: float(result.metrics[key])
                    for key in ("accuracy", "switch_errors", "mec_cost")
                    if key in result.metrics},
        "circuit_png": _circuit_png(result, out_dir),
    }


def build_stages(out_dir: Path) -> dict[str, Any]:
    return {name: stage_entry(run_stages(load_case(name), encoding), out_dir)
            for name, encoding in EXPORTS}


def write_stages_js(stages: dict[str, Any], path: Path) -> Path:
    """window.STAGES = {...}; scris atomic, în UTF-8, fără NaN sau Infinity."""
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    text = JS_PREFIX + json.dumps(stages, indent=1, ensure_ascii=False, allow_nan=False) + ";\n"
    tmp.write_text(text, encoding="utf-8")
    os.replace(tmp, path)
    return path


def export_stages(out_dir: Path = OUTPUT_DIR) -> Path:
    return write_stages_js(build_stages(out_dir), out_dir / "stages.js")


def _source_from_matrix(matrix: ReadMatrix) -> SourceData:
    """SourceData pentru export, cu pozițiile SNP-urilor luate din snp_ids."""
    positions = np.array([int(s) for s in matrix.snp_ids])
    starts = np.array([positions[np.flatnonzero(row != MISSING)[0]] for row in matrix.values])
    n_snps = len(positions)
    return SourceData(
        snp_positions=positions, ref_alleles=["A"] * n_snps, alt_alleles=["G"] * n_snps,
        read_ids=list(matrix.read_ids), read_starts=starts, values=matrix.values,
        technology=matrix.technology, special_labels=dict(matrix.special_labels),
        true_haplotypes=matrix.true_haplotypes,
    )


def export_inputs(out_dir: Path = OUTPUT_DIR) -> list[Path]:
    """Datele fiecărei instanțe în formatul de intrare al interfeței (results/inputs/)."""
    paths = []
    for name, _ in EXPORTS:
        path = out_dir / "inputs" / f"{name}.json"
        path.parent.mkdir(parents=True, exist_ok=True)
        save_reads_data_json(_source_from_matrix(load_case(name).matrix), path,
                             description=f"Instanța {name} din NeuroQ",
                             metadata={"instance": name})
        paths.append(path)
    return paths


def main() -> None:
    sys.stdout.reconfigure(encoding="utf-8")
    start = time.perf_counter()
    path = export_stages()
    inputs = export_inputs()
    print(f"scris {path} ({path.stat().st_size / 1024:.1f} KB) și {len(inputs)} fișiere de "
          f"intrare în {time.perf_counter() - start:.1f} s")


if __name__ == "__main__":
    main()
