"""Citirea și scrierea surselor în format JSON rar. Formatul e descris în docs/source_format.md.

Singurul loc unde indicii SNP din fișier (de la 1) se transformă în indicii din cod (de la 0).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Optional

import numpy as np

from neuroq.core import MISSING, SourceData

FORMAT_VERSION = 1
SNP_SPACING = 1000  # poziția implicită a SNP-ului k (de la 1) este k * 1000
DEFAULT_REF, DEFAULT_ALT = "A", "G"


def _read_json(path: Path) -> dict[str, Any]:
    text = Path(path).read_text(encoding="utf-8")
    if not text.strip():
        raise ValueError(f"{path}: fișierul este gol")
    try:
        data = json.loads(text)
    except json.JSONDecodeError as error:
        raise ValueError(f"{path}: JSON invalid ({error})") from error
    if not isinstance(data, dict):
        raise ValueError(f"{path}: rădăcina trebuie să fie un obiect JSON")
    return data


def _snp_row(snps: dict[str, int], n_snps: int, where: str) -> np.ndarray:
    """Transformă {"3": 1, "7": 0} (indici de la 1) într-un rând cu -1 pentru lipsă."""
    row = np.full(n_snps, MISSING)
    for key, value in snps.items():
        try:
            index = int(key)
        except ValueError:
            raise ValueError(f"{where}: cheia SNP {key!r} nu este un număr") from None
        if not 1 <= index <= n_snps:
            raise ValueError(f"{where}: SNP-ul {index} e în afara intervalului 1..{n_snps}")
        if value not in (0, 1):
            raise ValueError(f"{where}: SNP-ul {index} are valoarea {value!r}, nu 0 sau 1")
        row[index - 1] = value
    return row


def _snp_dict(row: np.ndarray) -> dict[str, int]:
    return {str(i + 1): int(v) for i, v in enumerate(row) if v != MISSING}


def _bits(text: str, n_snps: int, name: str) -> np.ndarray:
    if len(text) != n_snps or set(text) - {"0", "1"}:
        raise ValueError(f"truth.{name} trebuie să aibă {n_snps} caractere 0/1")
    return np.array([int(c) for c in text])


def _first_observed_start(row: np.ndarray, positions: np.ndarray) -> int:
    seen = np.flatnonzero(row != MISSING)
    return int(positions[seen[0]]) if len(seen) else 0


def _build(
    n_snps: int,
    rows: list[np.ndarray],
    read_ids: list[str],
    starts: list[Optional[int]],
    data: dict[str, Any],
) -> SourceData:
    if n_snps < 1:
        raise ValueError("n_snps trebuie să fie cel puțin 1")
    if not rows:
        raise ValueError("fișierul nu conține nicio citire")
    positions = np.array(data.get("snp_positions") or [(k + 1) * SNP_SPACING for k in range(n_snps)])
    if len(positions) != n_snps:
        raise ValueError(f"snp_positions trebuie să aibă {n_snps} valori")
    values = np.array(rows)
    read_starts = np.array([
        s if s is not None else _first_observed_start(row, positions) for s, row in zip(starts, rows)
    ])
    truth = data.get("truth")
    true_haplotypes = None
    if truth is not None:
        true_haplotypes = (_bits(truth["copy_a"], n_snps, "copy_a"),
                           _bits(truth["copy_b"], n_snps, "copy_b"))
    labels = {name: int(index) - 1 for name, index in (data.get("labels") or {}).items()}
    for name, column in labels.items():
        if not 0 <= column < n_snps:
            raise ValueError(f"eticheta {name} indică un SNP în afara intervalului 1..{n_snps}")
    return SourceData(
        snp_positions=positions,
        ref_alleles=list(data.get("ref_alleles") or [DEFAULT_REF] * n_snps),
        alt_alleles=list(data.get("alt_alleles") or [DEFAULT_ALT] * n_snps),
        read_ids=read_ids,
        read_starts=read_starts,
        values=values,
        technology=data.get("technology", "unknown"),
        special_labels=labels,
        true_haplotypes=true_haplotypes,
    )


def load_source_json(path: Path) -> SourceData:
    """Citește formatul NeuroQ (format_version 1) sau formatul interfeței (reads_data.json)."""
    data = _read_json(path)
    if "format_version" in data:
        if data["format_version"] != FORMAT_VERSION:
            raise ValueError(f"{path}: format_version {data['format_version']} nu e suportat")
        n_snps = int(data["n_snps"])
        reads = data.get("reads") or []
        rows = [_snp_row(r["snps"], n_snps, f"{path}, citirea {r.get('id', i + 1)}")
                for i, r in enumerate(reads)]
        ids = [str(r.get("id", f"r{i + 1}")) for i, r in enumerate(reads)]
        starts = [r.get("start") for r in reads]
        return _build(n_snps, rows, ids, starts, data)
    if "metadata" in data and "reads" in data:
        n_snps = int(data["metadata"]["num_snps"])
        reads = data["reads"] or []
        rows = [_snp_row(r, n_snps, f"{path}, citirea {i + 1}") for i, r in enumerate(reads)]
        ids = [f"r{i + 1}" for i in range(len(reads))]
        return _build(n_snps, rows, ids, [None] * len(reads), {"technology": "interface"})
    raise ValueError(f"{path}: lipsește format_version (sau metadata + reads pentru formatul interfeței)")


def save_source_json(source: SourceData, path: Path, meta: Optional[dict[str, Any]] = None) -> None:
    """Scrie sursa în formatul NeuroQ, versiunea 1."""
    n_snps = source.values.shape[1]
    payload: dict[str, Any] = {
        "format_version": FORMAT_VERSION,
        "n_snps": n_snps,
        "technology": source.technology,
        "snp_positions": [int(p) for p in source.snp_positions],
        "ref_alleles": list(source.ref_alleles),
        "alt_alleles": list(source.alt_alleles),
        "labels": {name: col + 1 for name, col in source.special_labels.items()},
        "reads": [
            {"id": rid, "start": int(start), "snps": _snp_dict(row)}
            for rid, start, row in zip(source.read_ids, source.read_starts, source.values)
        ],
    }
    if source.true_haplotypes is not None:
        copy_a, copy_b = source.true_haplotypes
        payload["truth"] = {"copy_a": "".join(map(str, copy_a)),
                            "copy_b": "".join(map(str, copy_b))}
    payload["meta"] = meta or {}
    Path(path).write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def save_reads_data_json(source: SourceData, path: Path, description: str = "") -> None:
    """Scrie doar citirile, în formatul interfeței (reads_data.json, SNP-uri de la 1)."""
    payload = {
        "metadata": {
            "num_snps": int(source.values.shape[1]),
            "num_reads": int(source.values.shape[0]),
            "description": description or f"Exportat din NeuroQ ({source.technology})",
        },
        "reads": [_snp_dict(row) for row in source.values],
    }
    Path(path).write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
