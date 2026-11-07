"""Validarea comună a unui fișier results.json (probă sau real)."""

import json
from pathlib import Path

from neuroq.core import ResultRow


def validate_results_file(path: Path, is_sample: bool) -> list[ResultRow]:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    assert data["is_sample"] is is_sample
    assert data["schema_version"] == "0.1"
    rows = []
    for raw in data["rows"]:
        row = ResultRow.from_dict(raw)
        assert row.to_dict() == raw

        assert all(isinstance(v, (int, float)) and not isinstance(v, bool)
                   for v in row.metrics.values())

        if row.method == "resources":
            # rând fără rezolvare: doar resurse logice și motivul pentru care nu s-a rulat
            assert row.encoding in {"position", "read"}
            assert row.resources is not None and "skipped_reason" in row.details
            assert {"n_variables", "n_couplings", "rzz_per_layer"} <= set(row.metrics)
            assert "haplotypes" not in row.details
            rows.append(row)
            continue

        haps = row.details["haplotypes"]
        assert len(haps["copy_a"]) == len(haps["copy_b"]) == row.instance["n_snps"]
        assert len(haps["read_assignment"]) == row.instance["n_reads"]
        assert set(haps["read_assignment"]) <= {-1, 0, 1}

        if "top_samples" in row.details:
            counts = [s["count"] for s in row.details["top_samples"]]
            assert counts == sorted(counts, reverse=True)
            for sample in row.details["top_samples"]:
                assert len(sample["bitstring"]) == row.resources.n_qubits
        rows.append(row)
    return rows
