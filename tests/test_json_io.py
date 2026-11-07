"""Teste pentru formatul JSON al surselor."""

import json
from pathlib import Path

import numpy as np
import pytest

from neuroq.core import SourceData
from neuroq.sources import SyntheticSource, load_source_json, save_reads_data_json, save_source_json

REPO = Path(__file__).resolve().parent.parent


def assert_same_source(a: SourceData, b: SourceData) -> None:
    assert np.array_equal(a.values, b.values)
    assert np.array_equal(a.snp_positions, b.snp_positions)
    assert np.array_equal(a.read_starts, b.read_starts)
    assert a.read_ids == b.read_ids and a.technology == b.technology
    assert a.ref_alleles == b.ref_alleles and a.alt_alleles == b.alt_alleles
    assert a.special_labels == b.special_labels
    for x, y in zip(a.true_haplotypes, b.true_haplotypes):
        assert np.array_equal(x, y)


def test_round_trip_is_identical(tmp_path):
    source = SyntheticSource.notebook_example().load()
    source.special_labels = {"V1": 1, "V2": 4}
    path = tmp_path / "source.json"
    save_source_json(source, path, meta={"seed": 0})
    assert_same_source(source, load_source_json(path))
    assert json.loads(path.read_text(encoding="utf-8"))["labels"] == {"V1": 2, "V2": 5}


def test_round_trip_through_interface_format(tmp_path):
    source = SyntheticSource.notebook_example().load()
    path = tmp_path / "reads_data.json"
    save_reads_data_json(source, path)
    loaded = load_source_json(path)
    assert np.array_equal(loaded.values, source.values)
    assert loaded.technology == "interface" and loaded.true_haplotypes is None


def test_reads_interface_file_from_main():
    loaded = load_source_json(REPO / "reads_data.json")
    assert loaded.values.shape == (17, 9)
    assert loaded.values[0].tolist() == [1, 0, 0, -1, -1, -1, -1, -1, -1]


def test_missing_snp_becomes_minus_one(tmp_path):
    path = tmp_path / "s.json"
    path.write_text(json.dumps({"format_version": 1, "n_snps": 4,
                                "reads": [{"snps": {"1": 0, "4": 1}}]}), encoding="utf-8")
    assert load_source_json(path).values.tolist() == [[0, -1, -1, 1]]


@pytest.mark.parametrize("content, message", [
    ("", "gol"),
    ('{"format_version": 1, "n_snps": 3, "reads": [{"snps": {"4": 1}}]}', "în afara intervalului"),
    ('{"format_version": 1, "n_snps": 3, "reads": [{"snps": {"0": 1}}]}', "în afara intervalului"),
    ('{"format_version": 1, "n_snps": 3, "reads": [{"snps": {"1": 2}}]}', "nu 0 sau 1"),
    ('{"format_version": 1, "n_snps": 3, "reads": []}', "nicio citire"),
    ('{"n_snps": 3}', "format_version"),
])
def test_bad_files_give_clear_errors(tmp_path, content, message):
    path = tmp_path / "bad.json"
    path.write_text(content, encoding="utf-8")
    with pytest.raises(ValueError, match=message):
        load_source_json(path)
