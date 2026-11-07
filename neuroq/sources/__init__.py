"""Aduce citirile unei regiuni: generator sintetic cu profiluri de tehnologie și format JSON."""

from neuroq.sources.json_io import load_source_json, save_reads_data_json, save_source_json
from neuroq.sources.synthetic import SyntheticSource

__all__ = ["SyntheticSource", "load_source_json", "save_reads_data_json", "save_source_json"]
