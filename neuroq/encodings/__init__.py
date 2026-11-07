"""Transformă ReadMatrix în IsingProblem: un qubit pe poziție sau un qubit pe citire."""

from neuroq.encodings.base import Encoding
from neuroq.encodings.position import PositionEncoding

__all__ = ["Encoding", "PositionEncoding"]
