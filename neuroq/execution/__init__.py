"""Rulează sau doar compilează circuite pe simulator ideal, cip IBM fals sau hardware real."""

from neuroq.execution.base import Executor
from neuroq.execution.bitorder import qiskit_to_project
from neuroq.execution.ideal import IdealSimulator

__all__ = ["Executor", "IdealSimulator", "qiskit_to_project"]
