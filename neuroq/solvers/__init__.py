"""Caută soluția unui IsingProblem: referință exactă, QAOA, RQAOA și baseline-ul de lipire."""

from neuroq.solvers.base import Solver
from neuroq.solvers.classical import ClassicalSolver
from neuroq.solvers.exact import ExactSolver
from neuroq.solvers.greedy import GreedyPhaser
from neuroq.solvers.qaoa import QAOASolver

__all__ = ["ClassicalSolver", "ExactSolver", "GreedyPhaser", "QAOASolver", "Solver"]
