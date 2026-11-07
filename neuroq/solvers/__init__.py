"""Caută soluția unui IsingProblem: referință exactă, QAOA, RQAOA și baseline-ul de lipire."""

from neuroq.solvers.base import Solver
from neuroq.solvers.exact import ExactSolver
from neuroq.solvers.qaoa import QAOASolver

__all__ = ["ExactSolver", "QAOASolver", "Solver"]
