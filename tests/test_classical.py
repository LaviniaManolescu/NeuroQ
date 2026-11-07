"""Teste pentru solverele clasice (pornind de la testul colegei, rescris în pytest)."""

import numpy as np
import pytest

from examples_data import X, Y, matrix_from, notebook_matrix, xy_matrix
from neuroq.encodings import PositionEncoding, ReadEncoding
from neuroq.encodings.position import build_evidence, pairwise_cost
from neuroq.solvers import ClassicalSolver, ExactSolver
from neuroq.solvers import _classical_algos as algos
from neuroq.sources.synthetic import NOTEBOOK_READS
from neuroq.transforms import fix_variable

ENCODINGS = [PositionEncoding(), ReadEncoding()]
NOTEBOOK_COPIES = {"01101", "10010"}


def copies(problem, bitstring: str) -> set[str]:
    haps = problem.decoder(bitstring)
    return {"".join(map(str, haps.copy_a)), "".join(map(str, haps.copy_b))}


def test_original_driver_on_notebook_evidence():
    a, d = build_evidence(NOTEBOOK_READS)
    x, c, _ = algos.solve_classical(a, d)
    # aserțiunea originală avea precedența greșită (`c==1 and ... or ...`)
    assert c == 1 and list(x) == [0, 1, 1, 0, 1]
    assert algos.cost(x, a, d) == pytest.approx(pairwise_cost(x, a, d))
    assert algos.exact_vectorized(a, d)[1] == 1


@pytest.mark.parametrize("encoding", ENCODINGS, ids=["position", "read"])
def test_best_reaches_exact_optimum_on_notebook(encoding):
    problem = encoding.encode(notebook_matrix())
    solution = ClassicalSolver("best").solve(problem)
    optimum = ExactSolver().solve(problem).energy
    # pe poziții optimul e 1; pe citiri aceeași împărțire costă 4 (alt Hamiltonian)
    assert solution.energy == pytest.approx(optimum)
    assert solution.energy == pytest.approx(1.0 if encoding.__class__ is PositionEncoding else 4.0)
    assert solution.energy == pytest.approx(problem.energy(solution.bitstring))
    assert copies(problem, solution.bitstring) == NOTEBOOK_COPIES


@pytest.mark.parametrize("encoding", ENCODINGS, ids=["position", "read"])
def test_greedy_tree_is_exact_on_clean_xy(encoding):
    problem = encoding.encode(xy_matrix())
    solution = ClassicalSolver("greedy_tree").solve(problem)
    assert solution.energy == pytest.approx(0.0)
    assert solution.metadata["components"] == 1
    assert copies(problem, solution.bitstring) == {X, Y}


@pytest.mark.parametrize("strategy", ["greedy_tree", "greedy_ls", "anneal", "best"])
def test_every_strategy_is_reproducible(strategy):
    problem = PositionEncoding().encode(notebook_matrix())
    first = ClassicalSolver(strategy, seed=5).solve(problem)
    second = ClassicalSolver(strategy, seed=5).solve(problem)
    assert first.bitstring == second.bitstring
    assert first.resources is None


def test_rejects_problems_with_local_fields():
    problem = fix_variable(PositionEncoding().encode(notebook_matrix()))
    with pytest.raises(ValueError, match="câmpuri h"):
        ClassicalSolver().solve(problem)


def test_rejects_unknown_strategy():
    with pytest.raises(ValueError):
        ClassicalSolver("quantum")


def _random_reads(rng: np.random.Generator) -> np.ndarray:
    """Ca în testul de stres al colegei: 3n citiri de 3–6 poziții, 5% erori."""
    n = int(rng.integers(8, 13))
    hidden = rng.integers(0, 2, n)
    reads = []
    for _ in range(3 * n):
        length = int(rng.integers(3, 7))
        start = int(rng.integers(0, n - length + 1))
        row = -np.ones(n, int)
        row[start:start + length] = hidden[start:start + length]
        if rng.random() < 0.5:
            row[start:start + length] = 1 - row[start:start + length]
        flips = rng.random(length) < 0.05
        row[start:start + length][flips] ^= 1
        reads.append(row)
    return np.array(reads)


def test_best_matches_exact_solver_on_20_random_instances():
    rng = np.random.default_rng(1)
    for t in range(20):
        problem = PositionEncoding().encode(matrix_from(_random_reads(rng)))
        found = ClassicalSolver("best", seed=t, restarts=10, steps=2000).solve(problem)
        assert found.energy == pytest.approx(ExactSolver().solve(problem).energy), f"instanța {t}"


@pytest.mark.parametrize("encoding", ENCODINGS, ids=["position", "read"])
@pytest.mark.parametrize("matrix", [notebook_matrix, xy_matrix], ids=["notebook", "xy"])
def test_exact_strategy_matches_exact_solver(encoding, matrix):
    problem = encoding.encode(matrix())
    found = ClassicalSolver("exact").solve(problem)
    assert found.energy == pytest.approx(ExactSolver().solve(problem).energy)
    assert found.metadata["exact_backend"] == "exact_vectorized"
