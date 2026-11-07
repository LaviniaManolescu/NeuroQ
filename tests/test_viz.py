"""Teste pentru figurile de demonstrație: fiecare funcție întoarce o figură pe notebook."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pytest  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402

from examples_data import notebook_matrix  # noqa: E402
from neuroq.encodings import PositionEncoding  # noqa: E402
from neuroq.encodings.position import build_evidence  # noqa: E402
from neuroq.solvers.qaoa import make_qaoa  # noqa: E402
from neuroq.sources.synthetic import NOTEBOOK_PAIR  # noqa: E402
from neuroq import viz  # noqa: E402

LABELS = [str(i + 1) for i in range(5)]


@pytest.fixture(autouse=True)
def close_figures():
    yield
    plt.close("all")


def problem():
    return PositionEncoding().encode(notebook_matrix())


def test_reads_evidence_and_ising_figures(tmp_path):
    m = notebook_matrix()
    a, d = build_evidence(m.values)
    p = problem()
    for fig in (viz.plot_reads_matrix(m.values, m.read_ids, m.snp_ids),
                viz.plot_evidence(a, d, LABELS),
                viz.plot_ising(p.J, p.offset, LABELS)):
        assert isinstance(fig, Figure)
    path = viz.save_figure(fig, tmp_path / "j.png")
    assert path.exists() and path.stat().st_size > 0


def test_graph_with_and_without_bitstring():
    p = problem()
    assert isinstance(viz.plot_graph(p.J, LABELS), Figure)
    fig = viz.plot_graph(p.J, LABELS, bitstring="01101", energy=p.energy("01101"))
    assert "energia = 1" in fig.axes[0].get_title()


def test_cost_and_top_samples():
    assert isinstance(viz.plot_qaoa_cost([[5.0, 3.0, 2.0], [6.0, 4.0]], optimum=1.0,
                                         best_start=0), Figure)
    top = [{"bitstring": "01101", "count": 400, "energy": 1.0},
           {"bitstring": "10010", "count": 380, "energy": 1.0},
           {"bitstring": "01100", "count": 90, "energy": 4.0}]
    assert isinstance(viz.plot_top_samples(top, optimum=1.0), Figure)


def test_circuit_small_draws_and_large_is_skipped():
    p = problem()
    qc, theta, _ = make_qaoa(p.J, p.h, p=1)
    bound = qc.assign_parameters([0.4, 0.3])
    assert isinstance(viz.plot_circuit(bound), Figure)
    big, _, _ = make_qaoa(np.triu(np.ones((14, 14)), k=1), np.zeros(14), p=1)
    assert viz.plot_circuit(big.assign_parameters([0.4, 0.3])) is None


def test_haplotypes_marks_differences():
    from neuroq.core import Haplotypes

    good = Haplotypes(np.array([1, 0, 0, 1, 0]), np.array([0, 1, 1, 0, 1]), np.zeros(9, int))
    fig = viz.plot_haplotypes(good, (NOTEBOOK_PAIR[0], NOTEBOOK_PAIR[1]))
    assert "identice" in fig.axes[0].get_title()  # oglinda contează ca potrivire
    bad = Haplotypes(np.array([0, 1, 1, 0, 0]), np.array([1, 0, 0, 1, 1]), np.zeros(9, int))
    fig = viz.plot_haplotypes(bad, (NOTEBOOK_PAIR[0], NOTEBOOK_PAIR[1]))
    assert "2 poziții diferite" in fig.axes[0].get_title()
