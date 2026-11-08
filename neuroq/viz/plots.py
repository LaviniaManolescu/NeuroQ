"""Figuri matplotlib pentru demonstrație. Fiecare funcție întoarce figura, fără să o afișeze.

Funcțiile primesc date gata calculate (matrici, șiruri, haplotipuri) și nu importă Qiskit:
`plot_circuit` doar apelează `draw("mpl")` pe circuitul primit.
"""

from __future__ import annotations

import logging
import math
from pathlib import Path
from typing import Any, Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.figure import Figure
from matplotlib.patches import Patch

from neuroq.core import Haplotypes

log = logging.getLogger(__name__)

DPI = 200
FONT_SIZE = 14
CELL_FONT = 12
MAX_CIRCUIT_QUBITS = 12
CIRCUIT_FOLD = 30  # câte coloane pe un rând, pentru circuitele adânci
NO_FOLD_DEPTH = 30  # sub această adâncime circuitul încape pe un singur rând
BASIC_GATES = {"h", "rz", "rzz", "rx", "measure", "barrier"}

COLOR_ZERO = "#6ACCE6"
COLOR_ONE = "#F7BD69"
COLOR_MISSING = "#EEF1F5"
COLOR_POSITIVE = "#C0392B"  # J > 0: respinge (vrea biți diferiți)
COLOR_NEGATIVE = "#2E6FD8"  # J < 0: atrage (vrea biți egali)
COLOR_VIOLATED = "#E0001B"

plt.rcParams.update({
    "font.size": FONT_SIZE,
    "axes.titlesize": FONT_SIZE + 2,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
})


def save_figure(fig: Figure, path: Path) -> Path:
    """Salvează figura la dpi 200, cu fundal alb."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=DPI, bbox_inches="tight", facecolor="white")
    return path


def _bit_text(bits: Sequence[int]) -> str:
    return "".join("-" if b == -1 else str(int(b)) for b in bits)


def _annotate_cells(ax, values: np.ndarray, fmt: str = "{:g}", skip=None) -> None:
    for (r, c), v in np.ndenumerate(values):
        if skip is not None and skip(v):
            continue
        ax.text(c, r, fmt.format(v), ha="center", va="center", fontsize=CELL_FONT)


def plot_reads_matrix(values: np.ndarray, read_ids: Sequence[str],
                      snp_ids: Sequence[str], title: str = "Citirile") -> Figure:
    """Citirile pe rânduri, SNP-urile pe coloane: 0, 1 sau gol (neobservat)."""
    n_reads, n_snps = values.shape
    fig, ax = plt.subplots(figsize=(max(5, 0.7 * n_snps + 2), max(3, 0.45 * n_reads + 1.5)))
    cmap = ListedColormap([COLOR_MISSING, COLOR_ZERO, COLOR_ONE])
    ax.imshow(values + 1, cmap=cmap, vmin=0, vmax=2, aspect="auto")
    _annotate_cells(ax, values, fmt="{:d}", skip=lambda v: v == -1)
    ax.set_xticks(range(n_snps), [str(i + 1) for i in range(n_snps)])
    ax.set_yticks(range(n_reads), list(read_ids))
    ax.set_xlabel("SNP")
    ax.set_title(f"{title}: {n_reads} citiri × {n_snps} SNP-uri")
    ax.legend(handles=[Patch(color=COLOR_ZERO, label="0 (referință)"),
                       Patch(color=COLOR_ONE, label="1 (alternativ)"),
                       Patch(color=COLOR_MISSING, label="neobservat")],
              loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3, frameon=False)
    fig.tight_layout()
    return fig


def _upper_heatmap(ax, matrix: np.ndarray, labels: Sequence[str], title: str, cmap: str,
                   vmin: float, vmax: float, fmt: str = "{:g}") -> None:
    n = len(matrix)
    shown = np.where(np.triu(np.ones((n, n), dtype=bool), k=1), matrix, np.nan)
    ax.imshow(shown, cmap=cmap, vmin=vmin, vmax=vmax)
    _annotate_cells(ax, shown, fmt=fmt, skip=lambda v: np.isnan(v) or v == 0)
    ax.set_xticks(range(n), labels, rotation=90 if n > 12 else 0)
    ax.set_yticks(range(n), labels)
    ax.set_title(title)


def plot_evidence(a: np.ndarray, d: np.ndarray, labels: Sequence[str]) -> Figure:
    """Contoarele a („la fel”) și d („opus”) alăturate, triunghiul de sus."""
    n = len(a)
    size = max(4.5, 0.5 * n + 2)
    fig, axes = plt.subplots(1, 2, figsize=(2 * size, size))
    top = max(1.0, float(max(a.max(), d.max())))
    _upper_heatmap(axes[0], a, labels, "a: voturi „la fel”", "Blues", 0, top)
    _upper_heatmap(axes[1], d, labels, "d: voturi „opus”", "Oranges", 0, top)
    fig.tight_layout()
    return fig


def plot_ising(J: np.ndarray, offset: float, labels: Sequence[str]) -> Figure:
    """Cuplajele J cu semn: roșu respinge (J > 0), albastru atrage (J < 0)."""
    n = len(J)
    size = max(5, 0.55 * n + 2)
    fig, ax = plt.subplots(figsize=(size + 1, size))
    limit = max(1.0, float(np.abs(J).max()))
    _upper_heatmap(ax, J, labels, f"J = (d − a) / 2,   constanta = {offset:g}", "RdBu_r",
                   -limit, limit)
    fig.colorbar(ax.images[0], ax=ax, shrink=0.8, label="J  (>0 respinge, <0 atrage)")
    fig.tight_layout()
    return fig


def _circle(n: int) -> np.ndarray:
    angles = np.pi / 2 - 2 * np.pi * np.arange(n) / n
    return np.column_stack([np.cos(angles), np.sin(angles)])


def plot_graph(J: np.ndarray, labels: Sequence[str], bitstring: Optional[str] = None,
               energy: Optional[float] = None, title: str = "Graful cuplajelor") -> Figure:
    """Graf circular: muchii cu valoarea J, culoare după semn, grosime după mărime.

    Cu `bitstring`, nodurile se colorează după bit, iar muchiile încălcate
    (J > 0 cu biți egali, J < 0 cu biți diferiți) devin roșii.
    """
    n = len(J)
    xy = _circle(n)
    fig, ax = plt.subplots(figsize=(8, 8))
    limit = max(1e-9, float(np.abs(J).max()))
    bits = [int(c) for c in bitstring] if bitstring is not None else None
    violated = 0
    for i in range(n):
        for k in range(i + 1, n):
            value = J[i, k]
            if value == 0:
                continue
            broken = bits is not None and ((value > 0) == (bits[i] == bits[k]))
            violated += int(broken)
            color = COLOR_VIOLATED if broken else (COLOR_POSITIVE if value > 0 else COLOR_NEGATIVE)
            style = "-" if value < 0 or broken else "--"
            ax.plot(*zip(xy[i], xy[k]), color=color, linestyle=style,
                    linewidth=1 + 4 * abs(value) / limit, alpha=0.9 if broken else 0.6, zorder=1)
            if n <= 12:
                mid = (xy[i] + xy[k]) / 2
                ax.text(*mid, f"{value:g}", fontsize=CELL_FONT - 2, ha="center", va="center",
                        bbox={"boxstyle": "round,pad=0.15", "fc": "white", "ec": "none"}, zorder=2)
    node_colors = ([COLOR_ONE if b else COLOR_ZERO for b in bits] if bits is not None
                   else ["#DCE1E8"] * n)
    ax.scatter(xy[:, 0], xy[:, 1], s=900, c=node_colors, edgecolors="#14243A", zorder=3)
    for i, (x, y) in enumerate(xy):
        ax.text(x, y, labels[i], ha="center", va="center", fontsize=CELL_FONT, zorder=4)
    parts = [title]
    if bitstring is not None:
        parts.append(f"șirul {bitstring}: {violated} muchii încălcate")
    if energy is not None:
        parts.append(f"energia = {energy:g}")
    ax.set_title("\n".join(parts))
    handles = [Patch(color=COLOR_NEGATIVE, label="J < 0 atrage"),
               Patch(color=COLOR_POSITIVE, label="J > 0 respinge")]
    if bits is not None:
        handles.append(Patch(color=COLOR_VIOLATED, label="încălcată"))
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.02), ncol=3,
              frameon=False)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.tight_layout()
    return fig


def plot_qaoa_cost(histories: Sequence[Sequence[float]], optimum: Optional[float] = None,
                   best_start: Optional[int] = None) -> Figure:
    """Valoarea medie a energiei la fiecare evaluare COBYLA, câte o curbă pe start."""
    fig, ax = plt.subplots(figsize=(9, 5))
    for s, history in enumerate(histories):
        label = f"start {s + 1}" + (" (ales)" if s == best_start else "")
        ax.plot(range(1, len(history) + 1), history, marker="o", markersize=3,
                linewidth=2.5 if s == best_start else 1.5, label=label)
    if optimum is not None:
        ax.axhline(optimum, color="#14243A", linestyle="--", label=f"optim exact = {optimum:g}")
    ax.set_xlabel("evaluare COBYLA")
    ax.set_ylabel("energia medie ⟨H⟩")
    ax.set_title("Antrenarea QAOA: costul pe iterații")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    return fig


def plot_top_samples(top: Sequence[dict[str, Any]], optimum: Optional[float] = None) -> Figure:
    """Cele mai probabile șiruri măsurate, cu energia fiecăruia scrisă pe bară."""
    top = list(top)[:10]
    fig, ax = plt.subplots(figsize=(10, 0.55 * len(top) + 2))
    names = [t["bitstring"] for t in top][::-1]
    counts = [t["count"] for t in top][::-1]
    energies = [t["energy"] for t in top][::-1]
    colors = ["#2F8A55" if optimum is not None and math.isclose(e, optimum) else COLOR_ZERO
              for e in energies]
    bars = ax.barh(names, counts, color=colors)
    for bar, energy in zip(bars, energies):
        ax.text(bar.get_width(), bar.get_y() + bar.get_height() / 2, f"  E = {energy:g}",
                va="center", fontsize=CELL_FONT)
    ax.set_xlabel("apariții")
    ax.set_title("Cele mai probabile 10 șiruri" + (" (verde = optim)" if optimum is not None else ""))
    ax.set_xlim(0, max(counts) * 1.3 if counts else 1)
    for label in ax.get_yticklabels():
        label.set_family("monospace")
    fig.tight_layout()
    return fig


def plot_confidence(
    samples: dict[str, int], energies: dict[str, float], optimum: float, chosen: str,
    p: Optional[int] = None,
) -> Figure:
    """Încrederea: ce procent din măsurători cade pe fiecare nivel de energie.

    `energies` dă energia fiecărui șir măsurat; nivelul optimului (cu tot cu oglinda lui, care
    are aceeași energie) e verde, iar nivelul șirului ales are contur.
    """
    total = sum(samples.values())
    by_energy: dict[float, int] = {}
    for bits, count in samples.items():
        by_energy[energies[bits]] = by_energy.get(energies[bits], 0) + count
    levels = sorted(by_energy)
    fig, ax = plt.subplots(figsize=(10, 4.2))
    colors = ["#2F8A55" if math.isclose(e, optimum) else COLOR_ZERO for e in levels]
    bars = ax.bar([f"{e:g}" for e in levels], [100 * by_energy[e] / total for e in levels],
                  color=colors)
    for bar, level in zip(bars, levels):
        if math.isclose(level, energies[chosen]):
            bar.set_edgecolor("black")
            bar.set_linewidth(2)
    hits = 100 * sum(c for b, c in samples.items() if math.isclose(energies[b], optimum)) / total
    ax.set_xlabel("energia șirului măsurat")
    ax.set_ylabel("% din măsurători")
    ax.set_title((f"p = {p}: " if p else "") + f"încrederea = {hits:.1f}% din măsurători pe "
                 f"energia optimă {optimum:g} (verde, oglinda inclusă); conturul = șirul ales")
    ax.title.set_fontsize(10)
    fig.tight_layout()
    return fig


def plot_circuit(circuit: Any, max_qubits: int = MAX_CIRCUIT_QUBITS,
                 title: Optional[str] = None) -> Optional[Figure]:
    """Desenează circuitul QAOA cu H, RZZ și RX vizibile; peste `max_qubits` doar loghează."""
    ops = dict(circuit.count_ops())
    if circuit.num_qubits > max_qubits:
        log.warning("circuit prea mare pentru desen: %d qubiți, porți %s",
                    circuit.num_qubits, ops)
        return None
    if set(ops) - BASIC_GATES:
        circuit = circuit.decompose(gates_to_decompose=list(set(ops) - BASIC_GATES))
    fold = -1 if circuit.depth() <= NO_FOLD_DEPTH else CIRCUIT_FOLD
    fig = circuit.draw("mpl", fold=fold, style={"fontsize": FONT_SIZE})
    fig.patch.set_facecolor("white")
    two_qubit = sum(1 for item in circuit.data if item.operation.num_qubits == 2)
    fig.suptitle(title or f"Circuitul QAOA: {circuit.num_qubits} qubiți, "
                 f"{two_qubit} porți cu doi qubiți, adâncime {circuit.depth()}",
                 fontsize=FONT_SIZE + 2)
    return fig


def _best_orientation(haps: Haplotypes, truth_a: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    a, b = np.asarray(haps.copy_a), np.asarray(haps.copy_b)
    if np.sum(b == truth_a) > np.sum(a == truth_a):
        return b, a
    return a, b


def plot_haplotypes(haps: Haplotypes, truth: tuple[np.ndarray, np.ndarray]) -> Figure:
    """Copiile decodate față de cele adevărate (aliniate cu oglinda), diferențele marcate."""
    truth_a, truth_b = (np.asarray(t) for t in truth)
    found_a, found_b = _best_orientation(haps, truth_a)
    rows = np.array([truth_a, found_a, truth_b, found_b])
    names = ["A adevărată", "A găsită", "B adevărată", "B găsită"]
    n = rows.shape[1]
    fig, ax = plt.subplots(figsize=(max(6, 0.7 * n + 3), 4))
    cmap = ListedColormap([COLOR_MISSING, COLOR_ZERO, COLOR_ONE])
    ax.imshow(rows + 1, cmap=cmap, vmin=0, vmax=2, aspect="auto")
    _annotate_cells(ax, rows, fmt="{:d}", skip=lambda v: v == -1)
    wrong = 0
    for row, truth_row in ((1, truth_a), (3, truth_b)):
        for c in np.flatnonzero(rows[row] != truth_row):
            wrong += 1
            ax.add_patch(plt.Rectangle((c - 0.5, row - 0.5), 1, 1, fill=False,
                                       edgecolor=COLOR_VIOLATED, linewidth=3))
    ax.set_yticks(range(4), names)
    ax.set_xticks(range(n), [str(i + 1) for i in range(n)])
    ax.set_xlabel("SNP")
    ax.axhline(1.5, color="#14243A", linewidth=2)
    if wrong == 0:
        verdict = "identice"
    elif wrong == 1:
        verdict = "1 poziție diferită (chenar roșu)"
    else:
        verdict = f"{wrong} poziții diferite (chenar roșu)"
    ax.set_title(f"Copiile găsite față de cele adevărate: {verdict}\n"
                 f"A = {_bit_text(found_a)}   B = {_bit_text(found_b)}")
    fig.tight_layout()
    return fig
