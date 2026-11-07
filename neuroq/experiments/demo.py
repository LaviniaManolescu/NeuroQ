"""Demonstrația live: rulează pipeline-ul și arată fiecare etapă în terminal și în figuri.

Rulare:
    python -m neuroq.experiments.demo --instance notebook
    python -m neuroq.experiments.demo --instance notebook --encoding citiri --fix --show
Figurile se salvează în results/figures/demo_<instanță>/. Fără --show se folosește backend-ul
Agg (doar salvare); cu --show figurile se deschid pe rând și demo-ul așteaptă Enter la final.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path
from typing import Optional

import matplotlib
import numpy as np

from neuroq.experiments.stages import (
    ENCODING_NAMES, INSTANCES, ROMANIAN, StagesResult, encoding_resources, load_case,
    resolve_encoding, run_stages,
)

DEFAULT_OUTPUT = Path("results") / "figures"


def _title(number: int, text: str) -> None:
    print(f"\n{'=' * 64}\n {number}. {text}\n{'=' * 64}", flush=True)


def _matrix_text(values: np.ndarray, row_names: list[str], col_names: list[str],
                 missing_dot: bool = False) -> str:
    """Tabel text; cu `missing_dot`, -1 (neobservat) apare ca punct."""
    width = max(4, *(len(c) for c in col_names)) + 1
    head = " " * 10 + "".join(c.rjust(width) for c in col_names)
    lines = [head]
    for name, row in zip(row_names, values):
        cells = "".join(("." if missing_dot and v == -1 else f"{v:g}").rjust(width) for v in row)
        lines.append(name[:9].ljust(10) + cells)
    return "\n".join(lines)


def _bits(copy) -> str:
    return "".join("-" if v == -1 else str(int(v)) for v in copy)


def live_cost(start: int, iteration: int, cost: float) -> None:
    """Callback pentru QAOA: o linie la fiecare evaluare COBYLA."""
    print(f"   start {start + 1}, evaluarea {iteration:3d}:  <H> = {cost:9.4f}", flush=True)


class FigureSink:
    """Salvează figurile numerotate și, cu --show, le deschide pe rând."""

    def __init__(self, directory: Path, suffix: str, show: bool) -> None:
        self.directory = directory
        self.suffix = suffix
        self.show = show
        self.saved: list[Path] = []

    def add(self, number: int, name: str, fig) -> None:
        from neuroq.viz import save_figure

        if fig is None:
            print("   (circuitul e prea mare pentru desen; vezi numerele de mai sus)")
            return
        path = save_figure(fig, self.directory / f"{number:02d}_{name}_{self.suffix}.png")
        self.saved.append(path)
        print(f"   figura: {path}")
        if self.show:
            import matplotlib.pyplot as plt

            plt.show(block=False)
            plt.pause(0.5)
        else:
            import matplotlib.pyplot as plt

            plt.close(fig)


def present(result: StagesResult, sink: FigureSink) -> None:
    """Tipărește etapele și desenează figurile pentru un rezultat deja calculat."""
    from neuroq import viz

    m = result.matrix
    kind = ROMANIAN[result.encoding]
    _title(1, f"Datele: {len(m.read_ids)} citiri × {len(m.snp_ids)} SNP-uri "
              f"(0 = referință, 1 = alternativ, . = neobservat)")
    print(_matrix_text(m.values, m.read_ids, [str(i + 1) for i in range(len(m.snp_ids))],
                       missing_dot=True))
    sink.add(1, "citiri", viz.plot_reads_matrix(m.values, m.read_ids, m.snp_ids))

    _title(2, f"Voturile pe perechi, codificarea pe {kind}: a = la fel, d = opus")
    print("a:\n" + _matrix_text(result.a, result.labels, result.labels))
    print("d:\n" + _matrix_text(result.d, result.labels, result.labels))
    sink.add(2, "a_d", viz.plot_evidence(result.a, result.d, result.labels))

    p = result.problem
    _title(3, f"Hamiltonianul Ising: J = (d - a) / 2, constanta = {p.offset:g}")
    print(_matrix_text(p.J, result.labels, result.labels))
    sink.add(3, "J", viz.plot_ising(p.J, p.offset, result.labels))

    _title(4, f"Graful: {p.n_qubits} variabile, {int((p.J != 0).sum())} cuplaje nenule")
    sink.add(4, "graf", viz.plot_graph(p.J, result.labels))

    qc = result.circuit
    two_qubit = sum(1 for item in qc.data if item.operation.num_qubits == 2)
    _title(5, f"Circuitul QAOA (p = 1): {qc.num_qubits} qubiți, {two_qubit} porți cu doi "
              f"qubiți, adâncime {qc.depth()}" + (" (un qubit fixat)" if result.fixed else ""))
    print(f"   porți: {dict(qc.count_ops())}")
    sink.add(5, "circuit", viz.plot_circuit(qc))

    history = result.cost_history
    _title(6, f"Antrenarea: {sum(len(h) for h in history)} evaluări COBYLA "
              f"(costul a fost afișat live mai sus)")
    print(f"   cel mai bun start: {result.best_start + 1}, <H> final = "
          f"{min(history[result.best_start]):.4f}, optim exact = {result.exact.energy:g}")
    sink.add(6, "cost_qaoa", viz.plot_qaoa_cost(history, result.exact.energy, result.best_start))

    _title(7, "Cele mai probabile 10 șiruri măsurate (2048 de măsurători)")
    for row in result.top_samples:
        mark = "  <- optim" if np.isclose(row["energy"], result.exact.energy) else ""
        print(f"   {row['bitstring']}   {row['count']:5d} apariții   E = {row['energy']:g}{mark}")
    sink.add(7, "top10", viz.plot_top_samples(result.top_samples, result.exact.energy))
    sink.add(8, "graf_solutie", viz.plot_graph(
        p.J, result.labels, bitstring=result.full_bitstring,
        energy=p.energy(result.full_bitstring), title="Soluția QAOA pe graf"))

    haps = result.haplotypes
    _title(8, "Copiile decodate față de cele adevărate")
    print(f"   găsite:    A = {_bits(haps.copy_a)}   B = {_bits(haps.copy_b)}")
    if m.true_haplotypes is not None:
        print(f"   adevărate: A = {_bits(m.true_haplotypes[0])}   B = {_bits(m.true_haplotypes[1])}")
        sink.add(9, "haplotipuri", viz.plot_haplotypes(haps, m.true_haplotypes))

    _title(9, "Rezultatul")
    for key in ("pair_cost", "accuracy", "switch_errors", "mec_cost", "optimum_frequency"):
        if key in result.metrics:
            print(f"   {key:18s} = {result.metrics[key]:g}")
    if result.other_encoding:
        other = result.other_encoding
        print(f"\n   Codificarea pe {ROMANIAN[other['encoding']]} nu se rulează aici: "
              f"{other['n_qubits']} qubiți, {other['n_couplings']} cuplaje nenule, "
              f"{other['rzz_per_layer']} porți rzz pe strat.")


def main(argv: Optional[list[str]] = None) -> list[Path]:
    sys.stdout.reconfigure(encoding="utf-8")
    sys.stderr.reconfigure(encoding="utf-8")
    parser = argparse.ArgumentParser(description="Demonstrația NeuroQ, etapă cu etapă.")
    parser.add_argument("--instance", choices=sorted(INSTANCES), default="notebook")
    parser.add_argument("--encoding", choices=["pozitii", "citiri"], default=None,
                        help="implicit: pozitii (notebook, illumina), citiri (nanopore)")
    parser.add_argument("--fix", action="store_true", help="fixează primul qubit (fix_variable)")
    parser.add_argument("--show", action="store_true", help="deschide figurile pe rând")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    if not args.show:
        matplotlib.use("Agg")

    start = time.perf_counter()
    try:
        kind = resolve_encoding(args.instance, args.encoding)
    except ValueError as error:
        if args.encoding in ENCODING_NAMES:
            other = encoding_resources(load_case(args.instance).matrix,
                                       ENCODING_NAMES[args.encoding])
            print(f"Codificarea pe {args.encoding} pentru {args.instance}: {other['n_qubits']} "
                  f"qubiți, {other['n_couplings']} cuplaje nenule, "
                  f"{other['rzz_per_layer']} porți rzz pe strat (nu se rulează).")
        parser.exit(2, f"eroare: {error}\n")
    print(f"NeuroQ, demonstrație: instanța {args.instance}, codificarea pe {ROMANIAN[kind]}"
          + (", cu fix_variable" if args.fix else ""))
    print("\nQAOA rulează acum (costul live, la fiecare evaluare):", flush=True)
    result = run_stages(load_case(args.instance), kind, args.fix, callback=live_cost)
    suffix = ROMANIAN[kind] + ("_fix" if args.fix else "")
    sink = FigureSink(args.out / f"demo_{args.instance}", suffix, args.show)
    present(result, sink)
    print(f"\nGata în {time.perf_counter() - start:.1f} s; {len(sink.saved)} figuri în "
          f"{sink.directory}", flush=True)
    if args.show:
        input("Apasă Enter ca să închizi figurile...")
    return sink.saved


if __name__ == "__main__":
    main()
