"""Demonstrația live: rulează pipeline-ul și arată fiecare etapă în terminal și în figuri.

Rulare:
    python -m neuroq.experiments.demo --instance notebook
    python -m neuroq.experiments.demo --instance notebook --encoding citiri --fix --show
    python -m neuroq.experiments.demo --instance nanopore --warm
    python -m neuroq.experiments.demo --instance all
    python -m neuroq.experiments.demo --instance illumina --p 3 [--fix] [--warm]
    python -m neuroq.experiments.demo --all-variants [--p 3]
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

    def __init__(self, directory: Path, suffix: str, show: bool, base: Optional[str] = None,
                 always_p: bool = False) -> None:
        self.directory = directory
        self.suffix = suffix
        self.base = base  # sufixul fără p; None = se deduce din sufix (cazul --p)
        self.always_p = always_p  # --all-variants: figurile pe p au mereu p în nume
        self.show = show
        self.saved: list[Path] = []

    def add(self, number: int, name: str, fig, suffix: Optional[str] = None) -> None:
        from neuroq.viz import save_figure

        if fig is None:
            print("   (circuitul e prea mare pentru desen; vezi numerele de mai sus)")
            return
        path = save_figure(fig, self.directory / f"{number:02d}_{name}_{suffix or self.suffix}.png")
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
    _title(5, f"Circuitul QAOA (p = {result.p}): {qc.num_qubits} qubiți, {two_qubit} porți cu doi "
              f"qubiți, adâncime {qc.depth()}" + (" (un qubit fixat)" if result.fixed else "")
              + (" (warm-start)" if result.warm else ""))
    print(f"   porți: {dict(qc.count_ops())}")
    sink.add(5, "circuit", viz.plot_circuit(qc))

    # fără --p (p = 1) numele figurilor rămân cele de dinainte; cu --p 3 fiecare p are figurile lui
    per_p = result.p > 1 or sink.always_p
    steps = result.layers if per_p else [result.qaoa]
    base = sink.base if sink.base is not None else sink.suffix.removesuffix(f"_p{result.p}")
    for k, step in enumerate(steps, start=1):
        tag = "_".join(part for part in (base, f"p{k}") if part) if per_p else None
        label = f" (p = {k})" if per_p else ""
        history = step.metadata["cost_history"]
        best = step.metadata["best_start"]
        _title(6, f"Antrenarea{label}: {sum(len(h) for h in history)} evaluări COBYLA "
                  f"(costul a fost afișat live mai sus)")
        print(f"   cel mai bun start: {best + 1}, <H> final = {min(history[best]):.4f}, "
              f"optim exact = {result.exact.energy:g}")
        sink.add(6, "cost_qaoa", viz.plot_qaoa_cost(history, result.exact.energy, best), tag)

        _title(7, f"Cele mai probabile 10 șiruri măsurate{label} (2048 de măsurători)")
        for row in step.metadata["top_samples"]:
            mark = "  <- optim" if np.isclose(row["energy"], result.exact.energy) else ""
            print(f"   {row['bitstring']}   {row['count']:5d} apariții   E = {row['energy']:g}{mark}")
        sink.add(7, "top10", viz.plot_top_samples(step.metadata["top_samples"], result.exact.energy),
                 tag)
        if per_p:
            energies = {b: result.solved.energy(b) for b in step.samples}
            sink.add(10, "incredere", viz.plot_confidence(
                step.samples, energies, result.exact.energy, step.bitstring, p=k), tag)
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


RUN_ALL = "all"


def run_instance(instance: str, kind: str, fix: bool, warm: bool, out: Path, show: bool,
                 p: int = 1):
    """Rulează demonstrația pe o instanță: tipărește etapele și salvează figurile."""
    start = time.perf_counter()
    print(f"NeuroQ, demonstrație: instanța {instance}, codificarea pe {ROMANIAN[kind]}"
          + (", cu fix_variable" if fix else "") + (", cu warm-start" if warm else "")
          + (f", p = 1..{p} (INTERP)" if p > 1 else ""))
    print("\nQAOA rulează acum (costul live, la fiecare evaluare):", flush=True)
    result = run_stages(load_case(instance), kind, fix, callback=live_cost, warm=warm, p=p)
    suffix = (ROMANIAN[kind] + ("_fix" if fix else "") + ("_warm" if warm else "")
              + (f"_p{p}" if p > 1 else ""))
    sink = FigureSink(out / f"demo_{instance}", suffix, show)
    present(result, sink)
    print(f"\nGata în {time.perf_counter() - start:.1f} s; {len(sink.saved)} figuri în "
          f"{sink.directory}", flush=True)
    return sink, result


def print_summary(runs: list, seconds: float) -> None:
    """Tabel cu toate instanțele rulate, la finalul lui --instance all."""
    print(f"\n{'#' * 64}\n Rezumat: {len(runs)} instanțe în {seconds:.1f} s\n{'#' * 64}")
    print(f"{'instanța':10s}{'codificarea':13s}{'qubiți':>7s}{'p':>3s}  {'șirul ales':14s}{'E':>6s}"
          f"{'acuratețe':>11s}{'QAOA':>8s}")
    for sink, result in runs:
        accuracy = result.metrics.get("accuracy")
        print(f"{result.instance:10s}{ROMANIAN[result.encoding]:13s}{result.solved.n_qubits:7d}{result.p:3d}  "
              f"{result.full_bitstring:14s}{result.problem.energy(result.full_bitstring):6g}"
              f"{'-' if accuracy is None else format(accuracy, '.0%'):>11s}"
              f"{result.qaoa_seconds:7.1f}s")
    print("\nHistograma măsurătorilor (aceleași date ca la etapa 5 din interfață):")
    for sink, result in runs:
        top10 = [path for path in sink.saved if "_top10_" in path.name][-1]  # ultimul p
        print(f"   {result.instance:10s}{top10}")


VARIANTS = {"qaoa": {}, "qaoa_fixed": {"fix": True}, "ws_qaoa": {"warm": True}}
COMPARISON_FILE = "comparatie_variante.png"
TABLE_FILE = "demo_variante.md"


def _mirror(bits: str) -> str:
    return bits.translate(str.maketrans("01", "10"))


def _optimum_frequency(result: StagesResult, step) -> float:
    shots = sum(step.samples.values())
    hits = sum(c for b, c in step.samples.items()
               if np.isclose(result.solved.energy(b), result.exact.energy))
    return hits / shots


def _status(result: StagesResult, full: str) -> str:
    if full == result.exact.bitstring:
        return "optimul"
    if full == _mirror(result.exact.bitstring):
        return "oglinda optimului"
    if np.isclose(result.problem.energy(full), result.exact.energy):
        return "alt optim (degenerat)"
    return "nu"


def variant_rows(result: StagesResult, variant: str) -> list[tuple]:
    """Un rând pe p: șirul ales (complet), energia lui, dacă e optimul, frecvența optimului."""
    from neuroq.experiments.stages import _with_fixed_bits

    rows = []
    for p, step in enumerate(result.layers, start=1):
        full = _with_fixed_bits(step.bitstring, result.solved.fixed_variables)
        rows.append((result.instance, variant, p, full, result.problem.energy(full),
                     _status(result, full), _optimum_frequency(result, step)))
    return rows


def table_text(rows: list[tuple], markdown: bool) -> str:
    head = ("instanța", "varianta", "p", "șirul ales", "energia", "este", "frecvența optimului")
    cells = [(i, v, str(p), s, f"{e:g}", st, f"{100 * f:.2f}%") for i, v, p, s, e, st, f in rows]
    if markdown:
        lines = ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
        return "\n".join(lines + ["| " + " | ".join(c) + " |" for c in cells]) + "\n"
    widths = [max(len(h), *(len(c[k]) for c in cells)) for k, h in enumerate(head)]
    fmt = "  ".join(f"{{:{w}s}}" for w in widths)
    return "\n".join([fmt.format(*head)] + [fmt.format(*c) for c in cells])


def run_all_variants(out: Path, p: int) -> tuple[list[Path], Path]:
    """Cele 3 instanțe × cele 3 variante (× p = 1..p): figuri separate pe variantă, o figură
    comparativă pe instanță, un tabel final."""
    import contextlib
    import io

    import matplotlib.pyplot as plt

    from neuroq import viz
    from neuroq.viz import save_figure

    start = time.perf_counter()
    saved: list[Path] = []
    rows: list[tuple] = []
    for instance in INSTANCES:
        kind = resolve_encoding(instance, None)
        columns, optimum = [], 0.0
        for variant, options in VARIANTS.items():
            began = time.perf_counter()
            result = run_stages(load_case(instance), kind, p=p, **options)
            sink = FigureSink(out / f"demo_{instance}" / variant, f"p{p}", False, base="",
                              always_p=True)
            with contextlib.redirect_stdout(io.StringIO()):  # detaliile etapelor: un rând mai jos
                present(result, sink)
            saved += sink.saved
            rows += variant_rows(result, variant)
            final, last = result.layers[-1], rows[-1]
            optimum = result.exact.energy
            columns.append({
                "title": f"{variant}\nales {last[3]}, E = {last[4]:g}\n"
                         f"optimul apare în {100 * last[6]:.1f}% din măsurători",
                "top": final.metadata["top_samples"], "histories": final.metadata["cost_history"],
                "best_start": final.metadata["best_start"]})
            print(f"{instance:9s} {variant:11s} p = 1..{p}: {len(sink.saved)} figuri în "
                  f"{sink.directory} ({time.perf_counter() - began:.1f} s)", flush=True)
        fig = viz.plot_variant_comparison(
            columns, optimum, f"{instance}: cele trei variante QAOA (p = {p})")
        path = save_figure(fig, out / f"demo_{instance}" / COMPARISON_FILE)
        saved.append(path)
        plt.close(fig)
        print(f"{instance:9s} figura comparativă: {path}", flush=True)
    print(f"\n{table_text(rows, markdown=False)}")
    table_path = out.parent / TABLE_FILE
    table_path.write_text(
        f"# Variantele QAOA pe cele trei instanțe (p = 1..{p})\n\n"
        + table_text(rows, markdown=True), encoding="utf-8")
    print(f"\nTabel: {table_path}; {len(saved)} figuri; {time.perf_counter() - start:.0f} s")
    return saved, table_path


def main(argv: Optional[list[str]] = None) -> list[Path]:
    sys.stdout.reconfigure(encoding="utf-8")
    sys.stderr.reconfigure(encoding="utf-8")
    parser = argparse.ArgumentParser(description="Demonstrația NeuroQ, etapă cu etapă.")
    parser.add_argument("--instance", choices=list(INSTANCES) + [RUN_ALL], default="notebook",
                        help="o instanță sau all (cele trei, una după alta, cu rezumat la final)")
    parser.add_argument("--encoding", choices=["pozitii", "citiri"], default=None,
                        help="implicit: pozitii (notebook, illumina), citiri (nanopore)")
    parser.add_argument("--fix", action="store_true", help="fixează primul qubit (fix_variable)")
    parser.add_argument("--warm", action="store_true",
                        help="warm-start QAOA, pornit din relaxarea continuă")
    parser.add_argument("--p", type=int, choices=[1, 2, 3], default=1,
                        help="straturile QAOA; la p > 1 antrenează p = 1..p în lanț (INTERP)")
    parser.add_argument("--all-variants", action="store_true",
                        help="cele 3 instanțe × qaoa, qaoa_fixed, ws_qaoa (× p = 1..--p), "
                             "figurile pe variantă în demo_<instanță>/<variantă>/")
    parser.add_argument("--show", action="store_true", help="deschide figurile pe rând")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    if not args.show:
        matplotlib.use("Agg")

    if args.all_variants:
        if args.fix or args.warm or args.encoding or args.instance != "notebook":
            print("--all-variants rulează toate instanțele și variantele; --instance, --encoding, "
                  "--fix și --warm se ignoră.\n")
        run_all_variants(args.out, args.p)
        return []
    if args.instance == RUN_ALL:
        if args.encoding:
            print("--encoding se ignoră la --instance all: fiecare instanță rulează pe codificarea ei "
                  "mică.\n")
        start = time.perf_counter()
        runs = []
        for name in INSTANCES:
            runs.append(run_instance(name, resolve_encoding(name, None), args.fix, args.warm,
                                     args.out, args.show, args.p))
            print()
        print_summary(runs, time.perf_counter() - start)
        saved = [path for sink, _ in runs for path in sink.saved]
    else:
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
        sink, _ = run_instance(args.instance, kind, args.fix, args.warm, args.out, args.show,
                               args.p)
        saved = sink.saved
    if args.show:
        input("Apasă Enter ca să închizi figurile...")
    return saved


if __name__ == "__main__":
    main()
