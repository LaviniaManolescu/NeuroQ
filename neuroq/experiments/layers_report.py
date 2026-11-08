"""Raportul QAOA pe p = 1, 2, 3 (antrenare în lanț, INTERP) pe cele trei instanțe.

Rulare:
    python -m neuroq.experiments.layers_report [--p-max 3] [--fix] [--warm]
Tipărește tabelul și îl scrie în results/layers_report.md. Nu modifică nimic altceva.
"""

from __future__ import annotations

import argparse
import itertools
import sys
import time
from pathlib import Path

from neuroq.encodings import PositionEncoding, ReadEncoding
from neuroq.execution import IdealSimulator, transpiled_depth
from neuroq.experiments.resource_comparison import SEED
from neuroq.experiments.stages import INSTANCES, ROMANIAN, load_case, resolve_encoding
from neuroq.solvers import ExactSolver
from neuroq.solvers.qaoa import train_layers, trained_circuit
from neuroq.transforms import fix_variable

OUTPUT = Path("results") / "layers_report.md"
HEADER = ("instanța", "p", "optim (freq)", "oglinda (freq)", "optim+oglindă", "cel mai probabil",
          "freq", "cost", "adâncime", "adâncime transpilată")


def _mirror(bits: str) -> str:
    return bits.translate(str.maketrans("01", "10"))


def _minimum_states(problem) -> list[str]:
    """Toate șirurile cu energia minimă (brute force; n ≤ 10)."""
    n = problem.n_qubits
    energies = {"".join(bits): problem.energy("".join(bits))
                for bits in itertools.product("01", repeat=n)}
    low = min(energies.values())
    return [b for b, e in energies.items() if abs(e - low) < 1e-9]


def rows_for(instance: str, p_max: int, fix: bool, warm: bool) -> list[tuple]:
    case = load_case(instance)
    kind = resolve_encoding(instance, None)
    problem = (PositionEncoding() if kind == "position" else ReadEncoding()).encode(case.matrix)
    solved = fix_variable(problem) if fix else problem
    optimum = ExactSolver().solve(solved)
    minima = _minimum_states(solved)
    layers = train_layers(solved, IdealSimulator(), p_max, seed=SEED, warm_start=warm)
    rows = []
    for p, solution in enumerate(layers, start=1):
        samples, shots = solution.samples, sum(solution.samples.values())
        mirror = _mirror(optimum.bitstring)
        top = max(samples, key=lambda b: (samples[b], b))
        circuit = trained_circuit(
            solved, solution.metadata["angles"], p,
            (solution.metadata["warm_start"] or {}).get("c"))
        rows.append((
            instance, p,
            f"{optimum.bitstring} {100 * samples.get(optimum.bitstring, 0) / shots:.2f}%",
            "-" if fix else f"{mirror} {100 * samples.get(mirror, 0) / shots:.2f}%",
            f"{100 * sum(samples.get(b, 0) for b in minima) / shots:.2f}% ({len(minima)} stări)",
            top, f"{100 * samples[top] / shots:.2f}%", f"{solved.energy(top):g}",
            circuit.depth(), transpiled_depth(circuit, SEED),
        ))
    return rows


def as_markdown(rows: list[tuple], note: str) -> str:
    lines = [f"# QAOA pe p = 1, 2, 3 (INTERP)\n\n{note}\n",
             "| " + " | ".join(HEADER) + " |", "|" + "---|" * len(HEADER)]
    lines += ["| " + " | ".join(str(c) for c in row) + " |" for row in rows]
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> None:
    sys.stdout.reconfigure(encoding="utf-8")
    parser = argparse.ArgumentParser(description="Raport QAOA pe straturi.")
    parser.add_argument("--p-max", type=int, default=3)
    parser.add_argument("--fix", action="store_true")
    parser.add_argument("--warm", action="store_true")
    parser.add_argument("--out", type=Path, default=OUTPUT)
    args = parser.parse_args(argv)
    start = time.perf_counter()
    rows: list[tuple] = []
    for instance in INSTANCES:
        print(f"{instance} ({ROMANIAN[resolve_encoding(instance, None)]}) ...", flush=True)
        rows += rows_for(instance, args.p_max, args.fix, args.warm)
    note = (f"Seed {SEED}, 2048 măsurători, COBYLA din 2 porniri (prima INTERP din p−1) pe "
            f"statevector; codificarea mică a fiecărei instanțe"
            + (", cu un qubit fixat" if args.fix else "") + (", warm-start" if args.warm else "")
            + ". Adâncimea transpilată: rz, sx, x, cx, nivel 3, fără hartă de qubiți.")
    table = as_markdown(rows, note)
    print("\n" + table)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(table, encoding="utf-8")
    print(f"scris în {args.out}, {time.perf_counter() - start:.0f} s")


if __name__ == "__main__":
    main()
