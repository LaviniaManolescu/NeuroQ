"""Experimentul principal: ambele codificări cu toți solverii pe trei instanțe separate.

Instanțe (fără ferestre și fără subeșantionare, nicio comparație între ele):
- exemplul din notebook (5 SNP-uri, 9 citiri);
- Illumina sintetic: 10 SNP-uri, 24 de citiri;
- Nanopore sintetic: 18 SNP-uri, 9 citiri.
Optimul de referință vine din ExactSolver până la 16 qubiți și din exact_vectorized
(ClassicalSolver, strategia "exact") până la 24. results.json și results.js se scriu doar
la final, prin fișiere temporare, ca o întrerupere să nu lase un fișier pe jumătate scris.
Rulare: python -m neuroq.experiments.resource_comparison
"""

from __future__ import annotations

import json
import os
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from neuroq.core import Haplotypes, IsingProblem, ReadMatrix, ResultRow, SelectionParams, Solution
from neuroq.encodings import Encoding, PositionEncoding, ReadEncoding
from neuroq.evaluation import accuracy, evaluate, mec_cost, switch_errors, top_samples
from neuroq.execution import IdealSimulator
from neuroq.experiments.notebook_example import haplotype_text, write_results
from neuroq.selection import WindowSelector
from neuroq.solvers import ClassicalSolver, ExactSolver, GreedyPhaser, QAOASolver
from neuroq.sources import SyntheticSource
from neuroq.sources.synthetic import ILLUMINA, NANOPORE, TechnologyProfile
from neuroq.transforms import fix_variable

SEED = 71
EXACT_SOLVER_MAX_QUBITS = 16
ILLUMINA_SNPS, ILLUMINA_READS = 10, 24
NANOPORE_SNPS, NANOPORE_READS = 18, 9
OUTPUT_DIR = Path("results")


@dataclass
class Case:
    """O matrice de citiri și descrierea ei pentru `instance`."""

    matrix: ReadMatrix
    instance: dict


def _whole_matrix(source: SyntheticSource) -> ReadMatrix:
    data = source.load()
    params = SelectionParams(first_snp=0, n_snps=data.values.shape[1], seed=SEED)
    return WindowSelector().select(data, params)


def _instance(case_id: str, matrix: ReadMatrix, profile: Optional[TechnologyProfile]) -> dict:
    return {
        "id": case_id,
        "source": "synthetic",
        "technology": matrix.technology,
        "profile": profile.name if profile else None,
        "first_snp": 0,
        "n_snps": len(matrix.snp_ids),
        "n_reads": len(matrix.read_ids),
        "error_rate": profile.error_rate if profile else None,
        "seed": SEED,
    }


def notebook_case() -> Case:
    matrix = _whole_matrix(SyntheticSource.notebook_example())
    return Case(matrix, _instance("notebook", matrix, None))


def illumina_case() -> Case:
    source = SyntheticSource.from_profile(ILLUMINA, ILLUMINA_SNPS, SEED, n_reads=ILLUMINA_READS)
    matrix = _whole_matrix(source)
    return Case(matrix, _instance("illumina", matrix, ILLUMINA))


def nanopore_case() -> Case:
    source = SyntheticSource.from_profile(NANOPORE, NANOPORE_SNPS, SEED, n_reads=NANOPORE_READS)
    matrix = _whole_matrix(source)
    return Case(matrix, _instance("nanopore", matrix, NANOPORE))


def _haplotype_details(haps: Haplotypes) -> dict:
    return {
        "copy_a": haplotype_text(haps.copy_a),
        "copy_b": haplotype_text(haps.copy_b),
        "read_assignment": [int(v) for v in haps.read_assignment],
    }


def _solver_row(
    case: Case, method: str, problem: IsingProblem, solution: Solution, optimum: float,
    transforms: list[str], extra_metrics: Optional[dict] = None,
    extra_details: Optional[dict] = None,
) -> ResultRow:
    haps = problem.decoder(solution.bitstring)
    metrics = evaluate(problem, solution, haps, case.matrix.true_haplotypes, optimum,
                       case.matrix.values)
    if "components" in solution.metadata:
        metrics["components"] = solution.metadata["components"]
    metrics.update(extra_metrics or {})
    details: dict = {"haplotypes": _haplotype_details(haps), **(extra_details or {})}
    if solution.samples:
        details["top_samples"] = top_samples(solution.samples)
    return ResultRow(instance=case.instance, method=method, encoding=problem.variable_kind,
                     transforms=transforms, metrics=metrics, resources=solution.resources,
                     details=details)


def _greedy_row(case: Case) -> ResultRow:
    haps = GreedyPhaser().phase(case.matrix, order_seed=SEED)
    metrics: dict = {"mec_cost": mec_cost(haps.copy_a, case.matrix.values)}
    if case.matrix.true_haplotypes is not None:
        truth_a = case.matrix.true_haplotypes[0]
        metrics["accuracy"] = accuracy(haps.copy_a, truth_a)
        metrics["switch_errors"] = switch_errors(haps.copy_a, truth_a)
    return ResultRow(instance=case.instance, method="greedy", encoding=None, transforms=[],
                     metrics=metrics, resources=None,
                     details={"haplotypes": _haplotype_details(haps)})


def _exact(problem: IsingProblem) -> tuple[Solution, str]:
    if problem.n_qubits <= EXACT_SOLVER_MAX_QUBITS:
        return ExactSolver(max_qubits=EXACT_SOLVER_MAX_QUBITS).solve(problem), "ExactSolver"
    return ClassicalSolver("exact").solve(problem), "exact_vectorized"


def _timed_qaoa(problem: IsingProblem, label: str) -> tuple[Solution, float]:
    start = time.perf_counter()
    solution = QAOASolver(IdealSimulator(), seed=SEED).solve(problem)
    seconds = time.perf_counter() - start
    print(f"  QAOA {label}: {problem.n_qubits} qubiți, {seconds:.1f} s", flush=True)
    return solution, seconds


def _encoding_rows(case: Case, encoding: Encoding) -> list[ResultRow]:
    problem = encoding.encode(case.matrix)
    label = f"{case.instance['id']}/{problem.variable_kind}"
    exact, backend = _exact(problem)
    fixed = fix_variable(problem)
    qaoa, qaoa_seconds = _timed_qaoa(problem, label)
    qaoa_fixed, fixed_seconds = _timed_qaoa(fixed, label + " fix")
    return [
        _solver_row(case, "exact", problem, exact, exact.energy, [],
                    extra_details={"exact_backend": backend}),
        _solver_row(case, "classical", problem, ClassicalSolver(seed=SEED).solve(problem),
                    exact.energy, []),
        _solver_row(case, "qaoa", problem, qaoa, exact.energy, [],
                    extra_metrics={"qaoa_seconds": round(qaoa_seconds, 2)}),
        _solver_row(case, "qaoa", fixed, qaoa_fixed, exact.energy, ["fix_variable"],
                    extra_metrics={"qaoa_seconds": round(fixed_seconds, 2)}),
    ]


def run_comparison(cases: list[Case]) -> list[ResultRow]:
    rows: list[ResultRow] = []
    for case in cases:
        print(f"{case.instance['id']}: {case.instance['n_snps']} SNP-uri, "
              f"{case.instance['n_reads']} citiri", flush=True)
        rows.append(_greedy_row(case))
        for encoding in (PositionEncoding(), ReadEncoding()):
            rows.extend(_encoding_rows(case, encoding))
    return rows


def write_outputs(rows: list[ResultRow], directory: Path) -> tuple[Path, Path]:
    """Scrie results.json și results.js (window.RESULTS = ...) atomic, prin fișiere .tmp."""
    directory.mkdir(parents=True, exist_ok=True)
    json_path, js_path = directory / "results.json", directory / "results.js"
    json_tmp, js_tmp = directory / "results.json.tmp", directory / "results.js.tmp"
    write_results(rows, json_tmp, note=(
        "Faza 2: notebook, Illumina (10 SNP-uri, 24 citiri) și Nanopore (18 SNP-uri, "
        "9 citiri), ambele codificări, toți solverii, simulator ideal."))
    payload = json.loads(json_tmp.read_text(encoding="utf-8"))
    js_tmp.write_text("window.RESULTS = " + json.dumps(payload, indent=2, ensure_ascii=False)
                      + ";\n", encoding="utf-8")
    os.replace(json_tmp, json_path)
    os.replace(js_tmp, js_path)
    return json_path, js_path


def main() -> None:
    # consola Windows redirecționată folosește cp1252, care nu are diacriticele românești
    sys.stdout.reconfigure(encoding="utf-8")
    start = time.perf_counter()
    rows = run_comparison([notebook_case(), illumina_case(), nanopore_case()])
    json_path, js_path = write_outputs(rows, OUTPUT_DIR)
    print(f"{len(rows)} rânduri în {time.perf_counter() - start:.1f} s: {json_path}, {js_path}",
          flush=True)


if __name__ == "__main__":
    main()
