"""Experimentul principal, partea mică: ambele codificări cu toți solverii pe ferestre mici.

Matrici: exemplul din notebook, o matrice Illumina și una Nanopore (profiluri sintetice).
Fiecare fereastră are cel mult 16 qubiți pe codificarea mai mare, ca să încapă în simulator.
Rulare: python -m neuroq.experiments.resource_comparison
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from neuroq.core import Haplotypes, IsingProblem, ReadMatrix, ResultRow, SelectionParams, Solution
from neuroq.encodings import Encoding, PositionEncoding, ReadEncoding
from neuroq.evaluation import accuracy, evaluate, mec_cost, switch_errors, top_samples
from neuroq.execution import IdealSimulator
from neuroq.experiments.notebook_example import haplotype_text, write_results
from neuroq.selection import Downsampler, WindowSelector
from neuroq.solvers import ClassicalSolver, ExactSolver, GreedyPhaser, QAOASolver
from neuroq.sources import SyntheticSource
from neuroq.sources.synthetic import ILLUMINA, NANOPORE, TechnologyProfile
from neuroq.transforms import fix_variable

SEED = 71
MAX_QUBITS = 16
MIN_WINDOW_SNPS = 4
SOURCE_SNPS = 30  # regiunea sintetică din care se taie fereastra
WINDOW_START = 10  # fereastra începe în interiorul regiunii, ca citirile să intre și din vecini
NANOPORE_WINDOW_SNPS = 12
NANOPORE_READS = 12
OUTPUT_DIR = Path("results")


@dataclass
class Case:
    """O matrice de citiri gata tăiată și descrierea ei pentru `instance`."""

    matrix: ReadMatrix
    instance: dict


def notebook_case() -> Case:
    data = SyntheticSource.notebook_example().load()
    params = SelectionParams(first_snp=0, n_snps=5, seed=SEED)
    matrix = WindowSelector().select(data, params)
    return Case(matrix, _instance("notebook-5snp", "notebook", params, matrix, None, None))


def illumina_case() -> Case:
    """Cea mai mică fereastră (>= 4 SNP-uri) cu <= 16 citiri; altfel subeșantionare la 16."""
    data = SyntheticSource.from_profile(ILLUMINA, SOURCE_SNPS, SEED).load()
    selector = WindowSelector()
    for n_snps in range(MIN_WINDOW_SNPS, MAX_QUBITS + 1):
        params = SelectionParams(first_snp=WINDOW_START, n_snps=n_snps, seed=SEED)
        matrix = selector.select(data, params)
        if len(matrix.read_ids) <= MAX_QUBITS:
            return Case(matrix, _instance("illumina-small", ILLUMINA.name, params, matrix,
                                          ILLUMINA, None))
    params = SelectionParams(first_snp=WINDOW_START, n_snps=MIN_WINDOW_SNPS,
                             target_reads=MAX_QUBITS, seed=SEED)
    matrix = selector.select(data, params)
    return Case(matrix, _instance("illumina-small", ILLUMINA.name, params, matrix,
                                  ILLUMINA, MAX_QUBITS))


def nanopore_case() -> Case:
    data = SyntheticSource.from_profile(NANOPORE, SOURCE_SNPS, SEED).load()
    params = SelectionParams(first_snp=WINDOW_START, n_snps=NANOPORE_WINDOW_SNPS,
                             target_reads=NANOPORE_READS, seed=SEED)
    matrix = WindowSelector().select(data, params)
    return Case(matrix, _instance("nanopore-small", NANOPORE.name, params, matrix,
                                  NANOPORE, NANOPORE_READS))


def _instance(
    case_id: str, technology: str, params: SelectionParams, matrix: ReadMatrix,
    profile: Optional[TechnologyProfile], downsampled_to: Optional[int],
) -> dict:
    return {
        "id": case_id,
        "source": "synthetic",
        "technology": technology,
        "profile": profile.name if profile else None,
        "first_snp": params.first_snp,
        "n_snps": params.n_snps,
        "n_reads": len(matrix.read_ids),
        "downsampled_to": downsampled_to,
        "error_rate": profile.error_rate if profile else None,
        "seed": SEED,
    }


def _haplotype_details(haps: Haplotypes) -> dict:
    return {
        "copy_a": haplotype_text(haps.copy_a),
        "copy_b": haplotype_text(haps.copy_b),
        "read_assignment": [int(v) for v in haps.read_assignment],
    }


def _solver_row(
    case: Case, method: str, problem: IsingProblem, solution: Solution,
    optimum: float, transforms: list[str],
) -> ResultRow:
    haps = problem.decoder(solution.bitstring)
    metrics = evaluate(problem, solution, haps, case.matrix.true_haplotypes, optimum,
                       case.matrix.values)
    if "components" in solution.metadata:
        metrics["components"] = solution.metadata["components"]
    details: dict = {"haplotypes": _haplotype_details(haps)}
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


def _encoding_rows(case: Case, encoding: Encoding) -> list[ResultRow]:
    problem = encoding.encode(case.matrix)
    if problem.n_qubits > MAX_QUBITS:
        raise ValueError(f"{case.instance['id']}: {problem.n_qubits} qubiți > {MAX_QUBITS}")
    exact = ExactSolver(max_qubits=MAX_QUBITS).solve(problem)
    fixed = fix_variable(problem)
    qaoa = QAOASolver(IdealSimulator(), seed=SEED)
    return [
        _solver_row(case, "exact", problem, exact, exact.energy, []),
        _solver_row(case, "classical", problem, ClassicalSolver(seed=SEED).solve(problem),
                    exact.energy, []),
        _solver_row(case, "qaoa", problem, qaoa.solve(problem), exact.energy, []),
        _solver_row(case, "qaoa", fixed, qaoa.solve(fixed), exact.energy, ["fix_variable"]),
    ]


def run_comparison(cases: list[Case]) -> list[ResultRow]:
    rows: list[ResultRow] = []
    for case in cases:
        rows.append(_greedy_row(case))
        for encoding in (PositionEncoding(), ReadEncoding()):
            rows.extend(_encoding_rows(case, encoding))
    return rows


def write_outputs(rows: list[ResultRow], directory: Path) -> tuple[Path, Path]:
    """Scrie results.json și results.js (window.RESULTS = ...) cu același conținut."""
    json_path = directory / "results.json"
    write_results(rows, json_path, note=(
        "Faza 2: ferestre mici (<= 16 qubiți), ambele codificări, toți solverii, simulator ideal."))
    js_path = directory / "results.js"
    payload = json.loads(json_path.read_text(encoding="utf-8"))
    js_path.write_text("window.RESULTS = " + json.dumps(payload, indent=2, ensure_ascii=False)
                       + ";\n", encoding="utf-8")
    return json_path, js_path


def main() -> None:
    start = time.perf_counter()
    rows = run_comparison([notebook_case(), illumina_case(), nanopore_case()])
    json_path, js_path = write_outputs(rows, OUTPUT_DIR)
    print(f"{len(rows)} rânduri în {time.perf_counter() - start:.1f} s: {json_path}, {js_path}")


if __name__ == "__main__":
    main()
