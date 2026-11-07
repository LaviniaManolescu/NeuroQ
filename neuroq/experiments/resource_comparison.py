"""Experimentul principal: toți solverii pe codificarea mică, un rând de resurse pe cea mare.

Instanțe (fără ferestre și fără subeșantionare, nicio comparație între ele):
- exemplul din notebook (5 SNP-uri, 9 citiri): ambele codificări, complet;
- Illumina sintetic (10 SNP-uri, 24 de citiri): complet pe poziții (10 qubiți); pe citiri
  (24 de qubiți) doar rândul de resurse, plus durata QAOA măsurată într-o rulare anterioară;
- Nanopore sintetic (18 SNP-uri, 9 citiri): complet pe citiri (9 qubiți); pe poziții
  (18 qubiți) doar rândul de resurse.
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
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from neuroq.core import (
    Haplotypes, IsingProblem, ReadMatrix, ResourceReport, ResultRow, SelectionParams, Solution,
)
from neuroq.encodings import Encoding, PositionEncoding, ReadEncoding
from neuroq.evaluation import accuracy, evaluate, mec_cost, switch_errors, top_samples
from neuroq.execution import IdealSimulator
from neuroq.experiments.notebook_example import haplotype_text, write_results
from neuroq.selection import WindowSelector
from neuroq.solvers import ClassicalSolver, ExactSolver, GreedyPhaser, QAOASolver
from neuroq.solvers.qaoa import make_qaoa
from neuroq.sources import SyntheticSource
from neuroq.sources.synthetic import ILLUMINA, NANOPORE, TechnologyProfile
from neuroq.transforms import fix_variable

SEED = 71
EXACT_SOLVER_MAX_QUBITS = 16
ILLUMINA_SNPS, ILLUMINA_READS = 10, 24
NANOPORE_SNPS, NANOPORE_READS = 18, 9
OUTPUT_DIR = Path("results")
ENCODINGS: dict[str, Encoding] = {"position": PositionEncoding(), "read": ReadEncoding()}
LOGICAL_BACKEND = "logical (necompilat)"

# Singurul punct măsurat pe Illumina, codificarea pe citiri (rulare anterioară, întreruptă).
ILLUMINA_READS_PREVIOUS_QAOA = {
    "method": "qaoa",
    "n_qubits": 24,
    "transforms": [],
    "qaoa_seconds": 1246,
    "status": "rulare anterioară întreruptă; un singur punct măsurat, fără cost sau acuratețe",
}


@dataclass
class Case:
    """O matrice de citiri, descrierea ei și ce se rulează pe fiecare codificare.

    `full` = codificările pe care rulează toți solverii; `resources_only` = codificările
    pentru care se scrie doar rândul de resurse, cu detaliile lui (motivul și altele).
    """

    matrix: ReadMatrix
    instance: dict
    full: tuple[str, ...] = ("position", "read")
    resources_only: dict[str, dict] = field(default_factory=dict)


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
    skipped = {
        "skipped_reason": ("codificarea pe citiri ar necesita 24 de qubiți; QAOA măsurat "
                           "într-o rulare anterioară a durat 1246 s"),
        "previous_measurement": ILLUMINA_READS_PREVIOUS_QAOA,
    }
    return Case(matrix, _instance("illumina", matrix, ILLUMINA),
                full=("position",), resources_only={"read": skipped})


def nanopore_case() -> Case:
    source = SyntheticSource.from_profile(NANOPORE, NANOPORE_SNPS, SEED, n_reads=NANOPORE_READS)
    matrix = _whole_matrix(source)
    skipped = {"skipped_reason": ("codificarea pe poziții ar necesita 18 qubiți; QAOA nu a "
                                  "fost rulat în acest experiment")}
    return Case(matrix, _instance("nanopore", matrix, NANOPORE),
                full=("read",), resources_only={"position": skipped})


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


def resources_row(case: Case, encoding: Encoding, extra_details: dict) -> ResultRow:
    """Rândul de resurse: variabile, J nenule și porți rzz ale unui strat QAOA (p=1)."""
    problem = encoding.encode(case.matrix)
    circuit, _, _ = make_qaoa(problem.J, problem.h, p=1)
    rzz = int(circuit.count_ops().get("rzz", 0))
    metrics = {
        "n_variables": problem.n_qubits,
        "n_couplings": int((problem.J != 0).sum()),
        "rzz_per_layer": rzz,
    }
    report = ResourceReport(n_qubits=problem.n_qubits, two_qubit_gates_before=rzz,
                            backend_name=LOGICAL_BACKEND)
    print(f"  resurse {case.instance['id']}/{problem.variable_kind}: {problem.n_qubits} "
          f"qubiți, {rzz} rzz", flush=True)
    return ResultRow(instance=case.instance, method="resources", encoding=problem.variable_kind,
                     transforms=[], metrics=metrics, resources=report,
                     details=dict(extra_details))


def run_comparison(cases: list[Case]) -> list[ResultRow]:
    rows: list[ResultRow] = []
    for case in cases:
        print(f"{case.instance['id']}: {case.instance['n_snps']} SNP-uri, "
              f"{case.instance['n_reads']} citiri", flush=True)
        rows.append(_greedy_row(case))
        for kind, encoding in ENCODINGS.items():
            if kind in case.full:
                rows.extend(_encoding_rows(case, encoding))
            elif kind in case.resources_only:
                rows.append(resources_row(case, encoding, case.resources_only[kind]))
    return rows


def write_outputs(rows: list[ResultRow], directory: Path) -> tuple[Path, Path]:
    """Scrie results.json și results.js (window.RESULTS = ...) atomic, prin fișiere .tmp."""
    directory.mkdir(parents=True, exist_ok=True)
    json_path, js_path = directory / "results.json", directory / "results.js"
    json_tmp, js_tmp = directory / "results.json.tmp", directory / "results.js.tmp"
    write_results(rows, json_tmp, note=(
        "Faza 2: notebook (ambele codificări), Illumina (10 SNP-uri, 24 citiri; solveri pe "
        "poziții, resurse pe citiri) și Nanopore (18 SNP-uri, 9 citiri; solveri pe citiri, "
        "resurse pe poziții), simulator ideal."))
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
