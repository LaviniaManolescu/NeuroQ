"""Runner-ul fazei 1: exemplul din notebook, codificare pe poziții, Exact și QAOA."""

from __future__ import annotations

import json
from pathlib import Path

from neuroq.core import IsingProblem, ReadMatrix, ResultRow, SelectionParams, Solution
from neuroq.encodings import PositionEncoding
from neuroq.evaluation import evaluate, top_samples
from neuroq.execution import IdealSimulator
from neuroq.selection import WindowSelector
from neuroq.solvers import ExactSolver, QAOASolver
from neuroq.sources import SyntheticSource

SCHEMA_VERSION = "0.1"
SEED = 71
DEFAULT_OUTPUT = Path("results") / "results.json"


def haplotype_text(copy) -> str:
    """Copia ca text: 0 și 1, iar un SNP nefazat (-1) devine '-'."""
    return "".join("-" if v == -1 else str(int(v)) for v in copy)


def _row(
    instance: dict, method: str, problem: IsingProblem, matrix: ReadMatrix,
    solution: Solution, optimum: float,
) -> ResultRow:
    haps = problem.decoder(solution.bitstring)
    details: dict = {
        "haplotypes": {
            "copy_a": haplotype_text(haps.copy_a),
            "copy_b": haplotype_text(haps.copy_b),
            "read_assignment": [int(v) for v in haps.read_assignment],
        }
    }
    if solution.samples:
        details["top_samples"] = top_samples(solution.samples)
    return ResultRow(
        instance=instance,
        method=method,
        encoding=problem.variable_kind,
        transforms=[],
        metrics=evaluate(problem, solution, haps, matrix.true_haplotypes, optimum, matrix.values),
        resources=solution.resources,
        details=details,
    )


def run_notebook_example() -> list[ResultRow]:
    params = SelectionParams(first_snp=0, n_snps=5, seed=SEED)
    data = SyntheticSource.notebook_example().load()
    matrix = WindowSelector().select(data, params)
    problem = PositionEncoding().encode(matrix)
    instance = {
        "id": "notebook-5snp-9reads",
        "source": "synthetic",
        "technology": matrix.technology,
        "first_snp": params.first_snp,
        "n_snps": params.n_snps,
        "n_reads": len(matrix.read_ids),
        "error_rate": None,
        "seed": SEED,
    }
    exact = ExactSolver().solve(problem)
    qaoa = QAOASolver(IdealSimulator(), seed=SEED).solve(problem)
    return [
        _row(instance, "exact", problem, matrix, exact, exact.energy),
        _row(instance, "qaoa", problem, matrix, qaoa, exact.energy),
    ]


def write_results(
    rows: list[ResultRow], path: Path, is_sample: bool = False,
    note: str = "Faza 1: exemplul din notebook, codificare pe poziții, simulator ideal.",
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema_version": SCHEMA_VERSION,
        "is_sample": is_sample,
        "note": note,
        "rows": [row.to_dict() for row in rows],
    }
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    write_results(run_notebook_example(), DEFAULT_OUTPUT)
    print(f"scris: {DEFAULT_OUTPUT}")
