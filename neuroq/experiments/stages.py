"""Etapele demonstrației, doar calcul: fără afișare, fără fișiere, fără matplotlib.

`run_stages` rulează pipeline-ul pe o instanță și întoarce tot ce trebuie desenat sau tipărit.
`demo.py` și notebooks/demo.ipynb o apelează și doar afișează rezultatul.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Any, Optional

import numpy as np

from neuroq.core import Haplotypes, IsingProblem, ReadMatrix, Solution
from neuroq.encodings import PositionEncoding, ReadEncoding
from neuroq.encodings.position import build_evidence
from neuroq.encodings.read import build_read_evidence
from neuroq.evaluation import evaluate
from neuroq.execution import IdealSimulator
from neuroq.experiments.resource_comparison import (
    SEED, Case, illumina_case, nanopore_case, notebook_case,
)
from neuroq.solvers import ExactSolver
from neuroq.solvers.qaoa import CostCallback, QAOASolver, make_qaoa, trained_circuit
from neuroq.transforms import fix_variable

INSTANCES = {"notebook": notebook_case, "illumina": illumina_case, "nanopore": nanopore_case}
ENCODING_NAMES = {"pozitii": "position", "citiri": "read", "position": "position", "read": "read"}
ROMANIAN = {"position": "pozitii", "read": "citiri"}
# Pe instanțele mari rulează doar codificarea mică; cealaltă primește doar numerele de resurse.
ALLOWED = {"notebook": ("position", "read"), "illumina": ("position",), "nanopore": ("read",)}


@dataclass
class StagesResult:
    """Tot ce produce demonstrația, etapă cu etapă."""

    instance: str
    encoding: str  # "position" sau "read"
    fixed: bool
    matrix: ReadMatrix
    labels: list[str]  # etichetele variabilelor (SNP-uri sau citiri)
    a: np.ndarray
    d: np.ndarray
    problem: IsingProblem  # problema completă (J, offset), înainte de fixare
    solved: IsingProblem  # problema dată lui QAOA (cu un qubit fixat, dacă `fixed`)
    exact: Solution
    qaoa: Solution
    circuit: Any  # circuitul QAOA antrenat, gata de desenat
    full_bitstring: str  # șirul QAOA cu bitul fixat reinserat, pe toate variabilele
    haplotypes: Haplotypes
    metrics: dict[str, float]
    other_encoding: Optional[dict[str, Any]]  # resursele codificării care nu se rulează
    qaoa_seconds: float = 0.0  # durata apelului QAOASolver.solve (antrenare + eșantionare)

    @property
    def cost_history(self) -> list[list[float]]:
        return self.qaoa.metadata["cost_history"]

    @property
    def best_start(self) -> int:
        return self.qaoa.metadata["best_start"]

    @property
    def top_samples(self) -> list[dict]:
        return self.qaoa.metadata["top_samples"]


def resolve_encoding(instance: str, encoding: Optional[str]) -> str:
    """Codificarea cerută (pozitii/citiri) sau cea implicită; refuză codificarea mare."""
    if instance not in INSTANCES:
        raise ValueError(f"instanță necunoscută: {instance!r}; alege din {sorted(INSTANCES)}")
    kind = ENCODING_NAMES.get(encoding or ROMANIAN[ALLOWED[instance][0]])
    if kind is None:
        raise ValueError(f"codificare necunoscută: {encoding!r}; alege pozitii sau citiri")
    if kind not in ALLOWED[instance]:
        raise ValueError(
            f"pe {instance} rulează doar codificarea pe {ROMANIAN[ALLOWED[instance][0]]}; "
            f"pentru {ROMANIAN[kind]} se afișează doar resursele")
    return kind


def encoding_resources(matrix: ReadMatrix, kind: str) -> dict[str, Any]:
    """Qubiții, J nenule și porțile rzz ale unui strat QAOA pentru o codificare."""
    problem = (PositionEncoding() if kind == "position" else ReadEncoding()).encode(matrix)
    circuit, _, _ = make_qaoa(problem.J, problem.h, p=1)
    return {
        "encoding": kind,
        "n_qubits": problem.n_qubits,
        "n_couplings": int((problem.J != 0).sum()),
        "rzz_per_layer": int(circuit.count_ops().get("rzz", 0)),
    }


def _labels(matrix: ReadMatrix, kind: str) -> list[str]:
    if kind == "position":
        return [str(i + 1) for i in range(len(matrix.snp_ids))]
    return [f"R{i + 1}" for i in range(len(matrix.read_ids))]


def _with_fixed_bits(bitstring: str, fixed: dict[int, int]) -> str:
    bits = list(bitstring)
    for index in sorted(fixed):
        bits.insert(index, str(fixed[index]))
    return "".join(bits)


def run_stages(
    case: Case, encoding: str, fix: bool = False, callback: Optional[CostCallback] = None,
) -> StagesResult:
    """Rulează toate etapele pe o instanță și o codificare; nu afișează nimic."""
    instance = case.instance["id"]
    kind = resolve_encoding(instance, encoding)
    matrix = case.matrix
    if kind == "position":
        a, d = build_evidence(matrix.values)
        problem = PositionEncoding().encode(matrix)
    else:
        a, d = build_read_evidence(matrix)
        problem = ReadEncoding().encode(matrix)
    solved = fix_variable(problem) if fix else problem
    exact = ExactSolver().solve(problem)
    start = time.perf_counter()
    qaoa = QAOASolver(IdealSimulator(), seed=SEED, callback=callback).solve(solved)
    qaoa_seconds = time.perf_counter() - start
    haps = solved.decoder(qaoa.bitstring)
    metrics = evaluate(solved, qaoa, haps, matrix.true_haplotypes, exact.energy, matrix.values)
    others = [k for k in ("position", "read") if k != kind and k not in ALLOWED[instance]]
    return StagesResult(
        instance=instance, encoding=kind, fixed=fix, matrix=matrix, labels=_labels(matrix, kind),
        a=a, d=d, problem=problem, solved=solved, exact=exact, qaoa=qaoa,
        circuit=trained_circuit(solved, qaoa.metadata["angles"], qaoa.metadata["p"]),
        full_bitstring=_with_fixed_bits(qaoa.bitstring, solved.fixed_variables),
        haplotypes=haps, metrics=metrics,
        other_encoding=encoding_resources(matrix, others[0]) if others else None,
        qaoa_seconds=qaoa_seconds,
    )


def load_case(instance: str) -> Case:
    if instance not in INSTANCES:
        raise ValueError(f"instanță necunoscută: {instance!r}; alege din {sorted(INSTANCES)}")
    return INSTANCES[instance]()
