"""Contractele dintre module: fișe de date cu câmpuri fixe.

Convenții comune:
- indicii sunt de la 0 (SNP-ul 1 din documentație este indicele 0);
- în matrice, 0 = litera de referință, 1 = litera alternativă, -1 = lipsă;
- un șir de biți are ordinea proiectului: primul caracter din stânga este
  qubitul 0 (poziția 1). Conversia din ordinea simulatorului se face în `execution`.

Aici nu se importă biblioteci cuantice, doar numpy și biblioteca standard.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Callable, Literal, Optional

import numpy as np

MISSING = -1
VariableKind = Literal["position", "read"]


def _check_matrix_values(values: np.ndarray) -> None:
    """Verifică că matricea are doar valorile -1, 0 și 1."""
    if values.ndim != 2:
        raise ValueError(f"matricea trebuie să fie 2D, are {values.ndim} dimensiuni")
    if not np.isin(values, (MISSING, 0, 1)).all():
        raise ValueError("matricea poate conține doar -1 (lipsă), 0 și 1")


@dataclass
class Haplotypes:
    """Cele două copii reconstruite și copia la care e pusă fiecare citire."""

    copy_a: np.ndarray
    copy_b: np.ndarray
    read_assignment: np.ndarray  # per citire: 0 = copia A, 1 = copia B, -1 = suspectă


@dataclass
class SourceData:
    """Tot ce aduce o sursă dintr-o regiune mare, înainte de selecție."""

    snp_positions: np.ndarray  # coordonate genomice ale SNP-urilor heterozigote
    ref_alleles: list[str]
    alt_alleles: list[str]
    read_ids: list[str]
    read_starts: np.ndarray  # coordonata de început a fiecărei citiri
    values: np.ndarray  # citiri × SNP-uri, cu 0, 1 și -1 pentru lipsă
    technology: str  # de exemplu "illumina" sau "nanopore"
    raw_letters: Optional[np.ndarray] = None  # citiri × SNP-uri, litere A/C/G/T sau ""
    special_labels: dict[str, int] = field(default_factory=dict)  # ex. {"V1": 3}
    true_haplotypes: Optional[tuple[np.ndarray, np.ndarray]] = None

    def __post_init__(self) -> None:
        _check_matrix_values(self.values)
        n_reads, n_snps = self.values.shape
        if len(self.snp_positions) != n_snps:
            raise ValueError("snp_positions nu are câte o valoare pe coloană")
        if len(self.ref_alleles) != n_snps or len(self.alt_alleles) != n_snps:
            raise ValueError("alelele nu au câte o valoare pe coloană")
        if len(self.read_ids) != n_reads or len(self.read_starts) != n_reads:
            raise ValueError("read_ids și read_starts nu au câte o valoare pe rând")


@dataclass
class SelectionParams:
    """O tăietură: fereastra de SNP-uri și, opțional, downsampling-ul."""

    first_snp: int  # indicele primului SNP din fereastră
    n_snps: int
    target_reads: Optional[int] = None  # None = fără downsampling
    seed: int = 0


@dataclass
class ReadMatrix:
    """Matricea citiri × SNP-uri; intrarea comună pentru codificări și baseline."""

    values: np.ndarray  # citiri × SNP-uri, cu 0, 1 și -1 pentru lipsă
    read_ids: list[str]
    snp_ids: list[str]
    technology: str
    special_labels: dict[str, int] = field(default_factory=dict)
    true_haplotypes: Optional[tuple[np.ndarray, np.ndarray]] = None

    def __post_init__(self) -> None:
        _check_matrix_values(self.values)
        n_reads, n_snps = self.values.shape
        if len(self.read_ids) != n_reads:
            raise ValueError("read_ids nu are câte o valoare pe rând")
        if len(self.snp_ids) != n_snps:
            raise ValueError("snp_ids nu are câte o valoare pe coloană")


@dataclass
class IsingProblem:
    """Hamiltonianul offset + Σ h·z + Σ J·z·z, cu z = +1 pentru bitul 0 și -1 pentru bitul 1.

    J pozitiv respinge (variabilele vor valori diferite), J negativ atrage.
    Căutăm energia minimă.
    """

    J: np.ndarray  # n × n, doar triunghiul de sus strict (i < j); restul zero
    h: Optional[np.ndarray] = None  # câmpuri locale; None înseamnă zero
    offset: float = 0.0
    variable_kind: VariableKind = "position"
    fixed_variables: dict[int, int] = field(default_factory=dict)  # index original → bit
    decoder: Optional[Callable[[str], Haplotypes]] = None

    def __post_init__(self) -> None:
        self.J = np.asarray(self.J, dtype=float)
        if self.J.ndim != 2 or self.J.shape[0] != self.J.shape[1]:
            raise ValueError("J trebuie să fie o matrice pătrată")
        if np.any(np.tril(self.J) != 0):
            raise ValueError("J poate avea valori doar deasupra diagonalei (i < j)")
        n = self.J.shape[0]
        self.h = np.zeros(n) if self.h is None else np.asarray(self.h, dtype=float)
        if self.h.shape != (n,):
            raise ValueError(f"h trebuie să aibă {n} valori")

    @property
    def n_qubits(self) -> int:
        return self.J.shape[0]

    def energy(self, bitstring: str) -> float:
        """Energia unui șir în ordinea proiectului: offset + Σ h·z + Σ J·z·z."""
        if len(bitstring) != self.n_qubits:
            raise ValueError(
                f"șirul are {len(bitstring)} biți, problema are {self.n_qubits} qubiți"
            )
        if set(bitstring) - {"0", "1"}:
            raise ValueError("șirul poate conține doar caracterele '0' și '1'")
        bits = np.fromiter((int(c) for c in bitstring), dtype=float, count=len(bitstring))
        z = 1.0 - 2.0 * bits
        return float(self.offset + self.h @ z + z @ self.J @ z)


@dataclass
class ResourceReport:
    """Resursele cuantice ale unei rulări sau compilări."""

    n_qubits: int
    two_qubit_gates_before: int  # înainte de compilare
    two_qubit_gates_after: Optional[int] = None  # după compilare pe cip
    depth: Optional[int] = None  # adâncimea după compilare
    backend_name: str = ""
    job_id: Optional[str] = None  # doar pentru hardware


@dataclass
class Solution:
    """Răspunsul unui solver pe un IsingProblem."""

    bitstring: str  # ordinea proiectului
    energy: float
    samples: Optional[dict[str, int]] = None  # șir → apariții, doar la metodele cuantice
    resources: Optional[ResourceReport] = None
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass
class ResultRow:
    """Un rând din results.json. Schema e descrisă în docs/results_schema.md."""

    instance: dict[str, Any]
    method: str  # "exact", "qaoa", "rqaoa", "greedy"
    encoding: Optional[str]  # "position", "read" sau None pentru greedy
    transforms: list[str] = field(default_factory=list)
    metrics: dict[str, float] = field(default_factory=dict)  # doar numere
    resources: Optional[ResourceReport] = None
    details: dict[str, Any] = field(default_factory=dict)  # ce nu e număr

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> ResultRow:
        data = dict(data)
        resources = data.pop("resources", None)
        return cls(
            **data,
            resources=ResourceReport(**resources) if resources is not None else None,
        )
