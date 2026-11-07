"""Metricile unei soluții: cost, acuratețe cu oglindă, switch errors, eșantioane."""

from __future__ import annotations

from typing import Optional

import numpy as np

from neuroq.core import Haplotypes, IsingProblem, Solution


def accuracy(copy_a: np.ndarray, truth_a: np.ndarray) -> float:
    """Fracția de SNP-uri corecte; oglinda (A↔B) contează ca soluție corectă.

    Un SNP nefazat (-1) contează ca greșeală în ambele orientări.
    """
    copy_a, truth_a = np.asarray(copy_a), np.asarray(truth_a)
    same = copy_a == truth_a
    mirror = copy_a == 1 - truth_a
    return float(max(same.mean(), mirror.mean()))


def switch_errors(copy_a: np.ndarray, truth_a: np.ndarray) -> int:
    """De câte ori se schimbă relația de fază între SNP-uri fazate consecutive (-1 se sare)."""
    copy_a, truth_a = np.asarray(copy_a), np.asarray(truth_a)
    phased = copy_a != -1
    same = (copy_a == truth_a)[phased]
    return int(np.count_nonzero(same[1:] != same[:-1]))


# după funcția mec_cost a colegei (algoritmii clasici de fazare), extinsă pentru SNP-uri nefazate
def mec_cost(copy_a: np.ndarray, reads: np.ndarray) -> int:
    """Minimum Error Correction: câte valori trebuie corectate ca fiecare citire să se
    potrivească cu copia A sau cu complementul ei. O valoare observată pe un SNP
    nefazat (-1) se numără mereu ca eroare.
    """
    copy_a = np.asarray(copy_a)
    phased = copy_a != -1
    total = 0
    for row in np.asarray(reads):
        seen = row >= 0
        known = seen & phased
        errors = int(np.count_nonzero(row[known] != copy_a[known]))
        total += min(errors, int(known.sum()) - errors) + int(np.count_nonzero(seen & ~phased))
    return total


def evaluate(
    problem: IsingProblem,
    solution: Solution,
    haplotypes: Haplotypes,
    truth: Optional[tuple[np.ndarray, np.ndarray]],
    optimum_energy: Optional[float],
    reads: Optional[np.ndarray] = None,
) -> dict[str, float]:
    """Metricile numerice ale unui rând din results.json; `reads` adaugă mec_cost."""
    metrics: dict[str, float] = {"pair_cost": problem.energy(solution.bitstring)}
    if reads is not None:
        metrics["mec_cost"] = mec_cost(haplotypes.copy_a, reads)
    if truth is not None:
        metrics["accuracy"] = accuracy(haplotypes.copy_a, truth[0])
        metrics["switch_errors"] = switch_errors(haplotypes.copy_a, truth[0])
    if solution.samples:
        shots = sum(solution.samples.values())
        energies = {bits: problem.energy(bits) for bits in solution.samples}
        metrics["mean_sample_cost"] = float(
            sum(energies[b] * c for b, c in solution.samples.items()) / shots
        )
        if optimum_energy is not None:
            hits = sum(c for b, c in solution.samples.items()
                       if np.isclose(energies[b], optimum_energy))
            metrics["optimum_frequency"] = hits / shots
    return metrics


def top_samples(samples: dict[str, int], k: int = 5) -> list[dict[str, object]]:
    """Cele mai frecvente k șiruri, descrescător după apariții."""
    ranked = sorted(samples.items(), key=lambda item: (-item[1], item[0]))[:k]
    return [{"bitstring": bits, "count": count} for bits, count in ranked]
