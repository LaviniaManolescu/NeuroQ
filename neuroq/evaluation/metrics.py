"""Metricile unei soluții: cost, acuratețe cu oglindă, switch errors, eșantioane."""

from __future__ import annotations

from typing import Optional

import numpy as np

from neuroq.core import Haplotypes, IsingProblem, Solution


def accuracy(copy_a: np.ndarray, truth_a: np.ndarray) -> float:
    """Fracția de SNP-uri corecte; oglinda (A↔B) contează ca soluție corectă."""
    same = np.asarray(copy_a) == np.asarray(truth_a)
    return float(max(same.mean(), 1.0 - same.mean()))


def switch_errors(copy_a: np.ndarray, truth_a: np.ndarray) -> int:
    """De câte ori se schimbă relația de fază (la fel / oglindă) între SNP-uri consecutive."""
    same = np.asarray(copy_a) == np.asarray(truth_a)
    return int(np.count_nonzero(same[1:] != same[:-1]))


def evaluate(
    problem: IsingProblem,
    solution: Solution,
    haplotypes: Haplotypes,
    truth: Optional[tuple[np.ndarray, np.ndarray]],
    optimum_energy: Optional[float],
) -> dict[str, float]:
    """Metricile numerice ale unui rând din results.json."""
    metrics: dict[str, float] = {"pair_cost": problem.energy(solution.bitstring)}
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
