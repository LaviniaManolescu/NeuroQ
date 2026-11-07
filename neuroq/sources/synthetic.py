"""Sursa sintetică: exemplul fix din notebook, un generator simplu și profiluri de tehnologie.

Valorile profilurilor vin din docs/arhitectura-detaliata.md („O precizare de onestitate”),
care rezumă valori publicate: Illumina are citiri de ~150 de perechi de baze și o eroare de
~0,1%; Nanopore are citiri de 1–100 kb și o eroare de ~5%; în genomul uman apare un SNP
heterozigot la 1.000–2.000 de baze. Profilul Illumina cu 2–3 SNP-uri pe citire este
deliberat optimist (o citire reală de 150 pb acoperă de obicei 0–1 SNP), ca în notebook.
"""

from __future__ import annotations

import math
import zlib
from dataclasses import dataclass
from typing import Optional

import numpy as np

from neuroq.core import MISSING, SourceData

# din notebook-ul organizatorilor: perechea ascunsă și cele 9 citiri (READS)
NOTEBOOK_PAIR = np.array([[0, 1, 1, 0, 1], [1, 0, 0, 1, 0]])
NOTEBOOK_READS = np.array([
    [0, 1, 1, -1, -1], [1, 0, 0, -1, -1],
    [-1, 1, 1, 0, -1], [-1, 0, 0, 1, -1],
    [-1, -1, 1, 0, 1], [-1, -1, 0, 1, 0],
    [0, -1, -1, 0, 1], [1, -1, -1, 1, 0],
    [0, 0, -1, -1, -1],
])

SNP_SPACING = 1000  # distanța, în litere, dintre două SNP-uri sintetice (modul simplu)
SNP_GAP_BP = (1000, 2000)  # distanța dintre două SNP-uri heterozigote, în profiluri
MAX_DRAWS_PER_READ = 50  # limită pentru citirile Nanopore care nu ating niciun SNP


@dataclass(frozen=True)
class TechnologyProfile:
    """Parametrii unei tehnologii de secvențiere, după valori publicate."""

    name: str
    error_rate: float
    coverage: float  # de câte ori e acoperit, în medie, fiecare SNP
    read_length_median_bp: float
    read_length_sigma: float = 0.0  # 0 = lungime fixă; altfel lognormală
    read_length_range_bp: tuple[float, float] = (1.0, math.inf)
    snps_per_read: Optional[tuple[int, int]] = None  # forțează 2–3 SNP-uri (Illumina optimist)

    def read_lengths(self, rng: np.random.Generator, size: int) -> np.ndarray:
        """Lungimile citirilor, în perechi de baze."""
        if self.read_length_sigma == 0:
            return np.full(size, self.read_length_median_bp)
        lengths = rng.lognormal(math.log(self.read_length_median_bp), self.read_length_sigma, size)
        return np.clip(lengths, *self.read_length_range_bp)


# Illumina: ~150 pb, ~0,1% eroare, acoperire ~30×; 2–3 SNP-uri pe citire (deliberat optimist).
ILLUMINA = TechnologyProfile(
    name="illumina", error_rate=0.001, coverage=30.0,
    read_length_median_bp=150.0, snps_per_read=(2, 3),
)
# Nanopore: 1–100 kb, mediană de câteva zeci de kb (lognormală), ~5% eroare, acoperire ~30×
# (puține citiri ca număr, pentru că fiecare e lungă).
NANOPORE = TechnologyProfile(
    name="nanopore", error_rate=0.05, coverage=30.0,
    read_length_median_bp=25_000.0, read_length_sigma=0.8,
    read_length_range_bp=(1_000.0, 100_000.0),
)
PROFILES = {p.name: p for p in (ILLUMINA, NANOPORE)}


class SyntheticSource:
    """Generează o pereche ascunsă de haplotipuri și citiri din ea, cu seed fix."""

    def __init__(
        self,
        n_snps: int = 5,
        n_reads: int = 9,
        snps_per_read: int = 3,
        error_rate: float = 0.0,
        seed: int = 0,
    ) -> None:
        if not 2 <= snps_per_read <= n_snps:
            raise ValueError("snps_per_read trebuie să fie între 2 și n_snps")
        if not 0.0 <= error_rate <= 1.0:
            raise ValueError("error_rate trebuie să fie între 0 și 1")
        self.n_snps = n_snps
        self.n_reads = n_reads
        self.snps_per_read = snps_per_read
        self.error_rate = error_rate
        self.seed = seed
        self._notebook = False
        self._profile: Optional[TechnologyProfile] = None
        self._target_reads: Optional[int] = None

    @classmethod
    def notebook_example(cls) -> SyntheticSource:
        """Exemplul fix din notebook: 5 SNP-uri, 9 citiri, A=01101, B=10010."""
        source = cls(n_snps=5, n_reads=9)
        source._notebook = True
        return source

    @classmethod
    def from_profile(
        cls, profile: TechnologyProfile, n_snps: int, seed: int, n_reads: Optional[int] = None
    ) -> SyntheticSource:
        """Citiri după un profil de tehnologie; adevărul depinde doar de (n_snps, seed).

        Fără `n_reads`, numărul de citiri rezultă din acoperirea profilului. Cu `n_reads`,
        se trag exact atâtea citiri din profil, fiecare cu cel puțin 2 SNP-uri.
        """
        if n_snps < 2:
            raise ValueError("n_snps trebuie să fie cel puțin 2")
        if n_reads is not None and n_reads < 1:
            raise ValueError("n_reads trebuie să fie cel puțin 1")
        source = cls(n_snps=n_snps, n_reads=2, snps_per_read=2, seed=seed)
        source._profile = profile
        source._target_reads = n_reads
        return source

    def load(self) -> SourceData:
        if self._profile is not None:
            return self._load_profile(self._profile)
        if self._notebook:
            values = NOTEBOOK_READS.copy()
            pair = NOTEBOOK_PAIR.copy()
            technology = "notebook"
            read_starts = np.array([self._first_covered(row) for row in values])
        else:
            values, pair, read_starts = self._generate()
            technology = "synthetic"
        n_reads, n_snps = values.shape
        return SourceData(
            snp_positions=np.arange(1, n_snps + 1) * SNP_SPACING,
            ref_alleles=["A"] * n_snps,
            alt_alleles=["G"] * n_snps,
            read_ids=[f"read{i + 1}" for i in range(n_reads)],
            read_starts=read_starts,
            values=values,
            technology=technology,
            true_haplotypes=(pair[0], pair[1]),
        )

    def _generate(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        rng = np.random.default_rng(self.seed)
        copy_a = rng.integers(0, 2, self.n_snps)
        pair = np.array([copy_a, 1 - copy_a])
        values = np.full((self.n_reads, self.n_snps), MISSING)
        starts = rng.integers(0, self.n_snps - self.snps_per_read + 1, self.n_reads)
        copies = rng.integers(0, 2, self.n_reads)
        for r, (start, copy) in enumerate(zip(starts, copies)):
            span = slice(start, start + self.snps_per_read)
            segment = pair[copy, span].copy()
            flips = rng.random(self.snps_per_read) < self.error_rate
            segment[flips] = 1 - segment[flips]
            values[r, span] = segment
        read_starts = (starts + 1) * SNP_SPACING
        return values, pair, read_starts

    def _truth_and_positions(self) -> tuple[np.ndarray, np.ndarray]:
        """Perechea adevărată și pozițiile SNP-urilor; depind doar de (n_snps, seed)."""
        rng = np.random.default_rng([self.seed, self.n_snps])
        copy_a = rng.integers(0, 2, self.n_snps)
        gaps = rng.integers(SNP_GAP_BP[0], SNP_GAP_BP[1] + 1, self.n_snps)
        return np.array([copy_a, 1 - copy_a]), np.cumsum(gaps)

    def _load_profile(self, profile: TechnologyProfile) -> SourceData:
        pair, positions = self._truth_and_positions()
        rng = np.random.default_rng([self.seed, self.n_snps, zlib.crc32(profile.name.encode())])
        if profile.snps_per_read is not None:
            spans = self._spans_by_snp_count(profile, positions, rng)
        else:
            spans = self._spans_by_length(profile, positions, rng)
        values = np.full((len(spans), self.n_snps), MISSING)
        starts = np.empty(len(spans), dtype=int)
        for r, (first, last, start_bp) in enumerate(spans):
            copy = rng.integers(0, 2)
            segment = pair[copy, first:last].copy()
            flips = rng.random(len(segment)) < profile.error_rate
            segment[flips] = 1 - segment[flips]
            values[r, first:last] = segment
            starts[r] = start_bp
        return SourceData(
            snp_positions=positions,
            ref_alleles=["A"] * self.n_snps,
            alt_alleles=["G"] * self.n_snps,
            read_ids=[f"{profile.name}{i + 1}" for i in range(len(spans))],
            read_starts=starts,
            values=values,
            technology=profile.name,
            true_haplotypes=(pair[0], pair[1]),
        )

    def _spans_by_snp_count(
        self, profile: TechnologyProfile, positions: np.ndarray, rng: np.random.Generator
    ) -> list[tuple[int, int, int]]:
        """Citiri care acoperă un număr fix de SNP-uri consecutive (Illumina optimist)."""
        low, high = profile.snps_per_read
        n_reads = self._target_reads or math.ceil(profile.coverage * self.n_snps / ((low + high) / 2))
        spans = []
        for _ in range(n_reads):
            k = int(rng.integers(low, high + 1))
            first = int(rng.integers(0, self.n_snps - k + 1))
            spans.append((first, first + k, int(positions[first])))
        return spans

    def _spans_by_length(
        self, profile: TechnologyProfile, positions: np.ndarray, rng: np.random.Generator
    ) -> list[tuple[int, int, int]]:
        """Citiri cu lungime în perechi de baze; acoperă SNP-urile din intervalul lor."""
        region = float(positions[-1] - positions[0] + 1)
        n_reads = self._target_reads or max(
            1, math.ceil(profile.coverage * region / profile.read_length_median_bp))
        min_snps = 1 if self._target_reads is None else 2  # cu n_reads fix, toate rămân în matrice
        spans: list[tuple[int, int, int]] = []
        draws = 0
        while len(spans) < n_reads and draws < MAX_DRAWS_PER_READ * n_reads:
            draws += 1
            length = float(profile.read_lengths(rng, 1)[0])
            start = rng.uniform(positions[0] - length, positions[-1])
            covered = np.flatnonzero((positions >= start) & (positions < start + length))
            if len(covered) >= min_snps:
                spans.append((int(covered[0]), int(covered[-1]) + 1, int(start)))
        if self._target_reads is not None and len(spans) < n_reads:
            raise ValueError(f"nu am putut trage {n_reads} citiri cu cel puțin 2 SNP-uri")
        return spans

    @staticmethod
    def _first_covered(row: np.ndarray) -> int:
        return int((np.flatnonzero(row != MISSING)[0] + 1) * SNP_SPACING)
