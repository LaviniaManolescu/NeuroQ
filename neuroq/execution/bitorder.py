"""Conversia ordinii biților: singurul loc din proiect unde se face.

Qiskit afișează qubitul 0 în dreapta; proiectul citește poziția 1 (qubitul 0) în stânga.
"""

from __future__ import annotations


def qiskit_to_project(key: str) -> str:
    """Transformă o cheie Qiskit (qubitul 0 în dreapta) în ordinea proiectului."""
    bits = key.replace(" ", "")
    if set(bits) - {"0", "1"}:
        raise ValueError(f"cheie de măsurare invalidă: {key!r}")
    return bits[::-1]
