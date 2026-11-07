"""Clasa de bază pentru executori."""

from __future__ import annotations

from abc import ABC, abstractmethod

from qiskit import QuantumCircuit

from neuroq.core import ResourceReport


def count_two_qubit_gates(circuit: QuantumCircuit) -> int:
    # din notebook-ul organizatorilor (funcția resources)
    return sum(item.operation.num_qubits == 2 for item in circuit.data)


class Executor(ABC):
    """Rulează un circuit măsurat și întoarce numărătorile în ordinea proiectului."""

    @abstractmethod
    def run(
        self, circuit: QuantumCircuit, shots: int, seed: int
    ) -> tuple[dict[str, int], ResourceReport]:
        """Circuitul trebuie să aibă deja măsurătorile (measure_all)."""
