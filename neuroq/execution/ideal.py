"""Simulator ideal, fără zgomot, pe Aer."""

from __future__ import annotations

from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator

from neuroq.core import ResourceReport
from neuroq.execution.base import Executor, count_two_qubit_gates
from neuroq.execution.bitorder import qiskit_to_project


class IdealSimulator(Executor):
    """AerSimulator fără zgomot; circuitul rulează așa cum este, fără compilare pe cip."""

    backend_name = "aer_simulator"

    def run(
        self, circuit: QuantumCircuit, shots: int, seed: int
    ) -> tuple[dict[str, int], ResourceReport]:
        # tiparul de rulare din notebook-ul organizatorilor
        raw = (
            AerSimulator(max_parallel_threads=1)
            .run(circuit, shots=shots, seed_simulator=seed)
            .result()
            .get_counts()
        )
        counts: dict[str, int] = {}
        for key, value in raw.items():
            bits = qiskit_to_project(key)
            counts[bits] = counts.get(bits, 0) + value
        two_qubit = count_two_qubit_gates(circuit)
        report = ResourceReport(
            n_qubits=circuit.num_qubits,
            two_qubit_gates_before=two_qubit,
            two_qubit_gates_after=two_qubit,
            depth=circuit.depth(),
            backend_name=self.backend_name,
        )
        return counts, report
