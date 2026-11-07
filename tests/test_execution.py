"""Teste pentru modulul execution."""

from qiskit import QuantumCircuit

from neuroq.execution import IdealSimulator, qiskit_to_project


def test_key_conversion_matches_notebook_assert():
    # notebook: decode_key('00101', 5) == [1, 0, 1, 0, 0]
    assert qiskit_to_project("00101") == "10100"


def test_x_on_qubit_zero_reads_as_first_position():
    qc = QuantumCircuit(3)
    qc.x(0)
    qc.measure_all()
    counts, report = IdealSimulator().run(qc, shots=100, seed=1)
    assert counts == {"100": 100}
    assert report.n_qubits == 3 and report.two_qubit_gates_before == 0
