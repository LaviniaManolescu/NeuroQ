"""Adâncimea unui circuit după transpilare într-un set de porți uzual pe cipuri IBM."""

from __future__ import annotations

from qiskit import QuantumCircuit
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

BASIS_GATES = ["rz", "sx", "x", "cx"]


def transpiled_depth(circuit: QuantumCircuit, seed: int = 71) -> int:
    """Adâncimea după transpilare în rz, sx, x, cx (nivel 3, fără hartă de qubiți, seed fix).

    Fără hardware: nu există hartă de cuplare, deci nu apar porți swap; numărul e un minim
    pentru un cip cu legături între oricare doi qubiți.
    """
    manager = generate_preset_pass_manager(
        optimization_level=3, basis_gates=BASIS_GATES, seed_transpiler=seed)
    return int(manager.run(circuit).depth())
