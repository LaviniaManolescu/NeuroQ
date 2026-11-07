"""Fixarea unui qubit: elimină oglinda și scade problema cu un qubit."""

from __future__ import annotations

import numpy as np

from neuroq.core import Haplotypes, IsingProblem


def _spin(value: int) -> float:
    if value not in (0, 1):
        raise ValueError("valoarea fixată trebuie să fie 0 sau 1")
    return 1.0 if value == 0 else -1.0


def _original_index(problem: IsingProblem, index: int) -> int:
    """Indicele din problema inițială, ținând cont de variabilele fixate înainte."""
    n_original = problem.n_qubits + len(problem.fixed_variables)
    remaining = [i for i in range(n_original) if i not in problem.fixed_variables]
    return remaining[index]


def fix_variable(problem: IsingProblem, index: int = 0, value: int = 0) -> IsingProblem:
    """Fixează qubitul `index` la `value`; cuplajele lui devin câmpuri locale h.

    Energia problemei noi pe un șir s este egală cu energia problemei vechi pe șirul
    cu bitul fixat reinserat la poziția `index`.
    """
    n = problem.n_qubits
    if not 0 <= index < n:
        raise ValueError(f"indicele {index} nu există într-o problemă cu {n} qubiți")
    s = _spin(value)

    full_J = problem.J + problem.J.T  # simetrică, ca să citim cuplajele cu q din ambele părți
    keep = [i for i in range(n) if i != index]
    J = np.triu(problem.J[np.ix_(keep, keep)], k=1)
    h = problem.h[keep] + s * full_J[index, keep]
    offset = problem.offset + s * problem.h[index]

    fixed = dict(problem.fixed_variables)
    fixed[_original_index(problem, index)] = value

    old_decoder = problem.decoder
    decoder = None
    if old_decoder is not None:
        def decoder(bitstring: str) -> Haplotypes:
            return old_decoder(bitstring[:index] + str(value) + bitstring[index:])

    return IsingProblem(
        J=J,
        h=h,
        offset=offset,
        variable_kind=problem.variable_kind,
        fixed_variables=fixed,
        decoder=decoder,
    )
