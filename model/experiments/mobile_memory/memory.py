"""Reversible write, conditional cage transfer, and selective record retention."""

import numpy as np

from .current import one_charge_read


RECORDS = (0, 1, 2, -2, -1)


def coherent_check():
    """Compare the full pointer circuit with its action on every input column."""
    charge_basis = np.eye(3, dtype=complex)
    pointer_basis = np.eye(len(RECORDS), dtype=complex)
    cage_basis = np.eye(2, dtype=complex)
    swap = cage_basis[:, ::-1]
    charge_columns = [
        one_charge_read(vector, (0, 1), 1) for vector in charge_basis.T
    ]
    space = len(charge_basis) * len(pointer_basis) ** 2 * len(cage_basis)
    input_size = len(charge_basis) * len(cage_basis)
    untouched = np.eye(len(pointer_basis) * len(cage_basis))
    write = np.zeros((space, space), complex)
    transfer = np.zeros((space, space), complex)
    copy = np.zeros((space, space), complex)
    expected_kept = np.zeros((space, input_size), complex)
    expected_erased = np.zeros((space, input_size), complex)
    blank = pointer_basis[:, :1]
    for shift, record in enumerate(RECORDS):
        sector = np.column_stack([columns.get(record, np.zeros(3)) for columns in charge_columns])
        column = pointer_basis[:, shift:shift + 1]
        control = column @ column.T
        translate = np.roll(pointer_basis, shift, axis=0)
        pulse = swap if record > 0 else cage_basis
        write += np.kron(sector, np.kron(translate, untouched))
        transfer += np.kron(charge_basis, np.kron(control, np.kron(pointer_basis, pulse)))
        copy += np.kron(charge_basis, np.kron(control, np.kron(translate, cage_basis)))
        expected_kept += np.kron(sector, np.kron(blank, np.kron(column, pulse)))
        expected_erased += np.kron(sector, np.kron(blank, np.kron(blank, pulse)))

    input_frame = np.kron(charge_basis, np.kron(blank, np.kron(blank, cage_basis)))
    inverse = write.conj().T
    kept = inverse @ transfer @ copy @ write @ input_frame
    erased = inverse @ transfer @ write @ input_frame
    residuals = {
        "write_unitarity": float(np.max(abs(inverse @ write - np.eye(space)))),
        "unwrite_with_public_copy": float(np.max(abs(kept - expected_kept))),
        "unwrite_without_public_copy": float(np.max(abs(erased - expected_erased))),
    }
    if max(residuals.values()) > 1e-12:
        raise ValueError("the conditional transfer failed the write and unwrite identity")
    return residuals
