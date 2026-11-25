"""Check the deposited words on all eight ternary ports, including leakage."""

from functools import lru_cache
from itertools import product

import numpy as np
import pytest
from scipy.linalg import expm

from .fixed_pointer import history


PORTS = (0, 1, 6, 10, 4, 5, 7, 11)
ENDPOINTS = tuple(product(range(3), repeat=2))
SWAP = np.eye(9)[:, [ENDPOINTS.index((b, a)) for a, b in ENDPOINTS]]
CONTACT = (np.eye(9) - SWAP) / 2 + np.diag([8 * (a - b) ** 2 for a, b in ENDPOINTS])
REFERENCE = np.diag([0, 8, 32])


@lru_cache(maxsize=None)
def gate(kind, duration):
    generator = {"contact": CONTACT, "reference": REFERENCE}[kind]
    return expm(-1j * np.pi * float(duration) * generator)


def apply_word(word, columns):
    for kind, a, b, duration in word:
        axes = [PORTS.index(a)]
        if kind == "contact":
            axes.append(PORTS.index(b))
        front = range(len(axes))
        operand = np.moveaxis(columns, axes, front)
        action = gate(kind, duration)
        updated = action @ operand.reshape(len(action), -1)
        columns = np.moveaxis(updated.reshape(operand.shape), front, axes)
    return columns


def port_loads(sign, source_address, records):
    loads = [1] * len(PORTS)
    loads[source_address] = 1 + sign
    loads[-1] = 1 - sign  # Balancing port stays coherent in all four input columns.
    for index, record in enumerate(records):
        start = 2 + 2 * index
        loads[start : start + 2] = (0, 2) if record == 0 else (2, 0)
    return tuple(loads)


def address_projectors(sign):
    flip = np.eye(2)[::-1]
    contact = 8.5 * np.eye(2) - flip / 2
    reference = np.diag([8 * (1 + sign) ** 2, 8])
    word = (
        expm(-1.5j * np.pi * contact)
        @ expm(-1j * np.pi * reference / 20)
        @ expm(-0.5j * np.pi * contact)
    )
    return [np.outer(row.conj(), row) for row in word]


@pytest.mark.parametrize("swapped", [False, True])
def test_full_ternary_history_matches_the_coherent_four_column_map(swapped):
    shape = (3,) * len(PORTS) + (4,)
    initial = np.zeros(shape, dtype=complex)
    ideal = np.zeros_like(initial)
    between_reads = np.eye(2)[::-1] if swapped else np.eye(2)
    for sign_index, sign in enumerate((-1, 1)):
        columns = slice(2 * sign_index, 2 * sign_index + 2)
        for address in (0, 1):
            initial[port_loads(sign, address, (1, 1)) + (2 * sign_index + address,)] = 1
        projectors = address_projectors(sign)
        for first, second in product((0, 1), repeat=2):
            branch = projectors[second] @ between_reads @ projectors[first]
            for address in (0, 1):
                key = port_loads(sign, address, (first, second)) + (columns,)
                ideal[key] = branch[address]

    actual = apply_word(history(swapped), initial)
    phase = np.exp(-1j * np.pi * 133 / 130)
    # Keep the same phase for both signs and all records. Off-code rows remain in the norm.
    error = (actual - phase * ideal).reshape(3 ** len(PORTS), 4)
    assert np.linalg.norm(error, ord=2) < 2e-12
