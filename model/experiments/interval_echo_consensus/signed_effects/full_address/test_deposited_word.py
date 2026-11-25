"""Replay the flat pulse deposit, without its router or its logical copy primitive."""

from fractions import Fraction
from functools import lru_cache
from itertools import groupby, product
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.executable_model import federation
from .. import collar_storage

ENDPOINTS = tuple(product(range(3), repeat=2))
SWAP = np.eye(9)[:, [ENDPOINTS.index((b, a)) for a, b in ENDPOINTS]]
CONTACT = (np.eye(9) - SWAP) / 2 + np.diag([8*(a-b)**2 for a, b in ENDPOINTS])
REFERENCE = np.diag([0., 8., 32.])
DEPOSIT = Path(__file__).with_name("results.json")


@lru_cache(maxsize=None)
def gate(kind, ratio):
    return expm(-1j * np.pi * float(Fraction(ratio)) *
                (CONTACT if kind == "contact" else REFERENCE))


@lru_cache(maxsize=None)
def propagate_local(word):
    """The observed output address and complex amplitude of every input column."""
    columns = np.eye(81, dtype=complex).reshape(3, 3, 3, 3, 81)
    for kind, a, b, ratio in word:
        axes = (a, b) if kind == "contact" else (a,)
        active = np.moveaxis(columns, axes, range(len(axes)))
        action = gate(kind, ratio)
        evolved = action @ active.reshape(len(action), -1)
        columns = np.moveaxis(evolved.reshape(active.shape), range(len(axes)), axes)
    columns = columns.reshape(81, 81)
    outputs = np.argmax(abs(columns), axis=0)
    return outputs, columns[outputs, range(81)]


def test_complete_flat_schedule_and_four_collar_deposits():
    archive = json.loads(DEPOSIT.read_text())
    board = federation(0)
    edges = set(board.seams)
    stations = [s["pointer"] for s in collar_storage.geometry()]
    cases = list(product((-1, 1), range(12), range(16)))
    source, reference = list(range(12)), 239
    loads = np.ones((len(cases), board.n), dtype=np.int8)
    for i, (sign, address, prior) in enumerate(cases):
        loads[i, source[address]] += sign
        loads[i, reference] -= sign
        for bit, (a, b) in enumerate(stations):
            loads[i, a], loads[i, b] = (0, 2) if prior & (1 << bit) else (2, 0)
    expected = loads.copy()
    for i, (_, address, prior) in enumerate(cases):
        for bit, (a, b) in enumerate(stations):
            expected[i, a], expected[i, b] = (0, 2) if (prior ^ address) & (1 << bit) else (2, 0)
    amplitudes = np.ones(len(cases), dtype=complex)
    def is_swap(step):
        return step[0] == "contact" and Fraction(step[-1]) == 1

    np.testing.assert_allclose(gate("contact", "1"), SWAP, atol=5e-14, rtol=0)
    calls = 0
    for swap_run, group in groupby(archive["native_pulse_word"], key=is_swap):
        steps = list(group)
        for kind, a, b, _ in steps:
            if kind == "contact":
                assert tuple(sorted((a, b))) in edges
            else:
                assert kind == "reference" and a == b and board.degree[a] == 5
        if swap_run:
            for _, a, b, _ in steps:
                loads[:, [a, b]] = loads[:, [b, a]]
        else:
            ports = sorted({p for _, a, b, _ in steps for p in (a, b)})
            assert len(ports) == 4
            action, phase = propagate_local(tuple((k, ports.index(a), ports.index(b), t)
                                                  for k, a, b, t in steps))
            columns = np.ravel_multi_index(loads[:, ports].T, (3,)*4)
            np.testing.assert_allclose(abs(phase[columns]), 1, atol=5e-13, rtol=0)
            loads[:, ports] = np.asarray(np.unravel_index(action[columns], (3,)*4)).T
            amplitudes *= phase[columns]
            calls += 1
    np.testing.assert_array_equal(loads, expected)
    expected_phase = np.exp(-1j*np.pi*133/260)**20
    np.testing.assert_allclose(amplitudes, expected_phase, atol=1e-12, rtol=0)
    assert calls == 20

    # Four different seams with the same certified local transfer. A common
    # deposit phase preserves all input amplitudes, including relative charge.
    _, basis, echo = collar_storage.gates(1.)
    pair = np.eye(9)[:, [ENDPOINTS.index((1, 1))]]
    images = (echo @ basis)[:, [ENDPOINTS.index(q) for q in ((2, 0), (0, 2))]]
    targets = np.kron(collar_storage.MEMORIES, pair)
    np.testing.assert_allclose(images, np.exp(-1j*np.pi/4)*targets, atol=8e-14, rtol=0)
