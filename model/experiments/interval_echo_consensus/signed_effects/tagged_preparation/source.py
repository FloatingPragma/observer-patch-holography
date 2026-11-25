"""Prepare the calibrated common density with an untouched material tag."""

from fractions import Fraction
from functools import lru_cache
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from ..full_address import measure as address
from ..sequential import require

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]
ORIGINAL = ROOT / "code/born_context_phase_lift/BORN_CONTEXT_WEB_PAYLOAD.v1.json"
TABLE = ROOT / "code/phase_operation_producer/PHASE_OPERATION_RECEIPT.v1.json"
ROUTE = (1, 82, 72, 73)
SITE = {"source": 73, "park": 76, "cross": 77, "host": 79}
PORTS = (0, 1, 76, 77, 79)


@lru_cache(maxsize=None)
def native_gate(kind, ratio):
    law = candidate.native_law
    if kind == "contact":
        h0, h1 = law.local_term_parts()
        full = expm(-1j * np.pi * ratio * (h0 + candidate.KAPPA * h1))[:, :9]
        require(np.linalg.norm(full[9:]) < 5e-12, "blank-memory leakage")
        return full[:9]
    require(kind == "reference", "unknown native operation")
    return np.diag(np.exp(-1j * np.pi * ratio * np.array([0, 8, 32])))


def prepare(tagged=True):
    census = json.loads(ORIGINAL.read_text())["context_web"]["diagonal_context"]
    count = census["realized_outcome_counts"]
    p = count[0] / sum(count)
    ratio = 2 * np.arccos(np.sqrt(p)) / np.pi
    states = np.zeros((3, 3, 3, 3, 3, 2), complex)
    for j, sign in enumerate((-1, 1)):
        states[1 + sign, 1, 2, 0, 1, j] = 1
    board = federation(0)
    require(1 in board.neighbours[0], "source rails lack their calibration contact")
    perm = np.arange(board.n)
    paths = list(zip(ROUTE, ROUTE[1:]))
    swap = np.eye(9)[:, np.arange(9).reshape(3, 3).T.ravel()]
    require(
        np.linalg.norm(native_gate("contact", 1) - swap, 2) < 3e-13,
        "route swap incorrect",
    )
    for a, b in paths:
        require(b in board.neighbours[a], "route crosses a missing seam")
        perm[[a, b]] = perm[[b, a]]
    require(
        perm[73] == 1 and all(perm[q] == q for q in PORTS if q != 1),
        "tag route disturbed a reserved load",
    )
    pulse_word = [("contact", 0, 1, ratio)]
    if tagged:
        for kind, a, b, t in address.local_word(SITE):
            if kind == "reference":
                require(board.degree[a] == 5, "unavailable degree reference")
            else:
                require(b in board.neighbours[a], "unavailable tag contact")
            pulse_word.append((kind, int(perm[a]), int(perm[b]), float(t)))
    for a, b in reversed(paths):
        perm[[a, b]] = perm[[b, a]]
    require(np.array_equal(perm, np.arange(board.n)), "tag route not restored")
    for kind, a, b, t in pulse_word:
        axes = (
            [PORTS.index(a)]
            if kind == "reference"
            else [PORTS.index(a), PORTS.index(b)]
        )
        front = range(len(axes))
        local = np.moveaxis(states, axes, front)
        gate = native_gate(kind, t)
        evolved = gate @ local.reshape(len(gate), -1)
        states = np.moveaxis(evolved.reshape(local.shape), front, axes)

    # Basis axes are charge sign, source address, tag bit. A balancing charge
    # at port 239 identifies the two input columns in the same total-load sector.
    result = np.zeros((2, 2, 2), complex)
    for j, sign in enumerate((-1, 1)):
        for rail in range(2):
            data = [1, 1]
            data[rail] += sign
            for tag, pair in enumerate(((2, 0), (0, 2))):
                result[j, rail, tag] = states[(*data, *pair, 1, j)]
                states[(*data, *pair, 1, j)] = 0
    require(np.linalg.norm(states) < 3e-12, "source, pointer or host leakage")
    require(np.linalg.norm(result[0] - result[1]) < 1e-12, "charge reference dephased")
    return result, float(ratio * np.pi)


def archived_effect(entry):
    def scalar(parts):
        real, radical, imag, imag_radical = [float(Fraction(n)) for n in parts]
        return real + np.sqrt(3) * radical + 1j * (imag + np.sqrt(3) * imag_radical)

    return np.array([[scalar(v) for v in row] for row in entry["effect"]])
