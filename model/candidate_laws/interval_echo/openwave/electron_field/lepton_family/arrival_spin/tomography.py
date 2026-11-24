"""Simulated local record tomography of the flux-heralded source-circulation channel."""

import argparse
from itertools import combinations
import json
import math
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[6]))
# The archive replay launches this script by filename from a temporary directory.
# ruff: noqa: E402
from model.candidate_laws.interval_echo.openwave.electron_field.lepton_family.arrival_spin import (
    circulation as circ,
)

SHOTS = 2**30
SEED = 191179
INPUTS = np.array([[1, 0], [0, 1], [1, 1], [1, -1], [1, 1j], [1, -1j]], dtype=complex)
INPUTS /= np.linalg.norm(INPUTS, axis=1)[:, None]


def sample(prob, rng):
    p = np.asarray(prob, float)
    if min(p) < -2e-12 or abs(sum(p) - 1) > 2e-12:
        raise ArithmeticError("The complete detector probabilities fail normalization")
    p = np.maximum(0.0, np.round(p, 12))
    return rng.multinomial(SHOTS, p / p.sum()).tolist()


def acquire(amplitude, pairs, rng):
    # One column per sharp exterior head word. Interference is local on the received rows.
    weights = np.sum(abs(amplitude) ** 2, axis=1)
    miss = 1 - sum(weights)
    diagonal = sample([*weights, miss], rng)
    records = []
    for i, j in pairs:
        a, b = amplitude[i], amplitude[j]
        other = sum(weights) - weights[i] - weights[j]
        settings = []
        for phase in (1.0, -1j):
            plus = np.vdot(a + phase * b, a + phase * b).real / 2
            minus = np.vdot(a - phase * b, a - phase * b).real / 2
            settings.append(sample([plus, minus, other, miss], rng))
        records.append(settings)
    return {"sharp_and_miss": diagonal, "pair_counts": records}


def invert(record, pairs):
    matrix = np.diag(np.array(record["sharp_and_miss"][:-1], dtype=complex) / SHOTS)
    for (i, j), (real, imag) in zip(pairs, record["pair_counts"], strict=True):
        # The analyzer with bra (1,-i) reads minus the imaginary part of rho_ij.
        matrix[i, j] = ((real[0] - real[1]) - 1j * (imag[0] - imag[1])) / (2 * SHOTS)
        matrix[j, i] = matrix[i, j].conjugate()
    return matrix


def study(board, route, sign, spur, rng):
    fixture = circ.tree.WrittenTree(board, route, sign)
    model = circ.Circulation(fixture, board, spur)
    source = np.zeros((len(model.states), 2))
    source[:2, 0], source[2:4, 1] = 2**-0.5, 2**-0.5
    eig, vec = np.linalg.eigh(model.hops_matrix + model.rings_matrix)
    at_receiver = (vec @ (np.exp(-8j * eig)[:, None] * (vec.T @ source)))[model.indices]
    mask = sum(
        bit
        for (a, b), bit in model.bits.items()
        if a in fixture.receiver and b in fixture.receiver
    )
    accepted = [model.states[i] for i in model.indices]
    labels = sorted({(p, w & mask) for p, w in accepted})
    exterior = sorted({w & ~mask for p, w in accepted})
    channel = np.zeros((len(labels), len(exterior), 2), complex)
    for (p, w), amplitudes in zip(accepted, at_receiver):
        channel[labels.index((p, w & mask)), exterior.index(w & ~mask)] = amplitudes
    pairs = list(combinations(range(len(labels)), 2))
    records = [acquire(channel @ inp, pairs, rng) for inp in INPUTS]
    zero, one, xplus, xminus, yplus, yminus = [invert(rec, pairs) for rec in records]
    off = (xplus - xminus + 1j * (yplus - yminus)) / 2
    # Transpose the input-reference factor of the Choi matrix. Keep preparation failures in the denominator.
    pt = np.block([[zero, off.conj().T], [off, one]]) / 2
    smallest = float(np.linalg.eigvalsh(pt)[0])
    m = len(labels)
    bins = 4 * len(INPUTS) * (m + 1 + 8 * len(pairs))
    eps = math.sqrt(math.log(2 * bins / 0.01) / (2 * SHOTS))
    # A categorical row has at most m+1 probabilities rounded to twelve decimal places.
    sampling_roundoff = (m + 3) * 1e-12
    bound = 2 * m * math.sqrt(2) * (eps + sampling_roundoff)
    return {
        "sign": sign,
        "spur": model.spur,
        "time": 8,
        "shots_per_setting": SHOTS,
        "local_records": labels,
        "input_order": ["0", "1", "X+", "X-", "Y+", "Y-"],
        "pair_order": pairs,
        "pair_outcomes": ["plus", "minus", "other heralded", "no signed flux"],
        "counts": records,
        "minimum_partial_transpose_eigenvalue": smallest,
        "simultaneous_eigenvalue_error": bound,
        "simultaneous_outcome_bins": bins,
        "certifies_entanglement": smallest + bound < 0,
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--output", type=Path, required=True)
    args = cli.parse_args()
    board, route = circ.tree.guide.longest_unpinned_route()
    rng = np.random.default_rng(SEED)
    data = [
        study(board, route, s, spur, rng) for spur in (False, True) for s in (-1, 1)
    ]
    if len({d["simultaneous_outcome_bins"] for d in data}) != 1:
        raise ArithmeticError("Four panels must have the same detector alphabet")
    with args.output.open("x") as out:
        json.dump(
            {"seed": SEED, "panels": data}, out, separators=(",", ":"), allow_nan=False
        )
        out.write("\n")
