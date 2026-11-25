"""Read the seam-direction energy symbol and test its interpretation as a frequency."""

import argparse
from functools import cache
import json
from math import pi
from pathlib import Path

import numpy as np
from scipy.linalg import eigh, expm

from model.candidate_laws.interval_echo.candidate import KAPPA, native_law
from model.executable_model import carrier
from model.experiments.interval_echo_consensus.refinement_cooling.dimer_band.measure import (
    gates,
    pointer,
)
from model.experiments.interval_echo_gauge.flag_measure.measure import path
from ..experiment import geometry

HERE = Path(__file__).resolve().parent
SHOTS = 2**34
DELTA = 1 / 32


def unit(v):
    return v / np.linalg.norm(v)


@cache
def chart():
    v = geometry()
    i, j = carrier.CELL_SEAMS[0]
    face = carrier.cell().faces()[0]
    directions = dict(
        port=v[0],
        edge=unit(v[i] + v[j]),
        face=unit(v[list(face)].sum(axis=0)),
        oblique=unit(np.array([1, 2, 3])),
    )
    return v / np.linalg.norm(v[i] - v[j]), directions


@cache
def preparation(direction, gradient, transverse):
    pos, axes = chart()
    axis = axes[direction]
    envelope = np.ones(12)
    if transverse:
        offset = next(x for x in pos if np.linalg.norm(np.cross(x, axis)) > 0.1)
        envelope = pos @ unit(np.cross(axis, offset))
    target = unit(envelope) * np.exp(1j * gradient * (pos @ axis))
    return prepare_target(target)


def prepare_target(target):
    sequence = list(path())
    psi = np.eye(12, dtype=complex)[sequence[0]]
    h = gates("native")[0][:2, :2]
    for offset, (i, j) in enumerate(zip(sequence, sequence[1:])):
        # Split the remaining population; subsequent private phase pulses
        # set the sign and phase after all material amplitudes are in place.
        remaining = np.linalg.norm(target[sequence[offset:]])
        duration = (
            2 * np.arccos(np.clip(abs(target[i]) / remaining, 0, 1))
            if remaining > 1e-14
            else 0
        )
        psi[[i, j]] = expm(-1j * duration * h) @ psi[[i, j]]
    duration = (np.angle(psi) - np.angle(target)) % (2 * pi) / (4 * KAPPA)
    psi *= np.exp(-1j * 4 * KAPPA * duration)
    np.testing.assert_allclose(psi, target, rtol=0, atol=3e-13)
    return psi


def local_block(scale):
    exchange, contrast = native_law.local_term_parts()
    full = scale * exchange + KAPPA * contrast
    ix = [native_law.SRC.index(x) for x in [(2, 0), (0, 2)]]
    np.testing.assert_allclose(np.delete(full[:, ix], ix, axis=0), 0, atol=1e-14)
    return full[np.ix_(ix, ix)]


@cache
def spectrum(which, scale):
    board = carrier.cell() if which == "cell" else carrier.federation(0)
    full = np.diag(4 * KAPPA * (6 - board.degree)).astype(float)
    block = local_block(scale)
    for i, j in board.seams:
        full[np.ix_([i, j], [i, j])] += block
    standard = np.ones(board.n) / np.sqrt(board.n)
    origin = float(standard @ full @ standard)
    np.testing.assert_allclose(full @ standard, origin * standard, atol=2e-12)
    return eigh(full - origin * np.eye(board.n))


def draw(rng, probs):
    # Use dyadic bins with a conservative distance-to-bin-boundary check.
    # Set the last weight by subtraction to keep all probabilities exact.
    raw = np.asarray(probs)
    np.testing.assert_allclose(raw.sum(), 1, atol=3e-12, rtol=0)
    grid = 2**25
    scaled = np.clip(raw[:-1], 0, 1) * grid
    if np.any(np.abs(scaled - np.floor(scaled) - 0.5) < 1e-7):
        raise ValueError("Uncertified probability bin boundary")
    counts = np.rint(scaled).astype(np.int64)
    p = np.r_[counts, grid - counts.sum()] / grid
    if min(p) < 0:
        raise ValueError("Negative residual bin")
    return rng.multinomial(SHOTS, p).tolist()


def quadratures(rng, energies, weights, holds, both):
    mean = np.exp(-1j * np.outer(holds, energies)) @ weights
    result = []
    for response in mean:
        # The analyser sign is reversed for sine, so a positive energy has
        # positive initial sine slope in the referenced clock.
        row = {"sin": draw(rng, [(1 - response.imag) / 2, (1 + response.imag) / 2])}
        if both:
            row["cos"] = draw(rng, [(1 + response.real) / 2, (1 - response.real) / 2])
        result.append(row)
    return result


def clock(rng, which, owner, psi, scale, holds, both=False):
    energies, vectors = spectrum(which, scale)
    amplitudes = vectors[12 * owner : 12 * owner + 12].T @ psi
    return quadratures(rng, energies, abs(amplitudes) ** 2, holds, both)


def assay():
    rng = np.random.default_rng(29151125)
    pointer()
    h, decoder, _ = gates("native")
    orientation = decoder[:2, :2]
    np.testing.assert_allclose(
        abs(orientation @ np.array([[1, 1], [1, -1]]) / np.sqrt(2)) ** 2,
        np.eye(2),
        atol=2e-13,
    )
    calibration = []
    for t in (pi / 12, pi / 3):
        calibration.append(draw(rng, abs(expm(-1j * t * h)[:2, 0]) ** 2))
    local, clocks, frequency = [], [], []
    boards = [("cell", 0), ("federation", 0), ("federation", 7)]
    for direction in chart()[1]:
        for gradient in (0, 1, 2, 3, 4):
            psi = preparation(direction, gradient, False)
            edges = []
            for i, j in carrier.CELL_SEAMS:
                signs = abs(orientation @ psi[[i, j]]) ** 2
                edges.append(draw(rng, np.r_[signs, 1 - signs.sum()]))
            local.append(dict(axis=direction, q=gradient, edges=edges))
            for which, owner in boards:
                for scale in (1, 3) if gradient in (0, 4) else (1,):
                    clocks.append(
                        dict(
                            board=which,
                            owner=owner,
                            axis=direction,
                            q=gradient,
                            scale=scale,
                            reads=clock(
                                rng, which, owner, psi, scale, [DELTA, 2 * DELTA]
                            ),
                        )
                    )
        for q in (0, 3):
            for transverse in (False, True) if direction == "port" else (False,):
                psi = preparation(direction, q, transverse)
                for which, owner in boards:
                    freq = clock(rng, which, owner, psi, 1, [1, 2, 3], both=True)
                    frequency.append(
                        dict(
                            board=which,
                            owner=owner,
                            axis=direction,
                            q=q,
                            transverse=transverse,
                            reads=freq,
                        )
                    )
    return dict(
        shots=SHOTS,
        calibration=calibration,
        local=local,
        energy=clocks,
        frequency=frequency,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    dest = HERE / "records.json"
    census = json.dumps(assay(), separators=(",", ":")) + "\n"
    if args.write:
        dest.write_text(census)
    elif census != dest.read_text():
        raise SystemExit("Seam symbol records differ")
