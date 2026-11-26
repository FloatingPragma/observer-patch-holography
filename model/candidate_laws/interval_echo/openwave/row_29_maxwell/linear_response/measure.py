"""Acquire edge and triangular-current outcomes on a freely released scalar material."""

import gzip
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from model.executable_model import carrier
from .. import experiment
from ..arrival_probe.measure import oriented_faces

HERE = Path(__file__).resolve().parent
SHOTS = 2**32
SEED = 290772
TIMES = (4, 16, 32, 64)
ANGLES = np.pi / np.array([16, 8, 4])
SOURCES = ("gateway_kick", "circular_kick", "circular_arm")
GAINS = (1.0, 0.5, 0.0)
EDGE = np.array(carrier.cell().seams)
FACE = np.array([f for f, _ in oriented_faces(carrier.cell(), experiment.geometry())])


def ports(psi):
    """Return +/- current, other local and external outcome probabilities."""
    a, b = psi[EDGE].T
    edge = np.column_stack((np.abs(a - 1j * b) ** 2 / 2, np.abs(a + 1j * b) ** 2 / 2))
    modes = np.fft.fft(psi[FACE], axis=1) / np.sqrt(3)
    pair = np.concatenate((edge, np.abs(modes[:, 1:]) ** 2))
    arrived = np.vdot(psi, psi).real
    return np.column_stack(
        (pair, arrived - pair.sum(axis=1), np.full(len(pair), 1 - arrived))
    )


def counts(prob, rng):
    if prob.min() < -2e-11 or np.abs(prob.sum(axis=-1) - 1).max() > 1e-11:
        raise ValueError("incomplete or negative detector probabilities")
    prob = np.maximum(prob, 0).round(13)
    prob /= prob.sum(axis=-1, keepdims=True)
    return rng.multinomial(SHOTS, prob)


def chart(level):
    board = carrier.federation(level)
    route, gate = experiment.path_and_gates(board)
    x, y, _ = experiment.frame(gate[0]).T
    sources = np.zeros((board.n, 2), dtype=complex)
    sources[gate[0], 0] = 1
    sources[:12, 1] = (x + 1j * y) / np.sqrt(2)
    h = experiment.candidate.generator(board)
    kinetic = board.laplacian() / 2
    residual = (
        h - kinetic - 6 * experiment.candidate.KAPPA * sparse.eye(board.n, format="csr")
    )
    if np.max(np.abs(residual.data), initial=0) > 1e-10:
        raise ValueError("candidate generator lost its one-unit reduction")
    return board, route, sources, kinetic


def replay():
    rng = np.random.default_rng(SEED)
    # owner, gain, wait, route position, source, angle, reversal, meter, exit
    raw = np.empty((2, 3, 4, 3, 3, 3, 2, 50, 4), dtype=np.int64)
    ground = np.empty((2, 3, 50, 4), dtype=np.int64)
    for level in (0, 1):
        board, route, sources, kinetic = chart(level)
        g = 1 / np.sqrt(board.n)
        for station in range(3):
            ground[level, station] = counts(ports(np.full(12, g)), rng)
        for law, gain in enumerate(GAINS):
            waves = sources.copy()
            previous = 0
            for ti, t in enumerate(TIMES):
                if gain:
                    waves = expm_multiply(-1j * gain * (t - previous) * kinetic, waves)
                previous = t
                for station, owner in enumerate(route[1:]):
                    point, circ = waves[12 * owner : 12 * owner + 12].T
                    for ai, theta in enumerate(ANGLES):
                        for si, sign in enumerate((1, -1)):
                            bank = raw[level, law, ti, station, :, ai, si]
                            kicks = (
                                g + g * np.expm1(1j * sign * theta) * point,
                                np.full(12, g),
                                np.cos(theta) * g + sign * np.sin(theta) * circ,
                            )
                            for i, psi in enumerate(kicks):
                                bank[i] = counts(ports(psi), rng)
    return {"flight": raw, "ground": ground}


def report(data):
    # There are fewer than 66,000 independent current expectations. Every
    # raw score lies within [-sqrt(3)/2, sqrt(3)/2] in full-exchange units.
    n = data["flight"].size // 4 + data["ground"].size // 4
    unit = np.sqrt(3 * np.log(2 * n / 0.01) / (2 * SHOTS))
    scale = np.r_[np.full(30, 0.5), np.full(20, np.sqrt(3) / 2)]
    flight = (data["flight"][..., 0] - data["flight"][..., 1]) * scale / SHOTS
    ground = (data["ground"][..., 0] - data["ground"][..., 1]) * scale / SHOTS
    odd = (flight[..., 0, :] - flight[..., 1, :]) / 2
    even = flight.mean(axis=-2) - ground[:, None, None, :, None, None, :]
    rows = []
    for li in range(2):
        for gi, gain in enumerate(GAINS):
            for station in (0, 1, 2):
                # Same reader counts yield the material current when calibrated
                # with this law's hopping gain. Zero gain is a stationary bank.
                o = odd[li, gi, :, station]
                e = even[li, gi, :, station]
                rows.append(
                    dict(
                        cells=(20, 80)[li],
                        gain=gain,
                        hops=station + 1,
                        max_odd_edge=float(np.max(np.abs(o[..., :30]))),
                        max_odd_face=float(np.max(np.abs(o[..., 30:]))),
                        max_even_face=float(np.max(np.abs(e[..., 30:]))),
                        max_circular_kick_odd_edge=float(
                            np.max(np.abs(o[:, 1, :, :30]))
                        ),
                    )
                )
    return {
        "raw_current_bound": float(unit),
        "original_launches_per_meter": SHOTS,
        "axes": [
            "board",
            "gain",
            "wait",
            "receiving_owner",
            "source",
            "angle",
            "sign",
            "edge_then_face",
            "outcome",
        ],
        "panels": rows,
    }


def main():
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = replay()
    out = json.dumps(report(data), indent=2, allow_nan=False) + "\n"
    if args.write:
        for key, a in data.items():
            with gzip.GzipFile(HERE / (key + ".npy.gz"), "wb", mtime=0) as dest:
                np.save(dest, a, allow_pickle=False)
        (HERE / "results.json").write_text(out)
    else:
        for key, a in data.items():
            with gzip.open(HERE / (key + ".npy.gz"), "rb") as src:
                np.testing.assert_array_equal(a, np.load(src, allow_pickle=False))
        assert (HERE / "results.json").read_text() == out


if __name__ == "__main__":
    main()
