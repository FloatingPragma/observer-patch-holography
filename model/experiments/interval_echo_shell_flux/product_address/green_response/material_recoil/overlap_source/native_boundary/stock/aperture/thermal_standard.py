"""Measure the blank-chart uniform double's insertion energy against a parked rail."""

import argparse
import json
from hashlib import sha256

import numpy as np

from . import thermal_pair as exp

HERE = exp.HERE


def standing_generator(board, raw, load):
    """Four times the six-slot standing-energy operator, including unused contacts."""
    matrix = np.zeros((board.n, board.n), int)
    for a, b in board.seams:
        for site, word in ((a, 3 * load), (b, load)):
            for outcome in np.flatnonzero(raw[:, word]):
                assert outcome in (load, 3 * load)
                matrix[a if outcome == 3 * load else b, site] += raw[outcome, word]
    # The standing six-slot boundary completion counts each unused contact's
    # ordinary squared-load energy. The isolated comparison rail has zero contacts.
    matrix += np.diag(4 * exp.motion.K * load**2 * (6 - board.degree))
    np.testing.assert_array_equal(matrix, matrix.T)
    return matrix


def measure(board, raw, load, seed):
    matrix = standing_generator(board, raw, load)
    row = matrix @ np.ones(board.n, dtype=int)
    assert np.all(row == row[0])
    energy = row[0] / 4
    time = 1 / 256
    high = np.max(np.abs(matrix).sum(axis=1)) / 4
    assert high * time < np.pi
    n = 2**32
    rng = np.random.default_rng(seed)
    counts = rng.binomial(
        n, (1 + np.array([np.cos(time * energy), np.sin(time * energy)])) / 2
    )
    z = complex(*(2 * counts / n - 1))
    err = np.arcsin(2 * np.sqrt(exp.contact.RISK / n) / abs(z))
    measured = (np.angle(z) + np.array([-err, err])) / time
    assert measured[0] <= energy <= measured[1]
    return {
        "ports": board.n,
        "uniform_energy": float(energy),
        "ticks": time,
        "energy_band": measured.tolist(),
        "counts": counts.tolist(),
        "shots": n,
        "spectral_bound": float(high),
    }


def run():
    board = exp.motion.cold.federation(0)
    raw = exp.contact.stock.kernels()[2]
    return {
        "double": measure(board, raw, 2, 677443),
        "single": measure(board, raw, 1, 677444),
        "ordinary_rail": "disconnected storage at the same conserved load",
        "sha256": {
            "thermal_energy.json": sha256(
                (HERE / "thermal_energy.json").read_bytes()
            ).hexdigest()
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    result = run()
    path = HERE / "thermal_standard.json"
    if args.write:
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        exp.contact.timing.lab.lab.cold.apparatus.check(
            json.loads(path.read_bytes()), result
        )
    print(result)
