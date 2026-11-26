"""Acquire the enrolled phase, flag, public-record, hand and raw-reset counts."""

import gzip
import json
from itertools import combinations
from pathlib import Path

import numpy as np

from . import apparatus as a
from ..signed_effects.sequential import GRID_BITS, grid_counts

HERE = Path(__file__).resolve().parent
N = 2**20
SEED = 177031
WAITS = (1, 2, 3, 5, 8)
MODES = ("signed", "two", "latched", "magnitude", "half_exchange")


def counts(probabilities, *ids, seed=SEED):
    p = np.asarray(probabilities).real.ravel()
    if p.min() < -1e-12 or abs(p.sum() - 1) > 1e-12:
        raise ArithmeticError("Incomplete probability accounting")
    # Integer draws on the 2**-40 grid give the same counts on every platform.
    rng = np.random.default_rng([seed, *ids])
    return grid_counts(rng, N, np.maximum(p, 0)).tolist()


def read_bank(bank, analyser, ids, seed=SEED):
    labels, probabilities = [], []
    for key, rho in bank.items():
        diagonal = np.diag(analyser.conj().T @ rho @ analyser).real
        labels.append(list(key))
        probabilities.extend(diagonal)
    bars = np.array(counts(probabilities, *ids, seed=seed)).reshape(-1, 3)
    return {"last_publications": labels, "terminal_ports": bars.tolist()}


def acquire():
    flags = {}
    for strength in (1.0, 0.5):
        settings = [
            a.flag_probability(theta, strength) for theta in (np.pi / 3, np.pi / 2)
        ]
        flags[str(strength)] = [
            counts([1 - p, p], 0, int(strength * 2), i) for i, p in enumerate(settings)
        ]
    hand = {
        str(n): [counts(p, 1, n, j) for j, p in enumerate(a.hand_probabilities(n))]
        for n in WAITS
    }
    panels, sharp = [], []
    for modeid, mode in enumerate(MODES):
        effects = (
            [a.triangle_blocks(edge) for edge in a.DUAL]
            if mode == "two"
            else [a.triangle_blocks(a.FINE, mode == "magnitude")]
        )
        strength = 0.5 if mode == "half_exchange" else 1.0
        p = [
            a.flag_probability(theta, strength)
            for theta in (np.pi / 3, np.pi / 2)[: len(effects)]
        ]
        for i in range(3):
            input_rho = np.diag(np.eye(3)[i])
            bank = a.history(input_rho, effects, p, 5, mode == "latched")
            sharp.append(
                dict(mode=mode, address=i, **read_bank(bank, np.eye(3), [2, modeid, i]))
            )
        for pairid, (i, j) in enumerate(combinations(range(3), 2)):
            for iphase in (0, 1):
                source = a.pair_word(i, j, iphase * np.pi / 2)[:, i]
                rho = np.outer(source, source.conj())
                for n in WAITS:
                    bank = a.history(rho, effects, p, n, mode == "latched")
                    for ophase in (0, 1):
                        read = a.pair_word(i, j, ophase * np.pi / 2)
                        receipt = read_bank(
                            bank, read, [3, modeid, pairid, iphase, n, ophase]
                        )
                        panels.append(
                            dict(
                                mode=mode,
                                pair=[i, j],
                                start_phase=iphase,
                                read_phase=ophase,
                                ticks=n,
                                **receipt,
                            )
                        )
    raw = a.raw_total_two()
    controls = [counts(row, 4, i) for i, row in enumerate([*raw.T, raw.mean(axis=1)])]
    free = [counts(abs(a.dimer(n * np.pi)[:, 0]) ** 2, 5, n) for n in (0, 1, 2)]
    return dict(
        shots=N,
        sampling_grid_bits=GRID_BITS,
        addresses=a.PORTS,
        flags=flags,
        hand=hand,
        sharp=sharp,
        panels=panels,
        raw_sharp_and_mixed=controls,
        free_returns=free,
    )


def report(data):
    flags = {k: [v[1] / N for v in vs] for k, vs in data["flags"].items()}
    phases = {
        n: np.arctan2(2 * banks[1][0] / N - 1, 2 * banks[0][0] / N - 1)
        for n, banks in data["hand"].items()
    }
    held = {}
    for row in data["panels"]:
        i, j = row["pair"]
        mode = row["mode"]
        p, q = flags["0.5" if mode == "half_exchange" else "1.0"]
        multiplier = 1 - p
        if mode == "two":
            multiplier = (1 - p if i == 0 else 1) * (1 - q if 1 in (i, j) else 1)
        if mode == "magnitude" and (i, j) == (0, 1):
            multiplier = 1
        n = row["ticks"]
        predicted = (
            1 + multiplier**n * int(row["start_phase"] == row["read_phase"])
        ) / 2
        observed = np.sum(row["terminal_ports"], axis=0)[i] / N
        if n > 1:
            held[mode] = max(held.get(mode, 0), abs(predicted - observed))
    return dict(
        flags=flags,
        clock_phases=phases,
        rates_per_tick=[-np.log1p(-p) for p in flags["1.0"]],
        rates_per_radian=[-np.log1p(-p) / phases["1"] for p in flags["1.0"]],
        max_held_probability_error=held,
    )


def plain(data):
    return json.dumps(data, separators=(",", ":"), allow_nan=False).encode() + b"\n"


if __name__ == "__main__":
    data = acquire()
    path = HERE / "records.json.gz"
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    out = parser.parse_args().output
    if out:
        out.write_bytes(gzip.compress(plain(data), mtime=0))
    elif gzip.decompress(path.read_bytes()) != plain(data):
        raise ArithmeticError("Count replay differs")
    print(json.dumps(report(data), indent=2))
