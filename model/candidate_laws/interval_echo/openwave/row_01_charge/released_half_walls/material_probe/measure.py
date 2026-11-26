"""Draw both native rail clicks and the independent sharp wall herald in each trial."""

import argparse
import hashlib
import json
from functools import cache
from itertools import product
from pathlib import Path

import numpy as np

from .. import ring
from .flight import (
    ANALYZERS,
    BOARDS,
    DWELLS,
    mount,
    probability,
    sensed_sum,
    verify_plates,
)

HERE = Path(__file__).resolve().parent
SHOTS = 2**24
SEED = 26011262
LAWS = {  # flight contrast, sensing contrast, distant apparatus
    "native": (8, 8, False),
    "distant": (8, 8, True),
    "uncoupled": (8, 0, False),
    "flight_off": (0, 8, False),
    "half": (4, 4, False),
}


def draw(key, p, size=SHOTS):
    values = np.array(p, dtype=float)
    if np.min(values) < -1e-14 or abs(sum(values) - 1) > 2e-11:
        raise ValueError("Invalid native outcome probabilities")
    values = values.clip(0)
    label = json.dumps([SEED, *key], separators=(",", ":")).encode("ascii")
    rng = np.random.default_rng(
        int.from_bytes(hashlib.sha256(label).digest()[:16], "big")
    )
    return rng.multinomial(size, values / values.sum()).tolist()


@cache
def source(n, amount, source_sign, flight, bonds):
    p = ring.released(n, amount, flight)[list(ring.HOLDS).index(3)]
    bits = ring.basis(n)[1]
    if source_sign < 0:
        bits = 1 - bits
    domain = bits + np.roll(bits, -1, axis=1) - 1
    gate = ((domain != 0).sum(axis=1) == 2) & np.all(
        domain[:, bonds] == source_sign, axis=1
    )
    return amount * bits.astype(int), gate, p / p.sum()


def collect():
    verify_plates()
    standard, measurements, setups = [], [], []
    for lvl, n in BOARDS:
        stations, remote = mount(lvl, n)
        setups.append(
            dict(
                level=lvl,
                length=n,
                contacts=[list(v) for v in stations],
                distant=[list(v) for v in remote],
                ring=ring.geometry(lvl, n),
            )
        )
        bonds = tuple(row[0] for row in stations)
        for k, side, x, y, s, t, read in product(
            (4, 8), range(2), range(3), range(3), (-1, 1), range(2), ANALYZERS
        ):
            j, a, b = stations[side]
            loads = np.zeros(n, dtype=int)
            loads[j], loads[(j + 1) % n] = x, y
            p = probability(sensed_sum(lvl, n, (a, b), loads, s, k), DWELLS[t], read)
            key = [lvl, n, k, side, x, y, s, t, read]
            standard.append([key, draw(["standard", *key], [1 - p, p])])
        for (name, (flight, k, distant)), amount, source_sign, probe_sign in product(
            LAWS.items(), (1, 2), (-1, 1), (-1, 1)
        ):
            words, gate, weights = source(n, amount, source_sign, flight, bonds)
            probes = remote if distant else [row[1:] for row in stations]
            # Aggregate the measured detuning strata, preserving all gate failures.
            detunings = [
                np.broadcast_to(
                    sensed_sum(lvl, n, tip, words, probe_sign, k), len(words)
                )
                for tip in probes
            ]
            categories, inverse = np.unique(
                np.column_stack((gate, *detunings)), axis=0, return_inverse=True
            )
            mass = np.bincount(inverse.ravel(), weights=weights)
            for t, read in product(range(2), ANALYZERS):
                bins = np.zeros((2, 2, 2))
                for (h, dl, dr), w in zip(categories, mass, strict=True):
                    pl, pr = (
                        probability(dl, DWELLS[t], read),
                        probability(dr, DWELLS[t], read),
                    )
                    bins[int(h)] += w * np.outer([1 - pl, pl], [1 - pr, pr])
                key = [lvl, n, name, amount, source_sign, probe_sign, t, read]
                measurements.append([key, draw(["held", *key], bins.ravel())])
    return dict(
        seed=SEED,
        shots=SHOTS,
        mounts=setups,
        standard=standard,
        measurements=measurements,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    raw = collect()
    if args.output:
        args.output.write_text(json.dumps(raw, separators=(",", ":")) + "\n")
    elif raw != json.loads((HERE / "counts.json").read_text()):
        raise SystemExit("Count replay differs")
    print(
        len(raw["standard"]),
        "calibration and",
        len(raw["measurements"]),
        "joint wall banks",
    )
