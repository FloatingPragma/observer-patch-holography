"""Acquire complete terminal event counts from two separate native winding stocks."""

import argparse
from itertools import combinations, product
from pathlib import Path

import numpy as np

from .sampling import integer_draws
from . import apparatus as lab

HERE = Path(__file__).parent
SEED = 11262731
SINGLE_SHOTS = 2**29
PAIR_SHOTS = 2**28
COORDS = np.array(list(combinations(range(12), 2)))
CYCLE = np.array([(1, 11), (0, 1), (1, 10)])
CENTERS = (1, 5)
CONTROLS = ((1, np.pi / 3), (0.5, np.pi / 2), (1, 0))


def singles(rng):
    rows = []
    for station, mark, sign in product(range(2), range(12), (-1, 1)):
        pmf = lab.banks((station,), (mark, mark), (sign, sign), 1, np.pi / 2, 1, COORDS)
        rows.append(integer_draws(rng, pmf, SINGLE_SHOTS))
    return np.array(rows).reshape(2, 12, 2, len(COORDS), 9)


def orders(rng, center, rate, time, quadrature, coords):
    anti = lab.geometry()[1]
    rows = []
    for order, signs, ends in product(
        ((0, 1), (1, 0)), product((-1, 1), repeat=2), product((0, 1), repeat=2)
    ):
        marks = [
            anti[q] if end else q for q, end in zip((0, center), ends, strict=True)
        ]
        pmf = lab.banks(order, marks, signs, rate, time, quadrature, coords)
        rows.append(integer_draws(rng, pmf, PAIR_SHOTS))
    return np.array(rows).reshape(2, 2, 2, 4, len(coords), 9)


def acquire():
    rng = np.random.default_rng(SEED)
    paths = lab.layout()["contacts"]
    standards = [integer_draws(rng, lab.standards(h), SINGLE_SHOTS) for h in (1, 0.5)]
    one = singles(rng)
    two = [orders(rng, q, 1, np.pi / 2, 0, COORDS) for q in CENTERS]
    held = [
        [orders(rng, 1, rate, time, phase, CYCLE) for phase in (0, 1)]
        for rate, time in CONTROLS
    ]
    return dict(
        standards=np.array(standards),
        single=one,
        double=np.array(two),
        held=np.array(held),
        coords=COORDS,
        cycle=CYCLE,
        path_lengths=np.array([len(p) - 1 for p in paths]),
        shots=np.array([SINGLE_SHOTS, PAIR_SHOTS]),
        seed=np.array(SEED),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = acquire()
    path = HERE / "counts.npz"
    if args.write:
        with path.open("xb") as out:
            np.savez_compressed(out, **data)
    else:
        with np.load(path) as saved:
            if set(saved.files) != set(data):
                raise ValueError("Different archive keys")
            for k in data:
                np.testing.assert_array_equal(saved[k], data[k])
    from .readback import decode
    import json

    print(json.dumps(decode(data), indent=2))


if __name__ == "__main__":
    main()
