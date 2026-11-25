"""Audit the common reference on every complete visible transition fibre."""

import argparse
import json
from fractions import Fraction

import numpy as np

from .measure import HERE, N, RADIUS, WORDS, forward_kernel, sweep_maps

EDGES = (3, 4, 5, 6)


def run():
    raw = sweep_maps(EDGES)
    forward = forward_kernel(EDGES)
    rng = np.random.default_rng(120074)
    rows = []
    for u in range(8):
        allow = np.array([i for i, w in enumerate(WORDS) if (w & 7) == u])
        ix = np.full(len(WORDS), -1, dtype=int)
        ix[allow] = np.arange(len(allow))
        table = ix[raw[:, allow]]
        if np.min(table) < 0:
            raise ValueError("Other visible ports changed")
        local = forward[np.ix_(allow, allow)]
        exact = np.array(
            [
                [
                    np.count_nonzero(table[:, i] == j) / len(table)
                    for j in range(len(allow))
                ]
                for i in range(len(allow))
            ]
        )
        np.testing.assert_allclose(local, exact, rtol=0, atol=4e-14)
        target = (np.asarray(WORDS)[allow] >> 3) & 1
        q = Fraction(int(target.sum()), len(target))
        if q != Fraction(4 - u.bit_count(), 5):
            raise ValueError("Bath fibre mismatch")
        for old in range(2):
            initial = int(np.flatnonzero(target == old)[0])
            forward_rate = float(np.linalg.matrix_power(local, 32)[initial] @ target)
            for port in (2, 3):
                # For port 2 transpose the adjacent rail labels before the bath service
                # and restore afterwards. u is in the coordinates of the three held rails.
                image = np.tile(np.arange(len(allow), dtype=np.uint8), (N, 1))
                hist = []
                while len(image):
                    image = np.take_along_axis(
                        image, table[rng.integers(len(table), size=len(image))], axis=1
                    )
                    done = (image == image[:, :1]).all(axis=1)
                    y = target[image[done, initial]]
                    hist.append(np.bincount(y, minlength=2).tolist())
                    image = image[~done]
                counts = np.sum(hist, axis=0)
                rows.append(
                    dict(
                        word=u,
                        old=old,
                        port=port,
                        expected=str(q),
                        outcome=counts.tolist(),
                        depths=hist,
                        ordinary_32=float(forward_rate),
                    )
                )
    # Exact probability of stopping on the first sweep and writing one on the exposed rail.
    # If a single sweep cannot coalesce, enumerate pairs of sweeps at the earliest such time.
    u = 7
    allow = np.array([i for i, w in enumerate(WORDS) if (w & 7) == u])
    collapse = raw[:, allow]
    depth1 = collapse[(collapse == collapse[:, :1]).all(axis=1)]
    if not len(depth1):
        raise ValueError("Expected a one-sweep witness")
    fraction = Fraction(sum((WORDS[i] >> 3) & 1 for i in depth1[:, 0]), len(depth1))
    return dict(
        shots=N,
        radius=RADIUS,
        rows=rows,
        first_sweep=dict(
            word=u,
            possible_tapes=len(raw),
            coalescing_tapes=len(depth1),
            conditional_one=str(fraction),
            reference_one="1/5",
        ),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    data = run()
    out = json.dumps(data, indent=2, allow_nan=False) + "\n"
    f = HERE / "transitions.json"
    if args.write:
        f.write_text(out)
    elif f.read_text() != out:
        raise ValueError("Transition replay mismatch")
    print("32 conditional read banks", data["first_sweep"])
