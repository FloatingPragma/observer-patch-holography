"""Can any single homogeneous conformal clock reproduce the paired pulse records?"""

import argparse
import json

import numpy as np

from .experiment import HERE, layout


def contrast(level, seed, dilute, reads, indices):
    first, second = reads[indices[0]], reads[indices[1]]
    final = reads[-1]
    _, _, _, stations, d0 = layout(level)
    misses = 0
    tested = 0
    witness = None
    for i, origin in enumerate(stations):
        targets = np.flatnonzero((final[i] >= 0) & (d0[i] > 0))
        codes = np.c_[d0[i, targets], first[i, targets]]
        keys, ids = np.unique(codes, axis=0, return_inverse=True)
        for g, (distance, stamp) in enumerate(keys):
            keep = targets[ids == g]
            seen = second[i, keep]
            assert seen.min() >= 0
            hist = np.bincount(seen - seen.min())
            covered = max(np.convolve(hist, np.ones(3, dtype=int), mode="full"))
            lost = len(keep) - int(covered)
            misses += lost
            tested += len(keep)
            span = int(seen.max() - seen.min())
            if witness is None or span > witness["span"]:
                witness = dict(
                    origin=origin,
                    initial_range=int(distance),
                    first=int(stamp),
                    destinations=[int(keep[seen.argmin()]), int(keep[seen.argmax()])],
                    second=[int(seen.min()), int(seen.max())],
                    span=span,
                )
    return dict(
        level=level,
        seed=seed,
        dilute=dilute,
        indices=indices,
        tested=tested,
        minimum_misses=misses,
        witness=witness,
    )


def run():
    cases = json.loads((HERE / "signals.json").read_text())["cases"]
    reports = []
    with np.load(HERE / "signals.npz") as raw:
        for p in cases:
            inputs = (p["level"], p["seed"], p["dilute"])
            reads = raw[f"{inputs[1]}_{int(inputs[2])}"]
            reports.extend(
                contrast(*inputs, reads, indices)
                for indices in ([0, 1], [1, 2], [5, 6])
            )
    return reports


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    path = HERE / "homogeneity.json"
    if args.write:
        path.write_text(json.dumps(result, indent=1) + "\n")
    else:
        assert result == json.loads(path.read_text())
        print("Verified the arrival-time obstruction without choosing a scale history.")


if __name__ == "__main__":
    main()
