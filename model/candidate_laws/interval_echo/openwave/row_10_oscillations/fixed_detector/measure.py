"""Check stationary detector response using only physical terminal counts."""

import argparse
import json
from itertools import combinations, permutations, product
from pathlib import Path

from ..joint_receivers.readback import interval
from ..portable_decoder.measure import MODES, OWNERS, PANELS, STARTS, acquire

HERE = Path(__file__).resolve().parent
TIMES = (320, 640, 960, 1280)
SEED = 26112691
SCORES = sorted({v for m in (0, 0.25, 0.5, 0.75, 1) for v in permutations((0, m, 1))})


def bands(counts):
    pairs = sum(counts[:3])
    return [interval(k, pairs) if pairs else (0.0, 1.0) for k in counts[:3]]


def score_range(box, weights):
    """Use sum p=1 as well as the component bounds. Missing cohorts span [0,1]."""
    low = max(
        sum(w * lo for w, (lo, hi) in zip(weights, box, strict=True)),
        1 - sum((1 - w) * hi for w, (lo, hi) in zip(weights, box, strict=True)),
    )
    high = min(
        sum(w * hi for w, (lo, hi) in zip(weights, box, strict=True)),
        1 - sum((1 - w) * lo for w, (lo, hi) in zip(weights, box, strict=True)),
    )
    return low, high


def strongest(first, second):
    candidates = []
    for weights in SCORES:
        a = [score_range(box, weights) for box in first]
        b = [score_range(box, weights) for box in second]
        left = [min(v[0] for v in a), max(v[1] for v in a)]
        right = [min(v[0] for v in b), max(v[1] for v in b)]
        gap = max(left[0] - right[1], right[0] - left[1])
        candidates.append(dict(gap=gap, weights=weights, ranges=[left, right]))
    return max(candidates, key=lambda row: row["gap"])


def decode(record):
    raw = record["banks"]
    banks = {tuple(row["key"]): row["counts"] for row in raw}
    expected = set(product(MODES, OWNERS, TIMES, range(3), PANELS))
    if (
        set(banks) != expected
        or len(raw) != len(expected)
        or record["starts"] != STARTS
        or any(
            len(v) != 7
            or sum(v) != STARTS
            or any(type(k) is not int or k < 0 for k in v)
            for v in banks.values()
        )
    ):
        raise ValueError("Missing or invalid original-start count bank")
    windows = {key: bands(v) for key, v in banks.items()}
    results = []
    for mode, owner, panel in product(MODES, OWNERS, PANELS):
        for first, second in combinations(TIMES, 2):
            by_time = [
                [windows[mode, owner, t, j, panel] for j in range(3)]
                for t in (first, second)
            ]
            results.append(
                dict(
                    mode=mode,
                    owner=owner,
                    panel=panel,
                    times=[first, second],
                    **strongest(*by_time),
                )
            )
    yields = [sum(v[:3]) * 1e12 / STARTS for k, v in banks.items() if k[0] == "native"]
    return dict(
        comparisons=results,
        native_pairs_per_trillion=[min(yields), max(yields)],
        resolved={
            mode: sum(r["mode"] == mode and r["gap"] > 0 for r in results)
            for mode in MODES
        },
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--acquire", type=Path)
    args = p.parse_args()
    if args.acquire:
        args.acquire.write_text(
            json.dumps(acquire(TIMES, SEED), separators=(",", ":"), allow_nan=False)
            + "\n"
        )
    else:
        print(
            json.dumps(
                decode(json.loads((HERE / "counts.json").read_text())),
                indent=2,
                allow_nan=False,
            )
        )
