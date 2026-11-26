"""Compare physical clock windows on single signed-daughter trial banks."""

import json
from itertools import combinations, product

from ..joint_receivers.readback import interval
from .measure import BASES, BENCH, HERE, MODES, OWNERS, SEED, STARTS, TIMES


def checked(rows, keys, length, size):
    banks = {tuple(r["key"]): r["counts"] for r in rows}
    if (
        len(rows) != len(banks)
        or set(banks) != set(keys)
        or any(
            len(c) != length
            or any(type(k) is not int or k < 0 for k in c)
            or sum(c) != size
            for c in banks.values()
        )
    ):
        raise ValueError("missing or invalid count bank")
    return banks


def bands(c):
    n = sum(c[:7])
    return [interval(k, n) if n else (0.0, 1.0) for k in c[:3]]


def separation(a, b):
    return any(
        lower > other_upper or other_lower > upper
        for (lower, upper), (other_lower, other_upper) in zip(
            bands(a), bands(b), strict=True
        )
    )


def read(data):
    if (data["starts"], data["calibration_starts"], data["seed"]) != (
        STARTS,
        BENCH,
        SEED,
    ):
        raise ValueError("different cohort or seed")
    b = checked(
        data["banks"], product(MODES, OWNERS, TIMES, BASES, range(3)), 12, STARTS
    )
    cal = checked(data["standards"], product((False, True), range(4)), 7, BENCH)
    resolved = []
    for owner, t, basis in product(OWNERS, TIMES, BASES):
        for a, c in combinations(range(3), 2):
            resolved.append(
                dict(
                    key=[owner, t, basis, a, c],
                    distinct=separation(
                        b["native", owner, t, basis, a], b["native", owner, t, basis, c]
                    ),
                )
            )
    controls = {}
    for mode in ("position", "one_contact", "contrast_four"):
        controls[mode] = sum(
            separation(b["native", *key], b[mode, *key])
            for key in product(OWNERS, TIMES, BASES, range(3))
        )
    totals = {}
    for mode in MODES:
        counts = [c for key, c in b.items() if key[0] == mode]
        totals[mode] = dict(
            accepted=sum(sum(c[:7]) for c in counts),
            capture_limits=[
                min(sum(c[:7]) for c in counts),
                max(sum(c[:7]) for c in counts),
            ],
            outside=sum(c[7] for c in counts),
            refused_swap=sum(c[8] for c in counts),
        )
    return dict(
        pairs=resolved,
        comparisons=controls,
        totals=totals,
        exposure=STARTS * 24,
        calibration_zero=[cal[q, 3] for q in (False, True)],
    )


if __name__ == "__main__":
    print(json.dumps(read(json.loads((HERE / "counts.json").read_text())), indent=2))
