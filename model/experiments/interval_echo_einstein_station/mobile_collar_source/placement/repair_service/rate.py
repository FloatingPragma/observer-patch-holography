"""Scope of the memory charge and the loss of a fixed patch at positive exchange rate."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json

import numpy as np
from model.candidate_laws.interval_echo import observations as obs
from . import compact, local
from .measure import HERE, compare


def derive():
    _, w, v = local.operators()
    identity = 64 * np.eye(len(w), dtype=np.int64)
    # Integer polynomial identities prove that norm equality implies a
    # fixed vector. No eigenvalue -1 is available in either sharp service.
    z = np.zeros_like(w)
    np.testing.assert_array_equal(w @ (w - identity) @ (2 * w + identity), z)
    np.testing.assert_array_equal(v @ (v - identity) @ (2 * v - identity), z)
    complemented = [
        (0 if m is None else local.MEM.index((1 - m[0], -m[1]))) * 4 + 3 - k
        for m in local.MEM
        for k in range(4)
    ]
    for gate in (w, v):
        np.testing.assert_array_equal(gate[np.ix_(complemented, complemented)], gate)
    for xy in range(4):
        assert v[xy, (xy % 2) * 2 + xy // 2] > 0
    # Measure the compact stationary source turnover by actual local
    # birth/return probabilities at each active edge; a fixed distribution
    # of each independent gate has the same expectation at every stage.
    b = obs.federation(0)
    rows = {}
    for den in (1, 2, 8, 32):
        row = compact.solve(b, "retained_exchange", Fraction(1, den))
        rows[str(den)] = row
    return dict(
        protocol_sha256=sha256((HERE / "RATE.md").read_bytes()).hexdigest(), rows=rows
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    result = derive()
    path = HERE / "rate.json"
    if args.write:
        path.write_text(json.dumps(result, indent=2) + "\n")
    else:
        compare(json.loads(path.read_text()), result)
    for d, row in result["rows"].items():
        print(
            d,
            row["first_cycle_contrast_below_quarter"],
            row["cycle256"]["near_minus_far"],
        )
