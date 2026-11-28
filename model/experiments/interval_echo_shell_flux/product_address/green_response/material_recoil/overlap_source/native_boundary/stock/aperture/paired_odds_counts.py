"""Test the equilibrium source with two local magnitude receipts per preparation."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import comb

import numpy as np
from scipy.optimize import brentq

from . import pair_contact as source

HERE = source.HERE
N = 1024
NULL = Fraction(59, 100)


def conditional_tail(table, ratio):
    """Exact lower tail after conditioning on both observed one-head totals."""
    (a, b), (c, k) = table
    row, col, n = c + k, b + k, a + b + c + k
    lo, hi = max(0, row + col - n), min(row, col)
    weights = [
        comb(row, j)
        * comb(n - row, col - j)
        * ratio.numerator ** (j - lo)
        * ratio.denominator ** (hi - j)
        for j in range(lo, hi + 1)
    ]
    return Fraction(sum(weights[: k - lo + 1]), sum(weights))


def read(banks):
    # Position in a log is its preparation identifier. A blank or
    # out-of-alphabet magnitude on either arm is recorded as a rejection.
    pairs = np.array([[int(x, 16) for x in s] for s in banks]).T
    magnitudes = source.stock.MAG[pairs]
    keep = ((pairs != 0) & (magnitudes != 1)).all(axis=1)
    cells = (magnitudes[keep] // 2) @ [2, 1]
    return np.bincount(cells, minlength=4).reshape(2, 2).tolist(), int((~keep).sum())


def bank(letters):
    return ["".join(format(x, "x") for x in side) for side in letters.T]


def run():
    previous = json.loads((HERE / "equilibrium_odds.json").read_bytes())
    floor = previous["results"]["full"]["any_inventory_from_filters"]["odds_floor"]
    assert float(NULL) < floor
    rng = np.random.default_rng(8164970)
    _, population = source.source_head_population()
    _, prior, _, _ = source.source_counts()
    p = prior / prior.sum()
    # Start on independent copies sampled from the 2200 source preparation.
    rows = []
    accepted = 0
    while accepted < N:
        pair = population[rng.integers(len(population))]
        rows.append(pair)
        if all(pair != 0) and all(source.stock.MAG[pair] != 1):
            accepted += 1
    arms = {"native": bank(np.array(rows))}
    active_a, active_b = p[1].sum(), p[:, 1].sum()
    codes = [
        int(np.flatnonzero((np.arange(11) != 0) & (source.stock.MAG == m))[0])
        for m in (0, 2)
    ]
    for tag, r in (("independent", 1.0), ("boundary", float(NULL))):
        x = brentq(
            lambda v: (1 - active_a - active_b + v) * v
            - r * (active_a - v) * (active_b - v),
            max(0, active_a + active_b - 1),
            min(active_a, active_b),
        )
        weights = [1 - active_a - active_b + x, active_b - x, active_a - x, x]
        draws = rng.choice(4, N, p=weights)
        arms[tag] = bank(np.array(codes)[np.column_stack([draws // 2, draws % 2])])
    out = {}
    for tag, logs in arms.items():
        table, reject = read(logs)
        assert sum(map(sum, table)) == N
        exact = conditional_tail(table, NULL)
        out[tag] = dict(
            local_head_logs=logs,
            table=table,
            rejects=reject,
            conditional_tail=float(exact),
            reject_at_one_percent=bool(exact <= Fraction(1, 100)),
        )
    return dict(
        sha256={
            p: sha256((HERE / p).read_bytes()).hexdigest()
            for p in ("PAIRED_ODDS_COUNTS.md", "equilibrium_odds.json")
        },
        null_minimum_odds=str(NULL),
        accepted_per_arm=N,
        results=out,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    opt = cli.parse_args()
    result = run()
    path = HERE / "paired_odds_counts.json"
    if opt.write:
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        assert json.loads(path.read_bytes()) == result
    for name, v in result["results"].items():
        print(name, v["table"], v["rejects"], v["conditional_tail"])
