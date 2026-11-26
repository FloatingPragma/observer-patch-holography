"""Predict the held ring measurements from independent head and local transport records."""

import argparse
import json
from itertools import product
from pathlib import Path

import numpy as np
from scipy.special import betaincinv, betainccinv

from .measure import (
    BASES,
    OWNERS,
    PANELS,
    SIZE,
    TIMES,
    analyzers,
    bank,
    sector_fields,
    source_weights,
)

FRESH = 26112611
TAIL = 0.01 / (2 * 16384)


def interval(k, n):
    if n <= 0 or not 0 <= k <= n:
        raise ValueError("Need a nonempty enrolled cohort and a valid count")
    return np.array(
        [
            betaincinv(k, n - k + 1, TAIL) if k else 0.0,
            betainccinv(k + 1, n - k, TAIL) if k < n else 1.0,
        ]
    )


def acquire():
    weights = source_weights()
    heads, stocks, held = [], [], []
    for j, w in enumerate(weights.T):
        heads.append(
            dict(
                key=[BASES[j // 3], j % 3],
                counts=bank([FRESH, "heads", j], [*w, 1 / 3 - w.sum(), 2 / 3]),
            )
        )
    for strength in (8, 4):
        field = sector_fields(strength, False)
        for owner in OWNERS:
            for i, t in enumerate(TIMES):
                for p, u in analyzers().items():
                    key = [strength, owner, t, p]
                    clicks = abs(field[owner][:, i, :] @ u.T) ** 2
                    for m, row in enumerate(clicks):
                        stocks.append(
                            dict(
                                key=[*key, m],
                                counts=bank(
                                    [FRESH, "stock", *key, m], [*row, 1 - row.sum()]
                                ),
                            )
                        )
                    for j, w in enumerate(weights.T):
                        final = w @ clicks
                        held.append(
                            dict(
                                key=[*key, j],
                                counts=bank(
                                    [FRESH, "held", *key, j],
                                    [
                                        *final,
                                        w.sum() - final.sum(),
                                        1 / 3 - w.sum(),
                                        2 / 3,
                                    ],
                                ),
                            )
                        )
    return dict(
        seed_key_prefix=FRESH,
        starts=SIZE,
        source_order=BASES,
        source_heads=heads,
        stock_clicks=stocks,
        held=held,
    )


def checked(rows, expected_width, n):
    if (
        len(rows) != expected_width
        or any(type(k) is not int or k < 0 for k in rows)
        or sum(rows) != n
    ):
        raise ValueError("Incomplete or nonintegral trial outcomes")
    return rows


def enrolled(records, expected, width, n):
    expected = set(expected)
    by_key = {tuple(r["key"]): checked(r["counts"], width, n) for r in records}
    if len(records) != len(expected) or set(by_key) != expected:
        raise ValueError("Missing, duplicated or unregistered cohort")
    return by_key


def decode(raw):
    n = raw["starts"]
    if (
        type(n) is not int
        or n != SIZE
        or raw["seed_key_prefix"] != FRESH
        or tuple(raw["source_order"]) != BASES
    ):
        raise ValueError("Unregistered exposure, seed or source order")
    sources = list(product(BASES, range(3)))
    rows = enrolled(raw["source_heads"], sources, 6, n)
    head = np.array([[interval(k, n) for k in rows[s][:4]] for s in sources])
    rows = enrolled(
        raw["stock_clicks"], product((8, 4), OWNERS, TIMES, PANELS, range(4)), 7, n
    )
    stock = {key: np.array([interval(k, n) for k in c[:6]]) for key, c in rows.items()}
    held = enrolled(raw["held"], product((8, 4), OWNERS, TIMES, PANELS, range(6)), 9, n)
    tally = {
        "compatible": 0,
        "held_bins": 0,
        "pooled_source_exclusions": 0,
        "pooled_source_bins": 0,
    }
    for strength in (8, 4):
        for owner in OWNERS:
            for t in TIMES:
                for p in PANELS:
                    key = (strength, owner, t, p)
                    detector = np.array([stock[*key, m] for m in range(4)])
                    for j in range(6):
                        counts = held[*key, j]
                        read = np.array([interval(k, n) for k in counts[:6]])

                        def excluded(source):
                            prediction = (source[:, None, :] * detector).sum(axis=0)
                            return (read[:, 0] > prediction[:, 1]) | (
                                read[:, 1] < prediction[:, 0]
                            )

                        tally["compatible"] += 6 - int(excluded(head[j]).sum())
                        tally["held_bins"] += 6
                        if j in (4, 5):
                            tally["pooled_source_exclusions"] += int(
                                excluded(head[3]).sum()
                            )
                            tally["pooled_source_bins"] += 6
    return tally


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    mode = cli.add_mutually_exclusive_group(required=True)
    mode.add_argument("--output", type=Path)
    mode.add_argument("--decode", type=Path)
    args = cli.parse_args()
    if args.output:
        args.output.write_text(json.dumps(acquire(), separators=(",", ":")) + "\n")
    else:
        print(json.dumps(decode(json.loads(args.decode.read_text())), indent=2))
