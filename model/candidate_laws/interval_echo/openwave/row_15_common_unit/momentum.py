"""Outer enclosures for every momentum share consistent with a target mass ratio."""

import gzip
import json

import numpy as np

from . import measure

PARTS = 100_000


def enclose(stocks):
    n0, n1 = measure.interval(stocks["neutral"]["curvature"])
    p0, p1 = measure.interval(stocks["positive"]["curvature"])
    v0, v1 = measure.interval(stocks["vacancy"]["curvature"])
    if min(n0, p0, v0) <= 0:
        raise ValueError("The enclosure assumes positive source curvatures")
    ticks = np.arange(PARTS + 1) / PARTS
    a, b = ticks[:-1], ticks[1:]
    # Bound each term separately throughout the closed bin. This can admit
    # extra f values along an edge and never discards a true positive split.
    lo = p0 * a * a + v0 * (1 - b) ** 2
    hi = p1 * b * b + v1 * (1 - a) ** 2
    keep = (lo <= measure.TARGET * n1 + 1e-11) & (hi >= measure.TARGET * n0 - 1e-11)
    transitions = np.flatnonzero(np.diff(np.r_[False, keep, False]))
    return [[int(a), int(b)] for a, b in transitions.reshape(-1, 2)]


def possible_splits():
    r = json.loads(gzip.decompress((measure.HERE / "counts.json.gz").read_bytes()))
    return {k: enclose(v["sources"]) for k, v in r["materials"].items()}


if __name__ == "__main__":
    print("Admissible closed cells in units of 1/100000:", possible_splits())
