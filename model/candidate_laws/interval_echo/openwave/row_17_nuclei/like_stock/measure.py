"""Whole-inventory port surveys and one-contact energy measurements."""

import argparse
import json
from collections import Counter
from itertools import combinations
from math import comb, log, sqrt
from pathlib import Path

import numpy as np

from model.candidate_laws.contrast_memory_spherical_v1 import law
from ...row_14_profile.measure import EDGES, native_geometry
from ...row_16_binding.signed_pair.material import Assembly

HERE = Path(__file__).resolve().parent
SHOTS = 2**24
SEED = 26112695
CONTRASTS = (8, 4, 0)


def shore_reads():
    echo = law.echo()
    result = []
    for n in range(3):
        col = law.SRC.index((n, 1))
        p = (abs(echo[9:, col].reshape(10, 9)) ** 2).sum(axis=1)
        signed = {law.LABELS[i][1] for i in np.flatnonzero(p)}
        if len(signed) != 1 or not np.isclose(p.sum(), 1, atol=1e-14, rtol=0):
            raise ValueError("Sharp stock lacks a unique signed shore read")
        result.append(signed.pop())
    return result


def survey(sign, n, reads, mask):
    census = Counter()
    candidates = []
    for support in combinations(range(12), n):
        stock = np.ones(12, dtype=int)
        stock[list(support)] += sign
        q = np.array(reads)[stock]
        if not np.array_equal(q, stock - 1):
            raise ValueError("Port standard changed")
        a, b = EDGES.T
        boundary = int(np.count_nonzero(q[a] - q[b]))
        ref = int(np.sum((q[mask] + 1) ** 2 - 1))
        census[boundary + ref, boundary] += 1
        candidates.append((boundary + ref, boundary, support))
    return census, candidates


def acquire():
    reads = shore_reads()
    mask, _ = native_geometry()
    rng = np.random.default_rng(SEED)
    rows = []
    for sign in (-1, 1):
        for n in range(7):
            census, candidates = survey(sign, n, reads, mask)
            for k in CONTRASTS:
                v, boundary, support = min(
                    candidates, key=lambda x: (k * x[0] + x[1] / 2, x[2])
                )
                material = Assembly(n if sign == 1 else 0, n if sign == -1 else 0, k)
                ids = np.flatnonzero(
                    np.all(material.w == sign * np.isin(np.arange(12), support), axis=1)
                )
                if len(ids) != 1:
                    raise ValueError("Selected word does not match the material")
                state = np.zeros(len(material.w))
                state[ids[0]] = 1.0
                raw = material.contact_bank(state, rng, SHOTS)
                rows.append(
                    dict(
                        sign=sign,
                        n=n,
                        k=k,
                        support=support,
                        min_potential=k * min(v for v, _ in census),
                        sharp_mean=k * v + boundary / 2,
                        census=[[a, b, c] for (a, b), c in sorted(census.items())],
                        meter=dict(levels=raw["levels"], counts=raw["counts"]),
                    )
                )
    return dict(seed=SEED, shots=SHOTS, signed_shore=reads, rows=rows)


def partitions(n, least=1):
    if n == 0:
        yield ()
    for i in range(least, n + 1):
        for rest in partitions(n - i, i):
            yield (i, *rest)


def decode(data):
    bounds, checks = {}, []
    if data["shots"] != SHOTS or data["signed_shore"] != [-1, 0, 1]:
        raise ValueError("Changed readout")
    for r in data["rows"]:
        sign, n, k = r["sign"], r["n"], r["k"]
        levels, counts = r["meter"]["levels"], r["meter"]["counts"]
        if sum(v[2] for v in r["census"]) != comb(12, n) or sum(counts) != SHOTS:
            raise ValueError("Missing stock or energy record")
        if len(levels) != len(counts) or any(
            type(c) is not int or c < 0 for c in counts
        ):
            raise ValueError("Invalid energy bins")
        observed = 39 * sum(x * y for x, y in zip(levels, counts, strict=True)) / SHOTS
        lo = min(0, k * ((1 + sign) ** 2 - 1))
        hi = max(k + 1, k * ((1 + sign) ** 2 - 1))
        error = 39 * (hi - lo) * sqrt(log(2 * 64 / 0.01) / (2 * SHOTS))
        checks.append(abs(observed - r["sharp_mean"]) <= error)
        lower = min(k * v for v, _, _ in r["census"])
        if lower != r["min_potential"] or (sign, n, k) in bounds:
            raise ValueError("Invalid potential census")
        bounds[sign, n, k] = (lower, observed + error)
    if len(bounds) != 42 or not all(checks):
        raise ValueError(
            "Incomplete binding survey or failed independent energy calibration"
        )
    rows = []
    for sign in (-1, 1):
        for k in CONTRASTS:
            for n in range(2, 7):
                parts = [p for p in partitions(n) if len(p) > 1]
                threshold = [
                    min(sum(bounds[sign, j, k][i] for j in p) for p in parts)
                    for i in (0, 1)
                ]
                lo, hi = bounds[sign, n, k]
                rows.append(
                    dict(
                        sign=sign,
                        k=k,
                        n=n,
                        ground=[lo, hi],
                        partitions=len(parts),
                        threshold=threshold,
                        separation=[threshold[0] - hi, threshold[1] - lo],
                    )
                )
    return dict(energy_checks=sum(checks), rows=rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--acquire", type=Path)
    args = parser.parse_args()
    if args.acquire:
        args.acquire.write_text(
            json.dumps(acquire(), separators=(",", ":"), allow_nan=False) + "\n"
        )
    else:
        print(
            json.dumps(decode(json.loads((HERE / "counts.json").read_text())), indent=2)
        )
