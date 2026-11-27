"""Exact ground brackets and prospective finite-clock comparison on two excesses."""

import argparse
import json
from fractions import Fraction
from itertools import product
from math import log2, pi
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.stats import beta

from ...row_16_binding.signed_pair.material import Assembly
from .measure import CONTRASTS, HERE, SHOTS

SEED = 26112696
STEP = 1 / 128


def sector(n, contrast):
    h = Assembly(n, 0, contrast).h.toarray()
    twice = np.rint(2 * h).astype(np.int64)
    assert np.array_equal(twice, 2 * h) and np.all(
        twice - np.diag(twice.diagonal()) <= 0
    )
    values, states = eigh(h)
    x = states[:, 0]
    if x.sum() < 0:
        x = -x
    if x.min() <= 0:
        raise ValueError("Numerical proposal lacks positive components")
    positive = [max(1, round(float(z) * 10**12)) for z in x]
    # Products are Python integers so the enclosure is independent of machine word size.
    hx = [sum(int(v) * w for v, w in zip(row, positive, strict=True)) for row in twice]
    ratios = [Fraction(v, 2 * w) for v, w in zip(hx, positive, strict=True)]
    return values, (min(ratios), max(ratios)), positive


def interval(k):
    p = 0.01 / (2 * 64)
    return [
        float(beta.ppf(p, k, SHOTS - k + 1)) if k else 0.0,
        float(beta.isf(p, k + 1, SHOTS - k)) if k < SHOTS else 1.0,
    ]


def gap(one, two):
    return [2 * one[0] - two[1], 2 * one[1] - two[0]]


def acquire():
    rng = np.random.default_rng(SEED)
    rows = []
    delays = None
    for contrast in CONTRASTS:
        a, one, x = sector(1, contrast)
        b, two, y = sector(2, contrast)
        separation = gap(one, two)
        if delays is None:
            size = float(max(Fraction(0), separation[0], -separation[1]))
            need = 2 * pi / STEP * 8 / size if size else 2**22
            delays = min(2**21, max(8192, 2 ** (1 + int(np.floor(log2(need))))))
        refs = [
            round(float(sum(br) / 2), 5) for br in (one, two, [2 * z for z in one])
        ] + [0.0]
        spectra = dict(
            one=a, two=b, separate=np.add.outer(a, a).ravel(), host=np.array([0.0])
        )
        banks = []
        for name, (ref_id, ref) in product(spectra, enumerate(refs)):
            z = (spectra[name] - ref) * STEP / (2 * pi)
            amp = np.sinc(delays * z) / np.sinc(z)
            prob = float(np.mean(amp**2))
            bright = int(rng.binomial(SHOTS, min(1.0, prob)))
            banks.append(
                dict(
                    material=name, reference_id=ref_id, counts=[bright, SHOTS - bright]
                )
            )
        rows.append(
            dict(
                k=contrast,
                one=[str(z) for z in one],
                two=[str(z) for z in two],
                positive_one=x,
                positive_two=y,
                refs=refs,
                banks=banks,
            )
        )
    return dict(seed=SEED, shots=SHOTS, step=STEP, delays=delays, rows=rows)


def decode(record):
    out = []
    for row in record["rows"]:
        one, two = ([Fraction(z) for z in row[key]] for key in ("one", "two"))
        low, high = gap(one, two)
        rates = {(r["material"], r["reference_id"]): r["counts"] for r in row["banks"]}
        if len(rates) != 16 or any(
            len(v) != 2 or min(v) < 0 or sum(v) != SHOTS for v in rates.values()
        ):
            raise ValueError("Incomplete clock banks")
        idx = 1 if low > 0 else 2
        signals = [rates[key, idx][0] for key in ("two", "separate")]
        out.append(
            dict(
                k=row["k"],
                one=[float(z) for z in one],
                two=[float(z) for z in two],
                separation=[float(low), float(high)],
                lower_reference=row["refs"][idx],
                joined_separate_bright=signals,
                limits=[interval(k) for k in signals],
            )
        )
    return dict(
        delays=record["delays"], first_zero=2 * pi / (STEP * record["delays"]), rows=out
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--acquire", type=Path)
    args = p.parse_args()
    if args.acquire:
        args.acquire.write_text(
            json.dumps(acquire(), separators=(",", ":"), allow_nan=False) + "\n"
        )
    else:
        print(
            json.dumps(
                decode(json.loads((HERE / "two_counts.json").read_text())), indent=2
            )
        )
