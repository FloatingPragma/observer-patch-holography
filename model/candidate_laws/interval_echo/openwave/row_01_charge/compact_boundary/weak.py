"""Second variation at the parked bridge and a separately simulated weak-opening read."""

import argparse
import json
import math
from functools import cache
from itertools import product
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

from . import measure as m

GAINS = (0.0, 1 / 16, 1 / 8)
TIMES = (2, 4)
STARTS = 2**32
SEED = 26112702


def enrollment():
    return product(m.SIGNS, (-1, 1), GAINS, TIMES)


@cache
def source(ss, ps):
    signs = tuple(sorted(s for s in (ss, ps) if s))
    idx, h0 = m.sector(signs, 8, False)
    _, h1 = m.sector(signs, 8, True)
    interaction = (h1 - h0) / 8
    # A blank bridge exchanges two unequal trits with amplitude -1/2.
    assert np.max(np.asarray(abs(interaction).sum(axis=1))) == 0.5
    first = m.initial(idx, ss, ps)
    return -1j * h0 / 8, -1j * interaction, first, m.outcome_bins(idx)


def evolve(op, initial):
    return expm_multiply(
        op, initial, start=0, stop=4, num=3, traceA=op.diagonal().sum()
    )[1:]


@cache
def quadratic(ss, ps):
    h0, v, first, bins = source(ss, ps)
    augmented = sparse.bmat(
        [[h0, None, None], [v, h0, None], [None, v, h0]], format="csr"
    )
    results = {}
    for t, row in zip(
        TIMES, evolve(augmented, np.r_[first, 0 * first, 0 * first]), strict=True
    ):
        a, b, c = np.split(row, 3)
        # With the bridge parked, one insertion changes the interior census.
        first_count = int(bool(ss))
        np.testing.assert_allclose(a[bins // 9 != first_count], 0, atol=1e-12, rtol=0)
        np.testing.assert_allclose(b[bins // 9 == first_count], 0, atol=1e-12, rtol=0)
        zeroth = np.bincount(bins % 9, weights=abs(a) ** 2, minlength=9)
        second = np.bincount(
            bins % 9, weights=abs(b) ** 2 + 2 * (a.conj() * c).real, minlength=9
        )
        results[t] = (zeroth, second)
    return results


@cache
def full(ss, ps, g):
    h0, v, first, bins = source(ss, ps)
    results = {}
    for t, row in zip(TIMES, evolve(h0 + g * v, first), strict=True):
        p = np.bincount(bins, weights=abs(row) ** 2, minlength=27)
        np.testing.assert_allclose(sum(p), 1, atol=2e-11, rtol=0)
        results[t] = p / sum(p)
    return results


def calibrate(g, ps):
    # Isolated physical port pair and its native host and signed stock.
    block = m.blocks(8)["blank"].astype(float) / 8
    diagonal = np.diag(np.diag(block))
    u = expm(-1j * np.pi * (diagonal + g * (block - diagonal)))
    a, b = (m.LAW.SRC.index(w) for w in ((1 + ps, 1), (1, 1 + ps)))
    return float(abs(u[b, a]) ** 2)


def acquire():
    rng = np.random.default_rng(SEED)
    banks = []
    for ss, ps, g, t in enrollment():
        p = full(ss, ps, g)[t]
        banks.append([ss, ps, g, t, rng.multinomial(STARTS, p).tolist()])
    cal = [
        [ps, g, int(rng.binomial(STARTS, calibrate(g, ps)))]
        for ps, g in product((-1, 1), GAINS)
    ]
    return dict(seed=SEED, starts=STARTS, banks=banks, calibration=cal)


def report(raw):
    if raw["seed"] != SEED or raw["starts"] != STARTS:
        raise ValueError("unexpected exposure or seed")
    eps = math.sqrt(math.log(2 * 4096 / 0.01) / (2 * STARTS))
    cal = {}
    for ps, g, hits in raw["calibration"]:
        if (ps, g) in cal or type(hits) is not int or not 0 <= hits <= STARTS:
            raise ValueError("invalid calibration")
        cal[ps, g] = abs(hits / STARTS - math.sin(math.pi * g / 2) ** 2) <= eps
    if cal.keys() != set(product((-1, 1), GAINS)) or not all(cal.values()):
        raise ValueError("incomplete or failed intensity standards")
    covered = rigorous = total = 0
    largest_residual = max_second = 0.0
    first_order = 0
    seen = set()
    for ss, ps, g, t, counts in raw["banks"]:
        key = (ss, ps, g, t)
        if (
            key in seen
            or len(counts) != 27
            or any(type(c) is not int or c < 0 for c in counts)
            or sum(counts) != STARTS
        ):
            raise ValueError("invalid or duplicate bank")
        seen.add(key)
        obs = np.array(counts).reshape(3, 9).sum(axis=0) / STARTS
        base, second = quadratic(ss, ps)[t]
        # For exact unitary evolution ||d^3 rho_g/dg^3||_1 <= (2t ||V||)^3.
        remainder = (abs(g) * t) ** 3 / 6
        diff = abs(obs - base - g**2 * second)
        if g:
            covered += int(sum(diff <= eps))
            rigorous += int(sum(diff <= eps + remainder + 2e-11))
            first_order += int(sum(abs(obs - base) > eps))
            total += 9
            largest_residual = max(largest_residual, float(max(diff)))
            max_second = max(max_second, float(max(abs(g**2 * second))))
    if seen != set(enrollment()):
        raise ValueError("missing weak aperture")
    return dict(
        marginal_half_width=eps,
        held_bins=total,
        quadratic_with_count_error=covered,
        quadratic_with_remainder=rigorous,
        resolved_against_parked_prediction=first_order,
        max_central_residual=largest_residual,
        largest_predicted_quadratic_shift=max_second,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--acquire", type=Path)
    args = p.parse_args()
    if args.acquire:
        record = acquire()
        args.acquire.write_text(json.dumps(record, separators=(",", ":")) + "\n")
    else:
        record = json.loads((m.HERE / "weak_counts.json").read_text())
    print(json.dumps(report(record), indent=2))
