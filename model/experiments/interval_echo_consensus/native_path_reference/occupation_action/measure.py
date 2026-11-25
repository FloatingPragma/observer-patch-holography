"""Native binary histories, their occupancy counts and the observed action per copy."""

import argparse
from fractions import Fraction
import json
from math import asin, comb, log, pi, sqrt
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.special import expit
from scipy.stats import beta

from model.candidate_laws.interval_echo.candidate import KAPPA, native_law
from model.executable_model import CELL_SEAMS
from ...refinement_cooling.dimer_band.measure import pointer

HERE = Path(__file__).resolve().parent
SHOTS = 262_144
SIZES = (1, 8, 16, 32, 64)
HALF = Fraction(1, 2)
THIRD = Fraction(1, 3)
HISTORICAL = (Fraction(107, 1431), Fraction(5, 508))


def rational_read(value):
    """Canonicalize roundoff only after checking the native pulse probability."""
    result = Fraction(float(value)).limit_denominator(1_000_000)
    if abs(float(result) - value) > 2e-13:
        raise ArithmeticError("read probability outside rational pulse certificate")
    return result


def native_kernel(flip_angles, exchange=True):
    edges = {frozenset(e) for e in CELL_SEAMS}
    if not {frozenset((2, 3)), frozenset((3, 4))} <= edges:
        raise ValueError("apparatus is off the local incidence")
    hopping, contrast = native_law.local_term_parts()
    h = (hopping if exchange else 0) + KAPPA * contrast
    columns = [native_law.SRC.index(pair) for pair in ((2, 0), (0, 2))]
    block = h[np.ix_(columns, columns)]
    np.testing.assert_allclose(np.delete(h[:, columns], columns, axis=0), 0, atol=2e-13)
    # The common diagonal drops out of every read. Check the generator that
    # gives the exact sin-squared half-angle certificate for both sharp starts.
    reduced = block - np.eye(2) * block[0, 0]
    expected = np.array([[0, -0.5], [-0.5, 0]]) if exchange else np.zeros((2, 2))
    np.testing.assert_allclose(reduced, expected, atol=2e-13)
    rows = []
    for start, angle in enumerate(flip_angles):
        amplitudes = expm(-1j * angle * reduced)[:, start]
        row = tuple(rational_read(p) for p in abs(amplitudes) ** 2)
        if sum(row) != 1:
            raise ArithmeticError("contact lost a pointer outcome")
        rows.append(row)
    return tuple(rows)


def count_law(n, m, p0, p1):
    """Convolution of the two initial rail strata with exact rational weights."""

    def row(size, p):
        return [comb(size, k) * p**k * (1 - p) ** (size - k) for k in range(size + 1)]

    law = [Fraction(0)] * (n + 1)
    for i, left in enumerate(row(n - m, p0)):
        for j, right in enumerate(row(m, p1)):
            law[i + j] += left * right
    return law


def log_mgf(x, p, p0, p1):
    def bernoulli(q):
        return np.logaddexp(np.log1p(-q), np.log(q) + p)

    return (1 - x) * bernoulli(p0) + x * bernoulli(p1)


def expected_fraction(x, p, p0, p1):
    def tilted(q):
        return expit(log(q / (1 - q)) + p)

    return (1 - x) * tilted(p0) + x * tilted(p1)


def log_rate(x, y, p0, p1):
    if y == 0:
        return -(1 - x) * log(1 - p0) - x * log(1 - p1)
    if y == 1:
        return -(1 - x) * log(p0) - x * log(p1)
    # The weighted mean lies between its two tilted Bernoulli components.
    target = log(y / (1 - y))
    logits = [log(q / (1 - q)) for q in (p0, p1)]
    lower, upper = target - max(logits), target - min(logits)
    for _ in range(60):
        mid = (lower + upper) / 2
        if expected_fraction(x, mid, p0, p1) < y:
            lower = mid
        else:
            upper = mid
    p = (lower + upper) / 2
    return p * y - log_mgf(x, p, p0, p1)


def chain_bank(rng, initial, kernel):
    # Binomial subdivision of each history prefix keeps every original launch.
    ones = int(rng.binomial(SHOTS, float(initial)))
    bank = [SHOTS - ones, ones]
    expected = [1 - initial, initial]
    for _ in range(4):
        next_bank, next_expected = [], []
        for word, (shots, weight) in enumerate(zip(bank, expected, strict=True)):
            row = kernel[word % 2]
            hits = int(rng.binomial(shots, float(row[1])))
            next_bank.extend((shots - hits, hits))
            next_expected.extend(weight * q for q in row)
        bank, expected = next_bank, next_expected
    return dict(counts_lexicographic=bank, expected=[str(p) for p in expected])


def fraction_bank(rng, n, m, kernel):
    p0, p1 = (r[1] for r in kernel)
    left = rng.binomial(n - m, float(p0), SHOTS)
    right = rng.binomial(m, float(p1), SHOTS)
    counts = np.bincount(left + right, minlength=n + 1).tolist()
    law = count_law(n, m, p0, p1)
    cdf_distance = np.max(
        abs(np.cumsum(counts) / SHOTS - np.cumsum([float(v) for v in law]))
    )
    return dict(
        copies=n,
        initial_ones=m,
        counts=counts,
        cdf_distance=round(float(cdf_distance), 12),
    )


def probability_interval(hits):
    # Each two-sided exact interval has total error .005, for .01 joint error.
    lo = 0 if hits == 0 else beta.ppf(0.0025, hits, SHOTS - hits + 1)
    hi = 1 if hits == SHOTS else beta.ppf(0.9975, hits + 1, SHOTS - hits)
    return lo, hi


def discriminator(counts, p0):
    # Residual g = bilinear corner cost minus the measured per-copy cost.
    n = len(counts) - 1

    def g(y, freq):
        linear = -(1 - y) * log(1 - p0) - y * log(p0)
        return linear + log(freq) / n

    quarters, halves = counts[n // 4], counts[n // 2]
    q_lo, q_hi = probability_interval(quarters)
    h_lo, h_hi = probability_interval(halves)
    raw = g(0.25, quarters / SHOTS) - 0.75 * g(0.5, halves / SHOTS)
    exact = (log(comb(n, n // 4)) - 0.75 * log(comb(n, n // 2))) / n
    return dict(
        measured=round(raw, 12),
        simultaneous_99=[
            round(g(0.25, q_lo) - 0.75 * g(0.5, h_hi), 12),
            round(g(0.25, q_hi) - 0.75 * g(0.5, h_lo), 12),
        ],
        exact_finite=round(exact, 12),
        limiting=round(log(32 / 27) / 4, 12),
        single_curvature=0,
    )


def contact_angle(flip):
    return 2 * asin(sqrt(float(flip)))


def run():
    rng = np.random.default_rng(620452)
    pointer()
    menus = {}
    for name, angles, exchange in (
        ("symmetric", [contact_angle(THIRD)] * 2, True),
        ("held_duration", [pi / 3] * 2, True),
        ("exchange_off", [contact_angle(THIRD)] * 2, False),
    ):
        kernel = native_kernel(angles, exchange)
        fair = native_kernel([pi / 2] * 2, exchange)[0][1]
        banks = [
            fraction_bank(rng, n, m, kernel)
            for n in SIZES
            for m in sorted({0, n // 2, n})
        ]
        menus[name] = dict(
            kernel=[[str(p) for p in row] for row in kernel],
            initial_one=str(fair),
            histories=chain_bank(rng, fair, kernel),
            fraction_banks=banks,
        )
    held = next(
        b
        for b in menus["symmetric"]["fraction_banks"]
        if b["copies"] == 64 and b["initial_ones"] == 0
    )
    control = int(rng.binomial(SHOTS, float(THIRD)))
    asym = native_kernel([contact_angle(q) for q in HISTORICAL])
    initial = native_kernel([contact_angle(Fraction(54356, 61511))] * 2)[0][1]
    return dict(
        shots_per_bank=SHOTS,
        menus=menus,
        held_test=discriminator(held["counts"], float(THIRD)),
        shared_draw_control=dict(
            copies=64, initial_ones=0, counts=[SHOTS - control] + [0] * 63 + [control]
        ),
        imported_chain=dict(
            kernel=[[str(p) for p in r] for r in asym],
            initial_one=str(initial),
            histories=chain_bank(rng, initial, asym),
        ),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    receipt = run()
    path = HERE / "results.json"
    if args.write:
        path.write_text(json.dumps(receipt, indent=2) + "\n")
    elif receipt != json.loads(path.read_text()):
        raise SystemExit("occupation action receipt differs")
    print(json.dumps(receipt["held_test"], indent=2))
