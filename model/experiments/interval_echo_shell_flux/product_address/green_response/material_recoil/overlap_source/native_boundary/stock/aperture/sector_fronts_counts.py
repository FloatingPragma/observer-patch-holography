"""First and second moments of the native register type counts at every reporting rate."""

from fractions import Fraction as Q

import numpy as np

from . import peer_population, sector_cadence
from .sector_retention import ReturnClock


def moments(bank):
    pairs = np.array(peer_population.fibers.BLANK_PAIRS)
    left, right = pairs.T
    width = len(bank[0][1])
    direction = np.zeros((width, len(pairs)), dtype=np.int64)
    direction[left, np.arange(len(left))] = -1
    direction[right, np.arange(len(right))] = 1
    tri = np.triu_indices(width)
    first = np.zeros(width, dtype=object)
    second = np.zeros(len(tri[0]), dtype=object)
    total_weight = sum(weight for weight, _ in bank)
    for weight, counts in bank:
        # The three family counts of blanks are a multivariate
        # hypergeometric allocation of B blanks in T slots.
        base = counts.copy()
        slots = base[left] + base[right]
        blanks = int(base[right].sum())
        total = int(slots.sum())
        den = max(total, 1)
        base[left], base[right] = slots, 0
        raw = den * base + blanks * (direction @ slots)
        counts_cov = (
            blanks
            * (total - blanks)
            * (
                direction
                @ (total * np.diag(slots) - np.outer(slots, slots))
                @ direction.T
            )
        )
        h = max(total - 1, 1)
        moment = (h * np.outer(raw, raw) + counts_cov)[tri]
        for vec, div, target in ((raw, den, first), (moment, h * den * den, second)):
            whole, rem = divmod(weight, div)
            obj = vec.astype(object)
            tails = rem * obj
            assert np.all(tails % div == 0)
            target += whole * obj + tails // div
    mean = np.array([Q(i, total_weight) for i in first], dtype=object)
    M = np.empty((width, width), dtype=object)
    M[tri] = [Q(i, total_weight) for i in second]
    M[(tri[1], tri[0])] = M[tri]
    cov = M - np.outer(mean, mean)
    assert np.all(cov.sum(0) == 0) and np.all(cov.diagonal() >= 0)
    return mean, cov


def group_readings(w, cat):
    # Last categories contain exactly one protected word each.
    return np.concatenate(
        ([w[cat == i].mean(0) for i in range(cat.max() + 1)], w[cat < 0])
    )


def cadence(subset, mats, flags, cat, counts, cov, size, bins):
    clock = ReturnClock((mats[2] @ mats[1] @ mats[0]).T.tocsr(), subset)
    power, eta = sector_cadence.mixing(clock, len(flags))
    quota = size // len(flags)
    weights = np.r_[np.array(counts, float), np.full(sum(cat < 0), quota)] / size
    gamma = np.array(cov, float) / size
    k = len(counts)
    # Center BEFORE summing the age series. In the limit of sparse
    # reporting the independent clocks forget the starting words.
    center = flags.mean(0)
    ticks = flags * clock.grid
    means = []
    errs = []
    error = np.zeros(3)
    for _ in range(2048):
        ticks, step = clock.advance(ticks)
        error += step
        means.append(group_readings(ticks / clock.grid - center, cat))
        errs.append(error.copy() + 2e-12)
    values = np.array(means)
    tail = (np.ptp(ticks / clock.grid, axis=0) + 2 * error + 1e-10) * power / (1 - eta)
    # Distinct-member covariance after independent reporting is p^2 X'MY.
    # The constant terms cancel because all rows and columns of M sum to zero.
    M = np.outer(weights, weights) - np.diag(weights)
    M[:k, :k] += gamma
    norm = float(abs(M).sum()) + 1e-10
    X, Y = values[:, :, 0], values[:, :, 1]
    XM = X @ M
    coeff = sum(np.convolve(col, other) for col, other in zip(XM.T, Y.T, strict=True))
    ex, ey, _ = np.asarray(errs).T
    A, B = abs(X).max(1), abs(Y).max(1)
    # The absolute entry sum of M bounds each bilinear coefficient. The
    # roundoff guard exceeds the dot, convolution and summation error for
    # these arrays (all dimensions <= 2048 and binary64 unit roundoff < 2^-52).
    guard = norm * (
        np.convolve(ex, B + ey)
        + np.convolve(A, ey)
        + 1e-10 * (1 + np.convolve(A + ex, B + ey))
    )
    out = []
    for row in bins:
        a, b = row["rate"]
        qhi = 1 - a
        cross = np.array(
            sector_cadence.evaluate(1 - b, qhi, coeff - guard, coeff + guard, 0.0)
        )
        tails = qhi ** len(values) * tail + 1e-13
        sX, sY = (sum(A + ex), sum(B + ey))
        if a > 0:
            sX, sY = min(sX, max(A + ex) / a), min(sY, max(B + ey) / a)
        omitted = norm * (tails[0] * sY + tails[1] * sX + tails[0] * tails[1]) + 1e-10
        cross += np.array([-omitted, omitted])
        gain = np.multiply.outer(row["rate"], cross)
        lo, hi = row["covariance_per_rate"]
        out.append(
            {
                "rate": [a, b],
                "per_rate": [lo + gain.min() - 1e-11, hi + gain.max() + 1e-11],
                "distinct_member_gain_per_squared_rate": cross.tolist(),
            }
        )
    return {
        "bins": out,
        "age_terms": len(values),
        "power": power,
        "contraction": eta,
        "tail_per_unit_reporting_rate": tail.tolist(),
        "absolute_count_covariance_sum_per_member": float(abs(gamma).sum()),
        "minimum_collective_gain_per_positive_rate": min(z["per_rate"][0] for z in out),
        "minimum_distinct_member_gain_per_squared_rate": min(
            z["distinct_member_gain_per_squared_rate"][0] for z in out
        ),
    }
