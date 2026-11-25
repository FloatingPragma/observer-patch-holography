"""Outward Clopper–Pearson enclosures from direct binomial tails."""

import numpy as np
from scipy.stats import binom


def binomial_limits(counts, launches, tail):
    values = np.asarray(counts)
    if (
        launches < 1
        or not 0 < tail < 0.5
        or not np.issubdtype(values.dtype, np.integer)
        or np.any(values < 0)
        or np.any(values > launches)
    ):
        raise ValueError("Invalid binomial counts or tail")
    unique, inverse = np.unique(values, return_inverse=True)
    frequency = unique / launches
    # Hoeffding brackets the endpoints; the checks bind the brackets to the
    # actual forward binomial tails, including near the ends of the alphabet.
    radius = np.sqrt(-np.log(tail) / (2 * launches))
    left = np.column_stack((np.maximum(0, frequency - radius), frequency))
    right = np.column_stack((frequency, np.minimum(1, frequency + radius)))
    active_low = unique > 0
    active_high = unique < launches

    def tails(p):
        return (
            binom.sf(unique[active_low] - 1, launches, p[active_low, 0]),
            binom.cdf(unique[active_high], launches, p[active_high, 1]),
        )

    low, high = tails(left)
    assert np.all(low <= tail) and np.all(high >= tail)
    low, high = tails(right)
    assert np.all(low >= tail) and np.all(high <= tail)
    for _ in range(64):
        midpoint = left + (right - left) / 2
        if not np.any((midpoint > left) & (midpoint < right)):
            break
        low, high = tails(midpoint)
        advance = np.zeros_like(left, dtype=bool)
        advance[active_low, 0] = low <= tail
        advance[active_high, 1] = high >= tail
        left = np.where(advance, midpoint, left)
        right = np.where(advance, right, midpoint)
    lower, upper = left[:, 0].copy(), right[:, 1].copy()
    lower[~active_low], upper[~active_high] = 0.0, 1.0
    assert np.all(binom.sf(unique[active_low] - 1, launches, lower[active_low]) <= tail)
    assert np.all(binom.cdf(unique[active_high], launches, upper[active_high]) <= tail)
    # Reconstruct input shape on both NumPy 1.x and 2.x.
    return lower[inverse].reshape(values.shape), upper[inverse].reshape(values.shape)
