"""Binomial count intervals checked against the forward tails, up to 2^52 launches."""

import numpy as np
from scipy.stats import beta, binom


def clopper_pearson(successes, launches, tail):
    """Equal-tail endpoints. Return the conservative side of each numerical bracket."""
    k = np.asarray(successes)
    if not (0 < launches <= 2**52 and launches == int(launches) and 0 < tail < 0.5):
        raise ValueError("Invalid launch count or one-sided tail")
    if np.any((k < 0) | (k > launches) | (k != np.floor(k))):
        raise ValueError("Invalid number of successes")
    k = k.astype(np.int64)
    center = k / launches
    result = []
    for upper in (False, True):
        if upper:
            estimate = beta.isf(tail, k + 1, np.maximum(launches - k, 1))
            forward = binom.cdf
            tail_count = k
            edge = 1.0
        else:
            estimate = beta.ppf(tail, np.maximum(k, 1), launches - k + 1)
            forward = binom.sf
            tail_count = k - 1
            edge = 0.0
        boundary = k == (launches if upper else 0)
        # A rounded quantile is only a starting guess. Even a finite, correctly
        # ordered guess can have the wrong tail by many orders of magnitude.
        estimate = np.where(boundary, edge, estimate)
        step = 1e-8 * abs(estimate - center)
        outer = np.clip(np.nextafter(estimate + (step if upper else -step), edge), 0, 1)
        tail_value = forward(tail_count, launches, outer)
        unsafe = ~boundary & (
            ~np.isfinite(tail_value)
            | (tail_value > tail)
            | (tail_value < 0.999 * tail)
            | ((outer < center) if upper else (outer > center))
        )
        out = np.array(outer, copy=True)
        if np.any(unsafe):
            counts = tail_count[unsafe]
            inside = center[unsafe]
            outside = np.full_like(inside, edge)
            for _ in range(80):
                mid = (inside + outside) / 2
                prob = forward(counts, launches, mid)
                if not np.all(np.isfinite(prob)):
                    raise ArithmeticError("Non-finite forward binomial tail")
                choose = prob <= tail
                outside = np.where(choose, mid, outside)
                inside = np.where(choose, inside, mid)
            out[unsafe] = outside
        result.append(out)
    return tuple(result)
