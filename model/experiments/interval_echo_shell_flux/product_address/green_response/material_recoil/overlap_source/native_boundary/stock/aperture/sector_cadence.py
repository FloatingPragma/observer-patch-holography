"""Certify the sign at every passive detection probability on the return process."""

import math

import numpy as np


def mixing(clock, n):
    # Lower-bound every entry of a power of the first-return kernel.
    # The exact dyadic products fit int64 since every lower row sum
    # is at most scale. Flooring preserves the entrywise lower bound.
    scale = 2**29
    first, error = clock.advance(np.eye(n, dtype=np.int64) * clock.grid)
    slop = np.ceil(error * clock.grid + 1).astype(np.int64)
    L = np.maximum(0, (first - slop) // (clock.grid // scale))
    assert np.all(L.sum(1) <= scale)
    power = 1
    while power < 256:
        L = (L @ L) // scale
        power *= 2
    overlap = int(L.min(0).sum())
    assert 0 < overlap <= scale
    return power, 1 - overlap / scale


def evaluate(qlo, qhi, low, high, remainder):
    # Interval Horner with directed binary64 rounding at each operation.
    lo, hi = -remainder, remainder
    for a, b in zip(low[::-1], high[::-1], strict=True):
        lo = math.nextafter(
            math.nextafter((qlo if lo >= 0 else qhi) * lo, -math.inf) + a, -math.inf
        )
        hi = math.nextafter(
            math.nextafter((qhi if hi >= 0 else qlo) * hi, math.inf) + b, math.inf
        )
    return lo, hi


def all_rates(series, tail, mixing_data):
    # q=1-p and Delta(p)/p=sum q^(k-1) Delta_k. Removing the factor
    # p permits a signed enclosure even for arbitrarily sparse reads.
    length = len(series) - 1
    m, eta = mixing_data
    rest = math.nextafter(tail * m / (1 - eta) + 1e-11, math.inf)
    low, high = series[1:].T

    def func(a, b):
        return evaluate(a, b, low, high, rest)

    left, right = 0.0, 1.0
    assert func(0, 0)[0] > 0

    def cover(lo, hi, sign):
        work = [(lo, hi)]
        accepted = []
        while work:
            a, b = work.pop()
            y = func(a, b)
            if y[0] > 0 if sign > 0 else y[1] < 0:
                accepted.append([a, b, y[0] if sign > 0 else y[1]])
            else:
                assert b - a > 2**-30
                mid = (a + b) / 2
                work.extend(((a, mid), (mid, b)))
        return accepted

    # Bound the response per unit reporting rate over bins fixed by the
    # clock alone. The sky calculation can look up a measured filling at
    # full precision and keep the remaining interpolation uncertainty.
    edges = sorted(
        {0.0, 1.0, *(2.0**-j * (1 + i / 8) for j in range(1, 25) for i in range(9))}
    )
    bins = [
        dict(rate=[a, b], per_rate=list(func(1 - b, 1 - a)))
        for a, b in zip(edges[:-1], edges[1:], strict=True)
    ]
    common = dict(
        power=m,
        oscillation_contraction_upper=eta,
        undiscounted_tail_bound=rest,
        sparse_response_per_unit_rate=list(func(1, 1)),
        rate_bins=bins,
    )
    if func(1, 1)[0] > 0:
        positive = cover(0, 1, 1)
        return dict(
            **common,
            positive_for_every_positive_rate=True,
            uniform_lower_per_unit_rate=min(row[2] for row in positive),
            positive_intervals=positive,
        )
    assert func(1, 1)[1] < 0
    # Pick the bracket numerically and then prove signs over both entire
    # complementary intervals and strict derivative sign inside the bracket.
    while right - left > 1 / 2**16:
        mid = (right + left) / 2
        lo, hi = func(mid, mid)
        if lo > 0:
            left = mid
        elif hi < 0:
            right = mid
        else:
            break
    pad = 2**-14
    left = max(0.0, left - pad)
    right = min(1.0, right + pad)
    pos = cover(0.0, left, 1)
    neg = cover(right, 1.0, -1)
    # Each full block has at most m coefficients and drops by eta.
    # After removing the leading powers the derivative tail has an
    # extra factor of the age, bounded by summing a geometric moment.
    der_tail = tail * m * ((length + m) / (1 - eta) + m * eta / (1 - eta) ** 2) + 1e-9
    factors = np.arange(1, length)
    der = evaluate(
        left, right, low[1:] * factors - 1e-10, high[1:] * factors + 1e-10, der_tail
    )
    assert der[1] < 0
    return dict(
        positive_rate_above=1 - left,
        negative_rate_below=1 - right,
        root_bracket=[1 - right, 1 - left],
        direction="positive above, negative below",
        unique_root=True,
        derivative_in_miss_probability=list(der),
        certified_blocks=[len(pos), len(neg)],
        **common,
    )
