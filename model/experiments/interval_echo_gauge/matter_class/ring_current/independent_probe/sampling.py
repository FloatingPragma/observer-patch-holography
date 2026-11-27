"""The enrolled conditional-binomial sampler with stable boundary comparisons."""

import math

import numpy as np
from scipy.stats import binom

from model.candidate_laws.interval_echo.openwave.row_30_collective.committed_heads.during_flight.array.measure import (
    GRID,
)


def resolving_cdf(k, n, w, mass):
    """Sum positive PMF terms about the mode using exact integer ratio differences."""
    mode = (n + 1) * w // mass
    radius = math.ceil(math.sqrt(n * math.log(2**91) / 2)) + 2
    lo = max(0, mode - radius)
    hi = min(n, mode + radius)
    if k < lo or k >= hi:
        raise ValueError("Quantile outside the verified central summation window")
    down = np.arange(mode, lo, -1, dtype=np.int64)
    up = np.arange(mode + 1, hi + 1, dtype=np.int64)
    left = np.cumprod(1 + (down * mass - (n + 1) * w) / ((n - down + 1) * w))
    right = np.cumprod(1 + ((n + 1) * w - up * mass) / (up * (mass - w)))
    total = 1 + left.sum() + right.sum()
    if k < mode:
        return left[mode - k - 1 :].sum() / total
    return (1 + left.sum() + right[: k - mode].sum()) / total


def quantiles(uniforms, n, w, mass):
    p = w / mass
    guessed = binom.ppf(uniforms, n, p).astype(np.int32)
    # At hundreds of millions of starts Boost's CDF/PPF can straddle an
    # integer threshold differently on two platforms. Resolve nearby CDF
    # endpoints without incomplete-beta evaluations or log-gamma cancellation.
    for shift in (-1, 0):
        k = guessed + shift
        check = (k >= 0) & (k < n) & (abs(binom.cdf(k, n, p) - uniforms) < 1e-7)
        for ix in np.flatnonzero(check):
            level = resolving_cdf(int(k[ix]), int(n[ix]), int(w[ix]), int(mass[ix]))
            if abs(level - uniforms[ix]) < 2e-12:
                raise ArithmeticError("Boundary too close for double-precision CDF")
            guessed[ix] = k[ix] + int(level < uniforms[ix])
    return guessed


def integer_draws(rng, pmf, shots):
    """Same grid, PCG64 words and conditional quantiles as the enrolled utility."""
    weights = np.rint(pmf * GRID).astype(np.int64)
    largest = np.argmax(weights, axis=-1)[..., None]
    remaining_weight = GRID - weights.sum(axis=-1, keepdims=True)
    np.put_along_axis(
        weights,
        largest,
        np.take_along_axis(weights, largest, axis=-1) + remaining_weight,
        axis=-1,
    )
    if np.min(weights) < 0:
        raise ArithmeticError("Signed sampling weight")
    out = np.empty(weights.shape, dtype=np.int32)
    remaining = np.full(weights.shape[:-1], shots, dtype=np.int32)
    mass = np.full(weights.shape[:-1], GRID, dtype=np.int64)
    for i in range(weights.shape[-1] - 1):
        w = weights[..., i]
        part = np.where(w == mass, remaining, 0)
        use = (w > 0) & (w < mass) & (remaining > 0)
        uniforms = (rng.integers(0, 1 << 52, size=np.count_nonzero(use)) + 0.5) / (
            1 << 52
        )
        if uniforms.size:
            part[use] = quantiles(uniforms, remaining[use], w[use], mass[use])
        out[..., i] = part
        remaining -= part
        mass -= w
    out[..., -1] = remaining
    return out
