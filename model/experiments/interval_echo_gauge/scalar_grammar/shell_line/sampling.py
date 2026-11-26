"""Replay the fixed PCG64 multinomial draws with a single-rounded binomial mode.

The BTPE proposal follows NumPy 2.3.5, numpy/random/src/distributions/distributions.c
(random_binomial_btpe), under its 3-clause BSD license in NUMPY-LICENSE.txt.
"""

from math import floor, log, sqrt

import numpy as np


def _btpe(rng, n, r):
    """Draw Binomial(n, r) for n*r > 30 and 0 < r <= 1/2."""
    q = 1 - r
    # (n+1) is an exactly represented integer throughout this experiment.
    # The product has the same single rounding as a fused n*r+r.
    center = (n + 1) * r
    mode = floor(center)
    mid = mode + 0.5
    width = floor(2.195 * sqrt(n * r * q) - 4.6 * q) + 0.5
    left, right = mid - width, mid + width
    height = 0.134 + 20.5 / (15.3 + mode)
    a = (center - left) / (center - left * r)
    left_rate = a * (1 + a / 2)
    a = (right - center) / (right * q)
    right_rate = a * (1 + a / 2)
    panel = width * (1 + 2 * height)
    left_tail = panel + height / left_rate
    total_area = left_tail + height / right_rate
    odds = r / q
    while True:
        u = rng.random() * total_area
        v = rng.random()
        if u <= width:
            return floor(mid - width * v + u)
        if u <= panel:
            x = left + (u - width) / height
            v = v * height + 1 - abs(mode - x + 0.5) / width
            if v > 1:
                continue
            y = floor(x)
        elif v == 0:
            continue
        elif u <= left_tail:
            y = floor(left + log(v) / left_rate)
            v *= (u - panel) * left_rate
        else:
            y = floor(right - log(v) / right_rate)
            v *= (u - left_tail) * right_rate
        if y < 0 or y > n:
            continue
        # BTPE acceptance is Pr(X=y)/Pr(X=mode). Evaluate this directly
        # by adjacent-binomial ratios for either tail on our finite banks.
        ratio = 1.0
        for k in range(min(y, mode) + 1, max(y, mode) + 1):
            step = (odds * (n + 1)) / k - odds
            ratio = ratio * step if y > mode else ratio / step
        if v <= ratio:
            return y


def multinomial(rng, trials, weights):
    """Resolve bins in historical C order, retaining empty and terminal bins."""
    remaining_weight = int(weights.sum())
    remaining_trials = trials
    counts = np.zeros_like(weights)
    for j, value in enumerate(weights[:-1]):
        if remaining_trials == 0:
            break
        weight = int(value)
        if weight == 0:
            continue
        p = weight / remaining_weight
        r = min(p, 1 - p)
        if remaining_trials * r <= 30:
            # NumPy uses inversion here, including the historical draw for p=1.
            chosen = rng.binomial(remaining_trials, p)
        else:
            chosen = _btpe(rng, remaining_trials, r)
            if p > 0.5:
                chosen = remaining_trials - chosen
        counts[j] = chosen
        remaining_trials -= chosen
        remaining_weight -= weight
    counts[-1] = remaining_trials
    return counts
