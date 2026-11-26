"""Uniform cut-spectrum bound on the independent contact-meter sensitivity box."""

from math import asin, pi, sin, sqrt

import numpy as np


def cut_enclosure(target, spacing, spacing_error, length=4096):
    """Bound the response at one fixed frequency in every pure-cut preparation."""
    m = 10  # nine active dimers, ten bonds
    js = np.linspace(7.996, 8.004, 33)
    hs = np.linspace(0.4997, 0.5003, 9)
    lower = np.full((2 * m, 2 * m), np.inf)
    upper = -lower
    for j in js:
        for h in hs:
            c = 2 * (j * np.eye(m) - h * np.eye(m, k=-1))
            e = np.linalg.svd(c, compute_uv=False)
            e = np.r_[e, -e]
            gaps = e[:, None] - e
            lower = np.minimum(lower, gaps)
            upper = np.maximum(upper, gaps)
    # Singular values are 2(|dJ|+|dh|)-Lipschitz. Each gap uses two.
    padding = 2 * ((js[1] - js[0]) + (hs[1] - hs[0])) + 1e-11
    endpoints = [
        f * t
        for f in (lower - padding, upper + padding)
        for t in (spacing - spacing_error, spacing + spacing_error)
    ]
    lo = target / 8 - np.max(endpoints, axis=0)
    hi = target / 8 - np.min(endpoints, axis=0)
    mid = (lo + hi) / 2
    half = (hi - lo) / 2

    def sin_floor(shift):
        # |sin| is 1-Lipschitz. Half-angle inputs span half/2 either way.
        return np.maximum(0, abs(np.sin((mid + shift) / 2)) - half / 2 - 1e-14)

    d = 2 * pi / length
    denominator = length * sin_floor(0) * sin_floor(-d) * sin_floor(d)
    numerator = sin(d / 2) ** 2 * np.minimum(
        1, abs(np.sin(length * mid / 2)) + length * half / 2 + 1e-12
    )
    bounds = np.ones_like(denominator)
    positive = denominator > 0
    bounds[positive] = np.minimum(
        1, numerator[positive] / denominator[positive] + 1e-12
    )
    min_k = (pi - asin(max(hs) / min(js))) / (m + 1)
    min_sine_norm = m / 2 - 1 / (2 * sin(min_k))
    max_entry_squared = 1 / (2 * min_sine_norm)
    return dict(
        contact_J=list(map(float, js[[0, -1]])),
        exchange_h=list(map(float, hs[[0, -1]])),
        grid=[len(js), len(hs)],
        mode_entry_squared_upper=float(max_entry_squared),
        frequency_padding=float(padding),
        cut_only_ceiling=float(
            sqrt(2 * m) * max_entry_squared * np.linalg.norm(bounds) + 1e-10
        ),
    )
