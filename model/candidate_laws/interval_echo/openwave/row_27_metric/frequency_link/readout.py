"""Fit modulation frequencies from counts and bound all aliases on the search interval."""

from functools import lru_cache

import numpy as np
from scipy.optimize import minimize_scalar

START = 3.0
T = START + np.arange(1024) / 16
D = np.array([2.0, 2.25])
N = 2**24
GRID = 2**32
F = np.linspace(0.4, 45, 4461)
L = np.log(2 * 1_000_000 / 0.01)


def errors(p):
    # Invert |y-p| <= sqrt(2 p L/N) + 2L/(3N); cover all claimed outcome means.
    a = L / N
    c = 2 * a / 3
    return c + a + np.sqrt(a * a + 2 * a * (p + c)) + 12 / GRID


@lru_cache(maxsize=1)
def moments():
    t = T[:512]
    h = t - t.mean()
    theta = np.outer(F, h)
    c = np.cos(theta)
    c -= c.mean(axis=1, keepdims=True)
    s = np.sin(theta)  # antisymmetric on the centered uniform grid
    return c, s, np.einsum("ij,ij->i", c, c), np.einsum("ij,ij->i", s, s)


def fit_at(f, h, y):
    design = np.column_stack((np.ones(len(h)), np.cos(f * h), np.sin(f * h)))
    coef = np.linalg.lstsq(design, y, rcond=None)[0]
    return coef, float(np.mean((design @ coef - y) ** 2))


def curve(f, coef, t):
    angle = f * (np.asarray(t) - T[:512].mean())
    return coef[0] + coef[1] * np.cos(angle) + coef[2] * np.sin(angle)


def read(counts):
    y = np.asarray(counts, dtype=float) / N
    y0 = y[:512]
    mean = y0.mean()
    centered = y0 - mean
    t = T[:512] - T[:512].mean()
    c, s, cc, ss = moments()
    power = (c @ centered) ** 2 / cc + (s @ centered) ** 2 / ss
    residual = np.sqrt(np.maximum(0, (centered @ centered - power) / 512))
    error = float(np.sqrt(np.mean(errors(y0) ** 2)))
    # For a globally nonnegative sinusoid B>=amp. At f>=0.4 and f<=45
    # the finite prefix mean is at least B-amp/(512*sin(0.4/(2*16))).
    amp_bound = (mean + error) / (1 - 1 / (512 * np.sin(0.4 / 32)))
    radius = error + amp_bound * np.sqrt(np.mean(t * t)) * 0.0050000001 + 1e-14
    ok = np.flatnonzero(residual <= radius)
    intervals = []
    for g in np.split(ok, np.flatnonzero(np.diff(ok) > 1) + 1):
        if len(g):
            intervals.append(
                [
                    max(0.4, float(F[g[0]] - 0.005000001)),
                    min(45, float(F[g[-1]] + 0.005000001)),
                ]
            )
    i = int(np.argmin(residual))
    fit = minimize_scalar(
        lambda f: fit_at(f, t, y0)[1],
        bounds=(F[max(0, i - 1)], F[min(len(F) - 1, i + 1)]),
        method="bounded",
        options={"xatol": 1e-11},
    )
    f = float(fit.x)
    coef, _ = fit_at(f, t, y0)
    tail = curve(f, coef, T[512:]) - y[512:]
    resolved = bool(np.std(y0) > error)
    return dict(
        frequency=f if resolved else None,
        coefficients=coef.tolist(),
        frequency_intervals=intervals,
        amplitude=float(np.hypot(*coef[1:])),
        prefix_error_bound=error,
        constant_excluded=resolved,
        held_rms=float(np.sqrt(np.mean(tail**2))),
        held_max=float(np.max(abs(tail))),
    )


def ratio(remote, local):
    a, b = remote["frequency_intervals"], local["frequency_intervals"]
    if (
        not (remote["constant_excluded"] and local["constant_excluded"])
        or len(a) != 1
        or len(b) != 1
    ):
        return None
    return dict(
        reading=remote["frequency"] / local["frequency"],
        allowed=[a[0][0] / b[0][1], a[0][1] / b[0][0]],
    )
