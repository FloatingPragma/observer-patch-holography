"""Dyadic full-H evolution and local energy projectors with rational error bounds."""

from fractions import Fraction as Q
from math import factorial, isqrt
import numpy as np


def certified_transfer(twice_h, exposure):
    """Enclose the 64-unit transfer on the nominal four-word apparatus."""
    assert exposure == 64 and twice_h.shape == (4, 4)
    # A scalar shift cancels in probabilities. Work in integer half units.
    scale = 2**75
    offset = (int(twice_h.diagonal().min()) + int(twice_h.diagonal().max())) // 2
    matrix = twice_h.astype(object) - offset * np.eye(4, dtype=int).astype(object)
    norm = max(sum(abs(int(z)) for z in row) for row in matrix)
    divisions = 512
    delay = Q(exposure, divisions)
    order = 64
    real = np.zeros((4, 4), dtype=object)
    imag = real.copy()
    power = np.eye(4, dtype=int).astype(object)
    for k in range(order + 1):
        coefficient = (-delay / 2) ** k / factorial(k)
        if k % 4 in (0, 1):
            (real if k % 4 == 0 else imag)[:] += coefficient * power
        else:
            (real if k % 4 == 2 else imag)[:] -= coefficient * power
        power = power @ matrix

    def scaled(z):
        return np.array([round(x * scale) for x in z.flat], dtype=object).reshape(4, 4)

    re, im = scaled(real), scaled(imag)
    # Frobenius norm of complex entry rounding is <4/scale.
    tail = (norm * delay / 2) ** (order + 1) / factorial(order + 1)
    err = tail + Q(4, scale)
    for _ in range(divisions.bit_length() - 1):
        cre = re @ re - im @ im
        cim = re @ im + im @ re
        re = (cre + scale // 2) // scale
        im = (cim + scale // 2) // scale
        err = err * (2 + err) + Q(4, scale)
        # Round the conservative operator-norm bound upwards at each step.
        err = Q(-(-err * scale // 1), scale)
    lower = Q(isqrt(1601 * scale**2), scale)
    b = np.array([[-40, -1], [-1, 40]], dtype=object)
    ident = np.eye(2, dtype=int).astype(object)
    high = (ident + b * (1 / lower)) * Q(1, 2)
    low = ident - high
    inp = np.kron(high, low)
    out = np.kron(low, high)
    approx = np.trace(out @ (re @ inp @ re.T + im @ inp @ im.T)) / scale**2
    delta = Q(1, 2 * scale) / lower
    joint = delta * (2 + delta)
    # The true four factors have operator norm one. Bound the change of
    # their trace by four times the telescoping product bound.
    bound = 4 * ((1 + max(err, joint)) ** 4 - 1)
    ticks = 10**12
    return {
        "exposure": exposure,
        "bits": 75,
        "substeps": divisions,
        "short_taylor_order": order,
        "probability_tick_denominator": ticks,
        "joint_probability_interval_ticks": [
            int((approx - bound) * ticks),
            int(-(-(approx + bound) * ticks // 1)),
        ],
    }
