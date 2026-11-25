"""A rational negative witness for the exact t=16 success map."""

from fractions import Fraction
from math import factorial

import numpy as np

from .measure import build


def witness():
    _, ix, terms = build(True)
    raw = sum(terms).astype(object)
    n = len(raw)
    grid = 2**72
    order = 20
    denominator = factorial(order) * 256**order
    power = np.eye(n, dtype=int).astype(object)
    re = np.zeros_like(power)
    im = np.zeros_like(power)
    for j in range(order + 1):
        factor = denominator // (factorial(j) * 256**j)
        re += (1, 0, -1, 0)[j % 4] * factor * power
        im += (0, -1, 0, 1)[j % 4] * factor * power
        power = raw @ power
    # exp(-i H/32) is expanded as exp(-i raw/256). Its integral Taylor
    # remainder uses unitary norm one along the real evolution parameter.
    bound = Fraction(max(sum(abs(x) for x in row) for row in raw), 256)
    error = bound ** (order + 1) / factorial(order + 1) + Fraction(2 * n, grid)
    re = re * grid // denominator
    im = im * grid // denominator
    for _ in range(9):
        re, im = (re @ re - im @ im) // grid, (re @ im + im @ re) // grid
        # Spectral norm product error. Each rounded real or imaginary
        # matrix entry differs from its exact product by less than 1/grid.
        assert error < 1
        error = 3 * error + Fraction(2 * n, grid)
    assert error < Fraction(1, 10**12)
    source = [ix[0, 0], ix[0, 7]]
    cr = np.zeros((4, 4), dtype=object)
    ci = np.zeros_like(cr)
    for a in (0, 7):
        rows = [ix[3, a + 8 * b] for b in (0, 7)]
        r = re[np.ix_(rows, source)].T.ravel()
        i = im[np.ix_(rows, source)].T.ravel()
        cr += np.outer(r, r) + np.outer(i, i)
        ci += np.outer(i, r) - np.outer(r, i)

    def transpose_reference(x):
        return x.reshape(2, 2, 2, 2).transpose(2, 1, 0, 3).reshape(4, 4)

    cr, ci = transpose_reference(cr), transpose_reference(ci)
    vr = np.array([255, -437, 572, -139], dtype=object)
    vi = np.array([0, -460, -89, 423], dtype=object)
    numerator = vr @ (cr @ vr - ci @ vi) + vi @ (cr @ vi + ci @ vr)
    expectation = Fraction(numerator, 2 * grid**2 * (vr @ vr + vi @ vi))
    # Projection and partial trace contract trace norm on Hermitian errors.
    # Reference transposition preserves Frobenius norm, bounded by trace norm.
    upper = expectation + 2 * error + error**2
    assert upper < Fraction(-81, 10000)
    return dict(
        time=16,
        witness_real=vr.tolist(),
        witness_imag=vi.tolist(),
        scale_bits=72,
        unitary_error_upper="1e-12",
        normalized_witness_on_original_starts_upper="-81/10000",
    )
