"""Coherent finite waiting experiment evaluated with one sparse recurrence."""

import numpy as np
from scipy import sparse
from scipy.special import jv, gammaln

LABELS = 32
STEP = 1 / 256
WAITS = (0, 0.25, 1)


def flights(h, enclosure, seeds, reference, repetitions=4, waits=WAITS):
    lo, hi = enclosure
    mid, rad = (lo + hi) / 2, (hi - lo) / 2
    scaled = (h - sparse.eye(h.shape[0]) * mid) / rad
    weights = np.ones(LABELS) / LABELS
    law = weights
    for _ in range(repetitions - 1):
        law = np.convolve(law, weights)
    delays = STEP * np.arange(len(law))
    n = 1
    x = rad * (max(waits) + delays[-1])
    # |J_k(x)| <= (x/2)^k/k!, by the integral representation.
    while n <= x or np.log(8) + n * np.log(x / 2) - gammaln(n + 1) > np.log(1e-15):
        n += 1
    orders = np.arange(n)
    coeff = []
    for t in waits:
        for selected in (False, True):
            s = t + delays if selected else np.array([t])
            w = law * np.exp(1j * reference * delays) if selected else np.ones(1)
            values = jv(orders[:, None], rad * s[None, :]) @ (np.exp(-1j * mid * s) * w)
            values *= (-1j) ** orders
            values[1:] *= 2
            coeff.append(values)
    coefficients = np.asarray(coeff)
    a = seeds.copy()
    b = scaled @ a
    out = coefficients[:, 0, None, None] * a + coefficients[:, 1, None, None] * b
    for k in range(2, n):
        a, b = b, 2 * (scaled @ b) - a
        for j in range(len(out)):
            out[j] += coefficients[j, k] * b
    return out, n
