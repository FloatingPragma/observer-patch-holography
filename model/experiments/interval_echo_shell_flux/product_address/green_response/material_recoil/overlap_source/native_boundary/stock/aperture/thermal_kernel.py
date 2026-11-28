"""Thermal diagonal and localized potential increments on the native pair chart."""

import math

import numpy as np
from numpy.polynomial import Chebyshev, Polynomial
from scipy import sparse
from scipy.special import iv

BETA = math.log(8 / 7)
DEGREE = 14


def bounds(h):
    diagonal = h.diagonal()
    radius = np.asarray(abs(h).sum(axis=1)).ravel() - abs(diagonal)
    return float(min(diagonal - radius)), float(max(diagonal + radius))


def baseline(h, beta=BETA, degree=DEGREE, *, energy=False):
    """Diagonal spectral moments at every address, bounded Chebyshev tail."""
    assert degree % 2 == 0
    lo, hi = bounds(h)
    center, scale = (hi + lo) / 2, (hi - lo) / 2
    matrix = (h - sparse.eye(h.shape[0], format="csr") * center) / scale
    t = beta * scale
    d = degree - int(energy)
    coeff = (
        2
        * np.exp(-beta * center)
        * iv(np.arange(d + 1), t)
        * (-1.0) ** np.arange(d + 1)
    )
    coeff[0] /= 2
    polynomial = Chebyshev(coeff).convert(kind=Polynomial)
    if energy:
        polynomial = polynomial * Polynomial([center, scale])
    c = polynomial.coef
    n = h.shape[0]
    out = np.zeros(n)
    # Consecutive Krylov inner products yield two diagonal moments at each layer.
    for start in range(0, n, 128):
        sites = np.arange(start, min(start + 128, n))
        col = sparse.csr_matrix(
            (np.ones(len(sites)), (sites, np.arange(len(sites)))), shape=(n, len(sites))
        )
        value = np.full(len(sites), c[0])
        for j in range(1, degree // 2 + 1):
            nxt = matrix @ col
            value += c[2 * j - 1] * np.asarray(nxt.multiply(col).sum(axis=0)).ravel()
            value += c[2 * j] * np.asarray(nxt.multiply(nxt).sum(axis=0)).ravel()
            col = nxt
        out[sites] = value
    m = d + 1
    tail = (
        2
        * np.exp(-beta * center + t * t / (4 * (m + 1)))
        * (t / 2) ** m
        / math.factorial(m)
        / (1 - t / (2 * (m + 1)))
    )
    # All rescaled row absolute sums are <= 1 by the Gershgorin enclosure.
    rounding = 128 * np.finfo(float).eps * (1 + abs(c).sum()) * (degree + 1) ** 2
    return out, float(tail * (max(abs(lo), abs(hi)) if energy else 1) + rounding)


def columns(matrix, site, order):
    v = np.zeros(matrix.shape[0])
    v[site] = 1
    powers = [v]
    for _ in range(order):
        powers.append(matrix @ powers[-1])
    return np.array(powers)


def single_updates(h, sites, strengths=(32, 33), beta=BETA, order=54, *, energy=False):
    # Telescope each matrix power about the exact one-site diagonal insertion.
    shift = float(max(h.diagonal()))
    base = sparse.eye(h.shape[0], format="csr") * shift - h
    assert min(base.data) >= 0
    scale = float(max(np.asarray(base.sum(axis=1)).ravel()) + max(strengths))
    base = base / scale
    t = beta * scale
    coeff = [math.exp(-beta * shift)]
    for k in range(1, order + 2):
        coeff.append(coeff[-1] * t / k)
    weights = np.array(
        [
            [coeff[i + j + 1] if i + j < order else 0 for j in range(order)]
            for i in range(order)
        ]
    )
    if energy:
        orders = np.add.outer(np.arange(order), np.arange(order)) + 1
        weights = weights * (shift - orders / beta)
    result = []
    for z in sites:
        powers = columns(base, z, order - 1)
        right = weights @ powers
        row = {}
        for u in strengths:
            mat = base.copy()
            mat[z, z] += u / scale
            left = columns(mat, z, order - 1)
            row[u] = u / scale * np.einsum("ij,ij->j", left, right)
        result.append(row)
    tail = coeff[order + 1] / (1 - t / (order + 2))
    rounding = 512 * np.finfo(float).eps * (order + 1) * math.exp(t - beta * shift)
    # Error on each positive increment from positive products and sums is relative,
    # plus the absolute truncation bound for the full nonnegative expansion.
    if energy:
        # Bound the differentiated exponential tail; all retained local heat
        # coefficients are positive. For this chart beta*shift < 1, so every
        # differentiated energy coefficient has the same (negative) sign.
        assert beta * shift < 1
        tail *= shift + (order + 1) / beta / (1 - t / (order + 1))
    return result, dict(tail=tail, rounding_relative=rounding, shift=shift, scale=scale)


def mixed_bounds(h, sites, distance, beta=BETA, depth=33):
    off = h.copy()
    off.setdiag(0)
    maxhop = max(np.asarray(abs(off).sum(axis=1)).ravel())
    maxdiagonal = max(-h.diagonal().min(), max(-h.diagonal()[sites] + depth))
    required = 2 * int(distance)
    # Dyson expansion by number of actual pair hops with diagonal dwell integrated.
    t = beta * maxhop
    top = (
        math.exp(beta * maxdiagonal)
        * t**required
        / math.factorial(required)
        / (1 - t / (required + 1))
    )
    # At the first mark one shortest round trip plus one insertion at each mark
    # is a positive word of degree 2d+2 in the shifted exponential.
    shift = max(h.diagonal())
    bottom = (
        math.exp(-beta * shift)
        * 32**2
        * 0.5**required
        * beta ** (required + 2)
        / math.factorial(required + 2)
    )
    return dict(
        distance=int(distance),
        first_mixed_order=required + 2,
        upper_each_diagonal=float(top),
        lower_at_first_site=bottom,
        inverse_temperature=beta,
        max_hop_rate=float(maxhop),
        max_diagonal_gain=float(maxdiagonal),
    )
