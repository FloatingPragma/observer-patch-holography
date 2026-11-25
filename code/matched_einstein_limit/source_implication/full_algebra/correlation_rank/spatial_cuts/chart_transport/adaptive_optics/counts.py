"""Finite binary mobile-interference records from a balanced source/standard switch."""

from fractions import Fraction
from math import isqrt

import numpy as np
from mpmath import iv
from scipy.linalg import qr
from scipy import sparse

from matched_einstein_limit.source_implication.directional_stress.translation.counts import (
    selection_norm,
)
from .. import apertures, mounts, transfer

EXPOSURE = 2**148
GRID = 2**192


def choose(rows):
    keys = [(a, i) for a in range(len(rows)) for i in (0, 1, 2, 4, 5, 6)]
    b = transfer.cuts.mid(np.array([rows[a][i] for a, i in keys]))
    norms = np.linalg.norm(b, axis=1)
    choices = []
    for power in (0.0, 1.0):
        p = qr((b / norms[:, None] ** power).T, mode="economic", pivoting=True)[2][:31]
        choices.append([keys[j] for j in sorted(p)])
    assert choices[0] != choices[1]
    return choices


def score_ranges():
    board, route, _, patches = apertures.apparatus()
    out = {}
    for a, patch in enumerate(patches):
        for q in (-1, 1):
            _, ops = mounts.cuts.native.op_table(board, route, patch, q, False)
            for i, name in enumerate(mounts.REGIONAL[:3]):
                op, unit = ops[name]
                own = op[::4, ::4]
                assert (op - sparse.kron(own, sparse.eye(4))).count_nonzero() == 0
                # These three detector operators are built from the blank
                # Hamiltonian and the partial shift, independently of head parity.
                c = Fraction.from_float(selection_norm(op, unit))
                assert c > 0
                out[a, 4 * int(q == 1) + i] = c
    return out


def coords(source):
    # The physical state is the supplied positive binary-rational block matrix
    # normalized to trace one. Keep the normalization and basis contractions directed.
    x = np.array(
        [
            [[iv.mpc(float(v.real), float(v.imag)) for v in row] for row in m]
            for m in source
        ],
        dtype=object,
    )
    trace = sum(x[p, i, i].real for p in range(2) for i in range(4))
    j = mounts.cuts.regional.cert.iso.basis(True)
    m = np.asarray(mounts.full_source.local_basis(), dtype=complex)
    result = []
    for axes in (j, m):
        for block in x:
            for axis in axes:
                result.append(
                    sum(
                        v.real * float(a.real) - v.imag * float(a.imag)
                        for v, a in zip(block.flat, axis.T.flat)
                    )
                    / trace
                )
    result.append(
        sum(x[0, i, i].real - x[1, i, i].real for i in range(4)) / (trace * iv.sqrt(8))
    )
    return np.array(result, dtype=object)


def random_binomial(n, prob, rng):
    """Exact rejection law for the specified positive dyadic success probability.

    A uniform-width envelope on geometric blocks dominates the binomial
    masses on both sides of its mode. Interval comparisons decide acceptance.
    """
    m = ((n + 1) * prob).__floor__()
    width = 1 + isqrt(n)
    if not width < m < n - width:
        raise ValueError("the near-fair detector needs both full shoulder blocks")
    right = Fraction(n - m - width, m + width + 1) * prob / (1 - prob)
    left = Fraction(m - width, n - m + width + 1) * (1 - prob) / prob
    assert width * (1 - right) >= 1 and width * (1 - left) >= 1
    accuracy = iv.dps
    try:
        while True:
            block = rng.randrange(3)
            if block == 2:
                while rng.getrandbits(1):
                    block += 1
            dist = block * width + rng.randrange(width)
            k = m + dist if rng.getrandbits(1) else m - 1 - dist
            if not 0 <= k <= n:
                continue
            bits = 128
            u = rng.getrandbits(bits)
            while True:
                iv.dps = max(100, int(bits * 0.31) + 65)
                ip = iv.mpf(prob.numerator) / prob.denominator
                lg = iv.loggamma
                loga = (
                    lg(m + 1)
                    - lg(k + 1)
                    + lg(n - m + 1)
                    - lg(n - k + 1)
                    + (k - m) * iv.log(ip / (1 - ip))
                    + max(0, block - 1) * iv.log(2)
                )
                p = iv.exp(loga)
                lo = iv.mpf(u) / (1 << bits)
                hi = iv.mpf(u + 1) / (1 << bits)
                if hi.b <= p.a:
                    return k
                if lo.a >= p.b:
                    break
                u = (u << 128) | rng.getrandbits(128)
                bits += 128
    finally:
        iv.dps = accuracy


def acquire(rows, keys, prices, source, rng):
    d = coords(source)
    values = []
    totals = []
    costs = []
    coupling = iv.mpf(0)
    for a, i in keys:
        c = prices[a, i]
        scale = iv.mpf(c.numerator) / c.denominator
        # Fair choice of unknown vs tracial stock with oppositely signed score:
        # the binary success probability is 1/2 + delta_mu/(4 C).
        prob = iv.mpf(".5") + sum(rows[a][i] * d) / (4 * scale)
        mid = (transfer.endpoint(prob, 0) + transfer.endpoint(prob, 1)) / 2
        dy = Fraction(round(mid * GRID), GRID)
        exact = iv.mpf(dy.numerator) / dy.denominator
        coupling += iv.sqrt(EXPOSURE / (2 * exact * (1 - exact))) * abs(prob - exact).b
        k = random_binomial(EXPOSURE, dy, rng)
        values.append(Fraction(2 * c * (2 * k - EXPOSURE), EXPOSURE))
        totals.append(k)
        costs.append(float(c))
    # Bernoulli relative entropy is bounded by chi-square. Pinsker then
    # bounds each count-bank TV, which can be summed over independent settings.
    assert coupling.b < iv.mpf("1e-18"), str(coupling)
    # Hoeffding with the actual score ranges, simultaneous for both endpoints,
    # all 31 channels and the three independently metered sources.
    rad = (
        iv.sqrt(sum(iv.mpf(v) ** 2 for v in costs))
        * 4
        * iv.sqrt(iv.log(2 * 31 * 2 * 3 / iv.mpf(".009")) / (2 * EXPOSURE))
    )
    return values, dict(
        success=totals,
        trials=EXPOSURE,
        score_ranges=costs,
        contrast_l2_error_upper=transfer.upper(rad),
        sampling_law_total_variation_upper=transfer.upper(coupling),
    )
