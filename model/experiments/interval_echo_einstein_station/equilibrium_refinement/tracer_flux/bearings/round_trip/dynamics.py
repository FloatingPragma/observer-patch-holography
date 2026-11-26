"""Native free evolution, uniform hitting-time bound and discounted two-leg first arrivals."""

from math import exp, factorial, sin, sqrt
import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.sparse import csr_matrix, eye
from scipy.sparse.csgraph import shortest_path
from .. import measure as bearings
from ... import measure as tracer


def setup():
    twice, basis, rim, outside = tracer.source(8)
    tag = np.array([a for a, _ in basis])
    _, edges = bearings.exchange_certificate(twice, tag)
    distances = shortest_path(
        csr_matrix(edges), directed=False, unweighted=True
    ).astype(int)
    h = np.array([0, 11, 14, 15])[distances]
    for target, row in enumerate(h):
        drift = (edges @ row) - 5 * row
        assert np.all(np.delete(drift, target) == -5)
    return twice, basis, tag, h, edges, outside


def wave(twice, t):
    n = len(twice)
    diag = twice.diagonal()
    centered = csr_matrix(twice) - int((max(diag) + min(diag)) // 2) * eye(
        n, format="csr"
    )
    x = centered * (t / 2)
    radius = float(np.max(np.asarray(abs(x).sum(1))))
    d = 24
    assert radius <= 2
    # Hermitian spectral Taylor remainder on the real interval plus sparse
    # row multiply-add errors propagated by the absolute row/column sum.
    eps = 2.0**-53
    gamma = (2 * np.max(np.diff(x.indptr)) + 6) * eps
    colerr = (
        radius ** (d + 1) / factorial(d + 1)
        + 4 * gamma * radius * exp(radius)
        + 4 * (d + 2) * eps * exp(radius)
    )
    value = np.eye(n, dtype=complex)
    term = np.eye(n)
    phases = (1, -1j, -1, 1j)
    for k in range(1, d + 1):
        term = (x @ term) / k
        value += phases[k % 4] * term
    return value, colerr


def certificate(u, err, tag, h, t, coherent=True):
    rate = 5 * sin(t / 2) ** 2
    n = len(tag)
    bound = 0.0
    sharp = 0.0
    for site in range(12):
        a = u[:, tag == site]
        for target in range(12):
            if target == site:
                continue
            delta = h[target, tag] - h[target, site]
            drift = a.conj().T @ (a * delta[:, None])
            drift += rate * np.eye(a.shape[1])
            sharp = max(sharp, float(max(abs(drift.diagonal()))))
            if coherent:
                bound = max(bound, float(max(abs(drift).sum(1))))
    eps = 2.0**-53
    m = sum(tag == 0)
    op = sqrt(m) * err
    # Absolute-value spectral norm of the n by m isometry is at most sqrt(m).
    # Hermitian dot errors are bounded in spectral norm via |A|^T |weights| |A|.
    rounding = 2 * (8 * n + 32) * eps * m * 15 * (1 + op) ** 2
    remainder = 15 * (2 * op + op**2) + rounding + 32 * m * eps
    # To certify the computed MAXIMUM ROW SUM, the spectral norm remainder can
    # be added after bounding the nominal spectral norm by its row-sum norm.
    return dict(
        conditional_drift_relative_upper=(bound + remainder) / rate,
        sharp_drift_relative_upper=(sharp + remainder) / rate,
        jump_rate=rate,
        column_amplitude_error_upper=err,
        column_probability_error_upper=3 * err + 32 * eps,
    )


def returns(p, tag, lam, one_step_error, hit_bound):
    """Failure to complete a,b,a before a fair geometric deadline; every start retained."""
    z = 1 - lam
    facilities = []
    one = np.zeros((len(tag), 12))
    residual = 0.0

    def certify(table, free):
        # Failure probability u satisfies u=lambda+(1-lambda)P u outside the target.
        nonlocal residual
        delta = table[free] - (lam + z * p[free] @ table)
        y = float(np.max(abs(table)))
        error = (
            float(np.max(abs(delta)))
            + (3 * len(tag) + 12) * np.finfo(float).eps * y
            + 4e-18
        )
        residual = max(residual, hit_bound * (error + one_step_error * y))

    for b in range(12):
        cut = tag == b
        free = ~cut
        system = np.eye(sum(free)) - z * p[np.ix_(free, free)]
        fac = cho_factor(system, check_finite=False)
        one[free, b] = cho_solve(fac, np.repeat(lam, sum(free)), check_finite=False)
        certify(one[:, b], free)
        facilities.append((fac, cut, free))
    first_error = residual
    lows = np.zeros((len(tag), 12))
    for b, (fac, cut, free) in enumerate(facilities):
        table = np.zeros_like(one)
        table[cut] = one[cut]
        rhs = lam + z * p[np.ix_(free, cut)] @ table[cut]
        table[free] = cho_solve(fac, rhs, check_finite=False)
        certify(table, free)
        for a in range(12):
            if a != b:
                lows[tag == a, b] = table[tag == a, a]
    return lows, first_error + residual
