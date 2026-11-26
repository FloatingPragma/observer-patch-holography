"""Free native four-mobile propagation and a certified projected Chebyshev calculation."""

from fractions import Fraction
from itertools import combinations
from math import factorial
import numpy as np
from scipy.sparse import csr_matrix, eye

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from .. import clock_ball, measure as regional


def hardware(tag_count):
    car = federation(0)
    ports, contacts, _, words = regional.apparatus(car, [0])
    radius = [
        len(path)
        for _, path in sorted(clock_ball.shortest_chains(ports, contacts, 0).items())
    ]
    basis = [
        (w, tags)
        for w in words
        for tags in combinations((a for a in ports if a not in w), tag_count)
    ]
    index = {state: i for i, state in enumerate(basis)}
    law = candidate.native_law
    seam = 2 * law.seam_term(candidate.KAPPA)
    for source, (a, b) in enumerate(law.SRC):
        expected = np.zeros(len(law.SRC))
        expected[source] = 2 * candidate.KAPPA * (a - b) ** 2
        if a != b:
            expected[source] += 1
            expected[law.SRC.index((b, a))] -= 1
        np.testing.assert_array_equal(seam[:, source], expected)
    data, rows, cols = [], [], []
    for i, (w, tags) in enumerate(basis):
        loads = np.zeros(12, dtype=int)
        loads[list(w)] = 2
        loads[list(tags)] = 1
        diag = 2 * candidate.KAPPA * sum(loads**2 * (6 - car.degree[:12]))
        for a, b in contacts:
            left, right = loads[a], loads[b]
            source = law.SRC.index((left, right))
            diag += seam[source, source]
            if left == right:
                continue
            swap = {a: b, b: a}
            state = tuple(
                tuple(sorted(swap.get(p, p) for p in side)) for side in (w, tags)
            )
            j = index[state]
            rows.append(j)
            cols.append(i)
            data.append(seam[law.SRC.index((right, left)), source])
        rows.append(i)
        cols.append(i)
        data.append(diag)
    twice_h = csr_matrix((data, (rows, cols)), shape=(len(basis), len(basis)))
    assert (twice_h != twice_h.T).nnz == 0
    assert (twice_h.data == np.rint(twice_h.data)).all()
    return twice_h, basis, radius


def bessel(z, k):
    """J_k at rational z; alternating-series remainder once terms decrease."""
    term = (z / 2) ** k / factorial(k)
    value = term
    j = 1
    while True:
        ratio = z * z / (4 * j * (j + k))
        term *= -ratio
        if ratio < 1 and abs(term) < Fraction(1, 10**21):
            return float(value)
        value += term
        j += 1


def propagators(twice_h, launches, times):
    """Enclose each probability column in l1, using the exact native Hermitian generator."""
    n = twice_h.shape[0]
    diag = twice_h.diagonal()
    center = int((diag.max() + diag.min()) // 2)
    centered = twice_h - center * eye(n, format="csr")
    radius = int(np.max(np.asarray(abs(centered).sum(axis=1))))
    scale = 1 << radius.bit_length()
    x = centered / scale
    z = [Fraction(t) * scale / 2 for t in times]
    d = int(max(z)) + 1
    while 4 * (max(z) / 2) ** (d + 1) / factorial(d + 1) > Fraction(1, 10**16):
        d += 1
    coeff = [
        [
            bessel(zi, k) * (1, -1j, -1, 1j)[k % 4] * (2 if k else 1)
            for k in range(d + 1)
        ]
        for zi in z
    ]
    prev = np.zeros((n, len(launches)))
    prev[launches, np.arange(len(launches))] = 1
    curr = x @ prev
    outputs = [cs[0] * prev + cs[1] * curr for cs in coeff]
    max_degree = int(np.max(np.diff(x.indptr)))
    # In 2*X*T - Tprev the sparse dot error is bounded in l2 by
    # gamma_(2m) || |X| ||_2 ||T||_2 plus one final rounded subtraction.
    # || |X| ||_2 <= 1, and assuming ||T||_2<=2 gives this per-step budget.
    eps = 2.0**-53
    local_error = (16 * max_degree + 16) * eps
    max_error = local_error * d * (d + 1) / 2
    assert max_error < 1
    for k in range(2, d + 1):
        prev, curr = curr, 2 * (x @ curr) - prev
        for u, cs in zip(outputs, coeff, strict=True):
            u += cs[k] * curr
    # ||U_j(X)|| <= j+1 propagates the local defect of the Chebyshev recurrence.
    # Each Bessel coefficient has truncation <1e-21 and rounding <=2 eps.
    # Floating scalar-vector products and accumulation of d+1 weighted terms.
    weight = max(sum(map(abs, cs)) for cs in coeff)
    amp_error = (
        1e-16 + weight * max_error + (d + 1) * 8 * eps + 8 * (d + 1) * eps * weight
    )
    probability_error = 2 * (2 * amp_error + amp_error**2 + 32 * eps)
    assert probability_error < 1e-7
    return [abs(u) ** 2 for u in outputs], dict(
        terms=d + 1,
        centered_doubled_scale=scale,
        max_row_terms=max_degree,
        column_probability_error_upper=float(np.nextafter(probability_error, np.inf)),
    )
