from math import factorial

import numpy as np
from scipy.linalg import expm
from scipy.sparse import csr_matrix

from .measure import ORDER, numerators, weights


def test_exact_numerator_and_two_port_order():
    # A nonuniform small graph gives different returns. Both transfer orders agree.
    adjacency = np.array([[0, 1, 1, 0], [1, 0, 1, 0], [1, 1, 0, 1], [0, 0, 1, 0]])
    degree = adjacency.sum(axis=1)
    times = (1, 4, 16)
    data = numerators(csr_matrix(adjacency), degree, [0, 3], [weights(t) for t in times])
    denom = (2**ORDER*factorial(ORDER))**2
    for t, row in zip(times, data, strict=True):
        h = (np.diag(degree)-adjacency)/2
        U = abs(expm(-1j*t*h))**2
        np.testing.assert_allclose([n/denom for n in row], [U[0, 0], U[3, 0], U[3, 3]],
                                   rtol=0, atol=8e-15)


def test_delay_certificate_accepts_a_common_nonlinear_clock():
    from fractions import Fraction as F
    from .common_delay import time_identity

    base = np.array([[F(1), F(0), F(1)], [F(-2), F(1), F(-1)], [F(4), F(5), F(6)]],
                    dtype=object)
    observed = base.copy()
    observed[2] += F(1, 7)*base[1]  # Substitute y = x + x²/7 in each channel.
    certificate = time_identity(base, observed)
    assert certificate["implied_squared_time_shift_coefficients"] == ["1/7"]*3
    assert certificate["residual_after_matching_near_return"] == ["0"]*3


def test_delay_coefficients_against_exact_two_port_oscillation():
    from fractions import Fraction as F
    from .common_delay import probability_polynomials

    adjacency = csr_matrix(np.array([[0, 1], [1, 0]]))
    coeffs = probability_polynomials(adjacency, np.array([1, 1]), [0, 1])
    # Return=(1+cos t)/2; transfer=(1-cos t)/2.
    for k, row in enumerate(coeffs[1:], start=1):
        actual = [F(x, 2**(2*k)*factorial(2*k)) for x in row]
        p = F((-1)**k, 2*factorial(2*k))
        assert actual == [p, -p, p]
