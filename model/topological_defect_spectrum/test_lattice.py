from fractions import Fraction as F
from itertools import product

import mpmath as mp
import numpy as np
import pytest
from scipy.linalg import eigvalsh
import sympy as s

from . import lattice


def test_full_original_weyl_box_phase_map_and_square():
    sites = tuple(product(range(-1, 2), repeat=3))
    phase = s.diag(*(lattice.phase(site) for site in sites for _ in (0, 1)))
    original = lattice.source_box(1, 2, 3, F(1, 2))
    transformed = lattice.transformed_box(1, 2, 3, F(1, 2))
    assert (phase.H*original*phase-transformed).applyfunc(s.simplify) == s.zeros(54)
    one = lattice.first_order(1, 2, 3, F(1, 2))
    q = one.H*one
    expected = sum((s.kronecker_product(*(q if j == axis else s.eye(3) for j in range(3)), s.eye(2))
                    for axis in range(3)), s.zeros(54))
    assert (transformed*transformed-expected).applyfunc(s.expand) == s.zeros(54)
    # Positivity of the three summands makes their joint kernel the tensor
    # product of the one-dimensional kernels, with both original spin rows.
    assert one.rank() == 2
    scalar_zero = one.nullspace()[0]
    kernel = phase*s.kronecker_product(scalar_zero, scalar_zero, scalar_zero, s.eye(2))
    assert kernel.rank() == 2
    assert (original*kernel).applyfunc(s.expand) == s.zeros(54, 2)


def test_jacobi_complete_finite_spectrum_matches_three_dimensional_source():
    one = np.asarray(lattice.jacobi(1)).astype(float)
    values = eigvalsh(one)
    np.testing.assert_allclose(values, [-np.sqrt(1.5), 0, np.sqrt(1.5)], atol=1e-15, rtol=0)
    predicted = sorted(sign*np.sqrt(sum(x*x for x in row))
                       for row in product(values, repeat=3) for sign in (-1, 1))
    actual = eigvalsh(np.asarray(lattice.source_box(1)).astype(complex))
    np.testing.assert_allclose(actual, predicted, atol=7e-15, rtol=0)
    for radius in (1, 2, 4):
        matrix = lattice.jacobi(radius)
        reverse = s.eye(matrix.rows)[:, list(reversed(range(matrix.rows)))]
        assert reverse*matrix*reverse == -matrix
        assert matrix.rank() == matrix.rows-1


def test_positive_series_certifies_bessel_normalization_and_recurrence():
    with mp.workdps(250):
        def number(value):
            return mp.mpf(value.numerator)/value.denominator
        for z in (F(1, 3), F(1), F(2)):
            for n in range(9):
                lo, hi = lattice.bessel_interval(n, z)
                value = mp.besseli(n, number(z))
                assert number(lo) <= value <= number(hi)
            partial = sum(mp.besseli(abs(n), number(z))**2 for n in range(-30, 31))
            assert abs(partial/mp.besseli(0, 2*number(z))-1) < mp.mpf('1e-65')
        for site in ((0, 0, 0), (1, -2, 3), (-4, 1, 2)):
            for axis in range(3):
                plus = tuple(n+int(j == axis) for j, n in enumerate(site))
                minus = tuple(n-int(j == axis) for j, n in enumerate(site))
                # Remove the exact source phase: each scalar lattice derivative
                # cancels its own onsite mass on the infinite Bessel vector.
                def envelope(point):
                    return mp.fprod(mp.besseli(abs(n), 1) for n in point)
                residual = (envelope(plus)-envelope(minus))/2+site[axis]*envelope(site)
                assert abs(residual) < mp.mpf('1e-95')


def test_lattice_domain_controls_and_conditional_scope():
    for values in ((0, 1, 1), (1, -1, 1), (1, 1, True)):
        with pytest.raises(ValueError):
            lattice.parameters(*values)
    for radius in (True, 0, -1, 1.):
        with pytest.raises(ValueError):
            lattice.jacobi(radius)
    with pytest.raises(ValueError):
        lattice.bessel_interval(0, 100, 1)
    result = lattice.certificate()
    assert result['full_exact_lattice_zero_mode_dimension'] == 2
    assert result['uniform_low_spectrum_continuum_rate_asserted'] is False
