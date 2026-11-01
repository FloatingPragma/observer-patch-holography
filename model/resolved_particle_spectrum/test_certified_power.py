from fractions import Fraction as F

import numpy as np
import pytest

from model.certified_reference import interval as iv
from model.resolved_particle_spectrum import certified_power as cp


def fractions(matrix):
    scale = 1 << matrix.bits
    return [[(F(matrix.real[i, j], scale), F(matrix.imag[i, j], scale))
             for j in range(matrix.dimension)] for i in range(matrix.dimension)]


def cmul(a, b):
    return a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0]


def product(a, b):
    d = len(a)
    return [[tuple(sum((cmul(a[i][k], b[k][j])[part] for k in range(d)), F())
                   for part in (0, 1)) for j in range(d)] for i in range(d)]


def exact_power(a, n):
    d = len(a)
    result = [[(F(int(i == j)), F()) for j in range(d)] for i in range(d)]
    for _ in range(n):
        result = product(result, a)
    return result


def row_error(actual, target):
    # Sum component absolute errors bounds complex row sum, independently
    # of the implementation's Frobenius error calculation.
    return max(sum(abs(x-y) for a, b in zip(ar, br, strict=True)
                   for x, y in zip(a, b, strict=True))
               for ar, br in zip(actual, target, strict=True))


def test_rounding_including_negative_ties_is_nearest():
    for numerator in range(-257, 258):
        q = F(numerator, 1024)
        rounded = cp._nearest_fraction(q, 8)
        assert abs(F(rounded, 256)-q) <= F(1, 512)


def test_interval_center_contains_exact_rational_unitary():
    real = [[F(3, 5), -F(4, 5)], [F(4, 5), F(3, 5)]]
    imag = [[0, 0], [0, 0]]
    base = cp.from_intervals(real, imag, bits=24)
    exact = [[(x, F()) for x in row] for row in real]
    assert row_error(fractions(base), exact) <= base.error
    assert base.error.denominator <= 1 << 24
    assert base.diagnostic_matrix().dtype == np.dtype('complex128')
    with pytest.raises(ValueError):
        base.real[0, 0] = 0


@pytest.mark.parametrize('exponent', [0, 1, 2, 3, 7, 31, 64, 129])
def test_rational_rotation_power_bound_against_independent_exact_arithmetic(exponent):
    re = [[F(3, 5), -F(4, 5)], [F(4, 5), F(3, 5)]]
    base = cp.from_intervals(re, [[0, 0], [0, 0]], bits=48)
    powered = cp.power(base, exponent)
    target = exact_power([[(x, F()) for x in row] for row in re], exponent)
    assert row_error(fractions(powered), target) <= powered.error
    assert powered.error.denominator <= 1 << 48
    expected = 0 if exponent == 0 else exponent.bit_count()+exponent.bit_length()-1
    assert powered.multiplications == expected


def test_millions_of_nontrivial_steps_are_logarithmic_and_exactly_checkable():
    # Cyclic Gaussian phases have an exact target at arbitrary exponent.
    base = cp.from_intervals([[0, 0], [0, 0]], [[0, 1], [1, 0]], bits=128)
    exponent = 7*(1 << 19)
    result = cp.power(base, exponent)
    assert np.array_equal(result.real, cp.identity(2, bits=128).real)
    assert not np.any(result.imag)
    assert result.multiplications == exponent.bit_count()+exponent.bit_length()-1
    assert result.multiplications < 30
    assert result.error < F(1, 10**28)


def test_complex_noncommuting_product_has_the_ordered_target():
    a = cp.from_intervals([[F(3, 5), 0], [0, F(3, 5)]],
                          [[F(4, 5), 0], [0, -F(4, 5)]], bits=48)
    b = cp.from_intervals([[F(3, 5), -F(4, 5)], [F(4, 5), F(3, 5)]],
                          [[0, 0], [0, 0]], bits=48)
    ab, ba = cp.multiply(a, b), cp.multiply(b, a)
    exact_a = [[(F(3, 5), F(4, 5)), (F(), F())], [(F(), F()), (F(3, 5), -F(4, 5))]]
    exact_b = [[(F(3, 5), F()), (-F(4, 5), F())], [(F(4, 5), F()), (F(3, 5), F())]]
    assert row_error(fractions(ab), product(exact_a, exact_b)) <= ab.error
    assert row_error(fractions(ba), product(exact_b, exact_a)) <= ba.error
    assert not np.array_equal(ab.imag, ba.imag)


def test_interval_kernel_and_explicit_endpoint_grids_both_enclose():
    sine, cosine = iv.sine_cosine(F(1, 7))
    real = [[cosine, -sine], [sine, cosine]]
    base = cp.from_intervals(real, [[0, 0], [0, 0]], bits=96)
    pairs = [[(F(x.lo, iv.SCALE), F(x.hi, iv.SCALE)) for x in row] for row in real]
    alternate = cp.from_intervals(pairs, [[0, 0], [0, 0]], bits=96)
    assert base.certificate() == alternate.certificate()
    for i in range(2):
        for j in range(2):
            center = F(base.real[i, j], 1 << base.bits)
            assert max(abs(center-pairs[i][j][0]), abs(center-pairs[i][j][1])) <= base.error


def test_exact_integer_expectation_normalization_and_complex_conjugation():
    base = cp.from_intervals([[F(3, 5), -F(4, 5)], [F(4, 5), F(3, 5)]],
                            [[0, 0], [0, 0]], bits=96)
    result = cp.expectation(base, [1, 0], [0, 1])
    assert result['integer_vector_norm_squared'] == 2
    assert abs(F(result['real_center'])-F(3, 5)) <= base.error
    assert abs(F(result['imag_center'])+F(4, 5)) <= base.error
    assert result == cp.expectation(base, [7, 0], [0, 7]) | {'integer_vector_norm_squared': 2}


@pytest.mark.parametrize('action', [
    lambda: cp.from_intervals([[.5]], [[0]]),
    lambda: cp.from_intervals([[0]], [[0]], bits=7),
    lambda: cp.from_intervals([[(F(1), F())]], [[0]]),
    lambda: cp.from_intervals([[1, 0]], [[0, 0]]),
    lambda: cp.CertifiedUnitary([[1.]], [[0]], 128, F()),
    lambda: cp.CertifiedUnitary([[1]], [[0]], 128, -F(1)),
    lambda: cp.power(cp.identity(1), -1),
    lambda: cp.power(cp.identity(1), True),
    lambda: cp.multiply(cp.identity(1), cp.identity(2)),
    lambda: cp.multiply(cp.identity(1), cp.identity(1, bits=64)),
    lambda: cp.expectation(cp.identity(1), [0]),
    lambda: cp.expectation(cp.identity(1), [.5]),
])
def test_unproved_numeric_or_malformed_inputs_are_rejected(action):
    with pytest.raises(ValueError):
        action()
