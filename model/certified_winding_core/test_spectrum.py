"""Independent integer, source-enclosure and full-spectrum countercontrols."""
from fractions import Fraction as F

import numpy as np
import pytest
import sympy as s
from scipy import sparse

from ..winding_sector_search import search
from . import intervals as iv, spectrum, core


def test_exact_limb_products_match_arbitrary_precision_including_negative_limbs():
    rng = np.random.default_rng(8192)
    a = rng.integers(-2**26, 2**26, (17, 31), dtype=np.int64)
    b = rng.integers(-2**26, 2**26, (31, 13), dtype=np.int64)
    actual = spectrum.exact_dot(a, b)
    expected = a.astype(object)@b.astype(object)
    assert np.array_equal(actual, expected)
    with pytest.raises(ValueError, match='domain'):
        spectrum.exact_dot(np.full((2, 2), 2**62, dtype=np.int64), np.ones((2, 2), dtype=np.int64))


@pytest.mark.parametrize('degree', [1, 2, 3, 4, -3])
def test_source_intervals_enclose_independent_high_precision_texture(degree):
    for point in ((1, 2, 3), (2, -1, -3), (0, 0, 1), (0, 0, -2), (0, 0, 0)):
        bounds = iv.direction(point, degree)
        x, y, z = point
        if x == y == 0:
            expected = [s.Integer(0), s.Integer(0), s.sign(degree)*s.sign(z)]
        else:
            a = (x+s.I*y)**abs(degree); b = (s.sqrt(x*x+y*y+z*z)-z)**abs(degree)
            den = s.expand_complex(a*s.conjugate(a))+b*b
            expected = [2*s.re(a)*b/den, 2*s.im(a)*b/den, (a*s.conjugate(a)-b*b)/den]
            expected = [s.sign(degree)*v for v in expected]
        for (lo, hi), value in zip(bounds, expected, strict=True):
            left = s.N(value-s.Rational(lo.numerator, lo.denominator), 90)
            right = s.N(s.Rational(hi.numerator, hi.denominator)-value, 90)
            assert left >= 0 and right >= 0


def test_full_source_matrix_and_rational_error_match_original():
    for degree in (1, 3, 4):
        hr, hi, bound = spectrum.source(1, degree, 1)
        actual = search.source_matrix(1, degree, 1).toarray()
        candidate = (hr.toarray()+1j*hi.toarray())/2**26
        assert np.linalg.norm(actual-candidate, 'fro')**2 <= float(bound)*(1+1e-6)


def test_bad_complete_basis_and_integer_overflow_are_rejected():
    factor = 1 << 20
    h = sparse.diags([factor, 2*factor], dtype=np.int64, format='csr')
    z = sparse.csr_matrix((2, 2), dtype=np.int64)
    e = np.array([factor, 2*factor], dtype=np.int64)
    q = factor*np.eye(2, dtype=np.int64); qi = np.zeros((2, 2), dtype=np.int64)
    good = spectrum.certify(h, z, e, q, qi, F(0), 20)
    assert good['ordered_eigenvalue_error_upper'] == '0'
    with pytest.raises(ArithmeticError, match='invertible'):
        spectrum.certify(h, z, e, np.zeros_like(q), qi, F(0), 20)
    with pytest.raises(ValueError, match='overflow'):
        spectrum.certify(h, z, e, np.full_like(q, np.iinfo(np.int64).min), qi, F(0), 20)


def test_weak_mass_kernel_exhaustion_and_complete_contact_bounds():
    case = spectrum.calculate(1, 3, 1)
    assert case['window_rank_exact'] == 6 and case['exact_kernel_dimension'] == 2
    assert F(case['nonzero_gap_lower']) > F(1, 20)
    failed = core.bounds(case, F(1, 65536))
    assert not failed['bound_proves_ground_and_charged_pair']
    admitted = core.bounds(case, F(1, 2**24))
    assert admitted['bound_proves_ground_and_charged_pair']
    assert admitted['exact_ground_multiplicity'] == admitted['each_pole_matrix_rank'] == 2
    assert admitted['full_Fock_dimension'] == 2**54
    assert F(admitted['charged_gap_interval'][0]) > 0
    case['degree'] = 2
    assert not core.bounds(case, F(1, 2**40))['bound_proves_ground_and_charged_pair']
