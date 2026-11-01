from fractions import Fraction as F
from math import isqrt

import numpy as np
import pytest

from . import spectrum as s


def test_exact_source_cyclic_rank_and_annihilation_match():
    result = s.exact_data()
    assert len(result['polynomial']) == 79
    assert result['polynomial'][0] == (1, 0)
    assert result['rank_prime'] == 97 and result['rank_root']**2 % 97 == 3
    assert 0 < result['hankel_determinant'] < 97
    _, a, b, rows = s._source()
    assert s._annihilates(result['polynomial'], rows, a.tolist(), b.tolist())
    changed = list(result['polynomial'])
    c, d = changed[-1]
    changed[-1] = (c+1, d)
    assert not s._annihilates(changed, rows, a.tolist(), b.tolist())


def test_wrong_modular_proposal_cannot_admit_a_polynomial(monkeypatch):
    original = s._bm
    def wrong(sequence, prime):
        polynomial, degree = original(sequence, prime)
        polynomial[-1] = (polynomial[-1]+1) % prime
        return polynomial, degree
    s.exact_data.cache_clear()
    monkeypatch.setattr(s, '_bm', wrong)
    with pytest.raises(ArithmeticError, match='exact source annihilation'):
        s.exact_data()
    s.exact_data.cache_clear()


def test_floating_eigenvalues_only_propose_exhaustive_exact_brackets(monkeypatch):
    coefficients = tuple(s._coefficient(x) for x in s.exact_data()['polynomial'])
    monkeypatch.setattr(np.linalg, 'eigvalsh', lambda _: np.zeros(560))
    with pytest.raises(ArithmeticError, match='exhaust78 exact sign changes'):
        s._brackets(coefficients)


def test_integer_interval_arithmetic_contains_exact_products():
    a, b = s._Interval(-3*s.SCALE, 2*s.SCALE), s._Interval(-5*s.SCALE, 7*s.SCALE)
    assert (a*b).fractions() == (F(-21), F(15))
    root = s._coefficient((0, 1))
    assert root.lo**2 < 3*s.SCALE**2 < root.hi**2
    negative = s._coefficient((2, -3))
    assert negative.lo == 2*s.SCALE-3*root.hi
    assert negative.hi == 2*s.SCALE-3*root.lo
    with pytest.raises(ArithmeticError, match='denominator'):
        s._ratio(a, s._Interval(-1, 1))


def _product_interval(a, b):
    values = [x*y for x in a for y in b]
    return min(values), max(values)


def test_all78_positive_weights_normalize_and_reproduce_exact_moments():
    result = s.certificate()
    rows = result['rows']
    assert len(rows) == result['exact_local_cyclic_dimension'] == 78
    weights = [tuple(map(F, row['normalized_weight_interval'])) for row in rows]
    gaps = [tuple(map(F, row['gap_above_ground_interval'])) for row in rows]
    assert all(0 < a <= b < 1 and b-a <= F(2, 10**24) for a, b in weights)
    assert all(a[1] < b[0] for a, b in zip(gaps, gaps[1:]))
    for power, exact in ((0, F(1)), (1, F(1, 2)), (2, F(131, 256))):
        terms = [_product_interval(w, (g[0]**power, g[1]**power)) for w, g in zip(weights, gaps)]
        lower, upper = sum(x[0] for x in terms), sum(x[1] for x in terms)
        assert lower <= exact <= upper
        assert upper-lower < F(1, 10**20)
    assert result['remaining_symmetry_component_lines'] == 76
    assert result['unresolved_local_probe_weight'] == '0'
    assert result['native_preparation_detector_or_finite_shot_resolution_executed'] is False
    assert result['infinite_volume_poles_or_relativistic_masses_identified'] is False


def test_known_symmetry_lines_are_exactly_one_tenth_without_coincident_residual():
    result = s.certificate()
    scale = 10**30
    r = isqrt(3*scale**2)
    lam = (F(r, 16*scale), F(r+1, 16*scale))
    for shift in (0, 1):
        target = tuple(x+shift for x in lam)
        matches = [row for row in result['rows']
                   if F(row['gap_above_ground_interval'][0]) <= target[1]
                   and F(row['gap_above_ground_interval'][1]) >= target[0]]
        assert len(matches) == 1
        low, high = map(F, matches[0]['normalized_weight_interval'])
        assert low <= F(1, 10) <= high
