"""Finite-source moments checked by combinatorics and response derivatives."""
import numpy as np
import sympy as s

from . import three_point


def _subset_moments(vectors, count):
    dimension = len(vectors)
    total = sum(vectors, s.zeros(3, 1))
    square = sum((v*v.T for v in vectors), s.zeros(3))
    mean = s.Rational(count, dimension)*total
    second = s.Rational(count, dimension)*square
    second += s.Rational(count*(count-1), dimension*(dimension-1))*(total*total.T-square)
    return mean, second


def test_moments_equal_independent_sampling_without_replacement_formula():
    _, _, basis, _, _ = three_point.source_data()
    core = three_point.local_preparation.CORE
    inside = [basis[i, :].T for i in core]
    outside = [basis[i, :].T for i in range(12) if i not in core]
    # P(c) = binomial(2,c) binomial(10,2-c) w_c from the native preparation.
    weights = [s.Rational(2, 225), s.Rational(1, 36), s.Rational(2, 45)]
    records = []
    for c in range(3):
        a, aa = _subset_moments(inside, c)
        b, bb = _subset_moments(outside, 2-c)
        records.append((s.binomial(2, c)*s.binomial(10, 2-c)*weights[c],
                        a+b, aa+bb+a*b.T+b*a.T))
    mean = three_point.clean(sum((p*m for p, m, _ in records), s.zeros(3, 1)))
    covariance = three_point.clean(sum((p*second for p, _, second in records), s.zeros(3))-mean*mean.T)
    for power in (1, 2):
        hmean = sum(p*c**power for c, (p, _, _) in enumerate(records))
        third = three_point.clean(sum((p*(c**power-hmean)*(second-m*mean.T-mean*m.T+mean*mean.T)
                                      for c, (p, m, second) in enumerate(records)), s.zeros(3)))
        actual = three_point.source_moments(power)
        assert actual['mean'] == mean
        assert actual['covariance'] == covariance
        assert actual['connected_third'] == third


def test_connected_tensor_is_actual_covariance_response():
    _, probabilities, _, records, core = three_point.source_data()
    p, q = np.asarray(probabilities, float), np.array([list(x) for x in records], float)
    for power in (1, 2):
        h = np.array(core, float)**power
        def covariance(theta):
            weights = p*np.exp(theta*h)
            weights /= weights.sum()
            mean = weights@q
            centered = q-mean
            return centered.T@(weights[:, None]*centered)
        step = 1e-4
        response = (covariance(step)-covariance(-step))/(2*step)
        np.testing.assert_allclose(response, np.asarray(three_point.source_moments(power)['connected_third'], float),
                                   rtol=2e-7, atol=2e-9)


def test_static_normalization_removes_rescaling_but_not_insertion_choice():
    first, second = (three_point.source_moments(i) for i in (1, 2))
    assert first['H_variance'] == s.Rational(644, 2025)
    assert second['H_variance'] == s.Rational(164, 225)
    # Positive scaling of H multiplies T and sqrt(VarH) by the same factor.
    for factor in (s.Rational(1, 7), s.Integer(11)):
        assert s.simplify(factor/s.sqrt(factor*factor)) == 1
    assert first['normalized_third'][1, 1].is_negative
    assert second['normalized_third'][1, 1].is_positive
    assert first['normalized_third'] != second['normalized_third']
