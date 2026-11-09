"""Original spectral weights, exact source moments and complete inverse scores."""
from fractions import Fraction as F
import json

import pytest
import sympy as s

from . import speed_concentration as sc


def moments(atoms):
    return (sum(w for x, w in atoms), sum(x*w for x, w in atoms),
            sum(w/x for x, w in atoms))


def test_uncut_identity_and_probability_bound_on_independent_atomic_laws():
    laws = [((F(1), F(1, 3)), (F(2), F(1, 4)), (F(1000), F(1, 10000))),
            ((F(1, 1000), F(1, 10000)), (F(3, 2), F(1, 2))),
            ((F(7, 3), F(1, 5)),)]
    for atoms in laws:
        w, h, chi = moments(atoms)
        for c in (F(1, 4), F(1), h/w):
            direct = sum(weight*(x-c)**2/x for x, weight in atoms)
            for eps in (F(1, 8), c, 2*c):
                result = sc.concentration((w, w), (h, h), (chi, chi), c, eps)
                outside = sum(weight for x, weight in atoms if abs(x-c) >= eps)
                assert result['weighted_residual_interval'] == [direct, direct]
                assert outside <= result['original_outside_probability_upper']
                tail = sum(x*weight for x, weight in atoms if x >= 2*c)
                assert tail <= result['ratio_first_moment_above_twice_speed_upper']
                assert result['original_inside_probability_lower'] <= w-outside


def test_tail_constant_is_sharp_at_upper_tube_boundary():
    c, eps, w = F(3, 2), F(2, 3), F(2, 5)
    x = c+eps
    row = sc.concentration((w, w), (w*x, w*x), (w/x, w/x), c, eps)
    assert row['original_outside_probability_upper'] == w
    assert row['original_inside_probability_lower'] == 0


def test_probability_weak_convergence_does_not_hide_first_energy_escape():
    n = 10000
    atoms = ((F(1), F(1)-F(1, n)), (F(n), F(1, n)))
    w, h, chi = moments(atoms)
    row = sc.concentration((w, w), (h, h), (chi, chi), 1, F(1, 4))
    assert row['weighted_residual_interval'][0] > F(9, 10)
    assert sum(x*weight for x, weight in atoms if x > 2) == 1


def test_positive_full_response_defect_allows_a_sharp_subbranch():
    atoms = ((F(1), F(1, 3)), (F(4), F(1, 6)))
    w, h, chi = moments(atoms)
    row = sc.concentration((w, w), (h, h), (chi, chi), h/w, F(1, 2))
    assert row['weighted_residual_interval'][0] > 0
    assert atoms[0][1] == F(1, 3)
    assert not row['actual_source_residual_vanishing_proved']


def test_complete_probability_normalization_is_retained():
    args = ((F(1, 4),)*2, (F(1, 2),)*2, (F(1, 8),)*2, 2, F(1, 10))
    first = sc.concentration(*args, normalization=1)
    diluted = sc.concentration(*args, normalization=3)
    assert first['original_inside_probability_lower'] == F(1, 4)
    assert diluted['original_inside_probability_lower'] == F(1, 36)
    assert diluted['weighted_residual_interval'] == [0, 0]
    assert not diluted['success_probability_divided_out']


def test_interval_errors_are_propagated_with_the_correct_signs():
    atoms = ((F(1), F(1, 4)), (F(3, 2), F(1, 5)))
    w, h, chi = moments(atoms)
    eps = F(1, 1000)
    exact = sc.concentration((w,)*2, (h,)*2, (chi,)*2, 1, F(1, 2))
    bounded = sc.concentration((w-eps, w+eps), (h-eps, h+eps),
                               (chi-eps, chi+eps), 1, F(1, 2))
    assert bounded['weighted_residual_interval'][1] == exact['weighted_residual_interval'][1]+4*eps
    assert bounded['original_inside_probability_lower'] <= exact['original_inside_probability_lower']
    assert bounded['original_outside_probability_upper'] >= exact['original_outside_probability_upper']


@pytest.mark.parametrize('args', [
    ((1, 1), (1, 1), (F(1, 2),)*2, 1, 1),
    ((2, 2), (2, 2), (2, 2), 1, 1),
    ((1, 0), (1, 1), (1, 1), 1, 1),
    ((1, 1), (0, 0), (1, 1), 1, 1),
    ((1, 1), (1, 1), (1, 1), 0, 1),
    ((1, 1), (1, 1), (1, 1), 1, 0),
])
def test_inconsistent_or_invalid_moments_are_rejected(args):
    with pytest.raises(ValueError):
        sc.concentration(*args)


@pytest.mark.parametrize('value', [True, 0.5, '1/2'])
def test_inexact_or_coerced_inputs_are_rejected(value):
    with pytest.raises(TypeError):
        sc.concentration((1, 1), (1, 1), (1, 1), value, 1)


def test_finite_gap_enclosure_against_an_independent_spectral_law():
    # Different energy and spatial labels; X is their ratio, not energy alone.
    events = ((F(1, 2), F(1, 4), F(1, 5)),
              (F(3), F(1, 2), F(2, 7)), (F(5), F(1), F(1, 11)))
    zeta, j, lam = F(1, 7), F(2), F(3)
    score = sum(w*zeta/(e/(j*k)+zeta) for e, k, w in events)/(lam*lam)
    exact = sum(w*j*k/e for e, k, w in events)
    row = sc.inverse_from_score((score, score), 0, lam, zeta, F(1, 2), j, (F(1, 4), 1))
    lo, hi = row['ratio_inverse_interval']
    assert lo <= exact <= hi
    assert row['ratio_gap_lower'] == F(1, 4)
    assert row['known_original_score_factor'] == zeta/(lam*lam)
    widened = sc.inverse_from_score((score, score), F(1, 1000), lam, zeta,
                                   F(1, 2), j, (F(1, 4), 1))
    assert widened['ratio_inverse_interval'][0] < lo
    assert widened['ratio_inverse_interval'][1] > hi


def test_joint_score_error_controls_energy_rounding_and_keeps_failures():
    zeta, j, sm, width = F(1, 10), F(2), F(1, 4), F(1, 1000)
    bound = sc.joint_score_error(zeta, j, sm, width, F(1, 50), F(1, 100))
    score = lambda e, k: zeta*j*k/(max(e, F(0))+zeta*j*k)
    for energy in (F(0), F(1, 100), F(2)):
        for k in (sm, F(1, 2), F(1)):
            assert abs(score(energy+width, k)-score(energy, k)) <= width/(zeta*j*sm)
    assert bound['absolute_complete_score_error_upper'] == F(1, 20)
    assert bound['all_failure_scores'] == 0
    assert bound['spatial_labels_require_source_recentered_chart']


def test_actual_periodic_source_detector_requires_transverse_label_recentring():
    # Four periodic sites have genuine translation labels, unlike the open cube.
    n, size = 4, 16
    integer, scale = sc.source.matrix(n, ((0, 1), (1, 2), (2, 3), (0, 3)), F(1), F(1, 4))
    marshall = s.diag(*(1 if sum((mask >> j) & 1 for j in (1, 3)) % 2 == 0 else -1
                       for mask in range(size)))
    h = marshall*s.Matrix(integer.tolist())*marshall/scale
    translation = s.zeros(size)
    for mask in range(size):
        shifted = sum(((mask >> j) & 1) << ((j+1) % n) for j in range(n))
        translation[shifted, mask] = 1
    global_z = s.diag(*(1 if mask.bit_count() % 2 == 0 else -1 for mask in range(size)))
    g = translation*global_z
    assert g*h == h*g and g**n == s.eye(size)
    for component, raw_character in (('z', -s.I), ('y', s.I)):
        packet = s.zeros(size)
        for j in range(n):
            spin = s.zeros(size)
            for mask in range(size):
                sign = 1-2*((mask >> j) & 1)
                if component == 'z':
                    spin[mask, mask] = s.Rational(sign, 2)
                else:
                    spin[mask ^ (1 << j), mask] = s.I*sign/2
            packet += s.I**j*spin
        assert (g*packet*g.H-raw_character*packet).applyfunc(s.simplify) == s.zeros(size)
    # The same physical k=pi/2 has opposite raw characters for S2 and S3.
    assert (s.I)/(-s.I) == -1


def test_exact_joint_gap_scale_and_regularizer_are_not_uniform_in_n():
    first = sc.joint_regularization(2, 3, F(1, 16))
    second = sc.joint_regularization(4, 3, F(1, 16))
    assert second['regularizer'] == first['regularizer']/8
    assert second['source_energy_regularizer_at_outer_symbol_over_J'] == first['source_energy_regularizer_at_outer_symbol_over_J']/32
    assert second['ratio_gap_lower'] == second['B_over_J']/(8*second['spatial_symbol_upper'])
    assert not second['uniform_scaling_limit_convergence_inferred']


def test_static_source_prefactors_and_interval_enlargement():
    row = sc.static_first_moment((F(2, 5),)*2, (F(1, 3),)*2, F(1, 32), (2, 2), (1, 1))
    assert row['ratio_first_moment_interval'] == [F(389, 960)]*2
    broad = sc.static_first_moment((F(1, 3), F(2, 5)), (F(1, 4), F(1, 3)),
                                  F(1, 32), (1, 2), (F(1, 2), 1))
    assert broad['ratio_first_moment_interval'][0] < row['ratio_first_moment_interval'][0]
    assert broad['ratio_first_moment_interval'][1] == row['ratio_first_moment_interval'][1]


def test_exact_actual_source_first_inverse_and_complete_six_records():
    row = sc.finite_source_control()
    assert row['first_moment'] == row['kappa']*2+row['B']*row['order']
    assert row['certified_inverse_upper'] == row['inverse_moment']
    assert len(row['records']) == 6
    assert sum(F(r['original_probability']) for r in row['records']) == 1
    assert sum(F(r['original_probability']) for r in row['records'] if r['preparation'] == 'success') == F(4, 5)
    assert any(F(r['original_probability']) == 0 for r in row['records'])
    assert not row['ratio_divisor_is_spatial_momentum']
    assert not row['additional_native_tape_executed']
    # The conservative B/8 gap cannot falsely certify equality at fixed zeta.
    universal = sc.inverse_from_score((row['complete_score'],)*2, 0, 1,
        row['regularizer'], row['B']/8, row['J'], (1, 1))
    assert universal['ratio_inverse_interval'][1] > row['inverse_moment']


def test_source_control_rejects_changed_source_matrix(monkeypatch):
    original = sc.source.matrix
    def changed(*args, **kwargs):
        matrix, scale = original(*args, **kwargs)
        matrix = matrix.copy()
        matrix[0, 0] += 1
        return matrix, scale
    monkeypatch.setattr(sc.source, 'matrix', changed)
    with pytest.raises(ArithmeticError, match='source frame'):
        sc.finite_source_control()


def test_same_source_moment_bounds_allow_corrected_sharp_and_two_speed_controls():
    m, kap = s.Rational(1, 4), s.Rational(1, 4)
    a, b = 12*kap/m, kap/m
    two = [(s.sqrt(a*bi), m*s.sqrt(bi/a)/4) for bi in (b/2, 3*b/2)]
    w = sum(weight for speed, weight in two)
    h = sum(speed*weight for speed, weight in two)
    chi = sum(weight/speed for speed, weight in two)
    assert s.simplify(h-kap/2) == 0 and s.simplify(chi-m*m/(24*kap)) == 0
    assert s.simplify(h*chi/w**2-1) == 7-4*s.sqrt(3)
    sharp_speed, sharp_w = s.sqrt(a*b), m*s.sqrt(b/a)/2
    assert s.simplify(sharp_speed*sharp_w-h) == 0
    assert s.simplify(sharp_w/sharp_speed-chi) == 0
    assert sharp_speed*sharp_w*(sharp_w/sharp_speed)/sharp_w**2-1 == 0
    assert kap-m*m > s.Rational(61, 14400)


def test_public_certificate_is_compact_and_fraction_detached():
    first = sc.certificate()
    expected = repr(sc.certificate())
    assert len(json.dumps(first, default=str)) < 15000
    first['finite_source_control']['B']._numerator = 900
    first['exact_source_speed_tube']['speed']._numerator = 901
    assert repr(sc.certificate()) == expected


@pytest.mark.parametrize('call', [
    lambda: sc.inverse_from_score((0, 2), 0, 1, 1, 1, 1, (1, 1)),
    lambda: sc.inverse_from_score((0, 1), 0, 1, 1, 1, 1, (0, 1)),
    lambda: sc.inverse_from_score((0, 1), 0, 1, 1, 0, 1, (1, 1)),
    lambda: sc.joint_score_error(1, 1, 1, 0, 2, 0),
    lambda: sc.joint_regularization(True, 1, 1),
    lambda: sc.joint_regularization(2, 1, 2),
    lambda: sc.static_first_moment((0, 1), (0, F(1, 2)), 1, (1, 1), (1, 1)),
])
def test_invalid_gap_score_or_source_parameters_fail(call):
    with pytest.raises((ValueError, TypeError)):
        call()
