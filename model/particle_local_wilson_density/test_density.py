from fractions import Fraction as F

import pytest
import sympy as s

from . import density as d


def test_full_native_kernel_identity_and_all_spectral_multiplicities():
    control = d.exact_controls()
    rough = control['noncommuting_native_symbol']
    assert rough['ordered_pair_identity_checks'] == 6
    assert rough['H_W_squared_spectral_ranks'] == {0: 16, 1: 28, 3: 16}
    assert rough['complete_auxiliary_symbol_probe_probabilities'] == {
        '0': s.Rational(4, 15), '1': s.Rational(7, 15), '3': s.Rational(4, 15)}
    assert rough['positive_plaquette_cost_sum'] == s.Rational(3, 5)
    assert rough['low_mode_count']**2 <= rough['universal_count_bound_squared']
    assert not control['point_coordinate_normalizable_preparation_claimed']


def test_native_representation_uses_actual_charge_rows_and_homomorphism():
    u = s.diag(s.eye(3), s.Matrix([[s.Rational(3, 5), 4*s.I/5],
                                 [4*s.I/5, s.Rational(3, 5)]]))
    v = s.diag(s.eye(3), s.Matrix([[s.Rational(5, 13), s.Rational(12, 13)],
                                 [-s.Rational(12, 13), s.Rational(5, 13)]]))
    assert d.clean(d.native_representation(u)*d.native_representation(v)) == d.native_representation(u*v)
    assert d.native_representation(s.eye(5)) == s.eye(15)


def test_positive_pair_split_on_noncommuting_native_rational_links():
    u = d.native_weak_link(0, s.Rational(3, 5), s.Rational(4, 5))
    v = d.native_weak_link(1, s.Rational(5, 13), s.Rational(12, 13))
    answer = d.decomposition((u, v, s.eye(15)))
    assert answer['identity_residual'] == s.zeros(60)
    assert answer['error'] == answer['error'].H
    assert answer['positive'] != s.zeros(60)
    assert answer['error'] != s.zeros(60)


def test_shift_commutator_is_unitary_transport_times_plaquette_defect():
    # Two genuine periodic sites, fifteen native rows each. The two inactive
    # directions retain their loops; an inhomogeneous weak link creates
    # curvature across the wrapping spatial shift.
    u = d.native_weak_link(0, s.Rational(3, 5), s.Rational(4, 5))
    v = d.native_weak_link(1, s.Rational(5, 13), s.Rational(12, 13))
    z = s.zeros(15)
    tx = s.BlockMatrix([[z, u], [s.eye(15), z]]).as_explicit()
    ty = s.diag(v, v.H)
    for x in (tx, tx.H):
        for y in (ty, ty.H):
            bracket = d.clean(d.comm(x, y))
            loop = d.clean(x*y*x.H*y.H)
            assert bracket == d.clean((loop-s.eye(30))*y*x)
            assert s.simplify(s.trace(bracket.H*bracket)
                              -s.trace((loop-s.eye(30)).H*(loop-s.eye(30)))) == 0
            assert bracket != s.zeros(30)


def test_same_density_plan_closes_local_weight_with_volume_free_constants():
    value = d.plan_bound()
    assert value['local_low_mode_trace_upper_squared'] <= F(5, 262144)
    assert value['simpler_local_low_mode_trace_upper_squared'] == F(5, 262144)
    assert value['bound_and_local_action_independent_of_total_volume']
    assert value['all_time_on_same_prepared_polynomial_field']
    assert value['original_probability_sum'] == 1
    assert not value['conditioned_outcomes_used']
    assert not value['interacting_Duhamel_record_bound_claimed']


def test_local_rank_cap_zero_density_and_region_scaling():
    zero = d.local_bound(0, cells=3, degree=128)
    assert zero['local_low_mode_trace_upper_squared'] == 0
    assert zero['one_body_local_rank'] == 180
    large = d.local_bound(6, cells=3, degree=0)
    assert large['local_low_mode_trace_upper_squared'] == 180**2
    assert large['auxiliary_uniform_local_mode_probe_low_probability_upper_squared'] == 1
    first = d.local_bound(F(1, 10**12), cells=1)
    triple = d.local_bound(F(1, 10**12), cells=3)
    assert triple['local_low_mode_trace_upper_squared'] == 9*first['local_low_mode_trace_upper_squared']
    assert triple['local_sign_minus_polynomial_HS_squared_upper']['rational_part'] == 3*first['local_sign_minus_polynomial_HS_squared_upper']['rational_part']


def test_invalid_or_inexact_domains_are_rejected():
    for value in (-1, 0.1, True):
        with pytest.raises(ValueError):
            d.local_bound(value)
    with pytest.raises(ValueError):
        d.local_bound(0, cells=0)
    with pytest.raises(ValueError):
        d.decomposition((s.eye(2), s.eye(2)))
    with pytest.raises(ValueError):
        d.decomposition((s.eye(2), 2*s.eye(2), s.eye(2)))
    with pytest.raises(ValueError):
        d.native_representation(s.diag(1, 1, 1, 1, -1))
