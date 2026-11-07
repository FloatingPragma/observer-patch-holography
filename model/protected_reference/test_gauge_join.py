from fractions import Fraction as F

import pytest
import sympy as s

from . import gauge_join as g


def test_rates_and_orientation_are_read_from_actual_native_columns():
    value = g.native_edge()
    assert value['forward_probability'] == F(25, 33)
    assert value['reverse_probability'] == F(8, 33)
    assert value['target_reference_weight'] == F(200, 7505)
    assert value['source_reference_weight'] == F(64, 7505)
    assert len(value['native_columns']) == 4
    assert all(row['orientation_flip'] == 1 for row in value['native_columns'])


@pytest.mark.parametrize('charge', [-4, -3, 1, 2, 6])
def test_same_carrier_reference_gauge_transport_and_modular_clock(charge):
    value = g.case_certificate(charge)
    assert value['full_dimension'] == value['active_dimension']+value['retained_inactive_dimension']
    assert value['same_state_modular_gap_multiples']['forward'] == -1
    assert value['complete_effect_identity']
    assert value['all_maps_preserve_protected_charge']


def test_failure_holds_are_required_and_carry_zero_actual_motion():
    maps = g.active_maps(6)
    size = maps['forward'].rows
    complete = maps['forward'].H*maps['forward']+maps['forward_hold'].H*maps['forward_hold']
    assert g.transport.link.zero(complete-s.SparseMatrix.eye(size))
    assert not g.transport.link.zero(maps['forward'].H*maps['forward']-s.SparseMatrix.eye(size))
    position = s.ImmutableSparseMatrix(size, size, {(i, i): i//(size//2) for i in range(size)})
    assert g.transport.link.zero(position*maps['forward_hold']-maps['forward_hold']*position)
    assert not g.transport.link.zero(position*maps['forward']-maps['forward']*position)


def test_success_applies_actual_crossed_spectator_orientation_flip():
    maps = g.active_maps(6)
    n = maps['internal_link_dimension']
    initial = s.SparseMatrix(4*n, 1, {(2*n, 0): 1})
    moved = maps['forward']*initial
    assert (moved.H*moved)[0] == 1
    assert all(n <= i < 2*n for i, _ in moved.todok())


def test_unbalanced_coin_changes_the_same_native_reference():
    charge = 6
    maps = g.active_maps(charge)
    state = g.active_reference(charge)
    incorrect = s.SparseMatrix.zeros(state.rows)
    for key in ('forward', 'forward_hold', 'reverse', 'reverse_hold'):
        k = maps[key]
        incorrect += k*state*k.H/2
    assert not g.transport.link.zero(incorrect-state)


def test_coherent_hopping_is_not_the_stationary_modular_generator():
    maps = g.active_maps(6)
    state = g.active_reference(6)
    hopping = maps['forward']+maps['reverse']
    assert not g.transport.link.zero(hopping*state-state*hopping)


def test_complete_direct_sum_preserves_active_inactive_cross_terms_as_specified():
    # Independent smallest partial-isometry control with an explicit third
    # inactive source configuration and a coherent memory bit. The same full
    # Kraus formula is verified on every input column, including cross terms.
    flip = s.Matrix([[0, 1], [1, 0]])
    transition = s.SparseMatrix(3, 3, {(0, 1): 1})
    t = s.kronecker_product(transition, flip)
    identity = s.eye(6)
    held_f, held_r = identity-t.H*t, identity-t*t.H
    alpha, beta = s.Rational(25, 33), s.Rational(8, 33)
    assert alpha*(t.H*t+held_f.H*held_f)+beta*(t*t.H+held_r.H*held_r) == identity
    cross = s.zeros(6); cross[0, 4] = 1
    out = alpha*(t*cross*t.H+held_f*cross*held_f.H)+beta*(t.H*cross*t+held_r*cross*held_r.H)
    assert out == alpha*cross
    # Treating the complement as a measured/discarded block would give zero.
    assert out != s.zeros(6)


def test_q_charge_is_not_assigned_individual_particle_statistics():
    value = g.certificate()
    assert not value['scope']['single_moving_internal_fiber_is_two_particle_colored_code']
    assert not value['scope']['all_native_seams_have_a_common_gauge_network']
    assert not value['scope']['stationary_state_is_unique_or_mixing_proved']


@pytest.mark.parametrize('charge', [-4, -3, 1, 2, 6])
def test_stationary_singlet_corner_keeps_reference_and_domain_explicit(charge):
    value = g.singlet_certificate(charge)
    assert value['state_is_annihilated_by_endpoint_Gauss_generators']
    assert value['target_probability'] == '25/33'
    assert value['source_probability'] == '8/33'
    assert value['two_orthonormal_singlets_form_reducing_corner']
    assert value['full_native_source_marginal_replaced_by_active_conditional_marginal']
    assert not value['literal_singlet_preparation_pulses_executed']
    assert value['one_unread_attempt_exact_preparation_in_each_sharp_charge_corner']


def test_one_attempt_corner_preparation_retains_untouched_reference():
    # Exact Choi-state control, including the two possible native Q signs.
    target = s.diag(s.Rational(25, 33), s.Rational(8, 33))
    bell = s.Matrix([1, 0, 0, 1])/s.sqrt(2)
    for q in (-1, 1):
        operators = (q*s.Matrix([[0, 1], [0, 0]]), s.diag(1, 0),
                     q*s.Matrix([[0, 0], [1, 0]]), s.diag(0, 1))
        probabilities = (s.Rational(25, 33),)*2+(s.Rational(8, 33),)*2
        out = s.zeros(4)
        for p, k in zip(probabilities, operators, strict=True):
            branch = s.kronecker_product(k, s.eye(2))*bell
            out += p*branch*branch.H
        assert out == s.kronecker_product(target, s.eye(2)/2)
