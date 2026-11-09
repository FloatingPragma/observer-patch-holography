from fractions import Fraction as F

import pytest
import sympy as s

from . import boundary as b


def test_literal_native_links_realize_complete_closed_spectrum():
    answer = b.exact_native_control()
    assert answer['exact_source_link_kernel_match']
    assert answer['native_spectral_control'] == {'1/16': 16, '1': 28, '7/2': 16}
    assert answer['all_sixty_columns_retained']
    assert answer['small_block_axis_derivative'] == 3
    assert answer['small_block_beta_derivative'] == 0
    assert answer['unit_electric_HS_squared'] == 9


def test_source_metric_and_all_projector_columns_are_retained():
    data = b.native_blocks()
    assert data['internal_generator_trace_metric_norm_squared'] == 8
    assert data['opposite_two_link_direction_norm_squared'] == 16
    assert [data[key].trace() for key in ('small', 'large', 'singlet')] == [16, 16, 28]
    assert data['small']+data['large']+data['singlet'] == s.eye(60)


def test_exact_tangent_satisfies_frechet_equations_at_multiple_gaps():
    for gap in (F(1, 4), F(1, 16), F(1, 256)):
        answer = b.symbol_family(gap)
        delta = s.Rational(gap.numerator, gap.denominator)
        assert answer['unit_electric_sign_derivative_norm'] == (1-delta)/(4*delta)
        assert answer['unit_electric_sign_derivative_HS_squared'] == (1-delta)**2/delta**2
        assert answer['low_mode_rank'] == 16


def test_complete_original_coordinate_records_and_retained_zero_labels():
    answer = b.pulse_records('1/4')
    assert len(answer['original_records']) == 12
    assert answer['complete_original_probability_sum'] == 1
    for row in answer['original_records']:
        assert row['original_probability'].is_nonnegative is True
    zero = [row for row in answer['original_records'] if row['original_probability'] == 0]
    assert len(zero) == 2 and all(row['block'] == 'singlet' for row in zero)
    assert answer['small_block_original_probability'] == s.Rational(4, 15)
    assert answer['small_block_conditional_flip_probability'] > s.Rational(1, 16)
    assert answer['coordinate_symbol_not_normalizable_preparation']


def test_vanishing_original_low_weight_does_not_bound_electric_cost():
    previous_weight, previous_cost = s.oo, 0
    for n in (4, 16, 256, 4096):
        answer = b.rare_packet_counter(s.Rational(1, n))
        assert answer['rare_symbol_low_mode_trace'] < previous_weight
        assert answer['rare_symbol_unit_electric_derivative_HS_squared'] > previous_cost
        assert answer['complete_probe_original_rare_flip_lower'] == s.Rational(1, 60*n)
        assert not answer['rare_branch_conditioning_used']
        assert not answer['counterfamily_is_the_existing_prepared_density']
        assert not answer['fixed_duration_interacting_response_impossible_claimed']
        previous_weight = answer['rare_symbol_low_mode_trace']
        previous_cost = answer['rare_symbol_unit_electric_derivative_HS_squared']


def test_inexact_or_invalid_gap_is_not_silently_accepted():
    for value in (0, -1, '1/2', '1', 0.25, True):
        with pytest.raises(ValueError):
            b.symbol_family(value)
    with pytest.raises(ValueError):
        b.rare_packet_counter('1/3')
