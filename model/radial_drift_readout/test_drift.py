"""Exact source, finite-time and complete-history regression controls."""
from fractions import Fraction as F

import pytest

from . import algebra, bounds, commands, controls


def test_actual_source_and_full_CAR_norms():
    value = algebra.certificate()
    assert value['actual_origin_block']['origin_coupling_Gram_at_a_v_one'] == [['3/2', '0'], ['0', '3/2']]
    assert value['exact_complete_CAR_star']['star_over_singular_value_spectrum'] == {
        '-2': 1, '-1': 4, '0': 6, '1': 4, '2': 1}
    row = bounds.local_drift(F(1))
    assert row['origin_parity_projector_commutator_norm_squared_exact'] == 6
    assert row['origin_sign_parity_commutator_norm_squared_exact'] == 24
    assert row['star_plus_center_contact_norm_squared_exact'] == 6+F(1, 16*65536**2)


@pytest.mark.parametrize('value', [F(0), F(1), F(6), F(1, 7), F(49, 9)])
def test_rational_square_root_encloses_without_floating_admission(value):
    lo, hi = bounds.sqrt_interval(value, 19)
    assert lo*lo <= value <= hi*hi
    assert hi-lo <= F(1, 2**19)


@pytest.mark.parametrize('value', [True, 0.1, -1])
def test_invalid_clock_and_scalar_values_rejected(value):
    with pytest.raises(ValueError):
        bounds.local_drift(value)


def test_velocity_spacing_scaling_and_volume_independent_domain():
    one = bounds.local_drift(F(1, 10), velocity=F(3), spacing=F(2), coupling=F(2), dimensions=1)
    three = bounds.local_drift(F(1, 10), velocity=F(3), spacing=F(2), coupling=F(2), dimensions=3)
    assert one['star_Fock_norm_squared_exact'] == F(9, 2)
    assert three['star_Fock_norm_squared_exact'] == F(27, 2)
    assert three['star_plus_center_contact_norm_squared_exact'] == F(55, 4)
    assert three['source_volume_factor'] == 1 and three['mass_texture_factor'] == 0
    assert three['timed_control_complete_diamond_error_upper'] == 2*three['timed_control_complete_half_trace_error_upper']


def test_positive_time_classifier_and_failed_timing_budget():
    witness = bounds.classifier_bound()
    assert 0 < witness['drift']['total_elapsed_Hamiltonian_time']
    assert witness['misclassification_strict_upper'] < F(1, 2500)
    assert witness['below_one_in_2500']
    assert not bounds.classifier_bound(F(1, 10000))['below_one_in_2500']
    assert not witness['record_law_is_exact_frozen_parity_mixture']


def test_all_timed_commands_preserve_the_source_word_and_finite_resources():
    normal = commands.timing_plan()
    faster = commands.timing_plan(F(1, 4*10**6))
    assert normal['native_pulses_total'] == 2432
    assert normal['fresh_logical_apparatus_bits'] == 736
    assert normal['integer_repair_attempts'] == 32
    assert normal['complete_round_word_sha256'] == faster['complete_round_word_sha256']
    assert normal['timed_command_sha256'] != faster['timed_command_sha256']
    duration = F(normal['positive_pulse_duration'])
    assert 2432*duration+32*F(normal['positive_nonpulse_time_per_round']) == F(1, 10**6)
    for left, right in zip(normal['one_round_commands'], faster['one_round_commands'], strict=True):
        assert left['twelve_exact_integrated_port_coefficients'] == right['twelve_exact_integrated_port_coefficients']
        for value, rate, fast in zip(left['twelve_exact_integrated_port_coefficients'],
                                     left['twelve_exact_port_rates'], right['twelve_exact_port_rates'], strict=True):
            assert F(value) == F(rate)*duration
            assert F(fast) == 4*F(rate)
    assert not normal['physical_port_bandwidth_or_clock_calibration_derived']


def test_zero_duration_is_not_an_admitted_finite_rate_schedule():
    with pytest.raises(ValueError):
        commands.timing_plan(F(0))


def test_independent_complete_Fock_interleaving_and_continuous_control():
    value = controls.certificate()
    assert value['interleaved_actual_source']['actual_record_probabilities_differ_from_frozen_mixture']
    assert value['interleaved_actual_source']['conditional_cross_parity_blocks_retained']
    assert value['simultaneous_continuous_control']['simultaneous_local_control_and_full_H_tested']
    assert value['simultaneous_continuous_control']['complete_coherent_pointer_output_retained']
