"""Focused exact-source, full-occupation and complete-record controls."""
from fractions import Fraction as F
from pathlib import Path
import sys

import pytest
import sympy as s

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from model.native_block_dirac import algebra, controller, reference, resources, source
from model.native_block_dirac.model import certificate, validate_pins


@pytest.fixture(scope='module')
def result():
    return certificate()


def test_full_constituent_gap_includes_both_ports_and_all_reference_sectors(result):
    b = result['complete_native_constituent']
    gaps = list(map(F, b['complete_sector_gap_lowers'].values()))
    assert len(gaps) == 3 and F(b['full_physical_constituent_gap_lower']) == min(gaps) > 0
    assert F(b['underlying_ec_current_Y']) == 1
    assert F(b['raw_full_kinematic_native_H_norm_upper']) > 480


def test_complete_CAR_encoding_and_singlet(result):
    E, P, aa = algebra.doublet_code()
    assert E[:, 3] == s.Matrix([0, 0, 0, 0, 0, 1/s.sqrt(2), 0, -1/s.sqrt(2), 0])
    assert P.rank() == 4
    assert aa[0]*aa[0].H+aa[0].H*aa[0] == P != s.eye(9)
    assert result['canonical_occupation_encoding']['canonical_CAR_pairs_checked'] == 16
    assert result['canonical_occupation_encoding']['complete_four_mode_Fock_code_dimension'] == 16


def test_both_frequency_spin_sums_and_noncollinear_boost(result):
    assert result['Dirac_frame_and_boost_algebra']['both_covariant_spin_sums_checked']
    D = s.Matrix(result['Dirac_frame_and_boost_algebra']['noncollinear_boost_Wigner_matrix'])
    assert (D-s.Matrix([[9, 1], [-1, 9]])/s.sqrt(82)).applyfunc(s.simplify) == s.zeros(2)


def test_actual_binary_cubic_and_all_record_reference_lift(result):
    g = result['finite_spatial_reference']
    assert g['actual_native_spin_lifts'] == 48 and g['exact_generator_products_checked'] == 96
    assert result['complete_reference_lift_control']['complete_lift_effect_sum_is_full_identity']


def test_driven_clock_noncommuting_control_and_raw_work_domain(result):
    assert result['noncommuting_pulse_clock_control']['noncommuting_intervention_control_passed']
    control = result['finite_member']['actual_source_clock_work_and_controller']
    assert control['total_elapsed_time'] == '1+5T/4'
    assert control['raw_constant_ground_energy_subtraction_is_not_used']
    norm = F(control['complete_raw_source_H_norm_upper'])
    assert F(control['net_bare_source_energy_change_absolute_upper_entire_schedule']) == 2*norm
    assert not control['autonomous_local_controller_proved']


def test_literal_complete_records_and_preparation_margin(result):
    member = result['finite_member']; rows = member['complete_sixteen_detector_records']
    assert len(rows) == 16
    assert s.simplify(sum(s.sympify(r['ideal_actual_source_probability']) for r in rows)) == 1
    for row in rows:
        low, high = map(F, row['conditional_native_preparation_interval'])
        assert 0 <= low <= high <= 1
        if row['first_down'] or row['second_down']: assert low == high == 0
    assert F(member['prepared_conditional_time_response_contrast_lower']) > F(8823, 100000)
    prep = member['native_preparation']
    assert F(prep['retained_purification_vector_error_upper']) < F(2336, 10**13)
    assert F(prep['all_preparation_success_probability_lower']) > F(9703, 10000)


def test_zero_filter_depth_is_deterministic_and_keeps_state_error():
    member = resources.finite_member(side=4, preparation_depth=0)
    prep = member['native_preparation']
    assert prep['complete_filter_records_retained'] == 0
    assert F(prep['all_preparation_success_probability_lower']) == 1
    assert F(prep['all_preparation_failure_probability_upper']) == 0
    assert F(prep['retained_purification_vector_error_upper']) > 0


def test_unconditional_amplification_retains_every_attempt_and_capacity(result):
    member = result['finite_member']; amp = member['optional_complete_preparation_amplification']
    assert amp['independent_full_attempts_retained'] == 6
    assert F(amp['all_attempts_failed_probability_upper']) < F(7, 10**10)
    assert F(amp['complete_record_TV_to_exact_ground_example_upper']) < F(1, 10**9)
    small = F(member['capacity']['all_preparation_sources_references_records_and_measurement_nats_interval'][1])
    large = F(amp['complete_source_preparation_reference_and_pointer_capacity_nats_interval'][0])
    assert 5*small < large < 6*small


def test_receipt_source_custody_and_no_empirical_mass_input(result):
    validate_pins(result)
    assert result['parent_native_join_receipt_sha256'] == source.PARENT_SHA
    assert not result['field_comparison']['mass_parameter_is_empirically_calibrated']
    broken = {**result, 'local_sha256': {'../escape': '0'*64}}
    with pytest.raises(ValueError, match='custody'):
        validate_pins(broken)


def test_immutable_research_scientific_comparison(result):
    comparison = result['independent_research_comparison']
    assert comparison['all_scientific_fields_equal'] and comparison['fresh_calculation_was_replayed']
    assert comparison['reference_sha256'] == 'df277440b874d997024a9446f4d84c2c086f572dccfbb283c39bc78861a87fec'
