"""Preserve the high-field band's source, grading and preparation boundaries."""
import copy
from fractions import Fraction as F

import pytest

from . import packet


@pytest.fixture(scope='module')
def owner():
    return packet.high_field.build()


def matched_owner(monkeypatch, value):
    """Alter producer and retained receipt together to test the consumer."""
    monkeypatch.setattr(packet.high_field, 'build', lambda: copy.deepcopy(value))
    monkeypatch.setattr(packet.strict_json, 'load', lambda _: copy.deepcopy(value))


def test_high_field_is_a_spectral_attachment_without_an_extra_object(owner, monkeypatch):
    matched_owner(monkeypatch, owner)
    value = packet.high_field_attachment()
    assert value['catalogue_entry'] == 'single_occupancy_spin_memory'
    assert value['additional_protected_objects_counted'] == 0
    assert value['additional_topological_classes_counted'] == 0
    assert value['identity_with_ordered_or_joint_small_field_limit_claimed'] is False
    assert value['physical_species'] == value['physical_masses'] == []
    rows = packet.rows()
    assert len(rows) == 11
    row = next(r for r in rows if r['id'] == value['catalogue_entry'])
    assert row['spectral_attachments'] == ['model/high_field_particle/PROOF.md']
    assert row['completeness'] == 'entire operator algebra of this code under the specified repair'


def test_complete_native_preparation_is_not_a_band_residue(owner, monkeypatch):
    matched_owner(monkeypatch, owner)
    value = packet.high_field_attachment()
    creation = value['complete_local_charge_creation']
    assert creation['successful_charge_by_sublattice'] == [-1, 1]
    assert creation['failure_and_out_of_code_records_retained'] is True
    assert creation['single_physical_fermion_under_grading_preserving_dictionary'] is False
    native = value['short_complete_native_preparation']
    assert native['actual_native_total_pulses'] == 167
    assert native['native_first_failure_slots'] == 608
    assert F(native['full_output_half_trace_error_from_exact_ground_instrument_upper']) < F(1, 100)
    assert native['ideal_exact_ground_success_probability_lower'] == '155233/155236'
    assert native['ideal_exact_ground_success_is_particle_band_weight'] is False
    assert native['infinite_volume_Yarotsky_threshold_certified_for_this_coupling'] is False


@pytest.mark.parametrize('section,key,value', [
    ('scope', 'physical_particle_assignments', ['electron']),
    ('scope', 'physical_mass_predictions', ['1']),
    ('scope', 'measured_mass_inputs_used', True),
    ('scope', 'finite_coupling_dispersion_remainder_certified', True),
    ('scope', 'all_source_regimes_excluded', True),
    ('scope', 'uniform_numeric_original_particle_residue_floor', '1/2'),
    ('scope', 'native_preparation_tape_executed_in_this_packet', 1),
    ('external_theorem_contract', 'explicit_high_field_admission_threshold', '1/400'),
    ('external_theorem_contract', 'band_existence', 'for every positive J/B'),
    ('external_theorem_contract', 'theorems_reproved_by_finite_matrix_tests', True),
    ('complete_local_charge_creation', 'success_renormalized', True),
    ('complete_local_charge_creation', 'single_physical_fermion_under_grading_preserving_dictionary', True),
    ('complete_local_charge_creation', 'abstract_full_Fock_completion_equals_literal_native_extension', True),
    ('complete_local_charge_creation', 'finite_coupling_success_is_a_pure_particle_band_projection', True),
    ('complete_local_charge_creation', 'success_renormalized', 0),
    ('short_complete_native_preparation', 'native_first_failure_slots', 0),
    ('short_complete_native_preparation', 'ideal_exact_ground_success_is_particle_band_weight', True),
    ('short_complete_native_preparation', 'infinite_volume_Yarotsky_threshold_certified_for_this_coupling', True),
    ('finite_ground_preparation_plan', 'native_tape_executed', True),
    ('exact_source_algebra', 'finite_periodic_double_edges_identified_with_bulk_unit_edges', True),
])
def test_matching_owner_and_receipt_cannot_promote_scope(owner, monkeypatch, section, key, value):
    forged = copy.deepcopy(owner)
    forged[section][key] = value
    matched_owner(monkeypatch, forged)
    with pytest.raises(ValueError, match='high-field'):
        packet.high_field_attachment()


def test_same_clock_speed_budget_rejects_matching_forgery(owner, monkeypatch):
    forged = copy.deepcopy(owner)
    forged['bounded_mass_speed_controls'][0]['physical_speed_upper'] = '1/1000'
    matched_owner(monkeypatch, forged)
    with pytest.raises(ValueError, match='same-clock speed bound'):
        packet.high_field_attachment()


def test_owner_receipt_mismatch_rejected(owner, monkeypatch):
    monkeypatch.setattr(packet.high_field, 'build', lambda: copy.deepcopy(owner))
    forged = copy.deepcopy(owner)
    forged['short_complete_native_preparation']['actual_native_total_pulses'] = 166
    monkeypatch.setattr(packet.strict_json, 'load', lambda _: forged)
    with pytest.raises(ValueError, match='owning source receipt differs'):
        packet.high_field_attachment()


def test_high_field_receipt_is_an_authenticated_parent():
    assert 'model/high_field_particle/results.json' in packet.PARENTS
