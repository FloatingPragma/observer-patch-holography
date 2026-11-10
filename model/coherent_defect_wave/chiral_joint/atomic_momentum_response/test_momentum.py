"""Actual source attachment, complete phase records and rejected substitutions."""
from itertools import product
import pytest
import sympy as s

from .model import momentum_probability, momentum_multiplicity, shift_spins, full_phase_probability
from .checks import native_translation_join, character_checks, complete_reference_records, rotation_translation_join
from .pump import full_local_operator_checks, actual_composite_translation_checks, complete_pump_records


def test_actual_native_translation_and_spin_restriction():
    row = native_translation_join()
    assert row['all_N2_spin_basis_actual_CAR_translation_controls'] == 768
    assert row['actual_occupied_PW_H_translation_controls'] == 3
    assert row['actual_all_current_translation_covariance_controls'] == 36
    assert not row['atomic_reference_is_full_PW_stationary_vacuum']
    for axis in range(3):
        delta = tuple(int(j == axis) for j in range(3))
        for mask in range(256):
            assert shift_spins(shift_spins(mask, 2, delta), 2, delta) == mask


def test_exact_cycle_character_and_all_momentum_multiplicities():
    row = character_checks()
    assert [setting['side'] for setting in row['executed_sizes']] == [2, 4, 6]
    assert momentum_probability(2, (0, 0, 0)) == s.Rational(23, 128)
    assert momentum_probability(2, (1, 0, 0)) == s.Rational(15, 128)
    assert sum(momentum_multiplicity(4, k) for k in product(range(4), repeat=3)) == 2**64
    assert full_phase_probability(2, (1, 1, 1)) == 0
    assert full_phase_probability(2, (0, 1, 0)) == 0
    with pytest.raises(ValueError): momentum_probability(3, (0, 0, 0))
    with pytest.raises(ValueError): momentum_probability(2, (2, 0, 0))
    with pytest.raises(ValueError): shift_spins(256, 2, (0, 0, 0))


def test_complete_AP_cubic_translation_semidirect_action():
    row = rotation_translation_join()
    assert row['exact_AP_one_fermion_rotation_translation_conjugations'] == 72
    assert row['exact_scalar_rotation_translation_conjugations'] == 72
    assert row['exact_oriented_link_label_rotation_translation_conjugations'] == 1728
    assert not row['interacting_spectral_cluster_or_gap_inferred_here']


def test_physical_purification_and_reference_character_control():
    row = complete_reference_records()
    assert row['local_purifying_singlet_norm_squared'] == '1'
    assert row['local_purifying_singlet_spin_checks'] == 3
    assert row['complete_active_purification_dimension'] == 65536
    assert len(row['source_only_complete_2N_phase_records']) == 64
    assert row['all_other_full_source_reference_phase_pairs']['count'] == 4088
    assert row['reference_character_control']['source_only_probabilities_unchanged']
    assert not row['source_and_joint_translation_labels_identified']


def test_full_anticommutator_and_physical_normalized_pump():
    row = full_local_operator_checks()
    assert row['actual_J_anticommutator_spectrum'] == {'0': 10, '1': 4, '2': 6}
    assert not row['J_squared_is_zero']
    assert row['uniform_pump_effect_bound'] == '0 <= B^dagger B <= I/2'
    assert row['three_exact_spin_intertwiners'] == 3
    assert row['total_even_parity_intertwiner']
    translated = actual_composite_translation_checks()
    assert translated['actual_full_carrier_neutral_operator_translation_controls'] == 960
    assert not translated['even_or_mixed_parity_grid_substitution_allowed']


def test_actual_ground_weight_and_every_retained_pump_outcome():
    row = complete_pump_records()
    weight = row['ground_response_ancestry']
    assert weight['explicit_per_setting_to_pump_factor'] == '4'
    lo, hi = map(s.Rational, weight['success_probability_interval'])
    assert 0 < lo < hi < s.Rational(1, 2)
    assert row['source_success_T_to_N'] == '-I'
    for setting in row['executed_sizes']:
        side = setting['side']; total = s.Integer(0)
        assert setting['positive_record_count'] == 2 * side**3
        assert setting['complete_record_count_including_zero_outcomes'] == (
            setting['positive_record_count'] + setting['all_other_outcomes']['count'])
        for record in setting['complete_positive_records']:
            p = s.Rational(record['background_probability']); total += p
            a, b = map(s.Rational, record['joint_success_probability_interval'])
            c, d = map(s.Rational, record['joint_failure_probability_interval'])
            assert a == p*lo and b == p*hi and a+d == b+c == p
            assert all(k % 2 for k in record['source_success_2N_phase_indices'])
            assert all(k % 2 == 0 for k in record['source_failure_2N_phase_indices'])
        assert total == 1


def test_custody_rejects_missing_input_and_changed_response_source(monkeypatch):
    from . import certificate as custody
    from . import pump
    with monkeypatch.context() as patch:
        patch.setattr(custody.action, 'SOURCE_FILES', (*custody.action.SOURCE_FILES, 'missing-atomic-momentum-seed.json'))
        with pytest.raises(FileNotFoundError, match='missing declared'):
            custody.input_hashes()
    with monkeypatch.context() as patch:
        patch.setattr(custody, 'import_closure', lambda seeds: set(seeds) | {'missing-atomic-momentum-import.py'})
        with pytest.raises(FileNotFoundError, match='missing atomic momentum import'):
            custody.input_hashes()
    real_loads = pump.json.loads
    def changed_response(text):
        value = real_loads(text)
        if 'actual_ground_composite_spectrum' in value:
            value['source']['coefficients']['up'] = '0'
        return value
    monkeypatch.setattr(pump.json, 'loads', changed_response)
    with pytest.raises(AssertionError): pump.success_weight()
