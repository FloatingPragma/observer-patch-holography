"""Native decoder, source algebra and complete cycle counterexamples."""
from dataclasses import replace
from fractions import Fraction as F
from . import source as s


def test_actual_E_effects_separate_all_H66_configurations():
    r = s.local_source_algebra()
    assert len({tuple(row) for row in r['actual_E_effect_signatures']}) == 66
    assert len(r['native_transport_spanning_tree']) == 65


def test_complete_native_cycle_does_not_identify_a_protected_charge():
    r = s.winding_control()
    assert r['cycle_probability_exact'] == '1/8'
    assert r['zero_chain_probability_exact'] == '1/4'
    assert r['one_repair_changes_prepared_zero_history_class_probability_exact'] == '1/2'
    assert sum(map(F, r['all_six_memory_class_probabilities'])) == 1
    assert not r['protected_Z6_source_charge_identified']


def test_decoder_rejects_forged_record_and_wrong_known_preparation():
    record = s.process.native.Record('repair', 'fixed', (0, 1), 'C', 0, 1, 1, 1)
    assert s.decode_positions((0,), (record,))['positions'] == (1,)
    assert s.decode_positions((1,), (record,))['status'] == 'incompatible_preparation'
    assert s.decode_positions((0,), (replace(record, current=-1),))['status'] == 'invalid_record'
    assert s.decode_positions((0,), (replace(record, raw=3),))['status'] == 'invalid_record'


def test_every_cut_and_separated_native_motion_are_complete():
    r = s.endpoint_motion(0)
    assert r['minimum_initial_endpoint_separation'] >= 3
    assert r['identification_failure_probability_exact'] == '0'
    assert sum(F(row[1]) for row in r['complete_decoded_position_distribution']) == 1
    assert sum(map(F, r['all_six_cut_current_mod6_probabilities'])) == 1
