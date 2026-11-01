"""Retained outputs, exact bounds, scope and cold source custody fail closed."""
from copy import deepcopy

import pytest

from tools import strict_json
from . import packet


@pytest.fixture(scope='module')
def replay():
    return packet.report()


def test_complete_retained_receipt_cold_replays(replay):
    retained = strict_json.load(packet.OUTPUT)
    packet.compare(retained, replay)
    packet.validate_accounting(retained)


@pytest.mark.parametrize('field,value', [('integrated_gate_accepted', True),
                                       ('physical_comparison_count_added', 1),
                                       ('physical_particle_identification', True),
                                       ('primitive_or_repair_law_changed', True),
                                       ('physical_hardware_execution_or_accuracy_proved', True)])
def test_physical_or_source_law_promotions_reject(replay, field, value):
    changed = deepcopy(replay); changed[field] = value
    with pytest.raises(ValueError):
        packet.compare(changed, replay)


def test_zero_export_slots_and_dark_modes_cannot_be_discarded(replay):
    changed = deepcopy(replay)
    changed['complete_native_Ramsey_history']['all40_native_repair_records'].pop()
    with pytest.raises(ValueError): packet.validate_accounting(changed)
    changed = deepcopy(replay)
    changed['complete_native_Ramsey_history']['pointer_settings'][0]['pointer_probabilities'].pop()
    with pytest.raises(ValueError): packet.validate_accounting(changed)
    changed = deepcopy(replay)
    changed['exact_source_accessible_spectrum']['all15_mode_pairs'].pop()
    with pytest.raises(ValueError): packet.validate_accounting(changed)


def test_small_per_leaf_errors_cannot_corrupt_source_abort_mass(replay):
    changed = deepcopy(replay)
    record = changed['complete_native_Ramsey_history']['common_prefix_source_records'][0]
    for row in record['all_first_export_failure_weights']: row[0] += 9e-11
    packet.compare(changed, replay)
    with pytest.raises(ValueError): packet.validate_accounting(changed)
    changed = deepcopy(replay)
    record = changed['complete_native_Ramsey_history']['common_prefix_source_records'][0]
    record['all_first_export_failure_weights'][record['pair_code_event_indices'][0]][3] = 1e-30
    with pytest.raises(ValueError): packet.validate_accounting(changed)


@pytest.mark.parametrize('bad', [-1e-30, float('nan'), float('inf'), True, 0, '0'])
def test_invalid_numerical_evidence_rejects(bad):
    with pytest.raises(ValueError):
        packet.compare({'probability': bad}, {'probability': 0.})


def test_exact_error_budget_cannot_be_replaced_by_a_float(replay):
    changed = deepcopy(replay)
    changed['exact_source_and_field_error_certificate']['complete_output_error_bounds'][0] = .008
    with pytest.raises(ValueError): packet.compare(changed, replay)
    with pytest.raises(ValueError):
        packet.compare({'host_numerical_admission': 3.1e-9}, {'host_numerical_admission': 3e-9})


def test_source_change_rejects_before_cold_worker(monkeypatch):
    changed = packet.source_bindings(); changed['model/dark_spectral_source/circuits.py'] = '0'*64
    monkeypatch.setattr(packet, 'source_bindings', lambda: changed)
    with pytest.raises(ValueError, match='fresh interpreter'): packet.report()


def test_transitive_source_manifest_binds_native_laws_templates_and_tests():
    paths = packet.source_bindings()
    for name in ('model/dark_spectral_source/circuits.py', 'model/dark_spectral_source/exact.py',
                 'model/dark_spectral_source/execute.py', 'model/dark_spectral_source/experiment.py',
                 'model/dark_spectral_source/test_source.py', 'model/dark_spectral_source/test_packet.py',
                 'model/dark_spectral_source/contract.json', 'model/dark_spectral_source/README.md',
                 'model/repair_spectral_sectors/results.json', 'model/certified_reference/results.json',
                 'model/certified_reference/source.py', 'model/certified_reference/certificate.py',
                 'model/observer/source_scalar/compiler.py', 'code/source_selection_model/response.json',
                 'paper/tex_fragments/DARK_SPECTRAL_SOURCE.tex'):
        assert name in paths
    assert 'model/dark_spectral_source/results.json' not in paths


def test_duplicate_keys_and_extra_claims_reject(replay):
    with pytest.raises(ValueError): strict_json.loads('{"done":false,"done":true}')
    with pytest.raises(ValueError): packet.compare({**replay, 'extra': True}, replay)
