"""Cold replay and refusal of unsupported timing or protection claims."""
from copy import deepcopy

import pytest

from tools import strict_json
from . import packet


@pytest.fixture(scope='module')
def replay():
    return packet.report()


def test_retained_receipt_matches_cold_replay(replay):
    packet.validate(strict_json.load(packet.OUTPUT), replay)


@pytest.mark.parametrize('key', [
    'physical_rate_latency_and_clock_contract_instantiated',
    'comparison_generator_substituted_for_actual_source',
    'interleaved_record_law_equals_frozen_parity_mixture',
    'postmeasurement_spectral_sector_protection_proved',
    'uncontrolled_field_repairs_included'])
def test_unproved_physical_or_protection_promotions_rejected(replay, key):
    value = deepcopy(replay)
    value['scope'][key] = True
    with pytest.raises(ValueError):
        packet.validate(value, replay)


def test_omitting_elapsed_drift_or_complete_native_controls_rejected(replay):
    value = deepcopy(replay)
    row = value['complete_drift_and_classifier']['positive_time_witness']
    row['drift']['timed_control_complete_half_trace_error_upper'] = '0'
    with pytest.raises(ValueError):
        packet.validate(value, replay)
    value = deepcopy(replay)
    value['timed_native_commands']['one_round_commands'].pop()
    with pytest.raises(ValueError):
        packet.validate(value, replay)


def test_proof_and_exact_physical_count_custody(replay):
    value = deepcopy(replay)
    value['source_bindings'].pop('model/radial_drift_readout/PROOF.md')
    with pytest.raises(ValueError):
        packet.validate(value, replay)
    value = deepcopy(replay)
    value['physical_comparison_count_added'] = False
    with pytest.raises(ValueError):
        packet.validate(value, replay)


def test_changed_source_after_import_rejected(monkeypatch):
    value = dict(packet._IMPORT_BINDINGS)
    value['model/radial_drift_readout/PROOF.md'] = '0'*64
    monkeypatch.setattr(packet, 'source_bindings', lambda: value)
    with pytest.raises(ValueError, match='after reused readout producer import'):
        packet._produce()
