"""Cold replay, custody and explicit scientific boundary countercontrols."""
from copy import deepcopy

import pytest

from tools import strict_json
from . import packet


@pytest.fixture(scope='module')
def replay():
    return packet.report()


def test_retained_receipt_matches_independent_cold_replay(replay):
    packet.validate(strict_json.load(packet.OUTPUT), replay)


@pytest.mark.parametrize('key', [
    'initial_parity_dephasing_inserted',
    'reused_mean_replaced_by_power_of_averaged_matrix',
    'ensemble_parity_probability_estimated_consistently_from_one_field',
    'postmeasurement_field_stays_in_initial_spectral_sector',
    'autonomous_Hamiltonian_and_repair_detector_derived',
    'continuum_interacting_QFT_or_Einstein_dynamics_derived'])
def test_invalid_scope_promotion_is_rejected(replay, key):
    altered = deepcopy(replay)
    altered['scope'][key] = True
    with pytest.raises(ValueError):
        packet.validate(altered, replay)


def test_quantum_coherence_and_complete_control_bound_cannot_be_deleted(replay):
    altered = deepcopy(replay)
    altered['history']['one_ordered_balanced_history']['cross_parity'] = '0'
    with pytest.raises(ValueError):
        packet.validate(altered, replay)
    altered = deepcopy(replay)
    altered['single_preparation_classifier']['full_32_round_half_trace_error_upper'] = '0'
    with pytest.raises(ValueError):
        packet.validate(altered, replay)


def test_proof_custody_and_exact_comparison_count_are_enforced(replay):
    altered = deepcopy(replay)
    altered['source_bindings'].pop('model/radial_reused_readout/PROOF.md')
    with pytest.raises(ValueError):
        packet.validate(altered, replay)
    altered = deepcopy(replay)
    altered['physical_comparison_count_added'] = False
    with pytest.raises(ValueError):
        packet.validate(altered, replay)


def test_changed_source_after_import_is_rejected(monkeypatch):
    altered = dict(packet._IMPORT_BINDINGS)
    altered['model/radial_reused_readout/PROOF.md'] = '0'*64
    monkeypatch.setattr(packet, 'source_bindings', lambda: altered)
    with pytest.raises(ValueError, match='after reused readout producer import'):
        packet._produce()
