"""Proof receipt scope, inherited custody and stale-import regressions."""
from copy import deepcopy

import pytest

from tools import strict_json
from . import packet


@pytest.fixture(scope='module')
def retained():
    return strict_json.load(packet.OUTPUT)


def test_receipt_preserves_original_source_and_broader_parent_status(retained):
    packet.validate(retained)
    assert retained['core']['controlled_27_site_cube']['ground_particle_number'] == 27
    assert retained['inherited_custody']['native_spin_observer']['historical_all_pass'] is False
    assert retained['all_pass'] is True


@pytest.mark.parametrize('field', ['source_selected_onsite_mass_or_texture_assumed',
                                 'cube_full_Fock_matrix_or_native_readout_executed',
                                 'filled_sea_native_preparation_executed',
                                 'detector_protected_against_interleaved_repair_or_Hamiltonian_drift',
                                 'volume_uniform_nonzero_charged_gap_proved',
                                 'completion_gate_accepted'])
def test_physical_or_execution_scope_promotions_rejected(retained, field):
    changed = deepcopy(retained)
    changed['scope'][field] = True
    with pytest.raises(ValueError, match='boundary'):
        packet.validate(changed)


def test_boolean_laundering_and_missing_failure_outputs_rejected(retained):
    for key, replacement in [('reference_columns_dephased', 0), ('native_first_failure_slots', 0),
                             ('quantum_reference_columns', 32), ('complete_logical_pointer_outcomes', 1)]:
        changed = deepcopy(retained)
        changed['native_execution'][key] = replacement
        with pytest.raises(ValueError, match='carrier'):
            packet.validate(changed)


def test_source_pins_and_after_import_change_are_rejected(retained, monkeypatch):
    changed = deepcopy(retained)
    changed['source_bindings'].pop(next(iter(changed['source_bindings'])))
    with pytest.raises(ValueError, match='manifest'):
        packet.validate(changed)
    altered = dict(packet._IMPORT_BINDINGS)
    altered['model/native_odd_box_core/core.py'] = '0'*64
    monkeypatch.setattr(packet, 'source_bindings', lambda: altered)
    with pytest.raises(ValueError, match='after.*import'):
        packet._produce()


def test_historical_parent_failure_cannot_be_promoted(retained):
    changed = deepcopy(retained)
    changed['inherited_custody']['native_spin_observer']['historical_all_pass'] = True
    with pytest.raises(ValueError, match='inherited_custody'):
        packet.validate(changed)


def test_negative_nonfinite_or_unbounded_host_diagnostic_rejected(retained):
    for number in (-1e-30, float('nan'), float('inf'), 3e-11):
        changed = deepcopy(retained)
        changed['native_execution']['complete_output_error_diagnostic'] = number
        with pytest.raises(ValueError, match='diagnostic'):
            packet.validate(changed)
