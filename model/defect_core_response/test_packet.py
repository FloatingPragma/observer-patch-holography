"""Receipt controls preserve source custody, pole weights and full histories."""
from copy import deepcopy

import pytest

from tools import strict_json
from . import packet


def stored():
    return strict_json.load(packet.OUTPUT)


def test_stored_receipt_matches_fresh_exact_checks():
    packet.validate(stored())


@pytest.mark.parametrize('damage', ['pole_weight', 'rank', 'quiet_probability', 'quiet_rate', 'archive', 'parent'])
def test_changed_physics_or_dropped_outputs_fail(damage):
    value = deepcopy(stored())
    if damage == 'pole_weight':
        value['response']['complete_source_poles']['unmodified_odd_cube'][
            'each_zero_mode_low_pole_weight_lower'] = '1/2'
    elif damage == 'rank':
        value['response']['complete_source_poles']['linear_texture_cube'][
            'each_addition_removal_pole_matrix_rank'] = 1
    elif damage == 'quiet_probability':
        value['response']['unit_rate_time_one_over_gap_lower_archive_controls'][
            'unmodified_odd_cube']['no_event_probability'] = '1'
    elif damage == 'archive':
        value['response']['unit_rate_time_one_over_gap_lower_archive_controls'][
            'linear_texture_cube']['every_event_and_abort_output_retained'] = False
    elif damage == 'quiet_rate':
        value['response']['declared_positive_rate_continuing_windows'][
            'unmodified_odd_cube']['rate_per_seam'] = '0'
    else:
        value['inherited_custody']['native_odd_box_core']['parent_expensive_receipt_reexecuted'] = True
    with pytest.raises(ValueError):
        packet.validate(value)


@pytest.mark.parametrize('flag', ['native_repair_protects_the_spectral_poles',
                                'physical_species_and_energy_units_identified',
                                'long_time_spectral_detector_executed',
                                'completion_gate_accepted'])
def test_unproved_physical_promotion_fails(flag):
    value = deepcopy(stored())
    value['scope'][flag] = True
    with pytest.raises(ValueError, match='scope'):
        packet.validate(value)


def test_missing_source_and_numeric_type_change_fail():
    value = deepcopy(stored())
    value['source_bindings'].pop(next(iter(value['source_bindings'])))
    with pytest.raises(ValueError, match='manifest'):
        packet.validate(value)
    value = deepcopy(stored())
    value['physical_comparison_count_added'] = False
    with pytest.raises(ValueError):
        packet.validate(value)


def test_cold_replay_rejects_mutated_in_process_response(monkeypatch):
    monkeypatch.setattr(packet.response, 'certificate', lambda: {'corrupted': True})
    # The worker computes from source; the caller detects its own inconsistent
    # implementation instead of accepting a counterfeit receipt.
    with pytest.raises(ValueError, match='response'):
        packet.report()
