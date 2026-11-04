from copy import deepcopy

import pytest
from tools import strict_json

from . import packet


@pytest.fixture(scope='module')
def replay():
    return packet.report()


def test_complete_cold_source_replay(replay):
    packet.validate(strict_json.load(packet.OUTPUT), replay)


@pytest.mark.parametrize('field', ['four_mode_slice_is_full_core_vacuum',
    'four_mode_slice_reduces_each_individual_record', 'same_interaction_is_strictly_site_local',
    'all_twelve_currents_are_conserved_by_fixed_Higgs_Hamiltonian',
    'physical_rest_masses_or_observed_species_derived', 'interacting_chiral_relativistic_continuum_derived',
    'interacting_electromagnetic_pixel_feedback_closed', 'integrated_gate_accepted'])
def test_unsupported_promotions_rejected(replay, field):
    changed = deepcopy(replay)
    changed['scope'][field] = True
    with pytest.raises(ValueError, match='complete replay'):
        packet.validate(changed, replay)


def test_missing_spatial_modes_or_record_outcomes_rejected(replay):
    for section, field in [('exact_algebra', 'localized_annihilators'),
                           ('complete_record_probe', 'fine_records_per_memory_input')]:
        changed = deepcopy(replay)
        changed[section][field] -= 1
        with pytest.raises(ValueError, match='complete replay'):
            packet.validate(changed, replay)


def test_modified_readout_energy_rejected(replay):
    changed = deepcopy(replay)
    changed['actual_source_controls']['charged_lepton_slices'][0]['addition_plus_removal_poles_and_weights'][0][0] = '0'
    with pytest.raises(ValueError, match='complete replay'):
        packet.validate(changed, replay)


def test_source_bindings_exclude_research_notes():
    assert not any(path.startswith('docs/') for path in packet.source_bindings())
