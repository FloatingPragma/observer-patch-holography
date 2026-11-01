from copy import deepcopy
import pytest
from . import packet


@pytest.fixture(scope='module')
def value():return packet._produce()


def test_cold_replay_matches_authenticated_direct_proof(value):
    assert packet.report()==value


def test_scope_cannot_be_extended_to_massive_interacting_or_other_ground(value):
    for key in ('nonzero_continuum_interaction_survives','massive_sector_supplied',
                'original_uncentered_ground_reidentified','unbounded_coupling_trajectory_decided',
                'physical_species_or_unit_dictionary_identified','Einstein_requirements_completed'):
        changed=deepcopy(value);changed['theorem'][key]=True
        with pytest.raises(ValueError):packet.validate(changed)


def test_both_charge_branches_and_cyclicity_are_indispensable(value):
    for key in ('both_vacuum_charge_branches_controlled','outside_initial_CAR_vacuum_cyclic'):
        changed=deepcopy(value);changed['theorem'][key]=False
        with pytest.raises(ValueError):packet.validate(changed)
    changed=deepcopy(value);changed['theorem']['massless_Weyl_field_count']=2
    with pytest.raises(ValueError):packet.validate(changed)
    changed=deepcopy(value);changed['bounds']['vacuum_duhamel_requires_uniform_input_mass_bound']=False
    with pytest.raises(ValueError):packet.validate(changed)


def test_energy_escape_cannot_be_changed_to_zero_equal_time_composite(value):
    changed=deepcopy(value);changed['theorem']['equal_time_composite_CAR_field_set_to_zero']=True
    with pytest.raises(ValueError):packet.validate(changed)
