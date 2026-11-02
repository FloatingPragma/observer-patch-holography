from copy import deepcopy
import pytest
from . import packet


@pytest.fixture(scope='module')
def value():
    return packet._produce()


def test_exact_receipt_accepts_and_discloses_inherited_custody(value):
    packet.validate(value)
    assert value['new_checks_pass'] is True
    assert value['prior_custody'][packet.PRIORS[0]]['historical_all_pass'] is False
    assert value['all_pass'] is False


@pytest.mark.parametrize('key,mutation', [
    ('physical_particle_species_identified', True),
    ('interacting_relativistic_continuum_constructed', True),
    ('first_spacing_counts_ground_multiplicity', False),
    ('infinite_GNS_gapless_particle_proved', True),
    ('large_carrier_native_pulse_execution_claimed', 0),
])
def test_unproved_physical_or_mistyped_scope_rejected(value, key, mutation):
    changed = deepcopy(value)
    changed['theorem'][key] = mutation
    with pytest.raises(ValueError):
        packet.validate(changed)


def test_missing_complete_record_outcomes_rejected(value):
    changed = deepcopy(value)
    changed['controls']['repair']['all_32_coin_raw_slots_retained'] = False
    with pytest.raises(ValueError):
        packet.validate(changed)


def test_historical_custody_cannot_be_promoted(value):
    changed = deepcopy(value)
    changed['all_pass'] = changed['dependency_custody_pass'] = True
    with pytest.raises(ValueError):
        packet.validate(changed)
