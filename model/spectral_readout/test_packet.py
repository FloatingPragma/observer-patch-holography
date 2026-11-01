from copy import deepcopy

import pytest

from . import packet


@pytest.fixture(scope='module')
def value():return packet._produce()


def test_complete_source_and_probability_law(value):
    packet.validate(value)
    assert [row['outcome_count'] for row in value['finite_source_readouts']]==[8,16,64]
    assert value['finite_source_readouts'][-1]['joint_error_probability_decimal']==['0.020209031433','0.020209031434']


def test_missing_outcome_or_unproved_native_resolution_is_rejected(value):
    changed=deepcopy(value);changed['finite_source_readouts'][0]['rows'].pop()
    with pytest.raises(ValueError):packet.validate(changed)
    changed=deepcopy(value);changed['native_first_harmonic']['all_higher_harmonics_natively_executed']=True
    with pytest.raises(ValueError):packet.validate(changed)


def test_uncomputed_state_and_mass_cannot_be_promoted(value):
    changed=deepcopy(value);changed['finite_volume_error_contract']['total_infinite_ground_bin_error']='0'
    with pytest.raises(ValueError):packet.validate(changed)
    changed=deepcopy(value);changed['spectral_measure']['relativistic_rest_mass_inferred']=True
    with pytest.raises(ValueError):packet.validate(changed)


def test_cold_worker_matches_in_process_certificate(value):
    assert packet.report()==value
