import copy
from fractions import Fraction as F
import numpy as np
import pytest
from . import packet, spectral


@pytest.fixture(scope='module')
def admitted():
    value, arrays = packet.produce()
    packet.validate(value, arrays)
    return value, arrays


def test_complete_true_source_joint_readout(admitted):
    value, arrays = admitted
    execution = value['execution']
    assert execution['pulses'] == 19228
    assert execution['native_failure_slots'] == 67982
    assert F(execution['complete_isometry_error']) < F(1, 10**10)
    assert execution['complete_output_distance_diagnostic'] < 1e-11
    assert arrays['success'].shape == (128, 1)
    assert arrays['failure_weights'].shape == (19228, 4)
    cube = value['cube_joint_bins']
    assert cube['all_source_excitation_energies_in_one_phase_period'] is True
    all_intervals = [bounds for row in cube['rows'] for bounds in row['energy_bin_probability_intervals']]
    assert sum(F(lo) for lo, _ in all_intervals) <= 1 <= sum(F(hi) for _, hi in all_intervals)
    assert all(0 < F(lo) <= F(hi) < 1 for lo, hi in all_intervals)


@pytest.mark.parametrize('field,value', [('native_full_cube_executed', True),
                                         ('sharp_scaling_shell_proved', True),
                                         ('ground_specific_dimer_split_error_zero', 1)])
def test_scope_and_type_promotions_rejected(admitted, field, value):
    original, arrays = admitted
    changed = copy.deepcopy(original)
    changed['theorem'][field] = value
    with pytest.raises(ValueError):
        packet.validate(changed, arrays)


def test_drop_native_complement_or_change_observer_probability_rejected(admitted):
    value, arrays = admitted
    with pytest.raises(ValueError):
        packet.validate(value, {'success': arrays['success']})
    altered = copy.deepcopy(value)
    altered['execution']['native_success_joint_probabilities'][0][0] += .01
    with pytest.raises(ValueError):
        packet.validate(altered, arrays)


def test_change_certified_bin_or_source_pin_rejected(admitted):
    value, arrays = admitted
    changed = copy.deepcopy(value)
    changed['cube_joint_bins']['rows'][0]['energy_bin_probability_intervals'][0] = ['0', '0']
    with pytest.raises(ValueError):
        packet.validate(changed, arrays)
    changed = copy.deepcopy(value)
    changed['source_bindings'][next(iter(changed['source_bindings']))] = '0' * 64
    with pytest.raises(ValueError):
        packet.validate(changed, arrays)


def test_source_snapshot_change_during_producer_rejected(monkeypatch):
    old = packet.bindings()
    calls = 0
    def changing():
        nonlocal calls
        calls += 1
        return old if calls == 1 else dict(old, forbidden='changed')
    monkeypatch.setattr(packet, 'bindings', changing)
    monkeypatch.setattr(packet.experiment, 'run', lambda: ({}, {}))
    with pytest.raises(ValueError, match='source changed'):
        packet.produce()


def test_warmed_cube_cache_rejects_float_clock(admitted):
    with pytest.raises(ValueError):
        spectral.cube_readout(4.0, F(1, 2))
    with pytest.raises(ValueError):
        spectral.cube_readout(4, .5)


def test_returned_cube_certificate_cannot_poison_cache(admitted):
    value = spectral.cube_readout()
    value['relativistic_mass_identified'] = True
    assert spectral.cube_readout()['relativistic_mass_identified'] is False
