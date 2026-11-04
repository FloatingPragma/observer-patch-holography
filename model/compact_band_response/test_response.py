from copy import deepcopy
from fractions import Fraction as F

import numpy as np
import pytest
from scipy import linalg as la

from . import response


def test_complete_source_projector_and_energy_derivatives():
    result = response.source_control()
    assert result['all_complement_modes_retained'] == 682
    assert float(result['projector_derivative_norm']) > .02
    assert F(result['rational_uniform_bounds']['spectral_projector_derivative_norm_upper']) == F(3, 29)
    assert result['hopping_energy_derivatives'][0].startswith('-0.0078')


def test_both_spectral_halves_are_required_for_complete_derivative():
    row = response.moving_band()
    energy, vector = la.eigh(row['h'])
    high = energy > 32
    q = vector[:, high]
    numerator = q.conj().T@row['T']@row['W']
    x = q@(numerator/(row['energies'][None, :]-energy[high, None]))
    positive_only = x@row['W'].conj().T+row['W']@x.conj().T
    assert la.norm(positive_only-row['dP']) > .02
    assert la.norm(row['C']@positive_only.conjugate()-positive_only@row['C']) > .02


def test_exact_force_detects_frozen_projector_deletion():
    control = response.force_control()
    assert control['moving_quartic_force_rank'] == 4
    assert control['moving_quartic_force_squared_trace'] == '4'
    assert control['freezing_the_localization_projector_omits_nonzero_force']


def test_thermal_readout_retains_work_and_spatial_splitting():
    result = response.thermal_control()
    work = float(result['mean_hopping_force'])
    energy = float(result['energy_derivative'])
    entropy = float(result['entropy_derivative'])
    assert work < -.002
    assert abs(entropy-64*energy) > .1
    assert abs(entropy-64*(energy-work)) < 1e-9
    assert not result['declared_reducing_slice_is_full_180_mode_Gibbs_state']


@pytest.mark.parametrize('parameter', [.499, 1.501, float('nan'), float('inf')])
def test_invalid_source_drive_is_rejected(parameter):
    with pytest.raises(ValueError):
        response.moving_band(parameter)


@pytest.mark.parametrize('mass', [7, 0, True, 8.5])
def test_gap_bound_requires_its_declared_domain(mass):
    with pytest.raises(ValueError):
        response.bounds(mass)


def test_receipt_rejects_scope_promotion_and_force_mutation():
    from . import packet
    replay = packet.produce()
    for field in ('conserved_relativistic_stress_constructed',
                  'regional_area_entropy_or_Einstein_stationarity_derived'):
        mutated = deepcopy(replay)
        mutated['scope'][field] = True
        with pytest.raises(ValueError):
            packet.validate(mutated, replay)
    mutated = deepcopy(replay)
    mutated['exact_missing_force_control']['moving_quartic_force_squared_trace'] = '0'
    with pytest.raises(ValueError):
        packet.validate(mutated, replay)
