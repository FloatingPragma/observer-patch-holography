"""Independent full-CAR, charge and source accounting checks."""
import math

import numpy as np
import pytest

from ..defect_field_bridge import joint
from . import balance, energy, packet, readout


@pytest.mark.parametrize('seed', [3, 1140, 2026])
def test_centered_quadratic_commutator_includes_all_scalar_terms(seed):
    rng = np.random.default_rng(seed)
    values = []
    for _ in range(2):
        h = rng.normal(size=(3, 3))+1j*rng.normal(size=(3, 3)); h += h.conj().T
        d = rng.normal(size=(3, 3))+1j*rng.normal(size=(3, 3)); d -= d.T
        values.append(energy.nambu(h, d))
    A, B = values
    a, b = [energy.quadratic_operator(x) for x in values]
    matrix = 1j*(A@B-B@A)
    assert np.linalg.norm(1j*(a@b-b@a)-energy.quadratic_operator(matrix)) < 1e-11
    assert abs(energy.normal_order_data(matrix)[2]) > .01


def test_same_small_oracle_exact_local_sum_derivative_and_zero_mutual_current():
    value = energy.small_energy_oracle()
    assert value['H'].shape == (16, 16)
    assert np.array_equal(value['H'], sum(value['h_sites']))
    assert np.array_equal(value['dH_dphi'], sum(value['dh_sites_dphi']))
    assert np.linalg.norm(value['J'][0, 1]) < 1e-15
    assert np.linalg.norm(value['H']@joint.charge()-joint.charge()@value['H']) < 1e-14


def test_three_site_all_sector_continuity_matches_independent_JW_matrices():
    data = energy.component_energy((27, 78), .37, sites=3, reference='common_static')
    H = energy.quadratic_operator(data['B'], data['reference_scalar'])
    local = [energy.quadratic_operator(B, data['local_reference_scalar']) for B in data['local_B']]
    assert np.linalg.norm(H-sum(local)) < 1e-12
    for x in range(3):
        incoming = np.zeros_like(H)
        for y in range(3):
            current = energy.quadratic_operator(energy.energy_current(data, y, x))
            assert np.linalg.norm(current-1j*(local[y]@local[x]-local[x]@local[y])) < 1e-12
            incoming += current
        assert np.linalg.norm(incoming-1j*(H@local[x]-local[x]@H)) < 1e-12


def test_reference_is_explicit_and_does_not_change_currents():
    bare = energy.component_energy((27, 78), reference='bare_propagation')
    shifted = energy.component_energy((27, 78), reference='common_static')
    assert bare['reference_scalar'] == 0
    assert shifted['reference_scalar'] > 0
    assert np.array_equal(bare['B'], shifted['B'])
    assert np.array_equal(energy.energy_current(bare, 0, 2), energy.energy_current(shifted, 0, 2))


def test_full_ring_current_has_exact_endpoint_operator_support_and_nonzero_readout():
    current = energy.ring_current_oracle()['operator']
    assert np.linalg.norm(current) > .2
    state, _ = packet.prepare_probe('current')
    assert abs(balance.expectation(current, state)-1/16) < 2e-9
    assert np.linalg.norm(current@joint.charge()-joint.charge()@current) < 1e-14


def test_individual_Pauli_pointer_is_not_a_neutral_energy_detector():
    # These terms cancel their charge-changing pieces only after summation.
    P = readout.pauli_matrix('IXXI')
    assert np.linalg.norm(P@joint.charge()-joint.charge()@P, 2) > 1


def test_every_grouped_source_instrument_preserves_charge_and_complete_effects():
    visited = set()
    for definition in readout.exact_quadratic_terms().values():
        for term in definition['terms']:
            key = (term['kind'], term['modes'], term['coefficient'])
            if key in visited:
                continue
            visited.add(key)
            pointer = readout.grouped_pointer(*key)
            assert pointer['unitary'].shape == (64, 64)
            assert np.linalg.norm(sum(k.conj().T@k for k in pointer['K'])-np.eye(16)) < 2e-7
            assert max(np.linalg.norm(k@joint.charge()-joint.charge()@k) for k in pointer['K']) < 2e-7
            reconstructed = sum(value*K.conj().T@K for value, K in zip(pointer['eigenvalues'], pointer['K'], strict=True))
            assert np.linalg.norm(reconstructed-pointer['operator']) < 2e-7
    assert len(visited) >= 10


@pytest.mark.parametrize('observable', ['H', 'h0', 'h1', 'drive_power_per_phase_rate', 'ring_J0to2'])
def test_grouped_measurements_reconstruct_full_observable_expectation(observable):
    value = energy.small_energy_oracle()
    matrices = {'H': value['H'], 'h0': value['h_sites'][0], 'h1': value['h_sites'][1],
                'drive_power_per_phase_rate': value['dH_dphi'], 'ring_J0to2': energy.ring_current_oracle()['operator']}
    state, _ = packet.prepare_probe('repair_sensitive')
    terms = readout.exact_quadratic_terms()[observable]
    measured = terms['constant']*np.linalg.norm(state)**2
    for term in terms['terms']:
        measured += readout.observe_group(term, state)['expectation']
    assert abs(measured-balance.expectation(matrices[observable], state)) < 2e-8


@pytest.mark.parametrize('policy,law,expected', [('coherent', 'C', .25), ('coherent', 'E', .5), ('flagged', 'C', .5), ('flagged', 'E', .5)])
def test_complete_repair_energy_transfer_is_not_assumed_zero(policy, law, expected):
    state, _ = packet.prepare_probe('repair_sensitive')
    K, _, _ = joint.repair(law, policy)
    output = np.concatenate([k@state for k in K], axis=1)
    H = energy.small_energy_oracle()['H']
    assert output.shape == (16, 128)
    change = balance.expectation(H, output)-balance.expectation(H, state)
    assert abs(change-expected) < 2e-8


@pytest.mark.parametrize('current', [-1, 0, 1])
def test_source_ordered_control_switching_work_and_full_record_balance(current):
    result = balance.protocol_balance(current=current)['report']
    assert len(result['complete_native_records']) == 32
    assert result['effective_duration'] == 1+current/4
    assert max(map(abs, result['segment_energy_drift'])) < 2e-8
    assert abs(sum(result['switching_work'])-result['reference_energy_change_before_repair']) < 2e-8
    assert abs(result['reference_energy_change_before_repair']+result['complete_repair_energy_transfer']
               -result['total_reference_energy_change']) < 2e-10


def test_invalid_energy_domains_rejected():
    with pytest.raises(ValueError):
        energy.component_energy((27, 78), reference='implicit')
    with pytest.raises(ValueError):
        energy.small_energy_oracle(float('nan'))
    with pytest.raises(ValueError):
        energy.nambu(np.zeros((2, 2)), np.eye(2))
