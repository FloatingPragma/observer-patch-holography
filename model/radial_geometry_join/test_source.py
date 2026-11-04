from fractions import Fraction as F

import numpy as np
import pytest

from . import source, bounds


def annihilator(mode, modes=6):
    value = np.zeros((2**modes, 2**modes))
    for mask in range(2**modes):
        if mask >> mode & 1:
            value[mask ^ (1 << mode), mask] = (-1)**(mask & ((1 << mode)-1)).bit_count()
    return value


def test_complete_reference_coupled_native_energy_instrument():
    result = source.quantum_execution()
    assert result['native_pulses'] == 40
    assert result['native_first_failure_slots'] == 146
    assert result['retained_quantum_reference_columns'] == 4
    assert result['success_postselected_or_renormalized'] is False
    assert F(result['nine_coin_attempt_trace_distance_upper']) < F(1, 10**9)


def test_full_Fock_work_and_charge_identities_with_cross_parity_coherences():
    modes = [annihilator(i) for i in range(6)]
    number = [a.T@a for a in modes]
    identity = np.eye(64)
    parity = (number[2]+number[3]-identity)@(number[2]+number[3]-identity)
    # Independent three-site radial chain, plus centered interaction.
    sigma = np.array([[0, 1], [1, 0]])
    one_body = np.zeros((6, 6), complex)
    for j, x in enumerate((-1, 0, 1)):
        if x:
            one_body[2*j:2*j+2, 2*j:2*j+2] = 4*(-1)**x*np.sign(x)*sigma
    star = np.zeros((64, 64), complex)
    for left in (0, 1):
        for a in (0, 1):
            for b in (0, 1):
                term = .5j*sigma[a, b]*modes[2*left+a].T@modes[2*(left+1)+b]
                star += term+term.conj().T
    h = star+sum(one_body[i, j]*modes[i].T@modes[j]
                 for i in range(6) for j in range(6))
    h += .1*sum((number[2*j]-.5*identity)@(number[2*j+1]-.5*identity) for j in range(3))
    ks = source.kraus(parity)
    dual = lambda a: sum(k.conj().T@a@k for k in ks)
    assert np.linalg.norm(dual(identity)-identity) < 1e-12
    assert np.linalg.norm(dual(sum(number))-sum(number)) < 1e-12
    assert np.linalg.norm(dual(parity)-parity) < 1e-12
    assert np.linalg.norm(dual(h)-h+18*star/25) < 1e-12
    assert np.linalg.norm(h@parity-parity@h) > 1
    rng = np.random.default_rng(431)
    vector = rng.normal(size=64)+1j*rng.normal(size=64)
    vector /= np.linalg.norm(vector)
    rho = np.outer(vector, vector.conj())
    output = sum(k@rho@k.conj().T for k in ks)
    dephased = parity@rho@parity+(identity-parity)@rho@(identity-parity)
    assert np.linalg.norm(output-(7*rho+18*dephased)/25) < 1e-12
    assert abs(np.trace(ks[1]@rho@ks[1].conj().T)-(.02+.96*np.trace(parity@rho))) < 1e-12
    effects = {}
    energies = {}
    for charge in (2, 3, 4):
        indices = [mask for mask in range(64) if mask.bit_count() == charge]
        energy, eigenvectors = np.linalg.eigh(h[np.ix_(indices, indices)])
        rank = 2 if charge == 3 else 1
        low = eigenvectors[:, :rank]
        effects[charge] = np.trace(low.conj().T@parity[np.ix_(indices, indices)]@low).real/rank
        energies[charge] = energy[0]
    assert abs(energies[2]-energies[4]) < 1e-12
    assert abs(effects[2]-effects[4]) < 1e-12
    assert min(effects[2], effects[4])-effects[3] > .5
    assert energies[2] > energies[3]


def test_same_radial_states_give_positive_geometry_contrast():
    from model.matter_geometry_response.geometry import enclosed_gram
    value = bounds.joined_bounds()
    e = value['ground_D0_expectation_interval'][1]
    assert e < F(1, 1000)
    assert value['actual_charged_gap_lower'] > 0
    assert value['ground_Gram_G01_minus_each_charged_Gram_G01_strict_lower'] > F(59, 1000)
    # Independent rational endpoint controls at the least-separated admitted
    # values; the uniform Lipschitz proof covers the entire input intervals.
    ground, _ = enclosed_gram(F(3, 25)*(1-2*e))
    charged, _ = enclosed_gram(F(3, 25)*(1-2*(1-e)))
    assert ground[0][1][0]-charged[0][1][1] > F(59, 1000)
    assert source.mean_matrix(e) != source.mean_matrix(1-e)


def test_observation_bound_retains_complete_control_error():
    value = bounds.observation_bound()
    assert value['fresh_attempts_per_preparation'] == 48000000
    assert value['inferred_ideal_Gram_contrast_strict_lower'] == F(1224999, 25000000)
    assert value['all_failure_records_retained_with_declared_binary_score_zero']


@pytest.mark.parametrize('value', [-1, 2, .5, True, F(-1), F(2)])
def test_reject_nonprobability_or_inexact_mean_parameters(value):
    with pytest.raises(ValueError):
        source.mean_matrix(value)


@pytest.mark.parametrize('wires', [(0, 0, 2), (False, 1, 2), (0, 1, -1)])
def test_reject_invalid_parity_wiring(wires):
    with pytest.raises(ValueError):
        source.parity_copy(*wires)


@pytest.mark.parametrize('operator', [np.diag([1., .5]),
    np.array([[0., 1.], [0., 0.]]), np.ma.array(np.eye(2), mask=[[0, 1], [0, 0]])])
def test_reject_nonprojector_or_masked_field_operator(operator):
    with pytest.raises(ValueError):
        source.kraus(operator)
