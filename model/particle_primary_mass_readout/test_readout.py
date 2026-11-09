"""Independent spectral coupling, complete-reference and interval controls."""
from copy import deepcopy
from fractions import Fraction as F

import numpy as np
import pytest
import sympy as s

from model.particle_common_contract import contract
from model.source_mass_shell import units
from . import readout as r


@pytest.fixture(scope='module')
def member():
    return contract.member(0)


def test_original_mixed_response_is_not_final_excitation():
    control = r.finite_control()
    rows = control['joint_energy_rows']
    active = [row for row in rows if row['weight']]
    assert len(active) == 1
    assert active[0]['transfer'] == 1 and active[0]['final_energy'] == 4
    assert active[0]['weight'] == F(1, 16)
    assert control['receipt']['exact_transport_cost'] == F(3, 16)
    bound = r.excitation_window((1, 1), F(1, 16), F(3, 16), 6)
    assert bound['original_final_excitation_weight_lower'] == F(1, 32)
    assert bound['inferred_final_excitation_interval'] == (1, 7)
    assert bound['rest_mass_identified'] is False


def test_complete_control_retains_the_original_idle_reference():
    data = r.finite_control()
    root = np.diag([np.sqrt(15/16), 1/4, 0])
    psi = root.reshape(9)
    outputs = [np.kron(np.array(k, dtype=complex), np.eye(3))@psi
               for k in (data['probe'], data['failure'])]
    probs = [np.vdot(x, x).real for x in outputs]
    assert np.allclose(probs, [1/16, 15/16], atol=1e-15)
    reference = sum(x.reshape(3, 3).T@x.reshape(3, 3).conj() for x in outputs)
    assert np.allclose(reference, root@root, atol=1e-15)
    h_joint = np.kron(np.diag([0, 3, 4]), np.eye(3))
    variance = np.vdot(psi, h_joint@h_joint@psi)-np.vdot(psi, h_joint@psi)**2
    assert abs(variance-135/256) < 1e-15


def test_transport_from_actual_rotated_projectors_including_negative_transfers():
    rotation = s.Matrix([[s.Rational(3, 5), -s.Rational(4, 5), 0],
                         [s.Rational(4, 5), s.Rational(3, 5), 0], [0, 0, 1]])
    energies = (-2, 1, 5)
    probabilities = (s.Rational(3, 4), s.Rational(1, 8), s.Rational(1, 8))
    h = rotation*s.diag(*energies)*rotation.H
    rho = rotation*s.diag(*probabilities)*rotation.H
    b = s.Matrix([[1, 2, -1], [0, -1, 1], [2, 0, 1]])/5
    assert sum(abs(x)**2 for x in b) < 1
    projectors = [rotation[:, i]*rotation[:, i].H for i in range(3)]
    masses = []
    transfer_first = final_first = cost = 0
    for i in range(3):
        for j in range(3):
            weight = s.trace(projectors[j]*b*projectors[i]*rho*projectors[i]*b.H)
            assert weight >= 0
            transfer = energies[j]-energies[i]
            final = energies[j]-energies[0]
            masses.append(weight)
            transfer_first += weight*transfer
            final_first += weight*final
            cost += weight*abs(final-transfer)
    epsilon = s.trace(rho*(h-energies[0]*s.eye(3)))
    assert sum(masses) == s.trace(rho*b.H*b)
    assert final_first-transfer_first == cost <= epsilon
    assert any(s.trace(projectors[j]*b*projectors[i]*rho*projectors[i]*b.H) > 0
               for i in range(3) for j in range(i))


def test_window_and_first_moment_bounds_keep_zero_response():
    result = r.transport_bound('1/16', 0, '1/16')
    assert result['coupling_transport_cost_upper'] == 0
    assert result['final_excitation_first_moment_upper'] == 0
    assert r.excitation_window((0, 0), 0, 1, 1, 0)['original_final_excitation_weight_lower'] == 0
    result = r.transport_bound('1/16', 1, '1/16')
    assert result['final_excitation_first_moment_upper'] == F(513, 16)
    assert result['success_probability_division'] is False


def test_energy_floor_is_only_needed_for_the_direct_absolute_route():
    assert r.absolute_to_excitation((4, 4))['final_excitation_interval'] is None
    bound = r.absolute_to_excitation((4, 5), ('-1/10', '1/10'))
    assert bound['final_excitation_interval'] == (F(39, 10), F(51, 10))
    shifted = r.absolute_to_excitation((14, 15), ('99/10', '101/10'))
    assert shifted['final_excitation_interval'] == bound['final_excitation_interval']
    with pytest.raises(ValueError):
        r.absolute_to_excitation((0, 1), (2, 3))


def test_transfer_momentum_needs_an_initial_momentum_allowance():
    result = r.momentum_from_transfer((0, 0), 1, '1/10')
    assert result['final_momentum_norm_interval'] == (0, 1)
    assert result['additional_original_weight_loss_upper'] == F(1, 10)
    assert result['initial_momentum_control_derived_from_cooling_alone'] is False


def test_conditional_shell_interval_encloses_all_monotone_corners():
    result = r.conditional_rest_interval((5, 6), (3, 4), ('9/10', '11/10'), '1/20')
    lo, hi = result['rest_energy_squared_interval']
    for energy in (5, 6):
        for momentum in (3, 4):
            for speed in (F(9, 10), F(11, 10)):
                for error in (-F(1, 20), F(1, 20)):
                    assert lo <= energy**2-speed**2*momentum**2+error <= hi
    mlo, mhi = result['rest_energy_interval']
    assert mlo*mlo <= lo and mhi*mhi >= hi
    exact = r.conditional_rest_interval((5, 5), (3, 3))
    assert exact['rest_energy_interval'] == (4, 4)
    assert exact['physical_momentum_or_shell_constructed'] is False
    with pytest.raises(ValueError):
        r.conditional_rest_interval((1, 1), (2, 2))


def test_primary_consumer_recomputes_source_and_keeps_missing_mass(member):
    data = r.common_contract_bound(member)
    r.verify_bound(member, data)
    payload = data['result']
    assert payload['transport_for_contraction']['coupling_transport_cost_upper'] == '1/16'
    assert payload['primary_rest_energy_interval'] is None
    assert payload['primary_physical_mass_GeV'] is None
    changed = deepcopy(data)
    changed['result']['primary_physical_mass_GeV'] = '1'
    with pytest.raises(ValueError):
        r.verify_bound(member, changed)


def test_supplied_primary_weight_has_explicit_binding_and_truth_boundary(member):
    binding = {'program': contract.PROGRAM, 'member_id': member['member_id'], **member['ids']}
    output = r.primary_excitation_bound(member, (1, 2), '1/2', 1,
        evidence_binding=binding, provenance=('conditional weight input for this test',))
    assert output['result']['original_final_excitation_weight_lower'] == '7/16'
    assert output['supplied_spectral_weight_truth_verified_by_arithmetic'] is False
    wrong = {**binding, 'clock_id': 'another source clock'}
    with pytest.raises(ValueError):
        r.primary_excitation_bound(member, (1, 2), '1/2', 1,
            evidence_binding=wrong, provenance=('test input',))


def test_declared_SI_units_are_available_without_empirical_clock_claim():
    by_length = r.declared_SI_comparison((4, 4), metres_per_source_coordinate=(1, 1))
    by_time = r.declared_SI_comparison((4, 4),
        seconds_per_source_time=(F(1, 299792458), F(1, 299792458)))
    assert by_length['energy_GeV_interval'] == by_time['energy_GeV_interval']
    doubled = r.declared_SI_comparison((4, 4), metres_per_source_coordinate=(2, 2))
    assert tuple(map(F, doubled['energy_GeV_interval'])) == tuple(
        F(x)/2 for x in by_length['energy_GeV_interval'])
    assert by_length['c_alone_fixes_absolute_scale'] is False
    assert by_length['mass_or_shell_identified_by_unit_declaration'] is False
    free = units.energy_conversion((4, 4), generator_id='test generator')
    assert free['energy_GeV_interval'] is None


def test_actual_finite_primary_band_transfers_without_probability_loss():
    from model.particle_primary_survival import survival
    member = contract.member(1)
    source_band = survival.selected_witness(member)
    result = r.primary_survival_excitation(member)
    assert result['inferred_final_excitation_interval'] == source_band['original_positive_band']
    assert result['original_final_excitation_weight_lower'] == source_band['original_positive_band_weight_lower']
    assert result['finite_filter_final_excitation_weight_lower'] == source_band['finite_apparatus_original_positive_band_weight_lower']
    assert result['additional_weight_loss'] == '0'
    assert result['positive_finite_excitation_weight_proved'] is True
    assert result['finite_filter_positive_excitation_weight_proved'] is True
    assert result['fixed_band_noncollapse_or_rest_mass_proved'] is False


@pytest.mark.parametrize('function,args', [
    (r.transport_bound, (-1,)), (r.transport_bound, (0.1,)),
    (r.excitation_window, ((1, 2), 2, 1, 1)),
    (r.excitation_window, ((1, 2), 0, 1, 0)),
    (r.momentum_from_transfer, ((0, 1), 1, 2)),
])
def test_invalid_exact_domains_are_rejected(function, args):
    with pytest.raises(ValueError):
        function(*args)
