import copy
import json
import subprocess
import sys
from fractions import Fraction as F

import numpy as np
import pytest
import sympy as s
from scipy.linalg import expm

from . import algebra, bounds, packet, polarization, readout, state_bounds


def test_exact_spin_algebra_and_actual_source_frame():
    assert len(algebra.local_certificate()['exact_two_site_axis_checks']) == 3
    assert all(row['full_selected_field_matrix_matches'] for row in algebra.source_frame_certificate())
    assert len(algebra.source_force_certificate()['source_force_rows']) == 3


@pytest.mark.parametrize('axis',range(3))
def test_full_three_site_rotation_and_remainder(axis):
    spin = [[np.array(a,dtype=complex) for a in site] for site in algebra.spins(3)]
    h = sum((sum(spin[x][a]@spin[y][a] for a in range(3))-np.eye(8)/4
             for x,y in ((0,1),(1,2))),np.zeros((8,8),complex))
    theta = (F(1,8),F(-1,11),F(2,15))
    generator = sum((float(theta[x])*spin[x][axis] for x in range(3)),np.zeros((8,8),complex))
    unitary = expm(1j*generator)
    first = 1j*(h@generator-generator@h)
    second = generator@(h@generator-generator@h)-(h@generator-generator@h)@generator
    p = F(1,7)
    remainder = (unitary.conj().T@h@unitary-h-first-second/2)/float(p)
    admitted = bounds.twist_remainder(p,[theta[0]-theta[1],theta[1]-theta[2]])
    assert np.linalg.norm(remainder,2) <= float(admitted)+2e-15
    # Independent trigonometric expression checks the current orientation.
    rotated = h.copy()
    for x,y in ((0,1),(1,2)):
        difference = float(theta[x]-theta[y])
        transverse = sum(spin[x][a]@spin[y][a] for a in range(3) if a != axis)
        current = sum((int(s.LeviCivita(axis,b,c))*spin[x][b]@spin[y][c]
                       for b in range(3) for c in range(3)),np.zeros((8,8),complex))
        rotated += (np.cos(difference)-1)*transverse-np.sin(difference)*current
    assert np.linalg.norm(rotated-unitary.conj().T@h@unitary,2) < 1e-14


def test_uniform_global_rotation_has_zero_exchange_remainder():
    for axis in range(3):
        spin = algebra.spins(2)
        h = sum((spin[0][a]*spin[1][a] for a in range(3)),s.zeros(4))-s.eye(4)/4
        assert algebra.comm(h,spin[0][axis]+spin[1][axis]) == s.zeros(4)
    assert bounds.twist_remainder(F(1,1000),[0,0]) == 0


def test_selection_field_is_not_removed_from_finite_operator_identity():
    spin = algebra.spins(2)
    selection = -(spin[0][0]-spin[1][0])/4
    generator = spin[0][1]+spin[1][1]
    assert algebra.comm(selection,generator) != s.zeros(4)


def test_exact_quadratic_scale_and_state_independent_remainders():
    for axis in ('P','Q'):
        coarse = bounds.compact_twist(F(1,8),1,1,axis)
        fine = bounds.compact_twist(F(1,16),1,1,axis)
        assert coarse['scaled_BCH_remainder_norm_upper'] == 4*fine['scaled_BCH_remainder_norm_upper']
        assert coarse['quadratic_coefficient_norm_upper'] == fine['quadratic_coefficient_norm_upper']
        assert coarse['bound_is_uniform_over_all_source_states_and_references']


def test_bond_double_commutator_bound_on_actual_noncommuting_path():
    spin = algebra.spins(3)
    edges = ((0,1),(1,2))
    h = sum((sum((spin[x][a]*spin[y][a] for a in range(3)),s.zeros(8))-s.eye(8)/4
             for x,y in edges),s.zeros(8))
    values = (F(1,3),F(-2,7))
    z = sum((s.Rational(value.numerator,value.denominator)*
             sum((spin[x][a]*spin[y][a] for a in (0,2)),s.zeros(8))
             for value,(x,y) in zip(values,edges)),s.zeros(8))
    double = algebra.comm(z,algebra.comm(h,z))
    exact = bounds.bond_kick(2,max(map(abs,values)),sum(map(abs,values)))
    assert double != s.zeros(8)
    assert np.linalg.norm(np.array(double,dtype=complex),2) <= float(exact['double_commutator_norm_over_J_upper'])
    assert bounds.bond_kick(6,F(1,1000),1)['local_overlap_count_upper'] == 198


def test_bond_phase_transport_has_a_vanishing_bound():
    coarse = bounds.phase_transport(F(1,4096),1,1,1,1)
    fine = bounds.phase_transport(F(1,8192),1,1,1,1)
    assert coarse['phase_transport_squared_norm_upper'] == coarse['uncapped_squared_norm_upper']
    assert 0 < fine['uncapped_squared_norm_upper'] < coarse['uncapped_squared_norm_upper']/2
    assert bounds.phase_transport(1,0,1,1,1)['phase_transport_squared_norm_upper'] == 0


def test_complete_local_reset_and_ground_density_bounds():
    row = algebra.reset_energy_certificate()
    assert row['complete_reset_Kraus_count'] == 4
    assert row['ground_stability_identity_exact']
    assert row['boundary_channel_change_spectrum_over_J'] == ['-1/4','3/4']
    for side in (1,2,4,16,64):
        bound = bounds.ground_energy_cube(side)
        assert bound['normalized_ground_energy_fluctuation_squared_over_J2_upper'] == F(18*(side-1),side**2)
        assert bound['finite_ground_density_above_thermodynamic_density_over_J_upper'] == F(3,side)
        assert not bound['ergodicity_assumed']
        assert not bound['individual_transverse_bond_limits_proved_scalar']


def test_complete_coefficient_readout_exact_effects_and_backaction():
    spin = algebra.spins(3)
    edges = ((0,1),(1,2))
    coefficients = (F(1,3),F(-1,5))
    plan = readout.plan(coefficients,2)
    effect = s.zeros(8)
    score_effect = s.zeros(8)
    for row in plan['complete_fine_records']:
        x,y = edges[row['edge_index']]
        axis = row['Pauli_axis']-1
        projector = (s.eye(8)/2+row['left_record']*spin[x][axis])*(s.eye(8)/2+row['right_record']*spin[y][axis])
        probability = s.Rational(row['setting_probability'])
        effect += probability*projector
        score_effect += probability*s.Rational(row['score'])*projector
    target = sum((s.Rational(value.numerator,value.denominator)*
                  sum((spin[x][a]*spin[y][a] for a in (0,2)),s.zeros(8))
                  for value,(x,y) in zip(coefficients,edges)),s.zeros(8))
    assert effect == s.eye(8)
    assert score_effect == target
    assert len(plan['complete_fine_records']) == 16
    assert not plan['score_is_spectral_measurement_of_entire_sum']
    # Tensoring the exact effect identities with I proves untouched-reference
    # preservation; this checks an entangled non-product input explicitly.
    vector = np.arange(1,17,dtype=float)+1j*np.arange(16,0,-1,dtype=float)
    vector /= np.linalg.norm(vector)
    total = np.vdot(vector,np.kron(np.array(effect,dtype=complex),np.eye(2))@vector).real
    assert abs(total-1) < 1e-14


def test_scalar_quadratic_control_does_not_force_one_frequency():
    row = algebra.hidden_force_control()
    assert row['observed_response_frequency_squares'] == ['1/2','3/2']
    assert row['both_observed_response_residues_are_nonzero']
    assert not row['model_is_the_ordered_spin_source']


def test_zero_coefficient_has_complete_identity_readout():
    for values in ((),(0,0)):
        row = readout.plan(values,2)
        assert row['complete_fine_records'] == [{'record':'identity','score':'0'}]
        assert row['independent_repetitions_for_Chebyshev_bound'] == 0
        assert not row['native_coherent_control_and_local_Pauli_pointer_services_required']


def test_actual_scalar_polarization_remainder_has_zero_ergodic_weight():
    row = polarization.exact_checks()
    assert row['gamma_axis2_pi_reverses_transverse_axis3_current']
    assert row['compact_test_derivative_integral'] == 0
    assert all(item['exact_weighted_difference_quotient_sum'] == '0'
               for item in row['exact_mesh_zero_sum_controls'])
    bound = polarization.compact_bounds(F(1,4096),1,1,1)
    assert bound['current_coefficient_maximum'] == F(1,2*4096**3)
    assert bound['difference_quotient_to_derivative_operator_error'] == F(27,4*4096)
    assert not bound['effective_ground_vector_convergence_modulus']


def test_actual_product_input_excludes_operator_norm_residual_limit():
    row = polarization.product_source_control()
    assert F(row['residual_expectation_limit']) == F(16384,99225)
    assert all(F(item['exact_actual_source_product_expectation']) > F(1,10)
               for item in row['exact_mesh_controls'])
    assert not row['source_residual_converges_in_operator_norm']


def test_selected_ground_boundary_second_moment_and_tail():
    for side in (1,2,3,4,16,64):
        row = state_bounds.regional_second_moment(side)
        assert row['selected_ground_equation_density_second_moment_over_J2_upper'] == F(36,side**2)
        assert row['best_density_second_moment_over_J2_upper'] <= F(36,side**2)
        assert row['volume_limit_precedes_zero_selection_field']
        threshold = F(12*side**2)
        assert state_bounds.excitation_tail(side,threshold) <= F(1,4)
    with pytest.raises(ValueError):
        state_bounds.excitation_tail(4,0)


def test_same_source_finite_ground_split_has_bounded_crossing_perturbation():
    from model.native_ordered_ground import source
    whole,scale = source.matrix(3,((0,1),(1,2)),F(1),F(1,4))
    inside,sinside = source.matrix(2,((0,1),),F(1),F(1,4))
    outside,soutside = source.matrix(1,(),F(1),F(1,4))
    h = whole/scale
    a = np.kron(np.eye(2),inside/sinside)
    b = np.kron(outside/soutside,np.eye(4))
    crossing = h-a-b
    assert np.array_equal(crossing@crossing,-crossing)
    assert np.array_equal(a@b,b@a)
    energies,vectors = np.linalg.eigh(h)
    ea,eb = np.linalg.eigvalsh(a)[0],np.linalg.eigvalsh(b)[0]
    vector = vectors[:,0]
    delta = a-ea*np.eye(8)
    assert -1-1e-14 <= energies[0]-ea-eb <= 1e-14
    assert np.linalg.norm(delta@vector) <= 1+1e-14
    assert np.linalg.norm((a+b-(ea+eb)*np.eye(8))@vector) <= 1+1e-14


def test_source_energy_phase_control_keeps_different_selected_phase_explicit():
    row = state_bounds.source_phase_control()
    assert row['source_Hamiltonian_and_all_its_energy_statistics_equal']
    assert row['one_site_half_trace_distance_strict_lower'] == '11/50'
    assert not row['both_states_are_the_same_plus_field_selected_state']


def test_actual_finite_gap_bound_does_not_prove_regional_purity():
    for side in (2,3,4):
        row = state_bounds.source_gap_control(side)
        assert row['regional_excited_subspace_probability_upper_from_boundary_and_gap'] == '1'
        assert not row['this_bound_certifies_a_pure_regional_ground']


def test_full_pauli_tomography_bound_uses_all_observables():
    assert state_bounds.regional_tomography_bound((F(1,3),F(1,4),0)) == F(25,576)
    assert state_bounds.regional_tomography_bound((2,2,2)) == 1
    assert state_bounds.regional_tomography_bound(()) == 0


@pytest.mark.parametrize('args',[(0,1,1,'P'),(F(1,2),0,1,'P'),(F(1,2),1,-1,'Q'),(1,1,1,'X')])
def test_invalid_scale_or_test_contract_is_rejected(args):
    with pytest.raises(ValueError):
        bounds.compact_twist(*args)


@pytest.fixture(scope='module')
def receipt():
    return packet.build()


@pytest.mark.parametrize('kind',['scope','outcome','bound','source'])
def test_scientific_payload_and_scope_mutations_fail(receipt,kind):
    value = copy.deepcopy(receipt)
    if kind == 'scope':
        value['scope']['sharp_dispersion_or_Lorentz_covariance_proved'] = True
    elif kind == 'outcome':
        value['complete_coefficient_readout']['complete_fine_records'].pop()
    elif kind == 'bound':
        value['compact_test_controls'][0]['scaled_BCH_remainder_norm_upper'] = '0'
    else:
        value['source_bindings'][next(iter(value['source_bindings']))] = '0'*64
    with pytest.raises(ValueError):
        packet.validate(value)


def test_actual_runtime_source_imports_are_bound():
    script='''
import json,sys
from pathlib import Path
from model.spin_twist_closure import packet
packet.build(); root=packet.ROOT.resolve(); local=set()
for module in tuple(sys.modules.values()):
    name=getattr(module,'__file__',None)
    if name:
        path=Path(name).resolve()
        if path.suffix=='.py' and path.is_relative_to(root):local.add(path.relative_to(root).as_posix())
assert 'model/native_ordered_ground/source.py' in local
print(json.dumps(sorted(local-set(packet.bindings()))))
'''
    result = subprocess.run([sys.executable,'-W','error','-c',script],cwd=packet.ROOT,
                            capture_output=True,text=True,check=True)
    assert json.loads(result.stdout) == []


def test_check_command_uses_the_requested_file(tmp_path,receipt):
    path = tmp_path/'requested.json'
    path.write_text(packet.canonical(receipt))
    args = [sys.executable,'-W','error','-m','model.spin_twist_closure','--check',str(path)]
    result = subprocess.run(args,cwd=packet.ROOT,capture_output=True,text=True)
    assert result.returncode == 0,result.stdout+result.stderr
    path.write_text('{"schema":"mutated"}')
    result = subprocess.run(args,cwd=packet.ROOT,capture_output=True,text=True)
    assert result.returncode != 0
