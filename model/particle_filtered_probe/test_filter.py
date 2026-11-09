"""Independent circuit, spectral and complete-reference checks."""
from copy import deepcopy
from fractions import Fraction as F

import numpy as np
import pytest
from scipy.linalg import expm
from scipy.integrate import quad
import sympy as s

from model.particle_common_contract import contract
from model.particle_physical_boundary_state import boundary
from model.spin_operational_position.bounds import generic_instrument_plan
from . import filter as f


@pytest.fixture(scope='module')
def control():
    return f.finite_control()


@pytest.fixture(scope='module')
def member():
    return contract.member(0)


def array(matrix):
    return np.array(matrix.evalf(18).tolist(), dtype=complex)


def test_declared_proof_and_runner_inputs_exist():
    assert all((f.ROOT/path).is_file() for path in f.SOURCE_INPUTS)


def test_positive_cell_integrals_are_actual_triangle_quadrature():
    for m in (2, 4, 10):
        q = f.quadrature('3/2', m)
        for row in q['rows']:
            lo, hi = row['interval']
            actual = quad(lambda x: (1-abs(x)/1.5)/1.5, float(lo), float(hi))[0]
            assert abs(actual-float(row['weight'])) < 1e-14
        assert sum(row['weight'] for row in q['rows']) == 1
    with pytest.raises(ValueError):
        f.quadrature(cells=3)


def test_actual_native_boundary_probe_is_a_complete_even_odd_unitary():
    native = boundary.native_degree_one_control()
    a = sum((s.Integer(-1)**q*p for q, p in native['projectors']), s.zeros(119))
    assert a == a.H and a*a == s.eye(119)
    assert set(a.diagonal()) == {-1, 1}
    assert sum(native['outcome_ranks'].values()) == 119


def test_full_lcu_matrix_executes_the_blank_success_and_all_failures(control):
    data = control['LCU']; w = array(data['unitary'])
    assert np.linalg.norm(w.conj().T@w-np.eye(12)) < 1e-13
    assert np.linalg.norm(w[:3, :3]-array(data['success'])) < 1e-14
    rho = array(control['rho'])
    probabilities = [np.trace(array(k)@rho@array(k).conj().T).real for k in data['kraus']]
    assert min(probabilities) >= -1e-14 and abs(sum(probabilities)-1) < 1e-14
    assert sum(probabilities[1:]) > 0
    # A classical mixture of the unitary nodes never writes a filter failure.
    classical_success = sum(float(row['weight'])*np.trace(array(u)@rho@array(u).conj().T).real
                            for row, u in zip(data['quadrature']['rows'], data['node_unitaries']))
    assert abs(classical_success-1) < 1e-14
    assert abs(classical_success-probabilities[0]) > 1e-5


def test_full_ramsey_statevector_includes_the_untouched_reference(control):
    rho = array(control['rho']); eigenvalues, eigenvectors = np.linalg.eigh(rho)
    root = (eigenvectors*np.sqrt(eigenvalues))@eigenvectors.conj().T
    psi = root.reshape(9)
    w = array(control['LCU']['unitary']); u = array(f.cayley(control['H'], '1', 16))
    echoed = np.kron(np.eye(4), u.conj().T)@w@np.kron(np.eye(4), u)
    initial = np.kron(np.eye(4)[0], psi)
    left = np.kron(w, np.eye(3))@initial
    right = np.kron(echoed, np.eye(3))@initial
    for imaginary in (False, True):
        phase = -1j if imaginary else 1
        joint = np.stack(((left+phase*right)/2, (left-phase*right)/2))
        # Layout: Ramsey bit, selector, source, reference.
        joint = joint.reshape(2, 4, 3, 3)
        probabilities = (np.abs(joint)**2).sum(axis=(2, 3))
        assert abs(probabilities.sum()-1) < 1e-13
        reference = np.einsum('brsi,brsj->ij', joint, joint.conj())
        assert np.linalg.norm(reference-rho.T) < 1e-13
        rows = f.ramsey_instrument(control['H'], control['LCU']['kraus'], imaginary=imaginary)
        for row in rows:
            exact_p = s.simplify(s.trace(control['rho']*row['K'].H*row['K']))
            assert abs(float(exact_p)-probabilities[row['Ramsey_bit'],row['filter_label']]) < 1e-13
        score = probabilities[0,0]-probabilities[1,0]
        emitted = control['receipt']['exact_original_ImF' if imaginary else 'exact_original_ReF']
        assert abs(score-float(F(emitted))) < 1e-13
        # Conditioning would recover F/F0, which is a different observable.
        conditional = score/probabilities[:,0].sum()
        assert abs(conditional-score) > 1e-5


def test_exact_projector_sum_and_continuous_time_integral(control):
    h, rho, a = map(array, (control['H'],control['rho'],control['A']))
    values, axes = np.linalg.eigh(h)
    r = axes.conj().T@rho@axes; aa = axes.conj().T@a@axes
    spectral_f0 = 0.; spectral_f = 0j
    for i, initial in enumerate(values):
        for j, final in enumerate(values):
            omega = final-initial
            original = r[i,i].real*abs(aa[j,i])**2
            multiplier = np.sinc(omega/32/np.pi)**4
            spectral_f0 += original*multiplier
            spectral_f += original*multiplier*np.exp(1j*omega)
    t = 1/16
    integral = np.zeros((3,3),complex)
    for i in range(3):
        for j in range(3):
            def entry(x):
                u = expm(-1j*x*h)
                return (1-abs(x)/t)/t*(u.conj().T@a@u)[i,j]
            integral[i,j] = quad(lambda x: entry(x).real,-t,t,points=[0],epsabs=1e-13)[0]+1j*quad(lambda x: entry(x).imag,-t,t,points=[0],epsabs=1e-13)[0]
    exact_f0 = np.trace(rho@integral.conj().T@integral).real
    evolution = expm(-1j*h)
    exact_f = np.trace(rho@integral.conj().T@evolution.conj().T@integral@evolution)
    assert abs(exact_f0-spectral_f0) < 1e-12
    assert abs(exact_f-spectral_f) < 1e-12
    receipt = control['receipt']; error = receipt['errors']
    assert abs(float(F(receipt['exact_original_F0']))-exact_f0) <= float(F(error['original_F0_error_to_continuous_triangle']))
    assert abs(float(F(receipt['exact_original_ReF']))-exact_f.real) <= float(F(error['original_ReF_or_ImF_error_to_continuous_triangle']))
    assert F(receipt['continuous_triangle_contrast_lower']) > F(1,4)
    band = receipt['positive_window_control']
    assert F(band['original_positive_band_weight_lower']) > F(1,16)
    actual_band_weight = r[0,0].real*abs(aa[1,0])**2*np.sinc(1/32/np.pi)**4
    assert float(F(band['original_positive_band_weight_lower'])) < actual_band_weight
    assert np.linalg.norm(array(control['LCU']['success'])-integral,2) <= float(F(error['filtered_probe_operator_error_upper']))


def test_one_sided_clock_and_original_energy_zero_are_not_silent_substitutions(control):
    a = control['A']; u = f.cayley(control['H'],'1',16)
    assert f.clean(u.H*a*u-u*a) != s.zeros(3)
    shift = s.Rational(7,3)
    # Exact Heisenberg evolution cancels scalar phases; a rational Cayley
    # approximation still has to be evaluated at the specified original H.
    shifted = f.cayley(control['H']+shift*s.eye(3),'1',16)
    assert f.clean(shifted.H*a*shifted-u.H*a*u) != s.zeros(3)


def test_primary_precision_schedule_and_ancestry_do_not_supply_survival(member):
    parameters = f.precision_parameters(member,response_error='1/1024')
    result = f.selected_plan(member,**parameters)
    assert F(result['errors']['original_contrast_error_to_continuous_triangle']) <= F(1,2048)
    assert result['clock_id'] == member['ids']['clock_id']
    assert result['reference_id'] == member['ids']['reference_id']
    assert not result['apparatus']['reference_actively_counter_evolved']
    assert not result['apparatus']['full_spatial_matrix_or_tape_executed']
    assert not result['primary_F0_or_Ramsey_contrast_computed']
    assert not result['finite_survival']['uniform_positive_survival_established']
    assert not result['calibration_or_learned_mass_inputs']
    receipt = f.bind_member(member,**parameters)
    damaged = deepcopy(receipt)
    damaged['result']['filter']['success_probability_renormalized'] = True
    damaged['result_sha256'] = contract.digest(damaged['result'])
    with pytest.raises(ValueError):
        f.verify_bound(member,damaged,**parameters)
    changed = deepcopy(member);changed['clock']['GeV_or_proper_time_conversion']='fit'
    with pytest.raises(ValueError):
        f.selected_plan(changed)


def test_factored_compiler_uses_the_inherited_all_input_bound():
    eps=F(1,1000000);old=generic_instrument_plan(2,8,eps)
    plan=f.factored_compiler_plan(2,8,eps)
    upper=34*2**plan['pulse_bound_factored_upper']['power_of_two']
    assert old['native_pulse_upper_bound'] <= upper
    assert F(29*upper,2**plan['port_fraction_bits']) <= eps/2
    q=plan['padded_unitary_bits']
    assert F(2*2**q*(q+2),2**plan['rank_one_tree_angle_bits']) <= eps/4


def test_contraction_dilation_retains_probe_failure_and_the_same_blank_block():
    h = s.diag(0,4); p = s.ones(2)/2
    result = f.contraction_instrument(h,p,'1/2',2,1)
    assert result['success'] == s.Matrix([[s.Rational(1,2),s.Rational(3,10)],
                                          [s.Rational(3,10),s.Rational(1,2)]])
    assert len(result['kraus']) == 4 and result['fine_labels'] == [(0,0),(0,1),(1,0),(1,1)]
    assert sum((k.H*k for k in result['kraus']),s.zeros(2)).applyfunc(s.simplify) == s.eye(2)
    assert len(f.ramsey_instrument(h,result['kraus'],'1/4',1)) == 8
    zero = f.contraction_instrument(h,s.zeros(2),'1/2',2,1)
    assert zero['success'] == s.zeros(2)
    assert sum((k.H*k for k in zero['kraus'][1:]),s.zeros(2)) == s.eye(2)


def test_filtered_lueders_effect_is_not_the_lcu_blank_effect():
    h = s.diag(0,4); p = s.ones(2)/2
    actual = f.filtered_lueders_instrument(h,p,'charge','1/2',2,1)
    e = actual['effect']; k0,k1 = actual['kraus']
    assert f.clean(k0.H*k0) == e
    assert f.clean(k0.H*k0+k1.H*k1) == s.eye(2)
    assert actual['LCU_blank_effect'] == e*e != e
    loop = f.filtered_lueders_instrument(h,2*p-s.eye(2),'plaquette','1/2',2,1)
    assert loop['effect'] == e
    # Execute the complete neutral dilation on a source/reference purification.
    rho = np.diag([.7,.3]); psi = np.sqrt(rho).reshape(4)
    state = np.kron(array(actual['unitary']),np.eye(2))@np.kron([1,0],psi)
    state = state.reshape(2,2,2)
    assert abs((abs(state)**2).sum()-1) < 1e-13
    assert np.linalg.norm(np.einsum('bsi,bsj->ij',state,state.conj())-rho.T) < 1e-13
    assert abs((abs(state[0])**2).sum()-np.trace(rho@array(e)).real) < 1e-13
    assert abs(np.trace(rho@array(e*e-e)).real) > .1
    history = f.compose_instruments(actual['kraus'],loop['kraus'])
    assert len(history) == 4
    assert f.clean(sum((k.H*k for k in history),s.zeros(2))) == s.eye(2)


def test_lueders_root_continuity_covers_complete_entangled_histories():
    delta = F(1,1024)
    exact = f.lueders_instrument(s.diag(0,1))
    perturbed = f.lueders_instrument(s.diag(s.Rational(1,1024),s.Rational(1023,1024)))
    psi = np.eye(2).reshape(4)/np.sqrt(2)
    def record_state(kraus):
        return np.array([np.outer(v,v.conj()) for k in kraus
                         for v in [np.kron(array(k),np.eye(2))@psi]])
    output_error = sum(np.linalg.norm(x,ord='nuc') for x in
                       record_state(exact['kraus'])-record_state(perturbed['kraus']))/2
    budget = f.lueders_error_budget(delta,history_length=2)
    assert output_error <= float(budget['one_complete_instrument_halftrace_error_upper'])
    first = f.compose_instruments(exact['kraus'],exact['kraus'])
    second = f.compose_instruments(perturbed['kraus'],perturbed['kraus'])
    history_error = sum(np.linalg.norm(x,ord='nuc') for x in record_state(first)-record_state(second))/2
    assert history_error <= float(budget['history_complete_halftrace_error_upper'])
    # The rational square-root estimate has no fixed absolute error floor.
    tiny = f.lueders_error_budget(F(1,2**400))
    assert tiny['one_complete_instrument_halftrace_error_upper'] < F(1,2**198)


def test_authenticated_field_menus_share_source_and_have_separate_lueders_plans(member):
    for kind in ('charge','plaquette'):
        descriptor = f.field_algebra.probe(member,kind=kind)
        parameters = f.precision_parameters(member,half_width=descriptor['half_width'])
        plan = f.selected_plan(member,probe_descriptor=descriptor,**parameters)
        assert plan['probe']['authenticated_menu'] == descriptor
        assert plan['filter']['Ramsey_outcomes_per_phase_setting'] == 4*parameters['cells']
        assert plan['finite_survival']['unfiltered_connected_variance_lower'] is None
        assert not plan['lueders_history_apparatus']['success_effect_equals_LCU_blank_effect']
        assert plan['lueders_history_apparatus']['positive_contraction_from_positive_quadrature']
        exact_parameters = {**parameters,'compiler_error':'0'}
        exact_plan = f.selected_plan(member,probe_descriptor=descriptor,**exact_parameters)
        assert not exact_plan['lueders_history_apparatus']['compiler']['compiled_tape_required']
        wrong = deepcopy(descriptor);wrong['clock_id'] = 'different'
        with pytest.raises(ValueError):
            f.selected_plan(member,probe_descriptor=wrong,**parameters)
        with pytest.raises(ValueError):
            f.selected_plan(member,probe_descriptor=descriptor,half_width='1/16')


def test_noncontractions_and_nonprojection_charge_menus_are_rejected():
    with pytest.raises(ValueError):
        f.hermitian_dilation(s.diag(1,2))
    with pytest.raises(ValueError):
        f.filtered_lueders_instrument(s.diag(0,1),s.eye(2)/2,'charge')
    with pytest.raises(ValueError):
        f.lueders_instrument(s.diag(-1,1))
    with pytest.raises(ValueError):
        f.compose_instruments((2*s.eye(2),),(s.eye(2)/2,))
    assert f.error_budget(0,0)['filtered_probe_operator_error_upper'] == 0


@pytest.mark.parametrize('bad', ['0','-1',True,0.25])
def test_exact_positive_filter_time_required(bad):
    with pytest.raises(ValueError):
        f.quadrature(bad)
