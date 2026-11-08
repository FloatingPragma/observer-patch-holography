import copy
from fractions import Fraction as F
import json
import subprocess
import sys
import numpy as np
import pytest
from scipy.linalg import expm
from . import algebra, asymptotic, charge_memory, finite, packet, zero_field


@pytest.fixture(scope='module')
def result():
    return packet.build()


def contains(interval, value, tolerance=0):
    return float(F(interval[0]))-tolerance <= value <= float(F(interval[1]))+tolerance


def test_schur_memory_by_independent_complex_matrix():
    # Complex source and probe prevent a real-conjugation shortcut.
    rng = np.random.default_rng(314159)
    raw = rng.normal(size=(7, 7))+1j*rng.normal(size=(7, 7))
    k = raw.conj().T@raw
    a = rng.normal(size=7)+1j*rng.normal(size=7)
    a /= np.linalg.norm(a)
    p = np.outer(a, a.conj())
    q = np.eye(7)-p
    mean = np.vdot(a, k@a).real
    r = q@k@a
    for eta in (.25, 1., 4.):
        g = np.vdot(a, np.linalg.solve(k+eta*np.eye(7), a)).real
        memory = np.vdot(r, np.linalg.solve(q@k@q+eta*np.eye(7), r)).real
        assert abs(memory-(mean+eta-1/g)) < 1e-12
        assert 0 < memory <= np.vdot(r, r).real/eta


def test_time_domain_memory_with_complex_data():
    from scipy.integrate import quad
    k = np.array([[2., 1j, .3], [-1j, 3., .5j], [.3, -.5j, 4.]])
    a = np.array([1., 0., 0.], complex)
    mean = np.vdot(a, k@a).real
    r = k[1:, 0]
    d = k[1:, 1:]
    c = lambda t: np.vdot(a, expm(-1j*t*k)@a)
    kernel = lambda t: np.vdot(r, expm(-1j*t*d)@r)
    for t in (.1, .7, 1.4):
        derivative = -1j*np.vdot(a, k@expm(-1j*t*k)@a)
        integral = quad(lambda s: (kernel(t-s)*c(s)).real, 0, t, epsabs=1e-12)[0]
        integral += 1j*quad(lambda s: (kernel(t-s)*c(s)).imag, 0, t, epsabs=1e-12)[0]
        assert abs(derivative+1j*mean*c(t)+integral) < 1e-11


def test_full_source_certificate_independent_diagonalization(result):
    _, _, raw, scale, _ = finite.verified_source()
    h = raw/scale
    e, v = np.linalg.eigh(h)
    z = np.array([1-2*(mask&1) for mask in range(256)])
    a = z*v[:, 0]
    gaps = e-e[0]
    weights = abs(v.conj().T@a)**2
    mean = float(weights@gaps)
    variance = float(weights@(gaps*gaps)-mean*mean)
    value = result['complete_source_memory']
    assert contains(value['mean_excitation_energy_interval'], mean)
    assert contains(value['complement_memory_strength_interval'], variance)
    for row in value['rows']:
        eta = float(F(row['eta']))
        g = float(weights@(1/(gaps+eta)))
        assert contains(row['normalized_resolvent_interval'], g)
        assert contains(row['complement_memory_interval'], mean+eta-1/g)
        assert row['all_spectral_weights_count'] == 256


def test_history_certificate_independent_full_unitary(result):
    _, _, raw, scale, _ = finite.verified_source()
    h = raw/scale
    _, v = np.linalg.eigh(h)
    g = v[:, 0]
    u = expm(1j*h/16)
    z0 = np.diag([1-2*(mask&1) for mask in range(256)])
    z1 = np.diag([1-2*((mask>>1)&1) for mask in range(256)])
    a = u@z0@u.conj().T
    c = np.vdot(g, a@z1@a@z1@g)
    row = result['complete_source_noncentral_history']
    assert contains(row['real_expectation_interval'], c.real)
    assert contains(row['imaginary_expectation_interval'], c.imag)
    assert contains(row['noncentral_history_variance_interval'], 1-abs(c)**2)
    assert F(row['noncentral_history_variance_interval'][0]) > F(6, 1000)
    for part in ('real', 'imaginary'):
        probabilities = row['complete_'+part+'_interferometer_probabilities']
        component = c.real if part == 'real' else c.imag
        assert contains(probabilities['plus'], (1+component)/2)
        assert contains(probabilities['minus'], (1-component)/2)


def test_polynomial_time_sign_and_unitary_remainder():
    h = np.array([[2, 1], [1, -1]], dtype=object)
    a = np.array([2, -3], dtype=object)
    b = np.array([1, 4], dtype=object)
    for sign in (-1, 1):
        r, i, divisor = finite._phase_polynomial_apply(h, a, b, 8, sign, 12)
        actual = np.array([complex(float(F(int(x), divisor)), float(F(int(y), divisor)))
                           for x, y in zip(r, i, strict=True)])
        expected = expm(sign*1j*np.array(h, float)/8)@(np.array(a, float)+1j*np.array(b, float))
        assert np.linalg.norm(actual-expected) < 1e-13


def test_exact_controls_and_readout_error():
    result = algebra.exact_controls()
    assert result['single_pole']['memory'] == 0
    assert result['two_pole']['memory'] == F(4, 9)
    error = algebra.readout_memory_error(F(1, 4), 28, F(1, 10**8), F(1, 10**8))
    assert 0 < error < F(1, 10000)
    with pytest.raises(ValueError):
        algebra.readout_memory_error(F(1, 4), 28, F(0), F(1, 100))


@pytest.mark.parametrize('energies,weights,eta', [([], [], 1), ([1], [2], 1),
    ([-1], [1], 1), ([1, 2], [2, -1], 1), ([1], [1], 0)])
def test_invalid_measure_rejected(energies, weights, eta):
    with pytest.raises(ValueError):
        algebra.spectral_memory(energies, weights, eta)


def test_same_generator_regularized_moments_without_first_moment_tightness():
    # A first moment escapes, but each fixed-filter observable in Eq11 converges.
    eta, lam = F(1, 4), F(1, 2)
    target_n = 1/(1+lam)**2
    errors = []
    for n in (100, 1000, 10000):
        energies, weights = [F(1), F(n)], [F(n-1, n), F(1, n)]
        norm = sum(w/(e+lam)**2 for e, w in zip(energies, weights))
        first = sum(w*e/(e+lam)**2 for e, w in zip(energies, weights))/norm
        second = sum(w*e*e/(e+lam)**2 for e, w in zip(energies, weights))/norm
        response = sum(w/((e+lam)**2*(e+eta)) for e, w in zip(energies, weights))/norm
        errors.append(abs(norm-target_n)+abs(first-1)+abs(second-1)+abs(response-1/(1+eta)))
    assert errors[2] < errors[1] < errors[0] < F(1, 10)


def test_actual_source_uniform_filter_probability_and_memory_bound():
    row = asymptotic.source_filter_bounds()
    assert row['Q_Weyl_amplitude']*row['comparison_P_Weyl_amplitude'] == 1
    assert row['unfiltered_norm_liminf_lower'] == F(25, 2304*1000**2)
    assert row['unfiltered_scaled_energy_limsup_upper'] == 12*1000**2
    assert 0 < row['eventual_filtered_norm_lower'] < row['filtered_norm_liminf_lower']
    assert row['eventual_filtered_memory_upper']*row['eventual_filtered_norm_lower'] == 1
    # Independent Jensen control on broad and sharply concentrated measures.
    for energies, weights in [([F(1), F(3)], [F(1, 8), F(3, 8)]),
                             ([F(1, 1000), F(1000)], [F(1, 2), F(1, 4)]),
                             ([F(2)], [F(1)])]:
        n = sum(weights)
        h = sum(e*w for e, w in zip(energies, weights))
        filtered = sum(w/(e+1)**2 for e, w in zip(energies, weights))
        assert filtered >= n**3/(h+n)**2


def test_two_speed_field_control_has_exact_source_moments_without_one_cone():
    row = asymptotic.two_speed_control()
    assert row['two_speeds_squared'] == ['6', '18']
    assert row['Q_first_energy_moment'] == '3/2'
    assert row['P_first_energy_coefficient_over_momentum_squared'] == '1/8'
    assert row['P_energy_over_momentum_ratio_defect'] == '7 - 4*sqrt(3)'
    assert row['one_common_Lorentz_speed'] is False
    assert row['is_a_state_or_scaling_limit_of_the_OPH_spin_source'] is False


def test_zero_field_same_source_and_full_response(result):
    h0, scale = zero_field.matrix()
    # Independent matrix from the diagonal and exchange terms, without field.
    expected = np.zeros((256, 256), dtype=np.int64)
    for mask in range(256):
        for x, y in finite.source.cube()[0]:
            if ((mask>>x)^(mask>>y))&1:
                expected[mask, mask] -= scale//2
                expected[mask^(1<<x)^(1<<y), mask] -= scale//2
    assert np.array_equal(h0, expected)
    e, v = np.linalg.eigh(h0/scale)
    a = v[:, 0]*np.array([1-2*(mask&1) for mask in range(256)])
    weights = abs(v.conj().T@a)**2
    gaps = e-e[0]
    mean = weights@gaps
    variance = weights@(gaps*gaps)-mean**2
    response = weights@(1/(gaps+.25))
    row = result['complete_same_box_zero_field_memory']
    assert contains(row['certificate']['ground_energy_interval'], e[0])
    assert contains(row['certificate']['gap_interval'], gaps[1])
    assert contains(row['memory_strength_interval'], variance)
    assert contains(row['memory_interval'], mean+.25-1/response)
    assert row['memory_explained_only_by_selection_field_torque'] is False
    assert F(row['memory_interval'][0]) > F(3, 10)


def test_zero_field_unique_ground_witness_rejects_deleted_level():
    with np.load(packet.ROOT/zero_field.WITNESS, allow_pickle=False) as archive:
        witness = {key: archive[key] for key in archive.files}
    witness['vectors'] = witness['vectors'][:, :-1]
    with pytest.raises(ValueError):
        zero_field.verify_witness(witness)


def test_actual_same_box_selection_field_memory_continuity(result):
    h0, scale = zero_field.matrix()
    energies, vectors = np.linalg.eigh(h0/scale)
    z = np.array([1-2*(mask&1) for mask in range(256)])
    a0 = z*vectors[:, 0]
    eta = F(1, 4)
    gap = F(result['complete_same_box_zero_field_memory']['certificate']['gap_interval'][0])
    gaps = energies-energies[0]
    weights = abs(vectors.conj().T@a0)**2
    mean0 = weights@gaps
    g0 = weights@(1/(gaps+float(eta)))
    sigma0 = mean0+float(eta)-1/g0
    for field in (F(1, 128), F(1, 1024)):
        hb, hbscale = finite.source.matrix(field=field)
        e, v = np.linalg.eigh(hb/hbscale)
        a = z*v[:, 0]
        gaps = e-e[0]
        weights = abs(v.conj().T@a)**2
        mean = weights@gaps
        response = weights@(1/(gaps+float(eta)))
        sigma = mean+float(eta)-1/response
        bound = zero_field.finite_field_continuity(gap, field, eta)
        assert abs(mean-mean0) <= float(F(bound['mean_difference_bound']))
        assert abs(response-g0) <= float(F(bound['resolvent_difference_bound']))
        assert abs(sigma-sigma0) <= float(F(bound['memory_difference_bound']))


def test_charge_resolved_memory_with_independent_spectral_projection(result):
    operator = charge_memory.charge.operators()
    h = operator['H']/operator['scale']
    e, v = np.linalg.eigh(h)
    ground = v[:, 0]
    created = operator['local_Z']@ground
    rows = result['complete_source_charge_resolved_memory']['channels']
    full_response = np.abs(v.conj().T@created)**2@(1/(e-e[0]+.25))
    reconstructed_response = 0.
    for row in rows:
        q = row['source_spin_charge']
        columns = operator['walsh'][:, np.array(operator['charges']) == q]/16
        image = columns@(columns.conj().T@created)
        probability = np.vdot(image, image).real
        image /= np.sqrt(probability)
        weights = abs(v.conj().T@image)**2
        gaps = e-e[0]
        mean = weights@gaps
        variance = weights@(gaps*gaps)-mean**2
        response = weights@(1/(gaps+.25))
        assert contains(row['original_channel_probability_interval'], probability)
        assert contains(row['conditional_mean_excitation_energy_interval'], mean)
        assert contains(row['conditional_memory_strength_interval'], variance)
        assert contains(row['conditional_memory_interval'], mean+.25-1/response)
        assert contains(row['original_complete_charge_and_energy_score_interval'], probability*.25*response)
        reconstructed_response += probability*response
    assert abs(reconstructed_response-full_response) < 1e-13
    assert result['complete_source_charge_resolved_memory']['physical_particle_count_assigned'] == 0


def test_fine_charge_label_grouping_does_not_prepare_coherent_branch(result):
    operator = charge_memory.charge.operators()
    h = operator['H']/operator['scale']
    e, v = np.linalg.eigh(h)
    ground = v[:, 0]
    created = operator['local_Z']@ground
    q = -1
    columns = operator['walsh'][:, np.array(operator['charges']) == q]/16
    amplitudes = columns.conj().T@created
    probability = np.vdot(amplitudes, amplitudes).real
    coherent = columns@amplitudes/np.sqrt(probability)
    fine_weights = abs(amplitudes)**2/probability
    coherent_energy = np.vdot(coherent, h@coherent).real-e[0]
    dephased_energy = sum(weight*np.vdot(column, h@column).real
                         for weight, column in zip(fine_weights, columns.T, strict=True))-e[0]
    assert abs(dephased_energy-coherent_energy) > .1
    assert np.sum(fine_weights**2) < .5
    row = result['complete_source_charge_resolved_memory']
    assert row['complete_creation_outcome_count'] == 2
    assert row['independent_fine_charge_audit_record_count'] == 256
    assert row['fine_X_audit_claimed_to_prepare_coherent_charge_branches'] is False


@pytest.mark.parametrize('order,regularizer,resolution', [
    (True, 1, 1), (F(0), 1, 1), (F(3, 4), 1, 1), (F(1, 1000), 0, 1),
    (F(1, 1000), 1, F(-1)), (.01, 1, 1)])
def test_uniform_estimate_rejects_unsupported_parameters(order, regularizer, resolution):
    with pytest.raises((ValueError, TypeError)):
        asymptotic.source_filter_bounds(order, regularizer, resolution)


@pytest.mark.parametrize('field', ['continuum_memory_numerically_evaluated',
    'noncentral_history_survives_field_scaling',
    'sharp_mass_shell_Lorentz_covariance_or_scattering_proved',
    'finite_gap_identified_as_physical_particle_mass'])
def test_scientific_promotion_rejected(result, field, monkeypatch):
    # Freeze the independently built expected value, not its comparison guard.
    expected = copy.deepcopy(result)
    monkeypatch.setattr(packet, 'build', lambda: expected)
    value = copy.deepcopy(result)
    value['scope'][field] = True
    with pytest.raises(ValueError):
        packet.validate(value)


def test_complete_probability_or_memory_mutation_rejected(result, monkeypatch):
    monkeypatch.setattr(packet, 'build', lambda: result)
    for mutate in ('probability', 'memory', 'level_count'):
        value = copy.deepcopy(result)
        if mutate == 'probability':
            del value['complete_source_noncentral_history']['complete_real_interferometer_probabilities']['minus']
        elif mutate == 'memory':
            value['complete_source_memory']['rows'][0]['complement_memory_interval'] = ['0', '0']
        else:
            value['complete_source_memory']['rows'][0]['all_spectral_weights_count'] = 255
        with pytest.raises(ValueError):
            packet.validate(value)


def test_fresh_owner_runtime_imports_all_bound():
    script = '''
import json, pathlib, sys
from model.spin_dynamical_memory import packet
value=packet.build()
root=packet.ROOT
missing=[]
for module in list(sys.modules.values()):
    raw=getattr(module,'__file__',None)
    if not raw: continue
    path=pathlib.Path(raw).resolve()
    if path.suffix=='.py' and path.is_relative_to(root):
        rel=path.relative_to(root).as_posix()
        if rel not in value['source_bindings']: missing.append(rel)
print(json.dumps(sorted(set(missing))))
'''
    process = subprocess.run([sys.executable, '-W', 'error', '-c', script],
                             cwd=packet.ROOT, check=True, capture_output=True, text=True)
    assert json.loads(process.stdout) == []


def test_response_closure_retains_original_scores_and_physical_boundaries(result, monkeypatch):
    expected = copy.deepcopy(result)
    monkeypatch.setattr(packet, 'build', lambda: expected)
    native = result['bounded_response_closure']['native_dimer']
    assert F(native['original_success_probability']) == F(16201, 48400)
    assert F(native['normalized_rank_one_defect']) == F(5184, 262472401)
    assert result['scope']['bounded_resolvent_reducing_criterion_proved'] is True
    for field in ('two_pair_U1_twists_alone_identify_a_wave_equation',
                  'actual_continuum_response_reducing_defect_proved_zero',
                  'arbitrary_Gram_response_projection_implemented'):
        changed = copy.deepcopy(result)
        changed['scope'][field] = True
        with pytest.raises(ValueError):
            packet.validate(changed)
    changed = copy.deepcopy(result)
    changed['bounded_response_closure']['native_dimer']['complete_record_probabilities'].pop()
    with pytest.raises(ValueError):
        packet.validate(changed)


def test_spectral_hierarchy_does_not_promote_Ritz_values_or_discard_original_failures(result,monkeypatch):
    expected=copy.deepcopy(result)
    monkeypatch.setattr(packet,'build',lambda:expected)
    row=result['source_spectral_window_hierarchy']['actual_source_joint_response']
    assert row['joint_block_rank']==3
    assert row['one_step_Ritz_is_an_exact_source_level'] is False
    for flag in ('uniform_growing_Gram_floor_proved','source_continuum_atom_inferred_from_finite_Ritz_values'):
        altered=copy.deepcopy(result)
        altered['scope'][flag]=True
        with pytest.raises(ValueError):packet.validate(altered)
    altered=copy.deepcopy(result)
    altered['source_spectral_window_hierarchy']['actual_source_joint_response']['filter_failure_probability']='0'
    with pytest.raises(ValueError):packet.validate(altered)
