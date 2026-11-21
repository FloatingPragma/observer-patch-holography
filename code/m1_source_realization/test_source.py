"""Mathematical, physical-domain, resource and hostile-input controls."""

import copy
import fnmatch
import itertools
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from m1_source_realization import check, circuits, model, topology, verify
from source_selection_model import verify_response


@pytest.fixture(scope='module')
def candidate():
    return model.candidate() | circuits.candidate() | {'topology': topology.candidate()}


def test_independent_full_replay(candidate):
    check.verify_evidence(candidate)


def test_committed_receipt(candidate):
    packet = verify.verify()
    check.same(packet['evidence'], candidate)


@pytest.mark.parametrize('roundoff', [-1.8656162355304537e-17, 1.8656162355304537e-17])
def test_coin_matrix_roundoff_does_not_relax_physical_scales(roundoff):
    path = 'evidence.coin.stages[2].matrix[0][0][1]'
    check.same(0., roundoff, path)
    check.same(roundoff, 0., path)
    with pytest.raises(ValueError):
        check.same(0., 1e-8, path)
    for field in ('phase_rate', 'code_event_time'):
        with pytest.raises(ValueError):
            check.same(0., 1e-14, 'evidence.resources.'+field)


def test_parent_source_ordered_response_and_charts():
    verify_response.verify(verify_response.strict_load(verify.ROOT/'code/source_selection_model/response.json'))


@pytest.mark.parametrize('d', [2, 3, 4, 5, 7])
def test_two_basis_projector_identity_and_general_channel(d):
    omega = np.eye(d).reshape(-1)/np.sqrt(d)
    pz = sum(np.outer(np.kron(z, z), np.kron(z, z)) for z in np.eye(d))
    px = sum(np.outer(np.kron(x.conj(), x), np.kron(x.conj(), x).conj()) for x in model.fourier(d).T)
    bell = np.outer(omega, omega)
    assert np.linalg.norm(pz@px-bell) < 1e-12
    assert np.linalg.eigvalsh(np.eye(d*d)-pz-px+bell).min() > -1e-12
    rng = np.random.default_rng(130+d)
    # Generic Stinespring isometry, not just the selected Weyl family.
    raw = rng.normal(size=(d*3, d))+1j*rng.normal(size=(d*3, d))
    v, _ = np.linalg.qr(raw)
    ks = list(v.reshape(3, d, d))
    j = model.choi(ks)
    fe = np.vdot(omega, j@omega).real
    assert fe >= model.classical_fidelity(ks, np.eye(d))+model.classical_fidelity(ks, model.fourier(d))-1-1e-12


@pytest.mark.parametrize('d', [2, 3, 4])
def test_a3_optimizer_beats_full_random_feasible_faces(d):
    rng = np.random.default_rng(720+d)
    ez, ex = .12, .17
    optimum = model.choi(model.selected_channel(d, ez, ex))
    def entropy(j):
        values = np.linalg.eigvalsh(j)
        return -sum(x*np.log(x) for x in values if x > 1e-14)
    h = entropy(optimum)
    omega = np.eye(d).reshape(-1)/np.sqrt(d)
    identity = np.outer(omega, omega)
    for _ in range(12):
        # Mix a generic channel with identity until both read constraints hold.
        raw = rng.normal(size=(d*d, d))+1j*rng.normal(size=(d*d, d))
        v, _ = np.linalg.qr(raw)
        ks = list(v.reshape(d, d, d))
        fz = model.classical_fidelity(ks, np.eye(d))
        fx = model.classical_fidelity(ks, model.fourier(d))
        weight = min(1, ez/(1-fz), ex/(1-fx))
        trial = weight*model.choi(ks)+(1-weight)*identity
        assert entropy(trial) <= h+1e-12
    dephase = model.choi([np.diag(np.eye(d)[j]) for j in range(d)])
    assert entropy(dephase) > entropy(identity)+.6
    # A changed reference really changes the optimization problem.
    biased = .8*identity+.2*np.eye(d*d)/(d*d)
    assert entropy(biased) < entropy(np.eye(d*d)/(d*d))


@pytest.mark.parametrize('d', [2, 4])
def test_code_transfer_preserves_entanglement_but_loses_full_state(d):
    ks = model.code_transfer_kraus(d)
    assert np.linalg.norm(sum(k.conj().T@k for k in ks)-np.eye(6)) < 1e-14
    bell = np.zeros((d, 6), complex)
    for j in range(d):
        bell[j, j] = 1/np.sqrt(d)
    rho = np.outer(bell.reshape(-1), bell.reshape(-1))
    out = model.act([np.kron(np.eye(d), k) for k in ks], rho)
    expected = np.zeros((d, d+1), complex)
    for j in range(d):
        expected[j, j] = 1/np.sqrt(d)
    assert np.linalg.norm(out-np.outer(expected.reshape(-1), expected.reshape(-1))) < 1e-14
    # Coherent +/- inputs outside the code have identical complete outputs.
    plus, minus = np.zeros(6), np.zeros(6)
    plus[0], plus[d] = 1/np.sqrt(2), 1/np.sqrt(2)
    minus[0], minus[d] = 1/np.sqrt(2), -1/np.sqrt(2)
    assert np.linalg.norm(model.act(ks, np.outer(plus, plus))-model.act(ks, np.outer(minus, minus))) < 1e-14
    failure = model.act(ks, np.diag(np.eye(6)[d]))
    assert failure[d, d] == 1 and np.trace(failure) == 1
    # Reset sender and overwrite receiver: the old receiver is fully traced.
    old = np.diag([.1, .2, .05, .15, .3, .2])
    full = np.kron(rho, old).reshape(d, 6, 6, d, 6, 6)
    sender = np.trace(full, axis1=2, axis2=5).reshape(6*d, 6*d)
    assert np.linalg.norm(sender-rho) < 1e-14
    decohered = model.act([np.kron(np.eye(d), k) for k in model.code_transfer_kraus(d, False)], rho)
    assert np.vdot(expected.reshape(-1), decohered@expected.reshape(-1)).real == pytest.approx(1/d)


def test_complete_response_cannot_include_buffer_or_sector_control():
    gs = [np.asarray(g, complex) for g in model.source_generators()]
    # Independent full-operator directions, not just different names for D.
    base = [np.kron(g, np.eye(2)) for g in gs]
    forbidden = np.kron(np.eye(6), np.diag([1j, -1j]))
    def rank(matrices):
        return np.linalg.matrix_rank(np.array([np.r_[m.real.ravel(), m.imag.ravel()] for m in matrices]))
    assert rank(base) == 12 and rank(base+[forbidden]) == 13
    base = [np.kron(np.eye(12), g) for g in gs]
    conditioned = np.kron(np.diag([1]+[0]*11), gs[0])
    assert rank(base+[conditioned]) == 13
    assert 4**2-1 > 12  # Logical SU4 is not the primitive carrier tangent.


def test_injective_readback_is_not_a_reversible_quantum_response():
    # Random Pauli readback with the complete basis/outcome transcript has
    # informationally complete classical probabilities. A common kernel
    # argument would be wrong; every Kraus branch still has deficient rank.
    vectors = [np.array([1, 0]), np.array([0, 1]),
               np.array([1, 1])/np.sqrt(2), np.array([1, -1])/np.sqrt(2),
               np.array([1, 1j])/np.sqrt(2), np.array([1, -1j])/np.sqrt(2)]
    ks = [np.outer(np.eye(6)[i], v.conj())/np.sqrt(3) for i, v in enumerate(vectors)]
    assert np.linalg.norm(sum(k.conj().T@k for k in ks)-np.eye(2)) < 1e-14
    assert all(np.linalg.matrix_rank(k) == 1 for k in ks)
    basis = [np.eye(2), np.array([[0, 1], [1, 0]]),
             np.array([[0, -1j], [1j, 0]]), np.diag([1, -1])]
    images = [model.act(ks, h).real.reshape(-1) for h in basis]
    assert np.linalg.matrix_rank(images) == 4
    # Even keeping every classical outcome and preparing the corresponding
    # eigenvector gives entanglement fidelity 1/2, not an inverse channel.
    measure_prepare = [np.outer(v, v.conj())/np.sqrt(3) for v in vectors]
    omega = np.eye(2).reshape(-1)/np.sqrt(2)
    j = model.choi(measure_prepare)
    assert np.vdot(omega, j@omega).real == pytest.approx(.5)


def test_nontrivial_final_qr_phase_cannot_be_discarded():
    # The reference coin happens to have determinant one and trivial final
    # phases. A second full unitary checks that the general compiler keeps
    # the final phase rather than depending on that accidental simplification.
    target = circuits.COIN@np.diag([1, 1, 1, np.exp(.37j)])
    diagonal, moves = model.decompose(target)
    actual, wrong = np.diag(diagonal), np.eye(4, dtype=complex)
    for i, k, g in moves:
        actual[[i, k], :] = g@actual[[i, k], :]
        wrong[[i, k], :] = g@wrong[[i, k], :]
    assert np.linalg.norm(actual-target) < 1e-13
    assert np.linalg.norm(wrong-target) > .3


@pytest.mark.parametrize('q,a,mass', [(2, .07, .7), (3, .013, 1.1)])
def test_every_spatial_matrix_entry_and_inverse(q, a, mass):
    actual = circuits.execute(q, a, mass)
    exact = check.closed_spatial(q, a, mass)
    assert np.linalg.norm(actual-exact) < 1e-12
    assert np.linalg.norm(actual.conj().T@actual-np.eye(len(actual))) < 1e-12


def test_single_particle_equivalence_does_not_imply_fermionic_flights():
    swap = np.eye(4)[[0, 2, 1, 3]]
    fermionic = swap@np.diag([1, 1, 1, -1])
    assert np.array_equal(swap[:, :3], fermionic[:, :3])
    assert np.linalg.norm(swap-fermionic) == 2
    psi = np.array([1, 0, 0, 1])/np.sqrt(2)
    assert abs(np.vdot(swap@psi, fermionic@psi)) < 1e-14


@pytest.mark.parametrize('d', [2, 3, 4])
def test_noisy_optimizer_time_refinement(d):
    def channel(t):
        return model.selected_channel(d, 0., (d-1)/d*(-np.expm1(-.7*t)))
    rng = np.random.default_rng(123)
    psi = rng.normal(size=d)+1j*rng.normal(size=d)
    psi /= np.linalg.norm(psi)
    rho = np.outer(psi, psi.conj())
    assert np.linalg.norm(model.act(channel(.2), model.act(channel(.3), rho))-model.act(channel(.5), rho)) < 1e-13
    repeated = rho.copy()
    for _ in range(40):
        repeated = model.act(channel(.5/40), repeated)
    assert np.linalg.norm(repeated-model.act(channel(.5), rho)) < 1e-13
    # Fixed error per subdivision is a different, increasingly destructive law.
    wrong = rho.copy()
    for _ in range(40):
        wrong = model.act(channel(.5), wrong)
    assert np.linalg.norm(wrong-repeated) > .1


def test_population_free_phase_bound_does_not_hide_vacuum_errors():
    for row in circuits.noise():
        assert row['half_trace_distance'] <= row['population_free_upper']+1e-12
    assert 1-(1-.02)**1000 > .999
    assert circuits.noise()[0]['vacuum_bitflip_failure'] < .08


def test_phase_noise_uv_cost_and_control_timing():
    kappa = np.arcsin(1/np.sqrt(6))/2
    for a in (1e-2, 1e-4, 1e-6):
        tau = a/np.sqrt(3)
        fixed = (-np.expm1(-2*.01))*kappa/tau
        shrinking = (-np.expm1(-2*a*a))*np.pi/tau
        assert fixed > .006/a and shrinking < 11*a
    r = circuits.resource_witness()
    assert 48*np.pi/r['drive_norm_bound']+128*r['code_event_time'] == pytest.approx(r['wall_tick']-r['flight_tick'])
    assert r['beat_frequency'] != pytest.approx(.08)
    assert r['preparation_time'] > 79 and r['mode_buffers'] > 10**30
    assert r['observable_swing_lower'] > .6752


def test_uncharged_identity_instruction_is_detected(monkeypatch):
    original = circuits.instructions
    monkeypatch.setattr(circuits, 'instructions', lambda a, m:
                        original(a, m)+[('phase', [0], np.array([[1.+0j]]))])
    # Identical quantum evolution still has a different resource/time budget.
    assert circuits.resource_witness()['pulses_per_cell_tick'] == 50
    with pytest.raises(ValueError):
        check.same(circuits.resource_witness(), check.expected_resources(), 'resources')


@pytest.mark.parametrize('theta,weight', [(-.4, 0.), (0., .37), (1.2, 1.)])
def test_full_readout_instrument_and_detector_effect(theta, weight):
    unitary = circuits.detector_unitary(theta)
    assert np.linalg.norm(unitary.conj().T@unitary-np.eye(17)) < 1e-13
    for mode in range(8):
        click, rest = circuits.destructive_read(mode, weight)
        ks = (click, *rest)
        assert np.linalg.norm(sum(k.conj().T@k for k in ks)-np.eye(17)) < 1e-14
        vacuum = np.zeros(17)
        vacuum[0] = 1
        assert np.linalg.norm(click@vacuum) == 0
    actual = next(r for r in circuits.readouts() if r['theta'] == theta and r['weight'] == weight)
    expected = next(r for r in check.expected_readouts() if r['theta'] == theta and r['weight'] == weight)
    check.same(actual, expected, 'readout')
    effect = check.matrix(actual['effect'], (17, 17))
    values = np.linalg.eigvalsh(effect)
    assert min(values) > -1e-13 and max(values) <= 1+1e-13
    assert sum(actual['first_click_probabilities'])+actual['no_click_trace'] == pytest.approx(1)


@pytest.mark.parametrize('field,bad', [('blank_events_per_processor', 1),
    ('seed_events', 1), ('departures_per_tree_node', 1), ('readout_buffer_events', 1),
    ('preparation_events_on_path', 672)])
def test_missing_serial_operations_fail_even_when_total_float_time_is_close(candidate, field, bad):
    p = copy.deepcopy(candidate)
    p['resources']['serial_schedule'][field] = bad
    with pytest.raises(ValueError):
        check.verify_evidence(p)


def test_missed_readout_mode_and_reversed_phase_fail(candidate):
    for field in ('first_click_probabilities', 'effect'):
        p = copy.deepcopy(candidate)
        if field == 'first_click_probabilities':
            p['readouts'][1][field][-1] = 0.
        else:
            p['readouts'][1][field][1][9][1] *= -1
        with pytest.raises(ValueError):
            check.verify_evidence(p)


def test_classical_source_loses_the_same_ramsey_readout_phase():
    psi = np.zeros(17, complex)
    psi[1] = psi[9] = 1/np.sqrt(2)
    probabilities = []
    for theta in (0., np.pi/2, np.pi):
        classical = circuits.classical_detector_effect(theta, 1.)
        assert np.linalg.norm(classical-np.diag([0.]+[.5]*16)) < 1e-13
        assert np.vdot(psi, classical@psi).real == pytest.approx(.5)
        transformed = circuits.detector_unitary(theta)@psi
        probabilities.append(float(sum(abs(transformed[1:9])**2)))
    assert probabilities == pytest.approx([1., .5, 0.])


@pytest.mark.parametrize('mutation', ['discard_no_click', 'classical_coherence', 'classical_probability'])
def test_readout_outcome_and_classical_ablation_cannot_be_faked(candidate, mutation):
    p = copy.deepcopy(candidate)
    row = p['readouts'][1]
    if mutation == 'discard_no_click':
        row['no_click_trace'] = 0.
    elif mutation == 'classical_coherence':
        row['classical_effect'] = copy.deepcopy(row['effect'])
    else:
        row['classical_probability'] = 1.
    with pytest.raises(ValueError):
        check.verify_evidence(p)


def test_basis_resets_are_explicit_cptp_interventions_not_ontic_optimizers():
    for level in (0, 1):
        ks = [np.outer(np.eye(6)[level], np.eye(6)[j]) for j in range(6)]
        assert np.array_equal(sum(k.T@k for k in ks), np.eye(6))
        for i, j in itertools.product(range(6), repeat=2):
            unit = np.zeros((6, 6))
            unit[i, j] = 1
            assert np.array_equal(model.act(ks, unit), (i == j)*np.diag(np.eye(6)[level]))
        assert all(np.linalg.matrix_rank(k) == 1 for k in ks)


@pytest.mark.parametrize('bad', [True, 1, 3, 5, 6, 7, 2., '2', None])
def test_unsupported_source_codes_cannot_bypass_properness(bad):
    with pytest.raises(ValueError):
        model.code_transfer_kraus(bad)


@pytest.mark.parametrize('d', [2, 3, 4])
def test_loose_a3_bounds_select_uniformity_instead_of_saturating(d):
    ks = model.selected_channel(d, 1., .99)
    assert np.linalg.norm(model.choi(ks)-np.eye(d*d)/(d*d)) < 1e-13


@pytest.mark.parametrize('bad', [-.01, 1.01, float('nan'), float('inf'), True, '0'])
def test_invalid_agreement_bounds_fail(bad):
    with pytest.raises(ValueError):
        model.selected_channel(2, bad, 0.)
    with pytest.raises(ValueError):
        model.selected_channel(2, 0., bad)


@pytest.mark.parametrize('groups', [1, 2, 4, 8])
def test_phase_generator_on_vacuum_plus_one_particle_has_no_population_factor(groups):
    modes = 2*groups
    size = modes+1
    rng = np.random.default_rng(610+groups)
    raw = rng.normal(size=(size, size))+1j*rng.normal(size=(size, size))
    rho = raw@raw.conj().T
    rho /= np.trace(rho)
    actual = np.zeros_like(rho)
    for g in range(groups):
        labels = np.zeros(size, int)
        labels[1+2*g:3+2*g] = (1, 2)
        actual += (labels[:, None] == labels[None, :])*rho-rho
    one = np.diag([0.]+[1.]*modes)
    diagonal = np.diag(np.diag(one@rho@one))
    grouped = np.zeros_like(rho)
    for g in range(groups):
        p = np.zeros_like(rho)
        p[1+2*g, 1+2*g] = p[2+2*g, 2+2*g] = 1
        grouped += p@rho@p
    assert np.linalg.norm(actual-(diagonal+grouped-one@rho-rho@one)) < 1e-13
    assert np.sum(abs(np.linalg.eigvalsh(actual)))/2 <= 2+1e-13


@pytest.mark.parametrize('mutation', [
    'control', 'phase', 'missing_coin_phase', 'dephasing_as_identity', 'drop_failure',
    'missing_modes', 'bad_flight', 'prep_hole', 'fixed_noise', 'retime', 'free_gates',
    'zero_noise', 'zero_time', 'empty_topology', 'duplicated_process_link', 'wrong_bridge',
    'bool_population', 'extra', 'nan'])
def test_hostile_evidence(candidate, mutation):
    p = copy.deepcopy(candidate)
    if mutation == 'control': p['controls']['rows'][0]['coefficients'][0] = ['0', '0']
    elif mutation == 'phase': p['entangler'][-1]['concurrence'] = 0.
    elif mutation == 'missing_coin_phase': p['coin']['diagonal'][0] = [-1., 0.]
    elif mutation == 'dephasing_as_identity': p['channels'][1]['x_fidelity'] = 1.
    elif mutation == 'drop_failure': p['instruments'][1]['matrix_unit_images'][35] = [[0., 0.]]*9
    elif mutation == 'missing_modes': p['spatial'][-1]['column_norms'].pop()
    elif mutation == 'bad_flight': p['spatial'][-1]['three_step_probe'][1][0] += .1
    elif mutation == 'prep_hole': p['preparations'][1]['amplitudes'][0][0][0] = [1., 0.]
    elif mutation == 'fixed_noise': p['noise'][0]['population_free_upper'] = .9
    elif mutation == 'retime': p['resources']['beat_frequency'] = .08
    elif mutation == 'free_gates': p['resources']['code_events_per_cell_tick'] = 96
    elif mutation == 'zero_noise': p['resources']['probability_noise_upper'] = 0.
    elif mutation == 'zero_time': p['resources']['code_event_time'] = 0.
    elif mutation == 'empty_topology': p['topology'] = []
    elif mutation == 'duplicated_process_link': p['topology'][1]['process_links'].append([0, 1])
    elif mutation == 'wrong_bridge': p['topology'][1]['coarsen'][0] = 2
    elif mutation == 'bool_population': p['spatial'][0]['q'] = True
    elif mutation == 'extra': p['declared_success'] = True
    elif mutation == 'nan': p['resources']['a'] = float('nan')
    with pytest.raises(ValueError):
        check.verify_evidence(p)


@pytest.mark.parametrize('raw', ['{}', '{"a":1,"a":2}', '{"a":NaN}', '{"a":Infinity}', '[]'])
def test_malformed_cli_input(tmp_path, raw):
    path = tmp_path/'bad.json'
    path.write_text(raw)
    result = subprocess.run([sys.executable, '-O', str(verify.HERE/'verify.py'), str(path)], capture_output=True)
    assert result.returncode != 0


def test_optimized_cli_accepts_reference_and_rejects_custody_and_physics(tmp_path):
    def run(path):
        return subprocess.run([sys.executable, '-O', str(verify.HERE/'verify.py'), str(path)], capture_output=True)
    assert run(verify.HERE/'receipt.json').returncode == 0
    original = verify.load(verify.HERE/'receipt.json')
    for kind in ('sources', 'claims', 'physics'):
        p = copy.deepcopy(original)
        if kind == 'physics':
            p['evidence']['resources']['phase_rate'] = 0.
        elif p[kind]:
            p[kind][next(iter(p[kind]))] = '0'*64
        else:
            p[kind]['OPH-UNREGISTERED'] = '0'*64
        path = tmp_path/f'{kind}.json'
        path.write_text(json.dumps(p))
        assert run(path).returncode != 0

