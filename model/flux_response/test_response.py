"""Independent Jordan-Wigner, gauge, finite-difference and degeneracy checks."""
import copy
import json
import math
from pathlib import Path

import numpy as np
import pytest
from scipy.linalg import expm

from .response import FluxRing, _json_value, _summary, report, verify_receipt, spectral_subspaces


def jw():
    operators = []
    for mode in range(8):
        value = np.ones((1, 1), complex)
        for j in reversed(range(8)):
            value = np.kron(value, np.array([[0, 1], [0, 0]]) if j == mode else
                            np.diag([1, -1]) if j < mode else np.eye(2))
        operators.append(value)
    return operators


@pytest.fixture(scope='module')
def independent_terms():
    a = jw()
    n = [op.conj().T@op for op in a]
    contact = sum(n[2*x]@n[2*x+1] for x in range(4))
    hops = [sum(1j*sign/16*a[2*x+s].conj().T@a[2*((x+1) % 4)+s]
                for s, sign in enumerate((1, -1))) for x in range(4)]
    return n, contact, hops


def independent_action(terms, flux, weights=(.25,)*4):
    n, contact, hops = terms
    return contact-sum(n)/2+sum(np.exp(1j*w*flux)*hop+
        np.exp(-1j*w*flux)*hop.conj().T for w, hop in zip(weights, hops, strict=True))


@pytest.fixture(scope='module')
def ring():
    return FluxRing(.4)


def test_same_complete_action_and_independent_derivatives(ring, independent_terms):
    assert np.linalg.norm(ring.generator-independent_action(independent_terms, .4), 2) < 2e-15
    from ..particle_response.dispersion import PeriodicSpectrum
    assert np.array_equal(FluxRing().generator, PeriodicSpectrum().generator)
    step = .001
    plus = independent_action(independent_terms, .4+step)
    minus = independent_action(independent_terms, .4-step)
    far_plus = independent_action(independent_terms, .4+2*step)
    far_minus = independent_action(independent_terms, .4-2*step)
    first = (far_minus-8*minus+8*plus-far_plus)/(12*step)
    second = (-far_minus+16*minus-30*ring.generator+16*plus-far_plus)/(12*step**2)
    assert np.linalg.norm(first-ring.current, 2) < 2e-13
    assert np.linalg.norm(second-ring.diamagnetic, 2) < 3e-10
    assert ring.vectors.shape == (256, 256)
    assert np.linalg.norm(ring.generator@ring.number-ring.number@ring.generator, 2) < 1e-14
    assert np.trace(ring.number) == 1024
    # Omitting the winding bond changes the operator: it is not the open dimer service.
    _, _, hops = independent_terms
    missing = independent_action(independent_terms, .4)-(
        np.exp(.1j)*hops[-1]+np.exp(-.1j)*hops[-1].conj().T)
    assert np.linalg.norm(missing-ring.generator, 2) > .1


def test_local_number_continuity_is_an_operator_identity_for_every_state(ring):
    for x, number in enumerate(ring.occupations):
        derivative = 1j*(ring.generator@number-number@ring.generator)
        flux = ring.bond_currents[(x-1) % 4]-ring.bond_currents[x]
        assert np.linalg.norm(derivative-flux, 2) < 2e-15
    response = ring.ground_response()
    assert np.ptp(response['stationary_bond_currents']) < 2e-15
    assert abs(response['stationary_bond_currents'][0]-response['flux_derivative']) < 2e-15


@pytest.mark.parametrize('weights', [(0., 0., 0., 1.), (.5, -.25, .25, .5)])
def test_only_total_loop_flux_matters_but_diamagnetic_piece_is_gauge_dependent(ring, weights):
    other = FluxRing(.4, weights=weights)
    transform = ring.gauge_to(other)
    assert np.linalg.norm(other.generator-transform@ring.generator@transform.conj().T, 2) < 1e-14
    assert np.max(abs(ring.energies-other.energies)) < 6e-15
    a, b = ring.ground_response(), other.ground_response()
    assert abs(a['flux_derivative']-b['flux_derivative']) < 1e-14
    assert abs(a['flux_curvature']-b['flux_curvature']) < 1e-14
    assert abs(a['diamagnetic_expectation']-b['diamagnetic_expectation']) > .003
    # A phi-dependent gauge changes individual derivative operators; a naive
    # covariant-current assertion would incorrectly discard that connection.
    assert np.linalg.norm(other.current-transform@ring.current@transform.conj().T, 2) > .05


def test_loop_flux_cannot_be_gauged_away_but_is_two_pi_periodic(ring):
    assert np.max(abs(ring.energies-FluxRing().energies)) > .001
    periodic = FluxRing(.4+2*math.pi)
    assert np.max(abs(ring.energies-periodic.energies)) < 6e-15
    assert abs(ring.ground_response()['flux_curvature']-
               periodic.ground_response()['flux_curvature']) < 2e-14


@pytest.mark.parametrize('flux', [0., .4, 1.2])
def test_ground_curvature_matches_independent_five_point_energy_derivative(flux, independent_terms):
    ring = FluxRing(flux); response = ring.ground_response()
    step = .01
    energies = [np.linalg.eigvalsh(independent_action(independent_terms, flux+j*step))[0]
                for j in (-2, -1, 0, 1, 2)]
    first = (energies[0]-8*energies[1]+8*energies[3]-energies[4])/(12*step)
    second = (-energies[0]+16*energies[1]-30*energies[2]+
              16*energies[3]-energies[4])/(12*step**2)
    assert abs(response['flux_derivative']-first) < 2e-11
    assert abs(response['flux_curvature']-second) < 1e-9
    # Both terms are necessary: deleting either must fail the same observable.
    assert abs(response['diamagnetic_expectation']-second) > .005
    assert abs(-response['paramagnetic_sum']-second) > .005


def test_all_current_weights_moments_and_time_domain_remain_complete(ring):
    response = ring.ground_response()
    ground = ring.vectors[:, 0]; excited = ring.current@ground
    assert response['complete_current_weights'].shape == (256,)
    assert response['complete_gaps'].shape == (256,)
    assert np.count_nonzero(response['complete_current_weights'] < 1e-24) > 150
    assert abs(response['current_zeroth_moment']-
               response['complete_current_weights'].sum()) < 1e-14
    for time in (-.4, 0., .7):
        propagator = expm(-1j*time*(ring.generator-ring.energies[0]*np.eye(256)))
        direct = np.vdot(excited, propagator@excited)
        spectral = response['complete_current_weights']@np.exp(-1j*time*response['complete_gaps'])
        assert abs(direct-spectral) < 1e-14
    assert response['numerical_moment_residual'] < 1e-14
    assert np.all(response['complete_current_weights'] >= 0)
    with pytest.raises(ValueError):
        response['complete_current_weights'][0] = 1


@pytest.mark.parametrize('flux', [.4, 1.2])
def test_free_control_has_closed_bloch_energy_and_curvature(flux):
    ring = FluxRing(flux, interaction=0.)
    response = ring.ground_response(); alpha = flux/4
    assert ring.ground_degeneracy == 1
    assert abs(response['ground_energy']+2*ring.speed*(math.cos(alpha)+math.sin(alpha))) < 1e-14
    assert abs(response['flux_derivative']-ring.speed*(math.sin(alpha)-math.cos(alpha))/2) < 1e-14
    assert abs(response['flux_curvature']-ring.speed*(math.cos(alpha)+math.sin(alpha))/8) < 1e-14
    assert response['paramagnetic_sum'] < 1e-25


def test_degenerate_free_and_atomic_controls_do_not_receive_nondegenerate_formula():
    free = FluxRing(interaction=0.)
    atomic = FluxRing(.4, speed=0.)
    for ring in (free, atomic):
        assert ring.ground_degeneracy == 16
        assert len(ring.energies) == 256
        with pytest.raises(ValueError, match='degeneracy'):
            ring.ground_response()
    assert np.count_nonzero(atomic.current) == 0
    assert np.count_nonzero(atomic.diamagnetic) == 0
    assert np.array_equal(atomic.generator, FluxRing(1.2, speed=0.).generator)


@pytest.mark.parametrize('kwargs', [dict(flux=True), dict(flux=float('nan')),
    dict(speed=-1), dict(interaction=-1), dict(weights=(1., 1., 0., 0.)),
    dict(weights=(0., 1.)), dict(weights=(0., 0., False, 1.)),
    dict(weights=(0., 0., float('inf'), 1.)), dict(flux=101)])
def test_invalid_inputs_fail_before_evaluation(kwargs):
    with pytest.raises(ValueError):
        FluxRing(**kwargs)


def test_gauge_equivalence_cannot_change_action_or_total_flux(ring):
    for other in (FluxRing(.5), FluxRing(.4, interaction=.9)):
        with pytest.raises(ValueError, match='same loop flux'):
            ring.gauge_to(other)


def test_retained_receipt_and_physical_scope():
    result = _json_value(report())
    assert verify_receipt(json.loads(Path(__file__).with_name('results.json').read_text()), result)
    assert result['integrated_completion_gate_accepted'] is False
    assert result['dynamical_photon_or_physical_QED_identified'] is False
    assert result['physical_clock_or_current_calibrated'] is False
    assert result['source_pulse_execution_added'] is False
    assert result['rigorous_numerical_error_bound'] is False


def test_degenerate_unitary_basis_rotations_preserve_the_complete_replay():
    from ..reference_observation.observation import frozen
    baseline = _json_value(report())
    rotated = copy.deepcopy(baseline)
    rng = np.random.default_rng(881)
    largest_weight_change = 0.
    for index, flux in enumerate((0., .4, 1.2)):
        ring = FluxRing(flux)
        old = ring.ground_response()['complete_current_weights']
        vectors = ring.vectors.copy()
        groups = spectral_subspaces(baseline['interacting_responses'][index])
        start = 0
        for group in groups:
            size = group['multiplicity']
            if size > 1:
                q, _ = np.linalg.qr(rng.normal(size=(size, size))+1j*rng.normal(size=(size, size)))
                vectors[:, start:start+size] = vectors[:, start:start+size]@q
            start += size
        ring.vectors = frozen(vectors)
        response = ring.ground_response()
        largest_weight_change = max(largest_weight_change,
            float(np.max(abs(old-response['complete_current_weights']))))
        rotated['interacting_responses'][index] = _json_value(_summary(ring))
        assert len(response['complete_current_weights']) == 256
        assert np.linalg.norm(ring.generator@vectors-vectors*ring.energies, 2) < 1e-13
    assert largest_weight_change > 1e-5  # far larger than any numerical replay budget
    assert verify_receipt(baseline, rotated)


@pytest.mark.parametrize('damage', ['cross_energy_weight', 'missing_state', 'energy',
                                    'negative_weight', 'nan_weight', 'scope', 'multiplicity'])
def test_basis_invariant_replay_rejects_changed_measure_or_missing_modes(damage):
    expected = json.loads(Path(__file__).with_name('results.json').read_text())
    changed = copy.deepcopy(expected)
    row = changed['interacting_responses'][0]
    if damage == 'cross_energy_weight':
        donor = int(np.argmax(row['complete_current_weights']))
        transfer = row['complete_current_weights'][donor]/2
        row['complete_current_weights'][donor] -= transfer
        row['complete_current_weights'][0] += transfer  # preserve total weight but change its energy
    elif damage == 'missing_state':
        row['complete_gaps'].pop(); row['complete_current_weights'].pop()
    elif damage == 'energy':
        row['complete_gaps'][-1] += .001
    elif damage == 'negative_weight':
        row['complete_current_weights'][1] = -1e-8
    elif damage == 'nan_weight':
        row['complete_current_weights'][1] = float('nan')
    elif damage == 'multiplicity':
        row['complete_gaps'][1] = row['complete_gaps'][0]
    else:
        changed['integrated_completion_gate_accepted'] = True
    with pytest.raises(ValueError):
        verify_receipt(expected, changed)


def test_cli_uses_the_same_basis_invariant_verifier(monkeypatch, tmp_path):
    import sys
    from . import response
    receipt = tmp_path/'receipt.json'
    packet = json.loads(Path(__file__).with_name('results.json').read_text())
    receipt.write_text(json.dumps(packet))
    calls = []
    def verifier(expected, actual):
        calls.append((expected, actual))
        raise ValueError('deliberate verifier rejection')
    monkeypatch.setattr(response, 'verify_receipt', verifier)
    monkeypatch.setattr(response, 'report', lambda: packet)
    monkeypatch.setattr(sys, 'argv', ['flux', '--check', str(receipt)])
    with pytest.raises(SystemExit, match='deliberate verifier rejection'):
        response.main()
    assert len(calls) == 1


@pytest.mark.parametrize('path,replacement', [
    (('integrated_completion_gate_accepted',), 0),
    (('dynamical_photon_or_physical_QED_identified',), 0.),
    (('complete_Fock_dimension',), 256.),
    (('atomic_control', 'ground_degeneracy'), 16.),
    (('interacting_responses', 0, 'retained_spectral_entries'), 256.),
    (('interacting_responses', 0, 'complete_gaps', 0), 0),
    (('interacting_responses', 0, 'complete_current_weights', 0), 0),
])
def test_strict_replay_rejects_equal_valued_type_substitutions(path, replacement):
    expected = json.loads(Path(__file__).with_name('results.json').read_text())
    changed = copy.deepcopy(expected)
    target = changed
    for key in path[:-1]:
        target = target[key]
    assert type(target[path[-1]]) is not type(replacement)
    assert target[path[-1]] == replacement
    target[path[-1]] = replacement
    with pytest.raises(ValueError, match='type'):
        verify_receipt(expected, changed)


@pytest.mark.parametrize('value', [float('nan'), float('inf'), -float('inf')])
def test_strict_replay_rejects_nonfinite_nonspectral_values_even_in_both_inputs(value):
    packet = json.loads(Path(__file__).with_name('results.json').read_text())
    packet['gauge_redistribution']['hamiltonian_residual'] = value
    with pytest.raises(ValueError, match='nonfinite'):
        verify_receipt(packet, packet)


@pytest.mark.parametrize('value', [float('nan'), float('inf'), -float('inf')])
def test_strict_serialization_rejects_nonfinite_diagnostics(value):
    with pytest.raises(ValueError, match='nonfinite'):
        _json_value({'diagnostic': value})


@pytest.mark.parametrize('location', ['root', 'nested'])
def test_strict_cli_rejects_duplicate_keys(monkeypatch, tmp_path, location):
    import sys
    from . import response
    packet = json.loads(Path(__file__).with_name('results.json').read_text())
    document = json.dumps(packet)
    if location == 'root':
        document = '{"integrated_completion_gate_accepted":true,'+document[1:]
    else:
        needle = '"ground_degeneracy": 16'
        assert needle in document
        document = document.replace(needle, '"ground_degeneracy": 1, '+needle, 1)
    # The permissive decoder would hide the first, contradictory value.
    assert json.loads(document) == packet
    receipt = tmp_path/'duplicate.json'
    receipt.write_text(document, encoding='utf-8')
    monkeypatch.setattr(response, 'report', lambda: packet)
    monkeypatch.setattr(sys, 'argv', ['flux', '--check', str(receipt)])
    with pytest.raises(SystemExit, match='duplicate'):
        response.main()
