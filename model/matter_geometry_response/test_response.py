"""Independent source, spectral, interval and scope controls."""
from copy import deepcopy
from fractions import Fraction as F
from functools import lru_cache
import math

import numpy as np
import pytest

from . import geometry, packet, source


def test_independent_charge_instrument_and_coherent_input():
    z = np.diag([1., -1.])
    k0, k1 = (4*np.eye(2)+3*z)/math.sqrt(50), 1j*(4*np.eye(2)-3*z)/math.sqrt(50)
    state = np.array([1., 2j])/math.sqrt(5)
    rho = np.outer(state, state.conj())
    assert np.allclose(k0.conj().T@k0+k1.conj().T@k1, np.eye(2))
    assert abs(np.trace(k1@rho@k1.conj().T)-float(source.record_probabilities(F(4, 5))[1])) < 1e-14
    output = k0@rho@k0.conj().T+k1@rho@k1.conj().T
    assert abs(output[1, 1]-rho[1, 1]) < 1e-14
    assert abs(output[0, 1]) > .01


def test_actual_emitted_quantum_words():
    result = source.quantum_execution()
    assert result['all_input_columns_executed'] is True
    assert result['fresh_fair_coin_words_per_attempt'] == 9


def test_controller_and_clock_fractions():
    for energy in (F(0), F(1, 2), F(1), F(2, 7)):
        weights = source.scheduler_weights(energy)
        z = F(24, 25)*(1-2*energy)
        assert weights[None] == F(1, 16)
        assert weights[(0, 1)] == (8+z)/256
        assert weights[(0, 5)] == (8-z)/256
        assert sum(weights.values()) == 1
        assert all(p > 0 for p in weights.values())


def test_all_native_histories_and_mean():
    result = source.native_execution()
    assert result['complete_native_histories_checked'] == 79872
    assert result['all_initial_load_mean_identities_checked_exactly']
    difference = result['correlated_probe_countercontrol']['actual_minus_product_mean']
    assert [F(x) for x in difference] == [F(-3, 1600), F(3, 800), F(0), F(0), F(0),
                                         F(-3, 1600), F(0), F(0), F(0), F(0), F(0), F(0)]


def test_exact_native_derivative_against_independent_resolvent():
    c = geometry.source()
    dp, dg = geometry.exact_derivative()
    l = np.array(geometry.laplacian(), float)
    d = np.zeros((12, 12))
    for edge, sign in (((0, 1), 1), ((0, 5), -1)):
        v = np.zeros(12)
        v[edge[0]], v[edge[1]] = 1., -1.
        d += sign*np.outer(v, v)
    assert np.allclose(d@d@d, 3*d)
    assert abs(np.linalg.norm(d, 2)-math.sqrt(3)) < 1e-14
    values, vectors = np.linalg.eigh(l)
    p = vectors[:, 1:4]@vectors[:, 1:4].T
    reduced = sum((np.outer(vectors[:, k], vectors[:, k])/(values[k]-(5-math.sqrt(5)))
                   for k in [0, *range(4, 12)]), np.zeros((12, 12)))
    independent = -reduced@d@p-p@d@reduced
    exact = np.array([[c.q5_float(x) for x in row] for row in dp])
    assert np.linalg.norm(independent-exact, 2) < 1e-14
    assert c.q5_float(dg[0][1]) > .3


@pytest.mark.parametrize('theta', [F(-3, 25), F(0), F(3, 25), F(1, 8)])
def test_integer_spectral_enclosure_against_independent_eigensolver(theta):
    gram, error = geometry.enclosed_gram(theta)
    values, vectors = np.linalg.eigh(np.array(geometry.laplacian(theta), float))
    p = vectors[:, 1:4]@vectors[:, 1:4].T
    exact_gram = p/np.sqrt(np.outer(np.diag(p), np.diag(p)))
    assert values[4]-values[3] > math.sqrt(5)+.5
    assert np.min(np.diag(p)) > 1/20
    assert error < F(1, 10**13)
    for i in range(12):
        for j in range(12):
            lo, hi = gram[i][j]
            assert float(lo)-1e-13 <= exact_gram[i, j] <= float(hi)+1e-13


def test_exact_finite_response_and_observation():
    result = geometry.certificate()
    assert F(result['finite_contrast_certified_lower']) > F(3, 40)
    assert result['numerical_eigenvectors_used_in_certificate'] is False
    assert F(packet.observation_certificate()['exp6_rational_lower']) > 400


@pytest.mark.parametrize('call', [lambda: source.selected_seam(True, 0),
    lambda: source.selected_seam(0, 256), lambda: source.record_probabilities(F(2)),
    lambda: geometry.laplacian(F(1, 7)), lambda: geometry.laplacian(0),
    lambda: geometry.enclosed_projector(F(0), precision=5)])
def test_invalid_interfaces_rejected(call):
    with pytest.raises(ValueError):
        call()


@lru_cache(maxsize=1)
def fresh_report():
    return packet.report()


def test_cold_receipt():
    value = fresh_report()
    retained = packet.strict_json.load(packet.OUTPUT)
    packet.validate(retained, value)


@pytest.mark.parametrize('key', ['fresh_probe_preparation_and_independence_are_inputs',
    'mean_iteration_for_arbitrary_correlated_probe_reuse',
    'joint_stress_entropy_area_or_Einstein_premises_instantiated',
    'source_state_map_to_spatial_fermion_core_or_defect_family',
    'physical_gate_accepted'])
def test_scope_mutations_rejected(key):
    value = deepcopy(fresh_report())
    value['scope'][key] = not value['scope'][key]
    with pytest.raises(ValueError):
        packet.validate(value, fresh_report())


def test_finite_contrast_mutation_rejected():
    value = deepcopy(fresh_report())
    value['geometry']['G01_energy0_minus_energy1_strict_lower'] = '1/5'
    with pytest.raises(ValueError):
        packet.validate(value, fresh_report())


def test_only_the_four_bounded_numeric_leakage_diagnostics_are_portable():
    replay = fresh_report()
    for resource in ('charged_endpoint_preparation_resources',
                     'interaction_and_scalar_preparation_resources',
                     'one_fair_coin_word_resources', 'pointer_resources'):
        value = deepcopy(replay)
        value['quantum_source'][resource]['maximum_numerical_valid_code_leakage'] = 1e-9
        packet.validate(value, replay)
        packet.validate(replay, value)
        value['quantum_source'][resource]['integrated_processor_norm_upper'] += 1e-12
        with pytest.raises(ValueError, match='fresh replay'):
            packet.validate(value, replay)
        value = deepcopy(replay)
        row = value['quantum_source'][resource]
        row['unregistered_leakage'] = row.pop('maximum_numerical_valid_code_leakage')
        with pytest.raises(ValueError, match='declared numerical'):
            packet.validate(value, replay)


def test_leakage_admission_rejects_bad_values_on_both_sides():
    replay = fresh_report()
    for invalid in (-1e-16, 1.0001e-9, float('nan'), float('inf'), 0, False, '0'):
        value = deepcopy(replay)
        value['quantum_source']['pointer_resources']['maximum_numerical_valid_code_leakage'] = invalid
        with pytest.raises(ValueError, match='admission'):
            packet.validate(value, replay)
        with pytest.raises(ValueError, match='admission'):
            packet.validate(replay, value)
