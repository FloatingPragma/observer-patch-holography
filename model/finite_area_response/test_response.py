"""Independent full quantum, complete-record and arithmetic response checks."""
import copy
from fractions import Fraction as F

import mpmath as mp
import numpy as np
import pytest
from scipy.linalg import block_diag, expm, expm_frechet
import sympy as s

from model.finite_area_response import response as r, packet
from model.native_equilibration import reciprocal_vacuum as source
from tools import strict_json


def test_source_action_and_independent_metric_are_exact():
    for matrix, geometry in zip(r.matrices(), r.GEOMETRIES, strict=True):
        actual = source.vacuum_repair.hamiltonian(2, 2, geometry)
        np.testing.assert_allclose(np.asarray(matrix, float), float(source.BETA)*actual,
                                   rtol=0, atol=1e-15)
    assert r.matrices()[0]*r.matrices()[1] != r.matrices()[1]*r.matrices()[0]
    metric = r.metric_certificate()
    assert metric['two_face_area'] == '2/r'
    assert metric['normalized_proper_volume'] == '1'
    assert metric['internal_edge_log_dimensions'] == ['0', '0']
    assert F(metric['branch_area_difference']) == F(40, 99)


def test_independent_LDL_spectrum_counts_and_controls():
    for matrix in r.matrices():
        ranges = r.spectral_enclosures(matrix)
        approximate = np.linalg.eigvalsh(np.asarray(matrix, float))
        assert len(ranges) == 8
        for i, (point, bound) in enumerate(zip(approximate, ranges, strict=True)):
            assert r.inertia(matrix, bound.lo) == i
            assert r.inertia(matrix, bound.hi) == i+1
            assert float(bound.lo)-2e-14 < point < float(bound.hi)+2e-14
    assert r.inertia(s.diag(-2, 1, 4), 0) == 1
    with pytest.raises(ValueError, match='zero LDL'):
        r.inertia(s.Matrix([[0, 1], [1, 0]]), 0)
    with pytest.raises(ValueError, match='symmetric'):
        r.inertia(s.Matrix([[1, 1], [0, 1]]), 0)


@pytest.mark.parametrize('point', [F(-3), F(-1, 7), F(0), F(1, 7), F(3)])
def test_outward_scalar_exponential_with_independent_high_precision(point):
    bound = r.exponential(point)
    with mp.workdps(100):
        actual = mp.exp(mp.mpf(point.numerator)/point.denominator)
        lo = mp.mpf(bound.lo.numerator)/bound.lo.denominator
        hi = mp.mpf(bound.hi.numerator)/bound.hi.denominator
        assert lo <= actual <= hi
        assert hi-lo < mp.mpf('1e-28')


def test_spectral_free_matrix_exponential_and_occupation_enclosures():
    with mp.workdps(90):
        for matrix in r.matrices():
            diagonal, _, _ = r.matrix_exponential_diagonal(matrix)
            exact = mp.matrix([[mp.mpf(int(x.p))/int(x.q) for x in row]
                               for row in matrix.tolist()])
            independent = mp.expm(3*mp.eye(8)-exact)
            for i, interval in enumerate(diagonal):
                lo = mp.mpf(interval.lo.numerator)/interval.lo.denominator
                hi = mp.mpf(interval.hi.numerator)/interval.hi.denominator
                assert lo <= independent[i, i] <= hi


def test_fixed_law_response_on_full_quantum_density_and_all_records():
    model = source.finite_model()
    p = model.probabilities.reshape(2, 8).sum(axis=1)
    area = np.repeat([20/9, 20/11], 8)
    rho = np.asarray(model.gibbs)
    reset = block_diag(p[0]*np.eye(8)/8, p[1]*np.eye(8)/8)
    epsilon = 1/3
    changed = (1-epsilon)*rho+epsilon*reset
    before_geometry = np.array([np.trace(changed[:8, :8]).real, np.trace(changed[8:, 8:]).real])
    np.testing.assert_allclose(before_geometry, p, rtol=0, atol=1e-15)
    final = source.apply_channel(model, changed)
    expected = r.certificate()['fixed_law_state_response']
    slope = np.dot(area, np.diag(final-rho).real)/epsilon
    assert float(F(expected['area_slope'][0]))-1e-13 < slope < float(F(expected['area_slope'][1]))+1e-13
    assert -0.000628 < slope < -0.000627
    for law in ('C', 'E'):
        records, weights, _, _ = source.records(model, law)
        assert len(records) == 512 and set(records[:, 7]) == {0, 1}
        transition = np.zeros((16, 16))
        for row, weight in zip(records, weights, strict=True):
            transition[row[6], row[0]] += weight
        np.testing.assert_allclose(transition.sum(axis=0), np.ones(16), rtol=0, atol=2e-15)
        probabilities = (1-epsilon)*model.probabilities+epsilon*np.repeat(p/8, 8)
        actual_area = area@(transition@probabilities-model.probabilities)/epsilon
        assert actual_area == pytest.approx(slope, abs=2e-14)
    injection = np.trace(model.hamiltonian@(reset-rho)).real
    assert 2.697 < injection < 2.698
    assert float(F(expected['mean_intervention_energy_injection_per_unit_epsilon'][0]))-1e-12 < injection
    assert injection < float(F(expected['mean_intervention_energy_injection_per_unit_epsilon'][1]))+1e-12


def test_complete_randomization_Kraus_family_preserves_geometry_and_trace():
    model = source.finite_model()
    rho = np.asarray(model.gibbs)
    output = np.zeros((16, 16), complex)
    completeness = np.zeros((8, 8))
    for i in range(8):
        for j in range(8):
            unit = np.zeros((8, 8))
            unit[i, j] = 1
            completeness += unit.T@unit/8
            operator = np.kron(np.eye(2), unit)
            output += operator@rho@operator.T/8
    np.testing.assert_array_equal(completeness, np.eye(8))
    p = model.probabilities.reshape(2, 8).sum(axis=1)
    np.testing.assert_allclose(output, block_diag(p[0]*np.eye(8)/8, p[1]*np.eye(8)/8), atol=2e-16, rtol=0)
    assert r.randomization_certificate()['replacement_Kraus_records'] == 64


def test_noncommuting_equilibrium_derivative_with_Frechet_and_secant_oracles():
    model = source.finite_model()
    beta = float(model.beta)
    occupation = np.diag([0]*4+[1]*4)
    perturbation = block_diag(occupation, occupation)
    area = np.diag(np.repeat([20/9, 20/11], 8))
    thermal, derivative = expm_frechet(-beta*model.hamiltonian, -beta*perturbation)
    partition = np.trace(thermal)
    tangent = derivative/partition-thermal*np.trace(derivative)/partition**2
    actual = np.trace(area@tangent).real
    report = r.equilibrium_response()
    assert float(report['response'].lo)-1e-14 <= actual <= float(report['response'].hi)+1e-14
    assert -0.000675 < actual < -0.000673
    wrong = -beta*perturbation@model.gibbs+beta*np.trace(perturbation@model.gibbs)*model.gibbs
    assert np.linalg.norm(wrong-wrong.conj().T) > 0.001
    assert np.linalg.norm(wrong-tangent) > 0.0005
    step = 1e-5
    states = []
    for sign in (1, -1):
        exp = expm(-beta*(model.hamiltonian+sign*step*perturbation))
        states.append(exp/np.trace(exp))
    measured = np.trace(area@(states[0]-states[1])).real/(2*step)
    assert measured == pytest.approx(actual, abs=3e-11)


def test_zero_response_controls_and_incomplete_spectrum_rejection():
    equal = ((r.Interval.point(1),)*8,)*2
    assert r.fixed_law_response(equal)['area'] == r.Interval.point(0)
    weights = tuple(r.spectral_weights(matrix)[1] for matrix in r.matrices())
    assert r.fixed_law_response(weights, (1, 1))['area'] == r.Interval.point(0)
    with pytest.raises(ValueError, match='complete'):
        r.fixed_law_response((weights[0][:-1], weights[1]))
    with pytest.raises(ValueError, match='faithful'):
        r.fixed_law_response(((r.Interval.point(0),)*8, weights[1]))
    with pytest.raises(ValueError, match='exact'):
        r.exponential(.5)
    with pytest.raises(ValueError, match='requires'):
        r.exponential(4)
    with pytest.raises(ValueError, match='division'):
        r.Interval.point(1)/r.Interval(-1, 1)


def test_cold_replay_and_scope_mutations():
    value = packet.report()
    assert packet.canonical(value) == packet.canonical(strict_json.loads(packet.OUTPUT.read_text()))
    for key in ('complete_interacting_QFT_or_Einstein_geometry_established',
                'same_reference_as_Weyl_particle_source'):
        bad = copy.deepcopy(value)
        bad['response'][key] = True
        with pytest.raises(ValueError, match='promotion'):
            packet.validate(bad)
    bad = copy.deepcopy(value)
    bad['response']['equilibrium_Hamiltonian_family_response']['requires_matching_bath_gaps_and_energy_readout_family'] = False
    with pytest.raises(ValueError, match='changed family'):
        packet.validate(bad)
