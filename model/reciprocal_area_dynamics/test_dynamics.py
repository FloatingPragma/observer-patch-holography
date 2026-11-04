"""Independent full-source trajectories, quantum entropies and path controls."""
import copy
from fractions import Fraction as F
from math import fsum

import mpmath as mp
import numpy as np
import pytest
from scipy.linalg import block_diag
import sympy as s

from . import bounds, dynamics, packet, records
from ..native_equilibration import reciprocal_vacuum as source
from tools import strict_json


@pytest.fixture(scope='module')
def model():
    return source.finite_model()


def entropy(state):
    eigenvalues = np.linalg.eigvalsh(state)
    assert eigenvalues.min() > -1e-13
    positive = eigenvalues[eigenvalues > 0]
    return float(-positive@np.log(positive))


def classical_entropy(probabilities):
    positive = probabilities[probabilities > 0]
    return float(-positive@np.log(positive))


def inside(value, encoded, tolerance=0):
    return float(F(encoded[0]))-tolerance <= value <= float(F(encoded[1]))+tolerance


def test_independent_high_precision_source_spectrum_and_exponentials():
    levels, weights, _, _ = bounds.source_data()
    with mp.workdps(160):
        def exact(value):
            value = F(value)
            return mp.mpf(value.numerator)/value.denominator
        for block, matrix in enumerate(bounds.prior.matrices()):
            numeric = mp.matrix([[exact(x) for x in row] for row in matrix.tolist()])
            values, _ = mp.eigsy(numeric)
            for index, value in enumerate(values):
                interval, weight = levels[8*block+index], weights[8*block+index]
                assert exact(interval.lo) <= value <= exact(interval.hi)
                assert exact(weight.lo) <= mp.exp(3-value) <= exact(weight.hi)
        for point in (F(-3), F(-1, 7), F(0), F(3)):
            interval = bounds.exp_bound(point)
            assert exact(interval.lo) <= mp.exp(exact(point)) <= exact(interval.hi)
    with pytest.raises(ValueError, match='domain'):
        bounds.exp_bound(4)
    with pytest.raises(ValueError, match='unresolved'):
        bounds.ONE/bounds.I(-1, 1)


def test_conditional_quantum_variation_and_full_energy_correction(model):
    pi = model.probabilities
    pstar = pi[8:].sum()
    tau = (model.gibbs[:8, :8]/(1-pstar), model.gibbs[8:, 8:]/pstar)
    energies = [float(np.trace(tau[g]@model.hamiltonian[8*g:8*g+8, 8*g:8*g+8]).real)
                for g in (0, 1)]
    beta, area = float(model.beta), np.repeat([20/9, 20/11], 8)
    eta = beta*(energies[1]-energies[0])/(20/9-20/11)
    assert .578 < eta < .579
    assert inside(eta, bounds.certificate()['conditional_thermal_area_coefficient'], 2e-14)
    constant = entropy(model.gibbs)+eta*(area@pi)
    for p in (.15, .4, pstar, .8):
        state = block_diag((1-p)*tau[0], p*tau[1])
        relative = (1-p)*np.log((1-p)/(1-pstar))+p*np.log(p/pstar)
        assert entropy(state)+eta*((1-p)*20/9+p*20/11) == pytest.approx(constant-relative, abs=2e-14)
    # A full coherent state carries the conditional energy term in (2).
    rng = np.random.default_rng(781)
    raw = rng.normal(size=(16, 16))+1j*rng.normal(size=(16, 16))
    state = raw@raw.conj().T+np.eye(16)
    state /= np.trace(state)
    p = float(np.trace(state[8:, 8:]).real)
    energy = float(np.trace(state@model.hamiltonian).real)
    expected_area = float(area@np.diag(state).real)
    excess = energy-((1-p)*energies[0]+p*energies[1])
    log_partition = -np.log(pi[0])-beta*model.energies[0]
    relative = -entropy(state)+beta*energy+log_partition
    assert entropy(state)+eta*expected_area-beta*excess == pytest.approx(constant-relative, abs=3e-14)
    assert abs(beta*excess) > .1


@pytest.mark.parametrize('law', ['C', 'E'])
def test_complete_source_records_and_exact_reversal(model, law):
    expected, reverse, _ = records.complete_labels(law)
    actual, weights, heat, _ = source.records(model, law)
    np.testing.assert_array_equal(actual, expected)
    assert len(reverse) == 512
    assert sum(actual[:, -1]) == 128
    assert sum(i == reverse[i] for i in range(512)) == 384
    for i, j in enumerate(reverse):
        assert reverse[j] == i
        assert actual[i, 0] == actual[j, 6]
        assert actual[i, 6] == actual[j, 0]
        assert heat[i] == -heat[j]
        assert np.log(weights[i]/weights[j]) == pytest.approx(-float(model.beta)*heat[i], abs=3e-15)


@pytest.mark.parametrize('law', ['C', 'E'])
def test_two_step_complete_path_fluctuation_and_heat_balance(model, law):
    rows, weights, heat, _ = source.records(model, law)
    _, reverse, _ = records.complete_labels(law)
    reverse = np.asarray(reverse)
    q = np.arange(1, 17, dtype=float)
    q /= q.sum()
    final = model.transition@model.transition@q
    starts = {i: np.flatnonzero(rows[:, 0] == i) for i in range(16)}
    first, second = [], []
    for i, row in enumerate(rows):
        for j in starts[row[6]]:
            first.append(i)
            second.append(j)
    first, second = np.asarray(first), np.asarray(second)
    assert len(first) == 16*32**2
    forward = q[rows[first, 0]]*weights[first]*weights[second]
    backward = final[rows[second, 6]]*weights[reverse[second]]*weights[reverse[first]]
    log_ratio = np.log(forward/backward)
    predicted = np.log(q[rows[first, 0]]/final[rows[second, 6]])-float(model.beta)*(heat[first]+heat[second])
    np.testing.assert_allclose(log_ratio, predicted, atol=4e-15, rtol=0)
    assert forward.sum() == pytest.approx(1, abs=5e-16)
    assert backward.sum() == pytest.approx(1, abs=5e-16)
    # Retain every path while avoiding a BLAS reduction-order error in the
    # sixteen-thousand-term normalization diagnostic.
    assert fsum(forward*np.exp(-log_ratio)) == pytest.approx(1, abs=5e-16)
    mean = float(forward@log_ratio)
    heat_total = float(forward@(heat[first]+heat[second]))
    assert heat_total == pytest.approx(model.energies@(final-q), abs=5e-15)
    assert mean == pytest.approx(classical_entropy(final)-classical_entropy(q)-float(model.beta)*heat_total, abs=5e-15)
    assert mean > .2


def test_memory_reconstructs_complete_source_area_and_energy(model):
    pi = model.probabilities
    p = np.array([pi[:8].sum(), pi[8:].sum()])
    inputs = [np.repeat(p/8, 8), np.arange(1, 17, dtype=float)/136, np.eye(16)[3]]
    source_rows, source_weights, _, _ = source.records(model)
    actual = np.zeros((16, 16))
    for row, probability in zip(source_rows, source_weights, strict=True):
        actual[row[6], row[0]] += probability
    for initial in inputs:
        result = dynamics.projected_trajectory(initial, 24, model=model)
        q = initial.copy()
        for n in range(25):
            assert result['geometry1_probability'][n] == pytest.approx(q[8:].sum(), abs=4e-15)
            assert result['energy'][n] == pytest.approx(model.energies@q, abs=5e-14)
            assert result['area'][n] == pytest.approx(np.repeat([20/9, 20/11], 8)@q, abs=8e-15)
            q = actual@q
    reset = dynamics.projected_trajectory(inputs[0], 8, model=model)
    assert reset['area'][1]-reset['area'][0] == pytest.approx(-.0006270566389034, abs=2e-15)
    assert reset['area'][7] < reset['area'][1]
    assert reset['memory_kernel'][0] > .0217
    # Dropping memory is a changed process, including on the conditional family.
    shifted = pi.copy()
    shifted[:8] *= .6/p[0]
    shifted[8:] *= .4/p[1]
    result = dynamics.projected_trajectory(shifted, 2, model=model)
    x = result['geometry1_probability']-p[1]
    assert abs(x[2]-result['instantaneous_coefficient']*x[1]) > .0005


def test_positive_memory_modes_and_rational_rank_certificate(model):
    result = bounds.certificate()
    assert result['minimal_scalar_memory_dimension'] == 14
    assert result['visible_nonconstant_relaxation_poles'] == 15
    assert all(F(lo) > 0 for lo, _ in result['orthogonal_memory_Gram_pivots'])
    last = result['orthogonal_memory_Gram_pivots'][-1]
    assert F(last[1])-F(last[0]) < F(11, 10**78)
    pi = model.probabilities
    root = np.sqrt(pi)
    p = pi[8:].sum()
    f = np.repeat([0., 1.], 8)-p
    unit = root
    normalized = root*f/np.sqrt(p*(1-p))
    symmetric = model.transition*root[None, :]/root[:, None]
    Q = np.eye(16)-np.outer(unit, unit)-np.outer(normalized, normalized)
    C = Q@symmetric@Q
    values, vectors = np.linalg.eigh((C+C.T)/2)
    assert np.max(np.abs(values[:2])) < 1e-15
    assert np.min(values[2:]) > .5
    assert np.max(values[2:]) < .9
    b = Q@symmetric@normalized
    weights = np.square(vectors.T@b)
    assert np.min(weights[2:]) > 1e-9
    for n, certified in enumerate(result['memory_kernel_first_eight']):
        assert inside(float(weights@values**n), certified, 2e-16)
    # Every exact minorization entry has a large numerical margin.
    assert np.min(model.transition@model.transition) > 1/131072
    assert np.linalg.eigvalsh(symmetric)[0] > .499999999999


def test_quantum_dephasing_and_same_process_entropy_area_heat_identity(model):
    rng = np.random.default_rng(134)
    raw = rng.normal(size=(16, 16))+1j*rng.normal(size=(16, 16))
    sigma = raw@raw.conj().T+np.eye(16)
    sigma /= np.trace(sigma)
    q = np.diag(model.frame.T@sigma@model.frame).real
    coherent_cost = classical_entropy(q)-entropy(sigma)
    assert coherent_cost > .1
    final = sigma.copy()
    for _ in range(3):
        final = source.apply_channel(model, final)
    qn = np.linalg.matrix_power(model.transition, 3)@q
    np.testing.assert_allclose(final, model.frame@np.diag(qn)@model.frame.T, rtol=0, atol=3e-16)
    beta, pi = float(model.beta), model.probabilities
    heat = float(np.trace(model.hamiltonian@(final-sigma)).real)
    production = entropy(final)-entropy(sigma)-beta*heat
    relative_loss = float(q@np.log(q/pi)-qn@np.log(qn/pi))
    assert production == pytest.approx(coherent_cost+relative_loss, abs=7e-15)
    pstar = pi[8:].sum()
    e0 = pi[:8]@model.energies[:8]/(1-pstar)
    e1 = pi[8:]@model.energies[8:]/pstar
    eta = beta*(e1-e0)/(20/9-20/11)
    delta_p = float(np.trace(final[8:, 8:]-sigma[8:, 8:]).real)
    delta_area = (20/11-20/9)*delta_p
    delta_excess = heat-(e1-e0)*delta_p
    assert entropy(final)-entropy(sigma)+eta*delta_area-beta*delta_excess == pytest.approx(production, abs=2e-15)
    assert production > 0


def test_flat_energy_control_removes_memory_and_invalid_inputs_refuse():
    model = source.finite_model((s.zeros(8), s.zeros(8)))
    initial = np.arange(1, 17, dtype=float)/136
    result = dynamics.projected_trajectory(initial, 8, model=model)
    assert np.max(abs(result['memory_kernel'])) < 1e-28
    assert np.max(abs(result['matter_forcing'])) < 1e-15
    assert result['instantaneous_coefficient'] == .5
    for state, count in ((np.ones(15)/15, 1), (np.ones(16), 1), (np.ones(16)/16, -1),
                         (np.ones(16)/16, True)):
        with pytest.raises(ValueError):
            dynamics.projected_trajectory(state, count, model=model)


def test_cold_receipt_replay_and_scope_controls():
    value = packet.report()
    assert packet.canonical(value) == packet.canonical(strict_json.loads(packet.OUTPUT.read_text()))
    for key in ('conditional_thermal_manifold_invariant_under_actual_law',
                'thermal_area_coefficient_is_universal_Newton_normalization',
                'regional_G4_or_Einstein_dynamics_proved'):
        wrong = copy.deepcopy(value)
        wrong['response'][key] = True
        with pytest.raises(ValueError, match='promotion'):
            packet.validate(wrong)
    wrong = copy.deepcopy(value)
    wrong['records']['laws']['E']['complete_records'] -= 1
    with pytest.raises(ValueError, match='branches'):
        packet.validate(wrong)
