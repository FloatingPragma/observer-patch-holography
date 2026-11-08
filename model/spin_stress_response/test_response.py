from fractions import Fraction as F

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.linalg import expm
import sympy as s

from model.native_spin_clock_join import energy
from . import response


def test_source_exact_marginal_and_nonzero_response():
    result = response.exact_response()
    assert result['uniform_regional_floor'] == '1/5'
    assert result['entropy_slope_strictly_positive'] is True
    assert result['modular_generator_equals_physical_local_energy'] is False
    assert result['physical_energy_response'] == '0'


def test_driven_operator_identity_with_independent_direct_commutator():
    data = energy.source(3, ((0, 1), (1, 2)))
    a = (s.Rational(2, 3), s.Rational(-1, 2), s.Rational(3, 2))
    rates = (s.Rational(1, 3), s.Rational(0), s.Rational(-1, 4))
    driven = response.weighted_source(3, ((0, 1), (1, 2)), a, rates)
    h = sum((x*e for x, e in zip(a, data['density'], strict=True)), s.zeros(8))
    assert driven['hamiltonian'] == h
    for x in range(3):
        source = rates[x]*data['density'][x]
        left = source+s.I*(h*driven['density'][x]-driven['density'][x]*h)
        divergence = sum((s.I*(driven['density'][x]*driven['density'][y]
                               - driven['density'][y]*driven['density'][x])
                          for y in range(3)), s.zeros(8))
        assert left+divergence == source
    # Omitting the second weight does not satisfy the actual driven action.
    wrong = sum((a[0]*data['currents'][0, y] for y in (1, 2)), s.zeros(8))
    assert s.I*(h*driven['density'][0]-driven['density'][0]*h)+wrong != s.zeros(8)


def test_all_graph_controls_include_nonzero_currents():
    rows = response.exact_currents()
    assert len(rows) == 3
    assert all(r['nonzero_current_pairs'] > 0 for r in rows)


def matrices():
    # Independent tensor construction and a nonstationary, noncommuting input.
    x = np.array([[0, 1], [1, 0]], complex)
    y = np.array([[0, -1j], [1j, 0]], complex)
    z = np.diag([1, -1]).astype(complex)
    pair = (sum(np.kron(p, p) for p in (x, y, z))-np.eye(4))/4
    h01 = np.kron(pair, np.eye(2))
    h12 = np.kron(np.eye(2), pair)
    h = h01+h12
    v = 0.7*h01-0.2*h12
    a = 1j*(h01@h12-h12@h01)
    raw = np.arange(1, 9)+1j*np.arange(8, 0, -1)
    psi = raw/np.linalg.norm(raw)
    return h, v, a, np.outer(psi, psi.conj())


def test_retarded_sign_and_uniform_second_order_error():
    h, v, a, rho = matrices()
    time = 0.7
    evolve = lambda p: expm(-1j*(h+p*v)*time)
    expectation = lambda p: np.trace(rho@evolve(p).conj().T@a@evolve(p)).real
    ah = evolve(0).conj().T@a@evolve(0)

    def integrand(t):
        u = expm(-1j*h*t)
        vs = u.conj().T@v@u
        return (1j*np.trace(rho@(vs@ah-ah@vs))).real

    slope = quad(integrand, 0, time, epsabs=1e-12)[0]
    finite_difference = (expectation(1e-5)-expectation(-1e-5))/(2e-5)
    assert abs(slope-finite_difference) < 1e-9
    assert abs(slope) > 1e-3
    for p in (-2.0, -0.3, 0.05, 1.3):
        remainder = abs(expectation(p)-expectation(0)-p*slope)
        assert remainder <= 2*p*p*np.linalg.norm(a, 2)*time*time*np.linalg.norm(v, 2)**2
        delta = evolve(p)@rho@evolve(p).conj().T-evolve(0)@rho@evolve(0).conj().T
        distance = np.abs(np.linalg.eigvalsh(delta)).sum()/2
        assert distance <= min(1, abs(p)*time*np.linalg.norm(v, 2))+1e-14


def matrix_log(a):
    values, vectors = np.linalg.eigh(a)
    return (vectors*np.log(values))@vectors.conj().T


def test_noncommuting_modular_identity_and_remainder():
    sigma = np.diag([0.1, 0.2, 0.3, 0.4]).astype(complex)
    k = -matrix_log(sigma)
    rng = np.random.default_rng(1145)
    for _ in range(12):
        raw = rng.normal(size=(4, 4))+1j*rng.normal(size=(4, 4))
        delta = raw+raw.conj().T
        delta -= np.trace(delta)*np.eye(4)/4
        delta *= 0.04/np.abs(np.linalg.eigvalsh(delta)).sum()
        rho = sigma+delta
        epsilon = np.abs(np.linalg.eigvalsh(delta)).sum()/2
        assert np.linalg.norm(sigma@delta-delta@sigma) > 0
        ds = -np.trace(rho@matrix_log(rho)).real+np.trace(sigma@matrix_log(sigma)).real
        relative = np.trace(rho@(matrix_log(rho)-matrix_log(sigma))).real
        modular = np.trace(delta@k).real
        assert abs(ds-(modular-relative)) < 1e-14
        assert 0 < relative <= epsilon**2/(0.1-epsilon)


def test_entropy_slope_and_energy_countercontrol_independently():
    time = np.pi/4
    p = lambda lam: (1+0.6*np.cos((1+lam)*time))/2
    entropy = lambda lam: -p(lam)*np.log(p(lam))-(1-p(lam))*np.log(1-p(lam))
    r = 3*np.sqrt(2)/10
    slope = 3*np.sqrt(2)*np.pi/80*np.log((1+r)/(1-r))
    assert abs((entropy(1e-5)-entropy(-1e-5))/2e-5-slope) < 1e-10
    assert slope > 0
    for lam in (-1/16, 1/16):
        modular = -(p(lam)-p(0))*np.log(p(0)/(1-p(0)))
        residual = modular-(entropy(lam)-entropy(0))
        epsilon = abs(p(lam)-p(0))
        assert 0 <= residual <= epsilon**2/(0.2-epsilon)


@pytest.mark.parametrize('floor,error', [(F(0), F(0)), (F(1, 5), F(1, 5)),
                                         (F(1, 5), F(-1)), (0.2, F(0))])
def test_first_law_rejects_unfaithful_or_unsupported_bound(floor, error):
    with pytest.raises(ValueError):
        response.first_law_bound(floor, error)


def test_rational_error_budget_and_source_error():
    row = response.response_bounds(F(3, 4), F(1), F(1), F(1, 16), F(1, 1000))
    assert row['ideal_halftrace'] == F(3, 64)
    assert row['linear_expectation_remainder'] == F(9, 2048)
    assert row['complete_halftrace_with_source_error'] == F(3, 64)+F(1, 1000)
    assert response.first_law_bound(F(1, 5), F(3, 64)) == F(45, 3136)
    assert response.modular_stability(F(1, 5), F(1, 100)) == {
        'actual_regional_floor': F(19, 100), 'modular_generator_error': F(1, 19)}


@pytest.mark.parametrize('bad', [True, 0.5, s.sqrt(2)])
def test_exact_weight_input_contract(bad):
    with pytest.raises(ValueError):
        response.weighted_source(2, ((0, 1),), (bad, 1), (0, 0))
