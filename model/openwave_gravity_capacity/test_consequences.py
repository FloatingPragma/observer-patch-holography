"""Independent source evolution, derivative, radial integral and rejection checks."""
from copy import deepcopy
from fractions import Fraction as F

import numpy as np
import pytest
from scipy.linalg import expm
import sympy as s

from model.native_equilibration import reciprocal_vacuum as native
from . import capacity, packet, response


@pytest.fixture(scope='module')
def full_source():
    return native.finite_model()


def midpoint(interval):
    return float((F(interval[0])+F(interval[1]))/2)


def test_complete_source_update_reproduces_all_metric_responses(full_source):
    result = response.source_response()
    p = full_source.probabilities
    reset = np.r_[np.full(8, p[:8].sum()/8), np.full(8, p[8:].sum()/8)]
    for epsilon in (0, .125, .5, 1):
        prepared = (1-epsilon)*p+epsilon*reset
        updated = full_source.transition@prepared
        assert abs(updated.sum()-1) < 2e-15
        for row in result['readouts'].values():
            labels = np.repeat([float(F(v)) for v in row['branch_values']], 8)
            measured = labels@(updated-prepared)
            assert abs(measured-epsilon*midpoint(row['response_per_unit_randomization'])) < 3e-15
        # Postselection on cross-geometry proposals discards rejection records.
        mutant = full_source.transition-np.diag(np.diag(full_source.transition))
        conditional = mutant@prepared
        conditional /= conditional.sum()
        if epsilon:
            assert abs(conditional[8:].sum()-updated[8:].sum()) > 1e-4


def test_metric_dictionary_and_source_proper_time_evolution(full_source):
    z = s.symbols('r', positive=True)
    row = response.metric_readouts(z)
    assert row['metric'].det() == -z*z
    assert s.simplify(row['metric'][0, 0]+row['metric'][2, 2]*z**3) == 0
    for h, r in zip(full_source.exact_hamiltonians, (F(9, 10), F(11, 10)), strict=True):
        matrix = np.asarray(h, float)
        psi = np.zeros(8, complex); psi[0] = 1/np.sqrt(2); psi[7] = 1j/np.sqrt(2)
        elapsed = .31
        actual = expm(-1j*matrix*elapsed)@psi
        proper = expm(-1j*(matrix/float(r))*(float(r)*elapsed))@psi
        np.testing.assert_allclose(actual, proper, atol=3e-15, rtol=0)
        wrong = expm(-1j*matrix*(float(r)*elapsed))@psi
        assert np.linalg.norm(actual-wrong) > .02


def test_metric_response_controls(full_source):
    flat = native.finite_model((full_source.exact_hamiltonians[0],)*2)
    p = flat.probabilities
    reset = np.r_[np.full(8, p[:8].sum()/8), np.full(8, p[8:].sum()/8)]
    changed = flat.transition@reset-reset
    assert abs(np.repeat([.9, 1.1], 8)@changed) < 3e-15
    with pytest.raises(ValueError):
        response.metric_readouts(-1)


def test_gravity_bound_from_independent_radial_integration():
    # A constant positive error source on the ball attains the exterior bound.
    r, radius, eps = s.symbols('r R epsilon', positive=True)
    q = s.symbols('q', nonnegative=True)
    enclosed = s.integrate(2*eps*q*q, (q, 0, radius))
    force = enclosed/r**2
    potential = s.integrate(enclosed/q**2, (q, r, s.oo))
    result = response.newton_error_bound(2, F(1, 1000), F(1, 2000), F(1, 4000), 1, 2)
    eps_pi = F(result['pi_times_Einstein_component_error'])
    assert F(force.subs({eps:s.Rational(eps_pi), radius:1, r:2})) == F(result['pi_times_Phi_force_error'])
    assert F(potential.subs({eps:s.Rational(eps_pi), radius:1, r:2})) == F(result['pi_times_Phi_potential_error'])
    assert F(result['pi_times_Phi_force_error']) == 4*F(result['pi_times_Psi_force_error'])
    # Dropping the three spatial residuals can miss the exact worst-case factor four.
    assert F(result['pi_times_Phi_force_error']) > eps_pi/F(24)
    with pytest.raises(ValueError):
        response.newton_error_bound(2, 1, 1, 1, 2, 1)
    with pytest.raises(ValueError):
        response.newton_error_bound(0, 1, 1, 1, 1, 2)


def test_capacity_w_matches_independent_continuity_derivative():
    a = s.symbols('a', positive=True)
    for exponent in (-3, -1, 0, 1, 3):
        rho = s.Rational(7, 5)/a**exponent
        w = -1-a*s.diff(rho, a)/(3*rho)
        assert s.simplify(w-(-1+s.Rational(exponent, 3))) == 0
        wrong = -1-s.Rational(exponent, 3)
        if exponent:
            assert s.simplify(w-wrong) != 0
    result = capacity.symbolic_certificate()
    assert result['simultaneous_N_kappa_rescaling_invisible_to_density']


def test_epoch_invariant_is_amplitude_free_and_sensitive_to_history():
    om, radiation, od = F(1, 4), F(1, 16), F(11, 16)
    for amplitude in (F(1), F(17), F(10**122)):
        normalized = []
        for a in (F(1, 2), F(1), F(2)):
            n = amplitude*a*a
            e2 = capacity.expansion_from_capacity(a, n/amplitude, om, radiation, od)
            remainder = e2-om/a**3-radiation/a**4
            assert n*remainder == od*amplitude
            normalized.append(e2)
            assert capacity.capacity_from_expansion(a, e2, om, radiation, od) == a*a
        if amplitude == 1:
            baseline = normalized
        else:
            assert baseline == normalized
    assert capacity.expansion_from_capacity(F(1, 2), 1, om, radiation, od) != baseline[0]
    with pytest.raises(ValueError):
        capacity.capacity_from_expansion(1, F(1, 4), om, radiation, od)
    with pytest.raises(ValueError):
        capacity.two_epoch_density_ratio(0, 2)
    with pytest.raises(ValueError, match='sum to one'):
        capacity.expansion_from_capacity(1, 1, om, radiation, 1)


def test_linearized_curvature_uses_spatial_residuals():
    checks = response.linearized_curvature()['checks']
    assert all(checks.values())
    # Static metric potentials with equal G00 can have different accelerations.
    x, y, z = s.symbols('x y z', real=True)
    psi, phi = x*x+y*y+z*z, 2*(x*x+y*y+z*z)
    lap = lambda f: sum(s.diff(f, q, q) for q in (x, y, z))
    assert 2*lap(psi) == 12
    assert 2*lap(phi-psi) == 12
    assert s.diff(phi, x) != s.diff(psi, x)


def test_scope_and_replay_reject_promotions():
    actual = packet.report()
    packet.validate(actual, actual)
    for key in ('native_source_supplies_Newton_limit', 'numeric_cosmological_Lambda_selected'):
        changed = deepcopy(actual)
        changed['scope'][key] = True
        with pytest.raises(ValueError):
            packet.validate(changed, actual)
