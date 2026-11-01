from fractions import Fraction as F
from math import asin, factorial, log, nextafter, pi

import numpy as np
import pytest
from scipy.linalg import expm

from . import bounds


def annihilators(modes):
    matrices = []
    for mode in range(modes):
        a = np.zeros((1 << modes, 1 << modes), complex)
        for state in range(1 << modes):
            if state & (1 << mode):
                a[state ^ (1 << mode), state] = (-1)**((state & ((1 << mode)-1)).bit_count())
        matrices.append(a)
    return matrices


def finite_control():
    a = annihilators(4)
    eye = np.eye(16)
    n = [x.conj().T@x for x in a]
    h = np.kron(np.array([[0., 1.], [1., 0.]]), np.eye(2))
    T = sum(h[i, j]*a[i].conj().T@a[j] for i in range(4) for j in range(4))
    V = sum(n[2*x]@n[2*x+1]-(n[2*x]+n[2*x+1])/2 for x in range(2))
    d = [(2*n[i ^ 1]-eye)@a[i] for i in range(4)]
    return a, n, d, T, V


def test_full_fock_normal_ordering_eliminates_one_excitation_both_branches():
    a, n, d, T, _ = finite_control()
    energies, vectors = np.linalg.eigh(T)
    ground = vectors[:, 0]
    costs = energies-energies[0]
    assert abs(energies[0]+2) < 1e-12
    assert sum(np.isclose(costs, 0)) == 1
    for i in range(4):
        assert abs(np.vdot(ground, n[i]@ground)-.5) < 1e-12
        for operator in (d[i], d[i].conj().T):
            weights = abs(vectors.conj().T@operator@ground)**2
            assert abs(sum(weights)-.5) < 1e-12
            assert sum(weights[~np.isclose(costs, 3)]) < 1e-25
        # Omitting the centering leaves an actual one-excitation component.
        uncentered = n[i ^ 1]@a[i]
        weights = abs(vectors.conj().T@uncentered@ground)**2
        assert abs(sum(weights[np.isclose(costs, 1)])-.125) < 1e-12


def test_nine_cross_contractions_exactly_cancel_degree_six():
    a = annihilators(4)
    rng = np.random.default_rng(173)
    h = rng.normal(size=(4, 4))+1j*rng.normal(size=(4, 4))
    h = (h+h.conj().T)/2
    w = expm(1j*.37*h)
    evolved = [sum(w[i, j]*a[j] for j in range(4)) for i in range(4)]
    A, B, C = evolved[1].conj().T, evolved[1], evolved[0]
    D, E, G = a[2].conj().T, a[3].conj().T, a[3]
    anti = lambda x, y: x@y+y@x
    # Keep all nine terms in their noncommuting order. Each cross pair is
    # independently checked to be a scalar CAR contraction.
    for x in (A, B, C):
        for y in (D, E, G):
            z = anti(x, y)
            assert np.max(abs(z-np.trace(z)/16*np.eye(16))) < 1e-13
    one = lambda x: anti(x, D)@E@G-D@anti(x, E)@G+D@E@anti(x, G)
    expanded = A@B@one(C)-A@one(B)@C+one(A)@B@C
    assert np.max(abs(anti(A@B@C, D@E@G)-expanded)) < 1e-13
    constants = bounds.contraction_constants()
    assert constants['cubic_cubic'] == 288
    assert constants['cubic_linear']+constants['linear_cubic'] == 24
    assert constants['total'] == 312


def test_small_centered_full_fock_response_respects_microtime_normalization():
    _, _, d, T, V = finite_control()
    u, s = F(1, 10**6), F(1, 100)
    rng = np.random.default_rng(811)
    f = rng.normal(size=4)+1j*rng.normal(size=4)
    f /= np.linalg.norm(f)
    D = sum(f[i].conjugate()*d[i] for i in range(4))
    values = []
    for H in (T, T+float(u)*V):
        _, vectors = np.linalg.eigh(H)
        g = vectors[:, 0]
        e = expm(1j*float(s)*H)
        evolved = e@D@e.conj().T
        values.append(np.vdot(g, (evolved@D.conj().T+D.conj().T@evolved)@g))
    bound = bounds.microtime_upper(u, s)
    assert 0 < float(bound) < .001
    assert abs(values[1]-values[0]) < float(bound)


def test_phase_volume_and_gaussian_constants_have_independent_normalization():
    for e in (.001, .03, .1, .3, .7, .999, 1.):
        assert (2*asin(e)/pi)**3 <= e**3+1e-15
    # Integral t^7/delta^2 exp(-t^2/(2delta^2)) dt = 2^3 3! delta^6.
    assert 32*2**3*factorial(3) == 1536
    assert bounds.free_low_energy_upper(F(1, 2)) == F(1, 2)
    assert bounds.free_low_energy_upper(1) == 1
    assert bounds.free_gaussian_upper(F(1, 4)) == F(3, 8)


def test_exact_exponential_and_relative_sqrt_majorants():
    assert log(2) > 2/3
    assert bounds.exp_upper(0) == 1
    assert bounds.exp_upper(F(2, 3)) == 2
    tiny = F(1, 2**3000)
    upper = bounds.sqrt_upper(tiny)
    assert upper*upper >= tiny
    assert upper < 2*F(1, 2**1500)
    x = F(11, 8)
    upper = bounds.sqrt_upper(x)
    assert upper*upper >= x
    assert upper*upper-x < F(1, 2**70)


def test_uniform_gaussian_mass_bound_decays_on_exact_refinements():
    results = []
    for denominator in (4, 8, 16):
        delta = F(1, denominator)
        exponent = int(F(243, 4)/(delta*delta))
        c = 1 << (2*exponent)
        result = bounds.gaussian_low_energy_upper(c, 1, 1, delta)
        results.append(result)
        assert result < 4*1536*delta**6
    assert results[0] > results[1] > results[2] > 0
    assert results[-1] < F(1, 5000)


def test_free_endpoint_and_asymptotic_delta_are_not_underflow_proofs():
    assert bounds.gaussian_low_energy_upper(4, 0, 1, F(1, 7)) == F(1, 128)
    assert bounds.gaussian_low_energy_upper(4, 0, 0, F(1, 7)) == 0
    assert bounds.asymptotic_delta_diagnostic(4, 0, 1) is None
    assert bounds.asymptotic_delta_diagnostic(4, 1, 40) == 10
    assert bounds.asymptotic_delta_diagnostic(1e300, 1e-300, 1) > 0
    assert bounds.asymptotic_delta_diagnostic(1e300, nextafter(1e300, 0), 1) > 0
    assert bounds.microtime_upper(1, 0) == 0
    cert = bounds.certificate()
    assert cert['floating_delta_is_diagnostic_only'] is True
    assert all(x['exact_comparison'] for x in cert['exact_refinement_examples'])


def test_vacuum_duhamel_preserves_cutoff_and_kinetic_derivative_terms():
    assert bounds.vacuum_duhamel_upper(2, 1, 1, 2, 16, F(1, 16)) == F(1, 2)
    assert bounds.vacuum_duhamel_upper(2, 0, 1, 2, 16, 1) == 0
    assert bounds.vacuum_duhamel_upper(0, 1, 1, 2, 16, 1) == 0
    assert bounds.vacuum_duhamel_upper(2, 1, 0, 0, 16, 1) == 0
    errors = [bounds.vacuum_duhamel_upper(1, 2, 1, 3, n, F(1, n**4)) for n in (8, 16, 32)]
    assert errors[0] > errors[1] > errors[2] > 0


@pytest.mark.parametrize('call', [
    lambda: bounds.occupation_defect_upper(-1),
    lambda: bounds.microtime_upper(.1, 1),
    lambda: bounds.microtime_upper(True, 1),
    lambda: bounds.microtime_upper(1, -1),
    lambda: bounds.free_low_energy_upper(-1),
    lambda: bounds.free_gaussian_upper(0),
    lambda: bounds.gaussian_low_energy_upper(0, 1, 1, 1),
    lambda: bounds.gaussian_low_energy_upper(1, -1, 1, 1),
    lambda: bounds.gaussian_low_energy_upper(1, 1, -1, 1),
    lambda: bounds.gaussian_low_energy_upper(1, 1, 1, 0),
    lambda: bounds.asymptotic_delta_diagnostic(1, 1, 1),
    lambda: bounds.asymptotic_delta_diagnostic(float('inf'), 1, 1),
    lambda: bounds.asymptotic_delta_diagnostic(1, float('nan'), 1),
    lambda: bounds.asymptotic_delta_diagnostic(True, 0, 1),
    lambda: bounds.asymptotic_delta_diagnostic(1e-300, 1e-310, 1e300),
    lambda: bounds.sqrt_upper(1, True),
    lambda: bounds.sqrt_upper(1, 257),
    lambda: bounds.exp_upper(1_000_000),
    lambda: bounds.vacuum_duhamel_upper(1, 1, 1, 1, 0, 1),
    lambda: bounds.vacuum_duhamel_upper(1, 1, 1, 1, 1, 2),
    lambda: bounds.vacuum_duhamel_upper(1, 1, 0, 1, 1, 1),
])
def test_invalid_or_unsupported_parameters_fail_closed(call):
    with pytest.raises(ValueError):
        call()
