"""Countertuning, all-valley and branch controls for the band-edge limit."""
from fractions import Fraction as F
import itertools
import math
import numpy as np
import pytest
from scipy.integrate import quad
from ..contact_scattering.onshell import angular_density
from . import limit as c


def test_exact_majorants():
    b = c.green_certificate()
    assert b['inside_W_remainder_coefficient'] == F(6568, 189) < 36
    assert b['outside_scaled_green_coefficient'] < 11
    assert b['inside_scaled_green_coefficient'] < 19
    assert b['edge_valleys'] == 8


@pytest.mark.parametrize('gamma,cap,eps', [
    (F(1), F(1), F(1, 1000)), (F(1, 100), F(100), F(1, 10)),
    (F(100), F(1, 100), F(2)), (F(3, 7), F(5, 11), F(1, 3))])
def test_schedule(gamma, cap, eps):
    cert = c.spacing_certificate(gamma, cap, eps)
    a = cert['spacing_upper']
    assert a*gamma <= 1 and a*a*cap <= F(1, 8)
    assert 836*a*cap <= 7*gamma
    assert cert['uniform_scattering_error_upper'] <= eps


@pytest.mark.parametrize('args', [(F(1), F(2), F(1)), (F(1, 10), F(1), F(1)),
                                  (F(0), F(1), F(1)), (0.001, F(1), F(1)),
                                  (F(1, 1000), F(0), F(1)), (True, 1, 1)])
def test_scattering_rejects_unproved_parameters(args):
    with pytest.raises(ValueError): c.scattering_certificate(*args)


def test_pole_boundary_and_no_gamma_zero_claim():
    p = c.pole_certificate(F(1, 256), 1)
    assert p['pole_offset_lower'] == F(1, 4)
    assert p['pole_offset_upper'] == 4
    for a, g in [(F(1, 255), 1), (F(1, 256), 0), (0.001, 1)]:
        with pytest.raises(ValueError): c.pole_certificate(a, g)


def test_source_contact_spin_resolvent():
    sigma = [np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]),
             np.array([[1, 0], [0, -1]])]
    singlet = np.array([0, 1, -1, 0])/math.sqrt(2)
    rng = np.random.default_rng(9473)
    for k in rng.uniform(-math.pi, math.pi, (12, 3)):
        speed = 1.7; energy = 9.0
        h = speed*sum(np.sin(k[j])*sigma[j] for j in range(3))
        free = np.kron(h, np.eye(2))-np.kron(np.eye(2), h)
        actual = np.vdot(singlet, np.linalg.solve(energy*np.eye(4)-free, singlet))
        expected = energy/(energy**2-4*speed**2*np.sum(np.sin(k)**2))
        assert abs(actual-expected) < 3e-16


def test_full_torus_green_identity_and_W0_bounds():
    n = 64
    q = np.arange(n)*2*math.pi/n
    C = (np.cos(q)**2)[:, None, None]+(np.cos(q)**2)[None, :, None]+(np.cos(q)**2)[None, None, :]
    eta = .4
    direct = np.mean(1/(C+eta))
    assert abs(c.W_diagnostic(eta).real-direct) < 1e-10
    assert 2/3 < c.W0_diagnostic() < 5
    a = .25; z = 2.; speed = math.sqrt(3)/a**2; E = 6/a**2+z
    source = np.mean(E/(E**2-4*speed**2*(3-C)))
    transformed = (a*a/2)*(1+a*a*z/6)*np.mean(1/(C+a*a*z+a**4*z*z/12))
    assert abs(source-transformed) < 1e-15


@pytest.mark.parametrize('eta', [.1, .01, .001, -.1, -.01, -.001])
def test_threshold_bound_and_all_valley_coefficient(eta):
    actual = c.W_diagnostic(eta)-c.W0_diagnostic()
    root = math.sqrt(eta) if eta>0 else 1j*math.sqrt(-eta)
    coefficient = 18 if eta>0 else 36
    assert abs(actual+2*root/math.pi) <= coefficient*abs(eta)
    # One-valley replacement has the wrong leading singularity.
    if abs(eta)==.001:
        assert abs(actual+root/(4*math.pi)) > .015


@pytest.mark.parametrize('t', [.001, .03, .125])
def test_imaginary_part_matches_independent_spherical_shell(t):
    # Convert the independent source contact density to spherical mean J.
    angular_mean = angular_density(2*math.sqrt(t), speed=1)*math.pi**2/t
    expected = -2*math.sqrt(t)*angular_mean/math.pi
    assert abs(c.W_diagnostic(-t).imag-expected) < 1e-11


@pytest.mark.parametrize('a', [F(1, 256), F(1, 1024), F(1, 4096)])
def test_scattering_and_pole_with_exact_bounds(a):
    cert = c.scattering_certificate(a, 1, 1)
    for s in [.001, .1, 1.]:
        out = c.scattering_diagnostic(float(a), 1., s)
        assert out['absolute_error'] < float(cert['uniform_scattering_error_upper'])
        assert abs(abs(out['source_phase'])-1) < 1e-14
        assert out['reflected_phase'] == out['source_phase'].conjugate()
    out = c.scattering_diagnostic(float(a), 1., 1.)
    assert abs(out['source_phase']-1) > 1
    assert abs(out['source_phase']-out['reflected_target']) > 1.99
    pole = c.pole_diagnostic(float(a), 1.)
    pcert = c.pole_certificate(a, 1)
    assert .25 < pole['offset'] < 4
    assert pole['absolute_error'] < float(pcert['pole_offset_error_upper'])


def test_without_countertuning_scattering_becomes_identity():
    errors=[]
    for a in [1/64, 1/256, 1/1024]:
        # U=1/a^2, separated from critical 2/(W0*a^2).
        D = 1-c.W0_diagnostic()/2-a*c.scaled_green_diagnostic(a, -1)
        errors.append(abs(D.conjugate()/D-1))
    assert errors[2] < errors[1] < errors[0] < .1


def test_all_eight_quadratic_edges_and_error_bound():
    p = (F(1, 3), F(-2, 5), F(3, 7)); a = F(1, 32)
    cert = c.dispersion_certificate(a, p)
    for signs in itertools.product([-1, 1], repeat=3):
        k = np.array(signs)*math.pi/2 + float(a)*np.array(p, dtype=float)
        deficit = (6-2*math.sqrt(3)*np.linalg.norm(np.sin(k)))/float(a*a)
        assert abs(deficit-float(cert['limit_energy'])) < float(cert['absolute_error_upper'])
    with pytest.raises(ValueError): c.dispersion_certificate(a, [1, 2])
    with pytest.raises(ValueError): c.dispersion_certificate(a, [1, 2, .5])


def test_principal_value_monomial_and_coefficient_majorant():
    R = .5; t = .08; root = math.sqrt(t)
    for n in range(1, 6):
        # Cauchy-weight integration independently evaluates the PV.
        actual = quad(lambda r: r**(2*n)/(r+root), 0, R,
                      weight='cauchy', wvar=root, epsabs=1e-12)[0]
        expected = sum(t**j*R**(2*(n-1-j)+1)/(2*(n-1-j)+1) for j in range(n))
        expected -= t**(n-.5)*math.atanh(root/R)
        assert abs(actual-expected) < 1e-12
