"""Independent symbol, coarea, series and stationary-scattering controls."""
from fractions import Fraction as F

import numpy as np
import pytest

from .onshell import (angular_density, pi_interval, scattering_eigenvalue,
                      sphere_coefficient, threshold_enclosure)
from .scattering import free_blocks
from ..bound_pair.fiber import SINGLET


def test_exact_coefficients_and_pi_enclosure():
    assert sphere_coefficient(0) == 1
    assert sphere_coefficient(1) == F(1, 2)
    assert sphere_coefficient(2) == F(11, 40)
    assert all(0 < sphere_coefficient(n) <= 1 for n in range(25))
    lo, hi = pi_interval()
    # Independent decimal truncation, not binary math.pi inside a tiny interval.
    known = F('3.14159265358979323846264338327950288419716939937510')
    assert lo < known+F(1, 10**50) and hi > known
    assert 0 < hi-lo < F(1, 10**43)


@pytest.mark.parametrize('energy',[F(1,64), F(1,16), F(1,8), F(3,16), -F(1,8)])
def test_exact_interval_against_independent_spherical_integral(energy):
    bound = threshold_enclosure(energy, degree=12)
    evaluated = angular_density(float(energy), order=80)
    # Account separately for floating quadrature roundoff at the tiny tails.
    assert float(bound['density_lower'])-3e-14 <= evaluated <= float(bound['density_upper'])+3e-14
    tighter = threshold_enclosure(energy, degree=24)
    assert tighter['density_lower'] >= bound['density_lower']
    assert tighter['density_upper'] <= bound['density_upper']
    assert bound['density_lower'] > 0
    assert bound['node_count'] == 8


def test_actual_four_spin_symbol_has_the_claimed_singlet_measure():
    k = np.array([[.23,.71,1.32], [.44,2.31,4.2], [2.51,.81,5.13]])
    c = .125
    blocks = free_blocks(k, np.zeros(3), c)
    for point, h in zip(k, blocks, strict=True):
        r = np.linalg.norm(np.sin(point))
        for z in (.17+.13j, -.31+.21j):
            actual = np.vdot(SINGLET, np.linalg.solve(z*np.eye(4)-h, SINGLET))
            assert abs(actual-z/(z*z-4*c*c*r*r)) < 2e-14
        energies, vectors = np.linalg.eigh(h)
        zero = np.abs(energies)<1e-12
        assert sum(zero) == 2
        assert np.linalg.norm(vectors[:,zero].conj().T@SINGLET)<1e-13
    # Each of eight nodes makes an equal contribution in the sine-coordinate map.
    base = np.array([.12,.24,.31])
    for mask in range(8):
        shifted = base+np.pi*np.array([(mask>>a)&1 for a in range(3)])
        assert abs(np.linalg.norm(np.sin(shifted))-np.linalg.norm(np.sin(base)))<1e-14


def test_rank_one_boundary_scattering_is_unitary_and_not_identity():
    # Independent stationary formula and optical theorem, with no assertion
    # that these freely chosen a,rho pairs are lattice boundary measurements.
    for a in (-2.,0.,1.,3.):
        for rho in (.001,.2,1.):
            for u in (0.,.1,1.):
                g = a-1j*np.pi*rho
                tau = u/(1-u*g)
                expected = 1-2j*np.pi*rho*tau
                actual = scattering_eigenvalue(a,rho,interaction=u)
                assert abs(actual-expected)<2e-15
                assert abs(abs(actual)-1)<3e-15
                assert abs(-tau.imag-np.pi*rho*abs(tau)**2)<2e-14
                if u:
                    assert abs(actual-1)>0
                else:
                    assert actual == 1


@pytest.mark.parametrize('bad',[0,F(1,4),F(-1,4),.125,True,float('nan')])
def test_threshold_endpoint_or_inexact_input_cannot_be_certified(bad):
    with pytest.raises(ValueError):
        threshold_enclosure(bad)


def test_small_energy_density_keeps_eight_nodes_and_quadratic_threshold():
    e,c = F(1,1024),F(1,8)
    result = threshold_enclosure(e,speed=c)
    pi_lo,pi_hi = pi_interval()
    leading_lo,leading_hi = e*e/(4*pi_hi*pi_hi*c**3),e*e/(4*pi_lo*pi_lo*c**3)
    assert result['density_lower'] > leading_lo
    assert result['density_upper'] < leading_hi*F(10001,10000)
    # Single-pocket deletion would remove seven eighths of the leading weight.
    assert result['density_lower'] > 7*leading_hi/8


def test_scattering_evaluation_cannot_emit_nonfinite_result():
    with pytest.raises(ValueError, match='floating evaluation domain'):
        scattering_eigenvalue(1e308,1e308,interaction=1e308)
    result = scattering_eigenvalue(1e300,1e300)
    assert np.isfinite(result) and abs(abs(result)-1)<1e-15


def test_retained_exact_threshold_receipt():
    import json
    from pathlib import Path
    from .onshell import report
    stored = json.loads(Path(__file__).with_name('onshell_results.json').read_text())
    current = report()
    for actual,expected in zip(current['rows'],stored['rows'],strict=True):
        assert actual['exact_enclosure'] == expected['exact_enclosure']
        assert abs(actual['independent_angular_quadrature']-expected['independent_angular_quadrature'])<2e-14
    assert current['physical_massless_particle_identified'] is False
