from fractions import Fraction as F

import mpmath as mp
import pytest
import sympy as sy

from . import diagnostics, dispersion, source


def mp_fraction(value):
    return mp.mpf(value.numerator)/value.denominator


def contains(interval, value):
    return mp_fraction(F(interval['lower_exact'])) <= value <= mp_fraction(F(interval['upper_exact']))


def test_original_centered_source_and_translation_are_exactly_bound():
    data = source.certificate()
    assert data['complete_Fock_dimension'] == 256
    assert data['finite_volume_selection_rule_retained']
    assert 'pi*Delta N_down' in data['relative_momentum_offset']


@pytest.mark.parametrize('order', [0, 1])
@pytest.mark.parametrize('z', [F(1), F(4), F(8), F(32), F(128)])
def test_integrated_taylor_bounds_against_independent_bessel(order, z):
    interval = dispersion.bessel_ratio(order, z)
    with mp.workdps(85):
        x = mp_fraction(z)
        actual = mp.besselk(order, x)/(mp.sqrt(mp.pi/(2*x))*mp.exp(-x))
        assert mp_fraction(interval.lower) <= actual <= mp_fraction(interval.upper)


@pytest.mark.parametrize('A,theta', [(4,F(1)), (8,F(1,2)), (16,F(-1)), (64,F(0))])
def test_exact_dispersion_intervals_and_uniform_shell_bound(A, theta):
    point = dispersion.point(A, theta)
    with mp.workdps(75):
        base = mp.sqrt(mp.pi/(2*A))*mp.exp(-A)
        energy, momentum = mp.mpf(0), mp.mpf(0)
        for m in range(1, 81, 2):
            energy += mp.besselk(1, m*A)*mp.cosh(m*mp_fraction(theta))/(m*base)
            momentum += mp.besselk(0, m*A)*mp.sinh(m*mp_fraction(theta))/(m*base)
        assert contains(point['E_over_M'], energy)
        assert contains(point['v_p_over_M'], momentum)
        assert contains(point['mass_shell_over_M_squared'], energy**2-momentum**2)
        bound = dispersion.uniform_bounds(A, 1)
        assert abs(energy**2-momentum**2-1) < mp_fraction(F(bound['absolute_mass_shell_residual_over_M_squared']['upper_exact']))


def test_uniform_limit_is_quantitative_and_tightens():
    values = [F(dispersion.uniform_bounds(A,1)['absolute_mass_shell_residual_over_M_squared']['upper_exact'])
              for A in [4,8,16,32,64,128]]
    assert all(a > b > 0 for a,b in zip(values,values[1:]))
    assert values[-1] < F(17,1000)


def test_fixed_source_is_not_mislabeled_relativistic():
    value = diagnostics.fixed_source()
    assert value['floating_diagnostic_only']
    assert value['rows'][-1]['E_squared_minus_v_squared_p_squared']-value['rows'][0]['E_squared_minus_v_squared_p_squared'] > .027


def test_anchor_and_nonunique_trajectories_are_symbolic_identities():
    A = sy.symbols('A', positive=True)
    M = sy.symbols('M', positive=True)
    v, A0 = sy.Rational(1,8), sy.pi/8
    M0 = sy.exp(-A0)/sy.pi
    a = 4*v*sy.exp(-A)/(M0*sy.sqrt(2*sy.pi*A))
    assert sy.simplify(a.subs(A,A0)) == 1
    assert sy.simplify((v/(2*a)).subs(A,A0)) == sy.Rational(1,16)
    assert sy.simplify((sy.pi*v/(a*A)).subs(A,A0)) == 1
    other = a*M0/M*(1+(M/M0-1)*sy.exp(-(A-A0)))
    assert sy.simplify(other.subs(A,A0)) == 1
    assert sy.limit((4/sy.pi)*sy.sqrt((sy.pi*v/(other*A))*(v/(2*other)))*sy.exp(-A), A,sy.oo) == M
    assert dispersion.anchored_mass()['trajectory_selected_by_unchanged_source'] is False


@pytest.mark.parametrize('A,theta', [(0,0),(1,1),(4,5),(4,True)])
def test_wrong_rapidity_domains_rejected(A,theta):
    with pytest.raises(ValueError): dispersion.point(A,theta)
