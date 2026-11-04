"""Independent spherical integrals, geometric controls and remainder checks."""
import itertools

import pytest
import sympy as s

from . import balls


@pytest.mark.parametrize('powers', [(0, 0, 0), (2, 0, 0), (2, 2, 0),
                                   (4, 0, 2), (2, 2, 2), (3, 0, 0)])
def test_moments_against_direct_spherical_coordinate_integrals(powers):
    theta, phi = s.symbols('theta phi', real=True)
    n = (s.sin(theta)*s.cos(phi), s.sin(theta)*s.sin(phi), s.cos(theta))
    value = s.prod(v**p for v, p in zip(n, powers, strict=True))*s.sin(theta)
    direct = s.integrate(s.integrate(value, (phi, 0, 2*s.pi)), (theta, 0, s.pi))
    assert s.simplify(direct-4*s.pi*balls.angular_average(powers)) == 0


def test_exact_first_variation_across_independent_anisotropic_metric_directions():
    x, y, z = balls.X
    scalars = [s.S.One, x+y*z, x*x+2*z*z, x*y*z+y**3]
    for i, j in itertools.combinations_with_replacement(range(3), 2):
        seed = s.zeros(3)
        seed[i, j] = seed[j, i] = 1
        h = balls.radial_metric(seed, 1+x/3)
        for f in scalars:
            result = balls.variation(h, f, s.Rational(2, 5), s.Rational(7, 3))
            assert result['remainder'] == 0
            assert result['fixed_generalized_volume_variation'] == result['leading_variation']


def test_independent_warped_metric_area_volume_and_origin_curvature():
    # dr²+r²(1+s*k*r²)dOmega² has exactly linear area/volume variations.
    r, parameter, k = s.symbols('r parameter k', positive=True)
    warp = r*s.sqrt(1+parameter*k*r*r)
    scalar_curvature = -4*s.diff(warp, r, 2)/warp+2*(1-s.diff(warp, r)**2)/warp**2
    assert s.simplify(s.limit(scalar_curvature, r, 0)+18*parameter*k) == 0
    h = balls.radial_metric(s.eye(3)/2)
    result = balls.variation(h, 0, s.Rational(1, 3), 2)
    assert result['delta_Guu'] == -9
    assert result['delta_area'] == 4*s.pi/s.Integer(81)
    assert result['delta_volume'] == 4*s.pi/s.Integer(1215)


def test_ordinary_volume_keeps_scalar_value_and_wrong_laplacian_coefficient():
    radius, Fstar = s.Rational(1, 3), s.Rational(5, 2)
    constant = balls.variation(s.zeros(3), 7, radius, Fstar)
    assert constant['fixed_generalized_volume_variation'] == 0
    assert constant['fixed_proper_volume_variation'] == 16*s.pi**2*radius**2*7
    quadratic = balls.variation(s.zeros(3), balls.X[0]**2, radius, Fstar)
    assert quadratic['fixed_proper_volume_variation'] == s.Rational(5, 2)*quadratic['leading_variation']
    assert quadratic['leading_variation'] != 0


def test_quartic_remainder_and_fourth_directional_majorant():
    x, y, z = balls.X
    radius = s.Rational(1, 2)
    result = balls.variation(s.zeros(3), x**4, radius, 1)
    assert result['remainder'] == 64*s.pi**2*radius**6/35
    assert balls.fourth_derivative_bound(x**4, radius) == 24
    assert abs(result['remainder']) < balls.remainder_bound(radius, 1, 0, 24)
    polynomial = x**6+2*x*x*y*y+z**4/3
    upper = balls.fourth_derivative_bound(polynomial, radius)
    parameter = s.Symbol('parameter', real=True)
    for position, direction in [((0, 0, 0), (1, 0, 0)),
                                ((s.Rational(1, 4), 0, 0), (s.Rational(3, 5), s.Rational(4, 5), 0)),
                                ((0, s.Rational(1, 4), s.Rational(1, 4)), (0, 0, 1))]:
        line = polynomial.subs({q: a+parameter*d for q, a, d in zip(balls.X, position, direction, strict=True)})
        derivative = s.diff(line, parameter, 4).subs(parameter, 0)
        assert abs(derivative) <= upper
    assert balls.certificate()['higher_degree_controls']


@pytest.mark.parametrize('damage', ['radius', 'Fstar', 'floating', 'nonradial', 'linear_metric'])
def test_first_variation_domain_is_enforced(damage):
    h, f, radius, Fstar = s.zeros(3), 1, 1, 1
    if damage == 'radius':
        radius = 0
    elif damage == 'Fstar':
        Fstar = -1
    elif damage == 'floating':
        f = .5
    elif damage == 'nonradial':
        h = balls.X[0]**2*s.eye(3)
    else:
        h = balls.X[0]*s.eye(3)
    with pytest.raises(ValueError):
        balls.variation(h, f, radius, Fstar)
