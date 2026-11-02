from fractions import Fraction as F
import numpy as np
import pytest
import sympy as s

from model.weyl_geometric_modular import entropy, geometry


def test_conformal_flow_and_exact_radial_coefficients():
    result = geometry.diamond_certificate()
    assert result['conformal_Killing_components'] == 16
    R, r, M, v = s.symbols('R r M v', positive=True)
    # The isotropic quadratic e=M*r^2/2 saturates the Hessian-norm remainder.
    exact = s.integrate((2*s.pi/v)*4*s.pi*r*r*(R*R-r*r)/(2*R)*(M*r*r/2), (r, 0, R))
    assert s.simplify(exact-4*s.pi**2*M*R**6/(35*v)) == 0
    assert geometry.matter_remainder_pi_squared(F(1, 10), 1, F(1, 8)) == F(1, 1093750)
    assert geometry.matter_remainder_pi_squared(1, 2, F(1, 8)) == 8*geometry.matter_remainder_pi_squared(1, 2)


@pytest.mark.parametrize('args', [(0, 1), (1, -1), (1, 1, 0), (True, 1), (1, 1.0)])
def test_invalid_exact_geometry_inputs(args):
    with pytest.raises(ValueError): geometry.matter_remainder_pi_squared(*args)


def test_flat_metric_does_not_acquire_a_nonzero_einstein_response():
    assert geometry.flat_geometry_countercontrol(F(2, 3), F(5, 7)) == F(10, 21)
    with pytest.raises(ValueError): geometry.flat_geometry_countercontrol(0, 1)


def test_binary_first_law_bound():
    assert entropy.certificate()['example_epsilon_1_10_radius_1_2'] == '1/150'
    for denominator in (2, 5, 10, 100):
        epsilon = F(1, denominator)
        r = F(3, 4)
        e = float(epsilon)
        actual = ((1+e)*np.log1p(e)+(1-e)*np.log1p(-e))/2
        assert 0 <= actual <= float(entropy.relative_entropy_bound(epsilon, r))
    for args in [(1, 1), (F(1,2), F(1,3)), (0.1, F(1,2)), (True, F(1,2))]:
        with pytest.raises(ValueError): entropy.relative_entropy_bound(*args)


def test_nondiagonal_quantum_preparation_satisfies_dpi_and_nonzero_first_law():
    X = np.array([[0., 1.], [1., 0.]])
    Z = np.diag([1., -1.]); I = np.eye(2)
    sigma = (I+.3*Z)/2; direction = (.2*X+.1*Z)/2
    assert np.linalg.eigvalsh(sigma+direction).min() > 0
    assert np.linalg.eigvalsh(sigma-direction).min() > 0
    assert np.linalg.norm(sigma@direction-direction@sigma) > 0
    log_sigma = np.diag(np.log(np.diag(sigma)))
    k = -np.trace(direction@log_sigma)
    assert abs(k) > .01
    S0 = -np.trace(sigma@log_sigma)
    for eps in (.25, .1, -.1, .01):
        rho = sigma+eps*direction
        lam = np.linalg.eigvalsh(rho)
        S = -np.sum(lam*np.log(lam))
        D = -S-np.trace(rho@log_sigma)
        binary = ((1+eps)*np.log1p(eps)+(1-eps)*np.log1p(-eps))/2
        assert -1e-15 <= D <= binary
        assert abs((S-S0)-(eps*k-D)) < 1e-14
        assert abs(S-S0-eps*k) <= eps*eps/(2*(1-.5**2))
