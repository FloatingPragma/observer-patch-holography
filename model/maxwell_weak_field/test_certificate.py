from fractions import Fraction as Q

import pytest
import sympy as s

from . import certificate as c


def test_vacuum_classification_and_restricted_magnetic_hessian():
    result = c.vacuum_certificate()
    assert result['electric_invariant_diagonal_dimension'] == 1
    assert result['closed_flux_dimension'] == 19
    assert result['magnetic_vacuum_classification']['rotation_covariant_local_gauge_dimension_exact'] == 1
    _, C, _ = c.local.carrier()
    # Constant-flux extensions are invisible, while a physical face perturbation is detected.
    assert C.T*s.ones(20, 20)*C == s.zeros(30)
    single_face = s.diag(1, *([0]*19))
    assert C.T*single_face*C != s.zeros(30)


def test_rotation_without_electric_diagonality_does_not_force_one_coefficient():
    _, C, rotations = c.local.carrier()
    electric = s.eye(30)+C.T*C
    assert all(U.T*electric*U == electric for U, _, _ in rotations)
    assert any(electric[i, j] != 0 for i in range(30) for j in range(30) if i != j)


def test_quartic_energy_force_and_constraint_identities():
    D, C, _ = c.local.carrier()
    epsilon, beta, lam = s.Rational(7, 3), s.Rational(11, 5), s.Rational(2, 7)
    a = s.symbols('a', real=True)
    for seed in range(1, 5):
        E = s.Matrix([s.Rational((i+seed)**2 % 13-6, 17) for i in range(30)])
        B = C*s.Matrix([s.Rational((i*seed+1) % 11-5, 19) for i in range(30)])
        energy, F = c.quartic_field(E, B, epsilon, beta, lam)
        # Independently differentiate the actual energy along its vector field.
        shifted_energy, _ = c.quartic_field(E+a*F[:30, :], B+a*F[30:, :], epsilon, beta, lam)
        assert s.diff(shifted_energy, a).subs(a, 0) == 0
        assert D.T*F[:30, :] == s.zeros(12, 1)
        assert sum(F[30:, :]) == 0
        G, K, _ = c.local.evolution(epsilon, beta)
        x = E.col_join(B)
        remainder = F-G*x
        norm_squared = (x.T*K*x)[0]
        force_squared = (remainder.T*K*remainder)[0]
        assert force_squared <= 6*lam**2/(epsilon*beta**3)*norm_squared**3
        assert norm_squared/2 <= energy <= norm_squared/2+lam*norm_squared**2/(4*beta**2)


def test_full_nonlinear_action_is_gauge_local_and_rotational():
    D, C, rotations = c.local.carrier()
    A = s.Matrix([s.Rational(i % 7-3, 11) for i in range(30)])
    chi = s.Matrix([s.Rational(i-5, 13) for i in range(12)])

    def action(potential):
        B = C*potential
        return sum(b*b/2+b**4/4 for b in B)

    assert action(A+D*chi) == action(A)
    assert all(action(U*A) == action(A) for U, _, _ in rotations)
    H = C.T*s.diag(*[1+3*b*b for b in C*A])*C
    assert c.local.is_local_gauge_hessian(H)
    assert H != C.T*C


def test_symbolic_nonlinearity_and_gauss_countercontrol():
    result = c.nonlinear_certificate()
    assert result['quartic_nonzero_cubic_correction_exact']
    counter = result['nonlinear_Gauss_counterexample']
    assert counter['Gauss_derivative'] == ['6', '-6', '0', '0', '-1', '1', '-1', '1', '0', '0', '0', '0']
    D, C, rotations = c.local.carrier()
    E, B = C.row(0).T, C[:, 0]
    F, skew = c.gyroscopic_field(E, B)
    assert D.T*E == s.zeros(12, 1)
    assert D.T*F[:30, :] != s.zeros(12, 1)
    # Reconstruct signed covariance of the entire skew operator, not only one output.
    for U, V, _ in rotations:
        _, rotated_skew = c.gyroscopic_field(U*E, V*B)
        assert rotated_skew == U*skew*U.T


def test_general_bound_exact_arithmetic_and_time_orientation():
    result = c.finite_time_bound(Q(1, 2), Q(3, 2), Q(1, 4), Q(1, 100), Q(-7, 3))
    assert result['energy_relative_radius'] == '1/24'
    assert result['energy_norm_squared_amplification'] == '25/23'
    assert result['error_norm_upper'] == str(Q(3, 4)*Q(25, 23)*Q(1, 10000)*Q(7, 3))
    assert result == c.finite_time_bound(Q(1, 2), Q(3, 2), Q(1, 4), Q(1, 100), Q(7, 3))
    assert c.finite_time_bound(0, 0, 1, Q(1, 2), 100)['error_norm_upper'] == '0'


@pytest.mark.parametrize('values', [
    (3, 1, 1, Q(1, 10), 1),  # quadratic lower bound vanishes
    (0, 1, 1, 1, 1),         # strict interior condition fails
    (1, 1, 1, Q(3, 4), 1),  # r<R alone is insufficient for confinement
    (-1, 1, 1, Q(1, 10), 1),
    (0.1, 1, 1, Q(1, 10), 1),
    (True, 1, 1, Q(1, 10), 1),
])
def test_general_bound_rejects_missing_or_rounded_hypotheses(values):
    with pytest.raises(ValueError):
        c.finite_time_bound(*values)


def test_quartic_global_bound_and_degenerate_boundary():
    bound = c.quartic_bound(Q(7, 3), Q(11, 5), Q(2, 7), Q(1, 100), Q(-3, 2))
    r, lam, beta = Q(1, 100), Q(2, 7), Q(11, 5)
    amplitude_squared = r*r+lam*r**4/(2*beta**2)
    assert Q(bound['trajectory_norm_squared_upper']) == amplitude_squared
    squared = 6*lam**2*Q(3, 2)**2*amplitude_squared**3/(Q(7, 3)*beta**3)
    assert Q(bound['error_norm_squared_upper']) == squared
    assert Q(bound['error_norm_upper'])**2 >= squared
    assert c.quartic_bound(1, 1, 0, 10, 100)['error_norm_upper'] == '0'
    with pytest.raises(ValueError):
        c.quartic_bound(1, 0, 1, Q(1, 100), 1)


def test_outward_root_handles_small_bounds_without_float_underflow():
    for value in (Q(0), Q(4), Q(2), Q(1, 10**80), Q(10**80, 7)):
        upper = c.sqrt_upper(value)
        assert upper**2 >= value
        if upper:
            assert (upper-Q(1, 2**48))**2 < value


def test_closed_flux_and_exact_input_admission():
    with pytest.raises(ValueError):
        c.quartic_field(s.zeros(30, 1), s.ones(20, 1))
    with pytest.raises(ValueError):
        c.quartic_field(s.zeros(30, 1), s.zeros(20, 1), electric=1.0)
    with pytest.raises(ValueError):
        c.quartic_field(s.I*s.ones(30, 1), s.zeros(20, 1))
