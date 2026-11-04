"""Independent curvature conventions and nonzero force-deletion controls."""
import pytest
import sympy as s

from . import source, ward
from .geometry import Geometry


def test_round_sphere_curvature_and_homogeneous_variational_equation():
    t, x, theta, phi = s.symbols('t x theta phi', real=True)
    geometry = Geometry((t, x, theta, phi), s.diag(-1, 1, 1, s.sin(theta)**2))
    assert s.trigsimp(geometry.scalar) == 2
    assert all(s.trigsimp(v) == 0 for v in geometry.divergence(geometry.einstein))
    result = ward.metric_variation_certificate()
    assert result['lapse_Euler_identity_exact'] and result['scale_Euler_identity_exact']


def test_curved_ward_keeps_improvement_gauge_and_yukawa_exchange():
    result = ward.curved_ward_certificate()
    assert result['exact_off_shell_residual'] == ['0']*4
    # For g=a(t)^2 eta and a=1+t², R=6a''/a³.
    assert s.Rational(result['curvature_at_point']) == s.Rational(2187, 250)
    terms = {key: s.Matrix(list(map(s.sympify, values)))
             for key, values in result['nonzero_force_terms'].items()}
    assert terms['total_divergence'] == terms['Higgs_EOM_term']+terms['nonabelian_gauge_force']+terms['Yukawa_exchange']
    for term in ('nonabelian_gauge_force', 'Yukawa_exchange', 'improvement_divergence'):
        assert terms[term] != s.zeros(4, 1)
    assert not result['field_equations_imposed']
    assert not result['scalar_sector_conserved_in_isolation']


def test_scalar_tensor_divergence_independently_from_coordinate_curvature():
    t, x, y, z = s.symbols('t x y z', real=True)
    geometry = Geometry((t, x, y, z), (1+t*t)**2*s.diag(-1, 1, 1, 1))
    F = 2+t+x*x
    divergence = geometry.divergence(2*geometry.curvature_response(F))
    expected = s.Matrix([-geometry.scalar*s.diff(F, q) for q in geometry.coordinates])
    assert all(s.cancel(a-b) == 0 for a, b in zip(divergence, expected, strict=True))
    assert divergence != s.zeros(4, 1)


def test_wald_and_volume_use_the_same_curvature_derivative():
    result = ward.wald_certificate()
    assert result['binormal_contraction'] == '-2*F'
    assert result['slice_curvature_contraction'] == '3*F/2'
    assert not result['entropy_of_a_source_state_identified']


def test_actual_higgs_action_and_two_sided_positive_variation():
    result = source.certificate()
    assert result['same_current_Higgs_equations_checked'] == 72
    assert result['family_Gram_commutator_nonzero']
    assert s.Rational(result['F_positive_lower']) > 2
    assert s.Rational(result['kernel_norm_squared_exact']) > 0
    x = s.Symbol('x', real=True)
    H, F, parameter = source.higgs_variation(x*x+1)
    assert s.expand(s.diff(F, parameter).subs(parameter, 0)) == x*x+1
    assert s.expand(F-(s.Rational(8, 3)-s.Rational(1, 6)*(H.H*H)[0])) == 0
    assert s.diff(F, parameter, 2) != 0  # first variation is not a finite affine F law


@pytest.mark.parametrize('parameter', [0, -1, .5, True])
def test_higgs_variation_does_not_divide_by_zero_or_accept_float_inputs(parameter):
    with pytest.raises(ValueError):
        source.higgs_variation(1, xi=parameter)
