from fractions import Fraction as F
import pytest
from . import bounds


def test_symbolic_modulus_never_exponentiates_large_powers():
    value=bounds.modulus(F(1,2**100),6)
    assert value['series_terms']>10**80
    assert value['mesh_is_strictly_positive'] is True
    assert value['mesh_float_value_supplied'] is False
    assert F(value['total_error_upper'])<F(value['epsilon'])


@pytest.mark.parametrize('epsilon,B',[(F(1),1),(F(1,2),6),(F(1,100),11)])
def test_exact_three_channel_budget(epsilon,B):
    v=bounds.modulus(epsilon,B)
    assert F(v['endpoint_error_upper'])==epsilon/4
    assert F(v['series_error_upper'])<=epsilon/4
    assert F(v['spatial_error_upper'])==epsilon/4
    n=v['series_terms'];d=F(v['delta'])
    assert n>=16/(epsilon*d)
    assert v['regularity']['K_12_exponent']==n*(2*n+1)


def test_regularity_recursion_matches_closed_form():
    k=2048
    for j in range(8):
        assert bounds.finite_K(j)==k
        k*=8192*12**(j+1)
    with pytest.raises(ValueError):bounds.finite_K(33)


def test_neumann_budget_and_monotone_shrinking():
    assert bounds.neumann_budget(bounds.RADIUS_RATIO)<F(1,4)
    assert bounds.normalized_error(bounds.RADIUS_RATIO/2)<bounds.normalized_error(bounds.RADIUS_RATIO)


@pytest.mark.parametrize('args',[(True,6),(0,6),(-1,6),(F(3,2),6),(F(1,2),False),(F(1,2),F(1,2)),(0.5,6)])
def test_invalid_modulus_inputs(args):
    with pytest.raises(ValueError):bounds.modulus(*args)


def test_fractional_boundary_constants():
    value=bounds.certificate()
    assert value['sharp_ball_Hs_to_Hs_over_12_constant']==2048
    assert value['positive_mesh_materialized_as_float'] is False
