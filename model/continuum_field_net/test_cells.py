from fractions import Fraction as F
import math

import pytest

from . import cells


def test_cell_projection_contraction_and_refinement_error():
    value=cells.certificate()
    assert value['continuum_norm_squared']=='16/15'
    assert value['continuum_derivative_norm_squared']=='8/3'
    errors=[F(row['projection_error_squared']) for row in value['refinements']]
    assert all(a>b>0 for a,b in zip(errors,errors[1:]))
    assert errors[-1]<F(1,100)


def test_cell_integrals_against_independent_midpoint_quadrature():
    spacing=F(1,4);values=cells.cell_integrals((1,0,-1),-1,1,spacing)
    for n,integral in values.items():
        lo=max(-1,float((n-F(1,2))*spacing));hi=min(1,float((n+F(1,2))*spacing))
        count=1000;dx=(hi-lo)/count
        approximate=sum(1-(lo+(j+.5)*dx)**2 for j in range(count))*dx
        assert abs(approximate-float(integral))<1e-8


def test_normalization_is_integral_squared_divided_by_cell_length():
    values=cells.cell_integrals((1,),F(-1,2),F(1,2),F(1,4))
    assert cells.pairing(values,values,F(1,4))==F(7,8)
    # Point sampling would miss the half-filled boundary cells and violate
    # the exact conditional-expectation identity.
    assert sum(x*x for x in values.values())!=cells.pairing(values,values,F(1,4))


@pytest.mark.parametrize('args',[((),-1,1,1),((1,),1,-1,1),((1,),-1,1,0),((1,),-1,1,.5)])
def test_invalid_cell_data_rejected(args):
    with pytest.raises(ValueError):cells.cell_integrals(*args)
