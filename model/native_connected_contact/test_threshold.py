from fractions import Fraction as F
from itertools import combinations

import pytest
import sympy as sp

from . import threshold as t
from .test_dressing import control,dense


def test_exact_primitive_projection_and_rayleigh_formula_on_full_control():
    p=control();pairs,g=dense(p);v=len(p.neighbours)
    f=sp.Matrix([1,1,-1,-1,0,0]);a=(f.T*f)[0];fourth=sum(x**4 for x in f)
    incidence=sp.Matrix([[int(i in pair) for pair in pairs] for i in range(v)])
    q=sp.eye(len(pairs))-incidence.T*(incidence*incidence.T).inv()*incidence
    b=sp.Matrix([f[i]*f[j] for i,j in pairs])
    u=[(-f[i]**2+a/sp.Rational(2*(v-1)))/(v-2) for i in range(v)]
    h=sp.Matrix([f[i]*f[j]-u[i]-u[j] for i,j in pairs])
    assert h==q*b and incidence*h==sp.zeros(v,1)
    norm=a*a/2-fourth/2-(fourth-a*a/sp.Rational(2*(v-1)))/(v-2)
    assert (h.T*h)[0]==norm
    edges=[(i,j) for i,row in enumerate(p.neighbours) for j in row if i<j]
    energy=sum((f[i]-f[j])**2*(a-f[i]**2-f[j]**2) for i,j in edges)
    correction=sum((f[i]**2-f[j]**2)**2 for i,j in edges)/sp.Rational(v-2)
    assert (h.T*g*h)[0]==energy-correction
    assert incidence*b != sp.zeros(v,1)  # Unprojected product is not the claimed trial.


def test_actual_exact_affine_trials_and_capacity_collapse():
    for level in range(3):
        row=t.trial(level)
        quotient=tuple(map(F,row['comparison_Rayleigh_quotient']))
        assert t.g.sign(quotient)>0 and t.g.at_most(quotient,t.upper(level))
        assert t.upper(level)<5
    assert t.upper(200)==F(550,240*4**200)
    assert t.upper(1)==t.upper(0)/4


def test_affine_symmetry_break_is_detected(monkeypatch):
    original=t.g.coordinates
    def wrong(level):
        rows,den=original(level);changed=list(rows);x=list(changed[0]);x[0]+=1;changed[0]=tuple(x)
        return tuple(changed),den
    monkeypatch.setattr(t.g,'coordinates',wrong)
    with pytest.raises(ArithmeticError,match='centered'):
        t.trial(0)


@pytest.mark.parametrize('value',[True,-1,0.0,'1'])
def test_exact_level_domain(value):
    with pytest.raises(ValueError):t.upper(value)
    with pytest.raises(ValueError):t.trial(value)
