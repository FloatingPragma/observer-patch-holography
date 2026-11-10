"""Exact independent projector and current coefficient controls."""
from fractions import Fraction as F
import sympy as s
from . import exact_current as current_audit

def test_independent_exact_current_radicals_and_positive_isolation():
    value=current_audit.calculate()
    expected=s.Rational(965,1456)+5*s.sqrt(6)/36-s.sqrt(14)/48+s.sqrt(21)/336+s.sqrt(91)/546+7*s.sqrt(546)/1872
    assert s.simplify(value['exact_sigma_squared_over_m_squared']-expected)==0
    low,high=map(F,value['outward_interval'])
    assert F('1.0435157377429135586350084518253')<low<=high<F('1.0435157377429135586350084518254')
    assert low>F(value['strict_sigma_over_m_lower'])**2
    assert value['complete_CRT_projector_ranks']==[8,16,8]
    assert value['rank']==16 and value['all_twelve_filled_vacuum_current_means']=='0'
    assert not value['binomial_sign_approximation_or_numerical_eigensolver_used']

def test_rational_interpolant_on_each_actual_positive_pair():
    for r in (3,7,13):
        z=s.Symbol('z');f=(3*r+s.Rational(1,4)-z)/(2*s.sqrt(r)*(r-s.Rational(1,4)))
        for sign in (-1,1):
            energy=s.sqrt(r)+s.Rational(sign,2)
            assert energy>0 and s.simplify(f.subs(z,energy**2)*energy)==1
