from fractions import Fraction as F
import sympy as s
import pytest
from .source import profile,energy,records,exact_gauss

def test_every_native_gauss_generator_on_both_charge_strings():
    assert exact_gauss('uc',3)==48
    assert exact_gauss('dc',2)==36

def test_boundary_extension_does_not_change_local_profile_or_energy():
    for q in(2,-1):
        a,b=profile(4,2,q),profile(9,2,q)
        assert a['source_charge_exact']==b['source_charge_exact'][:4]
        assert a['integer_flux_on_all_edges']==b['integer_flux_on_all_edges'][:3]
        assert a['complete_idle_energy_exact']==b['complete_idle_energy_exact']
        assert energy(q,3)>energy(q,2)>2
    with pytest.raises(ValueError):profile(3,3,2)

def test_complete_quadrature_records_and_work_for_all_preparations():
    for weight in(F(0),F(1,2),F(1)):
        row=records(2,2,s.pi/2,weight)
        for time in row['all_time_records']:
            ps=list(map(s.sympify,time['X_plus_X_minus_Y_plus_Y_minus_exact']))
            assert sum(ps)==1 and min(ps)>=0
        assert s.sympify(row['preparation_mean_work_exact'])+s.sympify(row['pointer_nonselective_mean_work_exact'])==s.Rational(energy(2,2))/2
    assert records(2,1,0,F(1,2))['all_time_records'][0]['X_plus_X_minus_Y_plus_Y_minus_exact']==['1/2','0','1/4','1/4']

def test_gauss_failure_is_a_preceding_test_not_an_extra_identity_effect():
    G=s.diag(1,1,0);I=s.eye(3)
    X=s.Matrix([[0,1,0],[1,0,0],[0,0,0]])
    Y=s.Matrix([[0,-s.I,0],[s.I,0,0],[0,0,0]])
    effects=[G*(I+sign*A)*G/4 for A in(X,Y)for sign in(1,-1)]+[I-G]
    assert sum(effects,s.zeros(3))==I
    assert sum([(I+sign*A)/4 for A in(X,Y)for sign in(1,-1)],s.zeros(3))+I-G!=I
