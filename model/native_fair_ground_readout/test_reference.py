from fractions import Fraction as F
from itertools import combinations
import math
import numpy as np
import pytest
from . import reference as r


@pytest.mark.parametrize('v,q,k',[(4,2,2),(7,3,4),(12,2,5),(20,9,7),(5,0,3),(5,5,3)])
def test_regional_ground_matches_exact_partial_trace(v,q,k):
    states=list(combinations(range(v),q));counts={}
    for s in states:
        local=sum(i<k for i in s);counts[local]=counts.get(local,0)+1
    assert r.dicke_weights(v,q,k)==tuple(F(counts.get(i,0),len(states)) for i in range(k+1))
    assert sum(r.dicke_weights(v,q,k))==1


@pytest.mark.parametrize('v,q,k,p',[(9,4,3,F(1,2)),(20,5,4,F(1,4)),(100,50,7,F(1,2)),(12,0,5,F()),(12,12,5,F(1))])
def test_capacity_bound_controls_complete_regional_trace_distance(v,q,k,p):
    finite=r.dicke_weights(v,q,k);limit=r.limiting_weights(k,p)
    distance=sum(abs(a-b) for a,b in zip(finite,limit,strict=True))/2
    assert distance<=r.regional_error(v,q,k,p)


def test_native_thermal_witness_is_exact_and_keeps_thermal_qft_possible():
    x=r.thermal_witness()
    assert x['weighted_transfer_commutator_norm_squared']=='37/86400'
    assert x['negative_energy_spectral_weight_lower']=='37/24000'
    assert x['same_state_positive_energy_vacuum_excluded'] is True
    assert x['thermal_QFT_excluded'] is False
    assert x['positive_underlying_H_contradicted'] is False


def test_reflected_spectral_weights_exactly_retain_both_signs():
    values=(F(1),F(3,4),F(1,2));z=sum(values)
    a=((0,1,2),(1,0,3),(2,3,0))
    positive=negative=F()
    for i in range(3):
        for j in range(3):
            weight=values[j]*a[i][j]**2/z
            if values[i]<values[j]:positive+=weight
            elif values[i]>values[j]:negative+=weight
    assert F(1,2)*positive<=negative<=2*positive
    assert positive>0 and negative>0


def test_thermal_reference_cannot_be_substituted_for_ground():
    values=(F(1),F(3,4),F(1,2));z=sum(values)
    assert 1-values[0]/z==F(5,9)
    assert r.cooling_error(3,F(1),2,0)==F(2,3)
    assert r.cooling_error(3,F(1),2,10)<r.cooling_error(3,F(1),2,2)


def test_regional_state_at_huge_capacity_uses_only_local_falling_factorials():
    v=10**122;q=v//3;k=12
    weights=r.dicke_weights(v,q,k)
    assert sum(weights)==1 and min(weights)>0
    assert len(weights)==13
    difference=sum(abs(a-b) for a,b in zip(weights,r.limiting_weights(k,F(1,3)),strict=True))/2
    assert difference<=r.regional_error(v,q,k,F(1,3))


def test_local_annihilator_survives_after_an_accessible_intervention():
    assert r.local_ground_witness(12,2)==F(5,66)
    v=10**122
    assert abs(r.local_ground_witness(v,v//3)-F(1,9))<F(1,10**120)
    assert r.local_ground_witness(12,0)==0


def test_global_dimension_free_cooling_length_is_sufficient():
    for v in range(2,9):
        length=r.cooling_length(v,F(1),v,F(1,1000))
        assert r.cooling_error(2**v,F(1),v,length)<F(1,1000)


def test_multitime_moment_budget_keeps_cross_terms():
    b=r.history_budget((1,2),(F(1,10),F(-1,5)),30,F(2,3),5)
    c1,d1=r.local_moments(1,5);c2,d2=r.local_moments(2,5)
    assert F(b['prefixes'][1]['second_energy_moment_norm_upper'])==d1+d2+2*c1*c2
    assert b['history_vector_norm_error_upper']==(F(1,10)*d1+F(1,5)*(d1+d2+2*c1*c2))/80


@pytest.mark.parametrize('v',[6,8,12])
def test_actual_ground_unitary_comparison_on_cycle(v):
    g=np.zeros((v,v))
    for i in range(v):
        j=(i+1)%v;g[i,i]+=1;g[j,j]+=1;g[i,j]-=1;g[j,i]-=1
    vals,u=np.linalg.eigh(g);vals=np.maximum(vals,0)
    alpha=1-2/v;h=-2*v*np.log(1-vals/(2*v))
    omega=np.ones(v)/np.sqrt(v);a=np.diag([1]+[0]*(v-1))
    vector=a@omega;t=.17
    actual=np.linalg.norm(u@((np.exp(-1j*t*h)-np.exp(-1j*t*vals))*(u.T@vector)))
    bound=abs(t)*np.linalg.norm(g@g@vector)/(4*v*alpha)
    assert 0<actual<=bound+1e-14
    comm=g@(g@a-a@g)-(g@a-a@g)@g
    assert np.linalg.norm(g@g@vector-comm@omega)<1e-12


@pytest.mark.parametrize('call',[lambda:r.dicke_weights(True,0,0),lambda:r.dicke_weights(2,3,1),
    lambda:r.limiting_weights(3,.5),lambda:r.history_budget((1,),(),30),
    lambda:r.history_budget((1,),(F(1),),0),lambda:r.cooling_error(2,F(),1,1)])
def test_premises_reject_invalid_or_inexact_inputs(call):
    with pytest.raises(ValueError):call()
