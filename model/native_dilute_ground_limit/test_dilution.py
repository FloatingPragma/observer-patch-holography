from fractions import Fraction as F
from math import comb
import numpy as np
import pytest
from scipy.linalg import expm
from . import dilution as d


@pytest.mark.parametrize('V,q,s',((12,2,2),(240,2,12),(240,0,12),(240,240,12),(20,19,3),(10,2,0)))
def test_exact_hypergeometric_vacuum_weight(V,q,s):
    expected=1-F(comb(V-s,q),comb(V,q))if q<=V-s else F(1)
    assert d.population(V,q,s)==expected


def test_huge_capacity_uses_only_fixed_local_factors():
    V=1<<400
    assert 0<d.population(V,2,12)<=F(24,V)
    result=d.history_bound(V,2,(2,2,4))
    assert result['complete_terminal_halftrace_upper']<F(1,2**120)


@pytest.mark.parametrize('args',((True,2,1),(12,True,1),(12,2,False),(0,0,0),(12,13,2),(12,2,13),(100,2,65)))
def test_malformed_or_unbounded_domains_rejected(args):
    with pytest.raises(ValueError):d.population(*args)


def test_complete_local_instrument_retains_zero_slots_and_rejects_charge_supply():
    I=((F(1),F()),(F(),F(1)));zero=((F(),F()),)*2
    p0=((F(1),F()),(F(),F()));p1=((F(),F()),(F(),F(1)))
    labels=('empty','occupied','zero')
    branches=(('empty',F(1),p0),('occupied',F(1),p1),('zero',F(1),zero))
    assert d.complete_number_preserving_instrument(labels,branches,1)==(F(1),F(),F())
    with pytest.raises(ValueError):d.complete_number_preserving_instrument(labels,branches[:-1],1)
    with pytest.raises(ValueError):d.complete_number_preserving_instrument(labels,branches[:-1]+(('zero',F(1),I),),1)
    X=((F(),F(1)),(F(1),F()))
    with pytest.raises(ValueError):d.complete_number_preserving_instrument(('imported_charge',),(('imported_charge',F(1),X),),1)


def test_source_H66_complete_endpoint_witness():
    row=d.native_endpoint_certificate()
    assert row['unequal_probability']==F(10,33)
    assert row['raw_slots']==25
    assert row['terminal_halftrace_sqrt_radicand']==F(85,363)
    assert row['public_native_resources']['attempts_executed']==1


def test_noncommuting_complete_controlled_history_has_joint_state_bound():
    # A one-particle cycle is an independent finite control, not the selected
    # federation. Retain quantum source and every instrument record.
    V=128;omega=np.ones(V)/np.sqrt(V)
    G=2*np.eye(V)-np.roll(np.eye(V),1,axis=0)-np.roll(np.eye(V),-1,axis=0)
    U=expm(-.7j*G)
    empty=np.eye(V);empty[0,0]=0;occupied=np.zeros((V,V));occupied[0,0]=1
    rotation=np.eye(V);rotation[:2,:2]=[[.6,.8],[-.8,.6]]
    assert np.linalg.norm(G@rotation-rotation@G)>1
    real_branches=[rotation@U@K@omega for K in (empty,occupied,np.zeros((V,V)))]
    ideal_branches=[omega,np.zeros(V),np.zeros(V)]
    distance=0
    for actual,wanted in zip(real_branches,ideal_branches,strict=True):
        diff=np.outer(actual,actual.conj())-np.outer(wanted,wanted.conj())
        distance+=np.linalg.eigvalsh(diff).__abs__().sum()/2
    bound=d.history_bound(V,1,(1,2))['complete_terminal_halftrace_upper']
    assert .1<distance<float(bound)<.5
    assert abs(sum(np.vdot(x,x).real for x in real_branches)-1)<1e-12


def test_exact_local_number_variance_and_connected_spectral_mass():
    for V,q in ((12,2),(240,2),(960,2),(240,120)):
        variance=F(q*(V-q),V*V)
        assert variance<=d.observable_bounds(V,q,1,1)['norm_one_variance_upper']
    assert d.observable_bounds(960,2,1,1)['norm_one_variance_upper']<d.observable_bounds(240,2,1,1)['norm_one_variance_upper']


def test_surviving_density_norm_and_support_routes_are_outside_dilution():
    for V in (16,32,64,128):
        # Finite density retains a norm-one local number variance.
        q=V//2
        assert F(q*(V-q),V*V)==F(1,4)
        assert d.population(V,q,1)==F(1,2)
        # Norm-sqrt(V) local fields evade the bounded-norm premise at q=2.
        assert V*F(2*(V-2),V*V)>1
        # Regions containing half the source evade bounded-support dilution.
        assert d.population(V,2,V//2)>F(3,4)


def test_source_density_ground_residue_and_eigenvalue_identity():
    result=d.density_closure_certificate()
    assert result['exact_density_intertwining_entries']==792
    assert result['exact_covariance_entries']==144
    assert result['mean_zero_unit_l2_density_residue']==F(5,33)


def test_exact_visibility_threshold_and_zero_sector():
    result=d.visibility_requirement(240,2,F(1),F(1,10))
    assert result['minimum_scaled_occupation']==F(1,40)
    assert result['minimum_integer_support']==3
    assert F(2*2,240)<F(1,40)<=F(3*2,240)
    assert d.visibility_requirement(240,2,F(1),F(1))['minimum_integer_support']==30
    empty=d.visibility_requirement(240,0,F(1),F(1,10))
    assert empty['positive_weight_impossible_in_empty_sector']and empty['minimum_integer_support']is None
    assert d.visibility_requirement(240,0,F(1),F())['minimum_integer_support']==0
    assert not result['sufficient_particle_identification_claimed']


@pytest.mark.parametrize('norm,weight',((True,F(1,2)),(F(),F()),(F(1),True),(F(1),F(2))))
def test_visibility_rejects_malformed_or_impossible_target(norm,weight):
    with pytest.raises(ValueError):d.visibility_requirement(240,2,norm,weight)
