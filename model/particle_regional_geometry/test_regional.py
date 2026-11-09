"""Actual source reduction, reference retention and modular-kernel controls."""
import pytest
import sympy as s

from model.particle_candidate_atlas.interacting_chiral import creator
from model.particle_regional_geometry import regional as r


@pytest.mark.parametrize('family',(0,1,2))
def test_actual_reduced_state_has_all_positive_support_and_pair_splitting(family):
    region=r.family_region(family)
    state=r.thermal.thermal_fixture()
    assert region['rho']==r.restrict_family(state['rho'],family)
    assert s.trace(region['rho'])==1
    assert region['ranks']==(1,6,2,13)
    assert sum(region['ranks'])==22
    assert all(x>0 for x in region['probabilities_per_state'])
    assert sum(n*p for n,p in zip(region['ranks'],region['probabilities_per_state'],strict=True))==1
    assert r.clean(region['rho']*(s.eye(64)-region['support']))==s.zeros(64)
    assert region['probabilities_per_state'][2]>region['probabilities_per_state'][3]
    assert region['global_number_commutator_vacuum_norm_squared']==16*(family+1)**2


def test_complete_signed_split_and_number_preserving_algebra_contract():
    state=r.thermal.thermal_fixture()
    labels=r.family_split(1)
    index=state['source']['basis'].index((1,(1<<0)|(1<<6)))
    assert labels[index][2]==-1
    region=r.family_region(1)
    assert r.inject_family(region['support'],1)==s.eye(307)
    a,b=region['probe'],region['number']+region['pair_projector']
    assert r.clean(r.inject_family(a*b,1)-r.inject_family(a,1)*r.inject_family(b,1))==s.zeros(307)
    assert s.simplify(s.trace(region['rho']*a)-s.trace(state['rho']*r.inject_family(a,1)))==0
    # Arbitrary odd CAR operators leave Q*=0 and do not give a homomorphism
    # on this fixed sector. Their complete ambient representation is larger.
    raising=creator(6,0)
    assert r.inject_family(raising,1)==s.zeros(307)
    assert r.inject_family(raising.H*raising,1)!=s.zeros(307)


def test_partial_trace_retains_an_entangled_untouched_reference():
    source=r.thermal.thermal_fixture()['source']
    i=source['basis'].index((1,(1<<0)|(1<<5)))
    j=source['basis'].index((1,(1<<1)|(1<<4)))
    vector=s.ImmutableSparseMatrix(614,1,{(2*i,0):1/s.sqrt(2),(2*j+1,0):s.I/s.sqrt(2)})
    reduced=r.restrict_family(vector*vector.H,0,2)
    local=s.ImmutableSparseMatrix(128,1,{(2*((1<<0)|(1<<5)),0):1/s.sqrt(2),
                                          (2*((1<<1)|(1<<4))+1,0):s.I/s.sqrt(2)})
    assert reduced==local*local.H
    assert s.trace(reduced)==1
    reference=s.Matrix(2,2,lambda a,b:sum(reduced[2*k+a,2*k+b] for k in range(64)))
    assert reference==s.eye(2)/2


def test_modular_pair_term_cannot_be_any_number_polynomial():
    region=r.family_region(0)
    a=region['probe']; n=region['number']; k=region['modular_H']
    assert r.clean(n*a-a*n)==s.zeros(64)
    assert r.clean(n*n*a-a*n*n)==s.zeros(64)
    assert r.clean(k*a-a*k)!=s.zeros(64)
    b,d=region['pair_vectors'][0],region['dark_vector']
    contrast=s.simplify((b.H*k*b)[0]-(d.H*k*d)[0])
    assert s.simplify(contrast-region['pair_modular_coefficient'])==0
    gauge=r.gauge_probe_certificate()
    assert gauge['source_gauge_generators']==12
    assert gauge['operator_commutators_checked']==36
    assert gauge['regional_pair_probe_is_gauge_invariant']


def test_cross_center_coherence_excludes_a_state_preserving_expectation():
    region=r.family_region(0)
    p0,_,p2=map(r.inject_family,region['occupation_projectors'])
    x=region['cross_center_coherence']
    assert r.clean(p0*x*p2+p2*x*p0-x)==s.zeros(307)
    assert region['cross_center_mean']==-s.Rational(131068,1077135)
    assert r.clean(p0*r.inject_family(region['probe'])*p2)==s.zeros(307)
    assert r.trace_expectation(x)==s.zeros(307)
    state=r.thermal.thermal_fixture()['rho']
    assert r.trace_expectation(state)!=state


def test_tracial_expectation_is_unital_retracts_and_preserves_trace():
    region=r.family_region(0)
    a=r.inject_family(region['probe']+region['number'])
    assert r.trace_expectation(s.eye(307))==s.eye(307)
    assert r.trace_expectation(a)==a
    state=r.thermal.thermal_fixture()['rho']
    averaged=r.trace_expectation(state)
    assert s.trace(averaged)==s.trace(state)==1
    assert r.trace_expectation(averaged)==averaged


def test_actual_entropy_tangent_lifts_to_a_faithful_same_source_state():
    value=r.entropy_response(0)
    assert value['declared_variation_size']==s.Rational(1,71809)
    assert value['perturbed_global_density_lower_bound']==s.Rational(3,4*71809)
    assert value['finite_relative_entropy_remainder_upper']==s.Rational(15,26*71809*511)
    delta=value['regional_tangent']; lifted=value['global_tangent']
    assert r.restrict_family(lifted)==delta
    for p in r.family_region(0)['occupation_projectors']:
        assert s.trace(delta*p)==0
    # Exact polynomial identity supplies the claimed global eigenvalue bound.
    assert r.clean(lifted*(lifted-s.eye(307)/4)*(lifted+s.eye(307)/26))==s.zeros(307)
    assert s.simplify(value['entropy_first_derivative']-s.log(s.Rational(13440,41873)))==0


def test_scalar_modular_flow_is_trivial_on_the_actual_retained_algebra():
    region=r.scalar_region()
    assert region['per_state_probabilities']==(s.Rational(262151,1077135),
                                              s.Rational(407492,1077135),s.Rational(407492,1077135))
    assert region['density_is_central_on_that_algebra']
    assert region['formal_number_modular_coefficient']<0
    assert region['global_H_NH_commutator_vacuum_norm_squared']==56
    for a in (s.diag(1,0,0),s.Matrix([[0,0,0],[0,0,1],[0,1,0]])):
        assert region['reduced_density']*a==a*region['reduced_density']


def test_domains_reject_wrong_source_carrier_or_family():
    with pytest.raises(ValueError,match='three declared'):
        r.family_split(3)
    with pytest.raises(ValueError,match='complete source'):
        r.restrict_family(s.eye(3)/3)
    with pytest.raises(ValueError,match='reference'):
        r.restrict_family(s.eye(307)/307,reference_dimension=0)
    with pytest.raises(ValueError,match='six-mode'):
        r.inject_family(s.eye(3))
