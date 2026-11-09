"""Actual source coefficients and independent finite information identities."""
from copy import deepcopy
from fractions import Fraction as F

import numpy as np
import pytest

from . import model as m


def test_native_schur_information_includes_global_center_coherence():
    r=m.native_control()
    assert abs(r['identity_residual'])<1e-13
    assert abs(r['regional_algebra_mutual_information']-np.log(5))<1e-13
    assert abs(r['representation_entropy']-np.log(108)/5)<1e-13
    assert r['center_coherence']>0.3
    assert not r['control_state_substituted_for_primary_source']


def test_arbitrary_multiplicity_correlations_and_cross_sector_coherence():
    rng=np.random.default_rng(3901)
    dims=[1,2];mult=[(2,2),(2,2)]
    X=rng.normal(size=(8,8))+1j*rng.normal(size=(8,8))
    rho=X@X.conj().T+np.eye(8);rho/=np.trace(rho)
    r=m.information_control(rho,dims,mult)
    assert abs(r['identity_residual'])<1e-12
    assert r['conditional_multiplicity_mutual_information']>0
    assert r['center_coherence']>0
    assert r['ambient_source_mutual_information']>r['regional_algebra_mutual_information']


def test_same_edge_term_does_not_fix_full_regional_information():
    rho=np.eye(8)/8
    bell=np.array([1,0,0,1])/np.sqrt(2)
    correlated=np.kron(np.eye(2)/2,np.outer(bell,bell))
    a=m.information_control(rho,[1,2],[(2,2),(2,2)])
    b=m.information_control(correlated,[1,2],[(2,2),(2,2)])
    assert abs(a['representation_entropy']-b['representation_entropy'])<1e-13
    assert abs(b['ambient_source_mutual_information']-a['ambient_source_mutual_information']-2*np.log(2))<1e-12


def test_source_matrix_element_uses_every_native_character_and_actual_normalization():
    blocks=m.cut.stress.channels.native_number_control()['blocks']
    assert sum(r['internal_dimension'] for r in blocks)==15
    charges={(sign*r['Q6'],-sign*r['Q6']) for r in blocks for sign in (-1,1)}
    assert len(charges)==10 and (0,0) not in charges
    assert F(1,7)/(2*15)==F(1,210)
    # Haar-orthogonal normalized characters each have amplitude1/30.
    A_blank=np.ones(10)/30
    assert abs(np.dot(A_blank,A_blank)-1/90)<1e-15


def test_strict_polynomial_pinching_work_and_coherence_bounds():
    E=np.array([-2.,-.3,.4,3.]);H=np.diag(E)
    # Q=4-H>0, with K'=2/(4-E) between1/3 and2.
    p=(4-E)**2;p/=sum(p);rho=np.diag(p);K=np.diag(-np.log(p))
    rng=np.random.default_rng(71);U,_=np.linalg.qr(rng.normal(size=(4,4)))
    P=[np.outer(U[:,i],U[:,i]) for i in range(4)]
    out=sum(Q@rho@Q for Q in P);X=out-rho
    W=np.trace(H@X);C=m.entropy(out)-m.entropy(rho)
    values,vectors=np.linalg.eigh(out);log_out=(vectors*np.log(values))@vectors.T
    reverse=np.trace(out@(log_out-K*(-1)))
    assert abs(C+reverse-np.trace(X@K))<1e-12
    assert W>0 and C>0 and W/3<=C+reverse<=2*W
    hs_sum=sum(np.linalg.norm(H@Q-Q@H,'fro')**2 for Q in P)
    assert W>=min(p)/6*hs_sum-1e-14
    assert C>=np.linalg.norm(X,'fro')**2/2-1e-14


def test_center_entropy_floor_needs_all_eleven_outcomes():
    c=1/32
    extremal=np.array([1-10*c]+[c]*10)
    H=-sum(extremal*np.log(extremal))
    assert H>=-10*c*np.log(c)
    rng=np.random.default_rng(331)
    for _ in range(30):
        weights=rng.dirichlet(np.ones(11))*(1-11*c)+c
        assert -sum(weights*np.log(weights))>=H-1e-13


@pytest.fixture(scope='module')
def actual():return m.bind_member(m.contract.member(0))


def test_actual_original_state_has_positive_factored_record_work_and_coherence(actual):
    r=actual['result']
    assert r['blank_to_each_crossing_native_character_H_matrix_element']=='-1/210'
    assert F(r['Hamiltonian_zero_center_commutator_HS_squared_lower'])==F(1,2205)
    c=r['source_density_floor'];assert c['negative_exponent']==582407133
    for key,power in [('center_coherence_lower_factored',2),('classical_center_record_mean_source_work_lower_factored',1)]:
        assert F(r[key]['rational_coefficient'])>0
        assert r[key]['negative_exponent']==power*c['negative_exponent']
    assert r['complete_classical_center_record_disturbs_global_source']
    assert not r['positive_capacity_uniform_information_bound_proved']
    assert not r['geometric_area_law_or_Einstein_equation_proved']


def test_binding_rejects_entropy_domain_or_physical_promotion(actual):
    m.verify_bound(m.contract.member(0),actual)
    altered=deepcopy(actual);altered['result']['geometric_area_law_or_Einstein_equation_proved']=True
    with pytest.raises(ValueError):m.verify_bound(m.contract.member(0),altered)
    with pytest.raises(ValueError):m.singlet_embedding([0],[(1,1)])
