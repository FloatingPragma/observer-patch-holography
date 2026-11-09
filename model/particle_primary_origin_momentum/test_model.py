"""Independent origin, quantum-history, locality and source-authentication controls."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product
import math

import numpy as np
import pytest
import sympy as s

from . import model as m


@pytest.mark.parametrize('level',[0,1,2])
def test_exact_cube_overlap_matches_literal_sets(level):
    g=m.origin_geometry(level);lo=g['collar_cells'];hi=g['source_side']-lo
    points=set(product(range(lo,hi),repeat=3))
    for shift in ((0,0,0),(1,-1,2),(g['origin_side'],0,0)):
        shifted={tuple(x+y for x,y in zip(p,shift)) for p in points}
        displacement=tuple(x*g['source_spacing'] for x in shift)
        bound=m.folner_bound(level,displacement)
        assert bound['origin_overlap_count']==len(points&shifted)
        assert bound['symmetric_difference_count']==len(points^shifted)
        assert bound['unit_norm_expectation_difference_upper']==F(len(points^shifted),len(points))


def test_collar_is_physically_divergent_and_fraction_vanishes_without_enumeration():
    values=[m.origin_geometry(n) for n in (6,12,18,24)]
    assert all(a['physical_collar']<b['physical_collar'] for a,b in zip(values,values[1:]))
    assert all(a['discarded_fraction']>b['discarded_fraction'] for a,b in zip(values,values[1:]))
    ratios=[m.folner_bound(n,('7/3','-2','1/5'))['unit_norm_expectation_difference_upper']
            for n in (6,12,18,24)]
    assert all(x>y for x,y in zip(ratios,ratios[1:]))
    assert ratios[-1]<F(1,1000000)
    assert m.origin_geometry(0)['origin_count']==8
    assert m.origin_geometry(24)['large_origin_register_allocated'] is False


def test_original_weight_restriction_is_sharp_and_has_no_success_division():
    level=1;g=m.origin_geometry(level);p=g['retained_fraction'];theta=1-p
    w=F(9,10);floor=m.collar_weight(level,w)
    # Saturate the outside probability at one; the remaining original mass
    # has exactly this conditional-origin probability, without conditioning on success.
    assert p*floor+theta==w
    assert m.collar_weight(level,F(1,100))==0
    assert m.collar_weight(0,F(3,100))==F(3,100)
    assert m.collar_weight(13,F(3,100))>F(69,10000)
    assert m.collar_weight(16,F(3,100))>F(185,10000)
    assert m.collar_weight(24,F(3,100))<F(3,100)
    with pytest.raises(ValueError):m.collar_weight(0,'2')


def test_same_origin_homomorphism_and_independent_origin_counterexample():
    result=m.same_origin_control()
    assert result['block_product_defect']==s.zeros(8)
    assert result['homomorphic_second_moment']==F(1,2)
    assert result['arithmetic_mean_second_moment']==F(1,4)
    assert result['same_origin']=={(0,0):F(1,2),(0,1):0,(1,0):0,(1,1):F(1,2)}
    assert all(value==F(1,4) for value in result['fresh_origins'].values())
    # The fixed input is not invariant under exchanging its two spatial sites.
    assert result['same_origin'][(1,1)]!=1


def test_complete_history_retains_entangled_reference_and_classical_origin():
    result=m.complete_reference_control()
    assert result['total_probability']==1
    assert result['origin_probabilities']==[F(1,2),F(1,2)]
    assert result['reference_marginal_difference']==s.zeros(4)
    assert len(result['records'])==8
    # Independent floating-point Kraus reconstruction of the entire cq output.
    eye=np.eye(4);Ns=(np.diag([0,1,0,1]),np.diag([0,0,1,1]))
    U=np.kron(np.array([[1,1],[1,-1]])/np.sqrt(2),np.eye(2))
    probabilities=np.array([.5,.25,.125,.125]);psi=np.zeros(16)
    psi[::5]=np.sqrt(probabilities);blocks=[]
    for origin,N in enumerate(Ns):
        effects=np.zeros((4,4))
        for r,t in product(range(2),repeat=2):
            ks=(N,eye-N);K=ks[t]@U@ks[r];effects+=K.T@K
            v=np.kron(K,eye)@psi/np.sqrt(2);blocks.append(np.outer(v,v))
        np.testing.assert_allclose(effects,eye,atol=1e-14)
    unrecorded=sum(blocks).reshape(4,4,4,4)
    ref=np.trace(unrecorded,axis1=0,axis2=2)
    np.testing.assert_allclose(ref,np.diag(probabilities),atol=1e-14)
    assert sum(np.trace(b) for b in blocks)==pytest.approx(1)


def test_complete_origin_error_is_average_not_origin_count_times_error():
    # Orthogonal origin blocks make half trace distance additive with their
    # classical input probabilities, including an arbitrary quantum reference.
    bell=np.array([1,0,0,1],dtype=complex)/np.sqrt(2)
    angles=[.1,.3,.8];joint=np.zeros((12,12),dtype=complex);individual=[]
    for i,t in enumerate(angles):
        u=np.diag([np.exp(1j*t),np.exp(-1j*t)])
        shifted=np.kron(u,np.eye(2))@bell
        diff=np.outer(shifted,shifted.conj())-np.outer(bell,bell.conj())
        joint[4*i:4*i+4,4*i:4*i+4]=diff/3
        individual.append(np.sum(np.abs(np.linalg.eigvalsh(diff)))/2)
    distance=np.sum(np.abs(np.linalg.eigvalsh(joint)))/2
    assert distance==pytest.approx(sum(individual)/3)
    assert distance<=max(individual)


def test_dyadic_bulk_bound_controls_the_actual_exponential_and_zero_time():
    member=m.contract.member(6)
    value=m.bulk_bound(member)
    assert value['all_origins_bulk_admitted']
    bits=value['operator_error']['negative_exponent']
    r=value['all_origins_support_distance_lower'];a=F(member['source']['spacing'])
    exact_log=math.log(float(value['Duhamel_prefactor']))-r*math.log(2)+2*float(
        value['weighted_rate_times_spacing_upper']*value['source_time_extent']/a)
    assert bits>0 and exact_log<=-bits*math.log(2)
    zero=m.bulk_bound(member,duration='0')
    assert zero['operator_error']['exact_zero']
    assert m.bulk_bound(m.contract.member(0))['all_origins_bulk_admitted'] is False


def test_rounded_regular_translation_defect_vanishes_with_spatial_smearing():
    rows=[m.spatial.translation_group_defect(n,('1/3','2/7','-1/9'),
        ('1/5','-1/7','2/9')) for n in (2,5,8)]
    assert all(row['finite_rounded_maps_are_exact_R3_group'] is False for row in rows)
    assert all(row['regular_quotient_maps_form_exact_R3_group'] is True for row in rows)
    assert rows[-1]['vanishing_upper']<rows[0]['vanishing_upper']


def test_invariant_cyclic_vector_does_not_force_nonzero_momentum():
    # Constant observable algebra C with its only state has a one-dimensional
    # GNS space. All translations and times act as identity.
    H=np.zeros((1,1));P=[np.zeros((1,1)) for _ in range(3)];Omega=np.ones(1)
    assert np.linalg.norm(H@Omega)==0
    assert all(np.linalg.norm(p@Omega)==0 for p in P)
    assert all(np.linalg.norm(H@p-p@H)==0 for p in P)


def test_source_consumer_regenerates_protocol_and_rejects_tampering():
    member=m.contract.member(0);receipt=m.bind_member(member);m.verify_bound(member,receipt)
    result=receipt['result']
    assert result['origin_original_weight_lower']=='3/100'
    assert result['observation_protocol']['complete_one_stage_record_count']==16
    assert result['relative_observation_canonical_H_and_P_strongly_commute']
    assert result['relative_observation_cyclic_vector_H_and_P_zero']
    assert not result['fixed_origin_selected_state_translation_invariant_proved']
    assert not result['nonzero_momentum_response_proved']
    bad=deepcopy(receipt);bad['result']['nonzero_momentum_response_proved']=True
    bad['result_sha256']=m.contract.digest(bad['result'])
    with pytest.raises(ValueError):m.verify_bound(member,bad)
    changed=deepcopy(member);changed['reference']['idle_H']=1
    with pytest.raises(ValueError):m.bind_member(changed)


def test_high_member_bound_preserves_positive_original_collar_weight():
    bound=m.bind_member(m.contract.member(13))['result']
    assert F(bound['origin_original_weight_lower'])>F(69,10000)
    assert bound['origin_finite_record_weight_lower']=='0'
    assert bound['same_source_bulk_bound']['operator_error']['negative_exponent']>0
    assert bound['observation_protocol']['same_complete_failure_law_as_aggregated_failure'] is False
