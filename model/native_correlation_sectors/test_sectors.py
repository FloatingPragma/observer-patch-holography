"""Independent native enumeration, finite linear algebra and raw-record controls."""
from fractions import Fraction as F
from itertools import combinations, product
from math import comb

import pytest
import sympy as sp

from . import sectors as s, records as r


def incidence(ports,degree):
    lo=tuple(combinations(range(ports),degree));hi=tuple(combinations(range(ports),degree+1))
    return sp.Matrix([[int(set(a)<=set(b)) for a in lo] for b in hi])


def native_counts(family,degree):
    states=tuple(combinations(range(len(family.weights)),degree));index={x:i for i,x in enumerate(states)}
    out=sp.zeros(len(states))
    for column,state in enumerate(states):
        for edge,coin in product(family.seams,(0,1)):
            image,_,_=r.native_step(state,edge,coin)
            out[index[image],column]+=1
    return states,out


@pytest.mark.parametrize('ports',(4,5,6))
def test_exhaustive_orthogonal_decomposition_and_primitive_norms(ports):
    u=[incidence(ports,k) for k in range(ports)]
    for degree in range(ports+1):
        du=u[degree].T*u[degree] if degree<ports else sp.zeros(comb(ports,degree))
        ud=u[degree-1]*u[degree-1].T if degree else sp.zeros(1)
        assert du-ud==(ports-2*degree)*sp.eye(comb(ports,degree))
    for q in range(ports+1):
        columns=[];blocks=[]
        for degree in range(min(q,ports-q)+1):
            w=(sp.eye(1) if degree==0 else sp.Matrix.hstack(*u[degree-1].T.nullspace()))
            assert w.cols==s.dimensions(ports,q)[degree]
            included=w
            for step in range(degree,q):included=u[step]*included/(step-degree+1)
            assert included.T*included==comb(ports-2*degree,q-degree)*(w.T*w)
            for prior in blocks:assert prior.T*included==sp.zeros(prior.cols,included.cols)
            blocks.append(included);columns.extend(included.columnspace())
        assert len(columns)==comb(ports,q)


@pytest.mark.parametrize('ports',(4,5,6))
def test_native_source_intertwines_all_degrees_and_preserves_primitives(ports):
    seams=tuple((i,i+1) for i in range(ports-1))
    family=s.source.Family((1,)*ports,seams,(1,)*len(seams))
    matrices=[]
    for q in range(ports+1):
        states,p=native_counts(family,q);matrices.append(p)
        for column,state in enumerate(states):
            actual=s.counts_step(family,{state:1},q)
            assert [actual.get(x,0) for x in states]==list(p[:,column])
    for q in range(ports):
        u=incidence(ports,q)
        assert matrices[q+1]*u==u*matrices[q]


@pytest.mark.parametrize('ports',(4,6,12))
def test_four_port_contrast_has_no_one_body_component(ports):
    h=s.contrast((0,1),(2,3),ports)
    assert not s.lower(h,ports,2)
    for q in range(ports+1):
        values=s.inclusion(h,ports,q,2) if q>=2 else {}
        if values:
            assert sum(values.values())==0
            for i in range(ports):assert sum(v for x,v in values.items() if i in x)==0
        actual=F(sum(v*v for v in values.values()),comb(ports,q))
        expected=4*s.norm_factor(ports,2,q) if 2<=q<=ports-2 else 0
        assert actual==expected


def test_complete_raw_record_algebra_includes_all_zero_slots():
    result=r.pair_certificate()
    assert result['two_E_raw_histories']==625
    assert result['live_two_E_histories']==36
    assert result['zero_two_E_histories']==589
    # Each map identity holds on every input/output matrix unit, not only
    # on a selected stationary sample.
    for a,b in product(r.raw_maps((0,1)),r.raw_maps((2,3))):
        mapping={i:b[2][j] for i,j in a[2].items() if j in b[2]}
        for i,j in product(mapping,repeat=2):
            fi=r.field({k for k in range(4) if i >> k & 1},(0,1),(2,3))
            fj=r.field({k for k in range(4) if mapping[j] >> k & 1},(0,1),(2,3))
            assert fi==a[3]*b[3] and fj==a[4]*b[4]


def test_naive_incoming_past_score_has_real_measurement_backaction():
    family=s.source.Family((1,)*4,((0,1),(1,2),(2,3)),(1,)*3)
    assert r.direct_probe(family,2,0,(0,1),(2,3))==F(2,3)
    assert r.direct_probe(family,2,0,(0,1),(2,3),reflected=False)==0


@pytest.mark.parametrize('bad',(True,1.0,-1))
def test_dimensions_reject_ambiguous_port_types(bad):
    with pytest.raises(ValueError):s.dimensions(bad,0)


@pytest.mark.parametrize('bad',(True,1.0,-1))
def test_degree_rejects_ambiguous_types(bad):
    with pytest.raises(ValueError):s.norm_factor(6,2,bad)


@pytest.mark.parametrize('edges',(((0,1),(1,2)),((1,0),(2,3)),((0,True),(2,3))))
def test_probe_requires_disjoint_canonical_seams(edges):
    with pytest.raises(ValueError):s.contrast(*edges,4)


def test_sparse_controls_fail_closed_on_degree_or_coefficient_changes():
    with pytest.raises(ValueError):s.lower({(0,1):1,(2,):1},4)
    with pytest.raises(ValueError):s.lower({(0,1):True},4)
    with pytest.raises(ValueError):s.lower({(0,1):1.0},4)
    with pytest.raises(ValueError):s.inclusion({(0,1):1},32,16,2)
    with pytest.raises(ValueError):s.norm_factor(4,2,1)
    with pytest.raises(ValueError):s.autocorrelation(0,True)
    with pytest.raises(ValueError):s.autocorrelation(0,13)


def test_biased_source_cannot_enter_the_fair_intertwiner():
    with pytest.raises(ValueError):s.counts_step(s.source.local_family(),{(0,1):1},2)


def test_full_probe_is_admitted_but_expanded_host_stop_is_not_a_success():
    program=r.probe_program(0,2,3)
    assert program['declared_resources']['attempts_upper']==7
    assert len(program['commands'])==5
    initial={(3,0,3,0):F(1)}
    stopped=s.source.interpreter.execute(program,initial,max_branches=1)
    assert stopped.status=='RESOURCE_STOP' and stopped.attempted_steps==0
    assert stopped.branches[0].density==initial
    assert stopped.branches[0].record.kind=='RESOURCE_STOP'
    with pytest.raises(ValueError):r.probe_program(0,2,s.source.interpreter.MAX_ATTEMPTS)


def test_raw_map_evaluator_rejects_malformed_or_oversized_inputs():
    for ports,masks in ((True,None),(100,None),(4,(True,)),(4,(16,)),(4,(1,1))):
        with pytest.raises(ValueError):r.raw_maps((0,1),ports,masks)
