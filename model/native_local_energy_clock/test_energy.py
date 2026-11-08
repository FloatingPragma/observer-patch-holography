from fractions import Fraction as F
from itertools import combinations, product

import numpy as np
import pytest
from scipy.linalg import eigh

from . import energy as e
from model.observer.quantum_mixing import seams
from model.federated_observer import source


@pytest.mark.parametrize('q',range(5))
def test_exact_complete_sector_moments_and_currents(q):
    edges=tuple(combinations(range(4),2))
    row=e.exact_energy_certificate(4,edges,q)
    assert row['all_seam_continuity_identities_exact']
    assert row['repair_exchange_energy_change_on_P_reference'] >= 0
    assert row['repair_energy_conservation'] == (q in (0,4))


@pytest.mark.parametrize('level',(0,1))
def test_source_federation_counts(level):
    family=source.federation(level)
    m,delta,a=e.graph_data(family.ports,tuple(sorted(family.seams)))
    assert (family.ports,m,a)==(240*4**level,630*4**level,5400*4**level)
    assert delta==6


def test_growing_density_bound_is_quantitative_without_allocating_global_state():
    values=[]
    for level in (6,12,24,100):
        size=4**level
        row=e.local_clock_bound(240*size,630*size,5400*size,80*size,6,1,F(1))
        assert row['canonical_variance']<=F(38*6+6,9)*630*size
        values.append(row['density_clock_dynamics_upper'])
    assert all(a>b for a,b in zip(values,values[1:]))
    assert values[-1] < F(1,10**24)


def test_h66_exact_same_reference_repair_heating():
    row=e.exact_energy_certificate(12,tuple(sorted(seams())),2)
    assert row['dimension']==66
    assert row['trace_P_unnormalized']==56
    assert row['repair_exchange_energy_change_on_P_reference']==F(37,840)
    assert not row['repair_energy_conservation']
    from model.native_energy_measurement import spectral
    states,_,g=e.matrices(12,tuple(sorted(seams())),2)
    index={state:i for i,state in enumerate(states)}
    order=[index[tuple(state)] for state in spectral.synthesis.tomo.compiled.quantum.STATES]
    assert np.array_equal(60*np.eye(66,dtype=np.int64)-g[np.ix_(order,order)],
                          np.array(spectral.source_counts(),dtype=np.int64))


def test_log_energy_heating_and_complete_detector_error_budget():
    states,_,g=e.matrices(12,tuple(sorted(seams())),2)
    phi_g=e.source_numerator(g,states,tuple(sorted(seams())))/120
    lam,v=eigh(g.astype(float));h=(v*(-60*np.log1p(-lam/60)))@v.T
    change=float(np.trace(h@(g-phi_g))/(60*56))
    assert change>float(F(37,840))
    from model.native_energy_measurement import spectral
    intervals=spectral.energy_labels()
    assert all(0<=lo<=hi<=30 and hi-lo<=F(1,10**30) for lo,hi in intervals)
    budget=e.h66_detector_heating_bound()
    assert budget['complete_measured_score_increase_lower']>F(34,1000)
    assert e.h66_detector_heating_bound(preparation_error=F(1,100))['complete_measured_score_increase_lower']<0


def test_local_frequency_lemma_and_dynamics_on_actual_h66():
    edges=tuple(sorted(seams()))
    states,_,g=e.matrices(12,edges,2); m=len(edges)
    lam,v=eigh(g.astype(float))
    p=1-lam/(2*m); rho=p/p.sum()
    assert p.min()>=2/3-1e-14
    h=-2*m*np.log(p)
    mu=float(np.dot(rho,lam)); var=float(np.dot(rho,(lam-mu)**2))
    a=1/(1-mu/(2*m))
    x=np.diag([state[0] for state in states]).astype(complex)
    b=g@x-x@g; b2=g@b-b@g
    bound=(np.linalg.norm(b,2)*np.sqrt(var)+np.linalg.norm(b2,2)/2)/(2*m*(2/3)**2)
    spectral=v.T@x@v
    d=(h[:,None]-h[None,:]-a*(lam[:,None]-lam[None,:]))*spectral
    weighted=lambda z:np.sqrt(np.sum(abs(z)**2*rho[None,:]))
    assert weighted(d)<=bound+1e-13
    for t in (-2.,-.1,0.,.5,3.):
        actual=(np.exp(1j*t*(h[:,None]-h[None,:]))-
                np.exp(1j*t*a*(lam[:,None]-lam[None,:])))*spectral
        assert weighted(actual)<=abs(t)*bound+1e-13


def test_fair_source_complete_channel_self_adjoint():
    edges=tuple(combinations(range(4),2)); states,_,_=e.matrices(4,edges,2)
    d=len(states); columns=[]
    for i in range(d):
        for j in range(d):
            unit=np.zeros((d,d),dtype=np.int64);unit[i,j]=1
            columns.append(e.source_numerator(unit,states,edges).reshape(-1))
    superoperator=np.column_stack(columns)
    assert np.array_equal(superoperator,superoperator.T)


@pytest.mark.parametrize('bad',(True,1.2,'3'))
def test_exact_input_types_refused(bad):
    with pytest.raises(ValueError): e.moments(12,30,240,bad)


def test_exact_matrix_allocation_ceiling_is_preserved():
    with pytest.raises(ValueError):e.matrices(25,((0,1),),0)
    with pytest.raises(ValueError):e.matrices(12,((0,1),),6)
    with pytest.raises(ValueError):e.matrices(4,((0,1),),5)


def test_exact_square_root_enclosures():
    for x in (F(0),F(1),F(2),F(123456789,987654321),F(2)**400):
        y=e.sqrt_upper(x)
        assert y*y>=x
        assert (y-F(1,2**40))**2<x or not x


def test_binary_triangle_floor_and_dimension_three_negative_control():
    def triangle(dimension):
        states=tuple(product(range(dimension),repeat=3));index={x:i for i,x in enumerate(states)}
        out=np.zeros((len(states),len(states)),dtype=np.int64)
        for a,b in combinations(range(3),2):
            out+=np.eye(len(states),dtype=np.int64)
            for j,x in enumerate(states):
                y=list(x);y[a],y[b]=y[b],y[a]
                out[index[tuple(y)],j]-=1
        return states,out
    _,binary=triangle(2)
    assert not np.any(binary@(3*np.eye(8,dtype=np.int64)-binary))
    states,ternary=triangle(3)
    sign=np.array([0 if len(set(x))!=3 else (-1)**sum(x[i]>x[j] for i,j in combinations(range(3),2))
                   for x in states],dtype=np.int64)
    assert np.any(sign)
    assert np.array_equal(ternary@sign,6*sign)


def _complete_local_menus():
    edges=tuple(combinations(range(4),2));states,local,g=e.matrices(4,edges,2)
    identity=np.eye(len(states),dtype=object)
    zero=np.zeros_like(identity)
    n0=np.diag([state[0] for state in states]).astype(object)
    n2=np.diag([state[2] for state in states]).astype(object)
    symmetric=identity-local[0].astype(object)*F(1,2)
    pairs=((n0,identity-n0),(symmetric,identity-symmetric),(n2,identity-n2))
    labels=('first','second','reserved_zero')
    as_tuple=lambda value:tuple(tuple(x for x in row) for row in value.tolist())
    menus=tuple(tuple(zip(labels,(as_tuple(a),as_tuple(b),as_tuple(zero)))) for a,b in pairs)
    return states,g,labels,menus


def test_complete_history_law_noncommuting_local_menus():
    states,g,labels,menus=_complete_local_menus()
    kraus=tuple(e.complete_rational_instrument(labels,menu) for menu in menus)
    assert np.any(kraus[0][0]@kraus[1][0]-kraus[1][0]@kraus[0][0])
    lam,v=eigh(g.astype(float));m=6
    p=1-lam/(2*m);rho=(v*(p/p.sum()))@v.T
    h=-2*m*np.log(p)
    mu=float(np.dot(p/p.sum(),lam));a=1/(1-mu/(2*m))
    times=(F(1,10),F(3,10),F(7,10))
    def laws(energies):
        output={():rho.astype(complex)}
        for duration,menu in zip(times,kraus,strict=True):
            unitary=(v*np.exp(-1j*float(duration)*energies))@v.T
            following={}
            for history,density in output.items():
                evolved=unitary@density@unitary.conj().T
                for label,k in zip(labels,menu,strict=True):
                    matrix=np.asarray(k,dtype=float)
                    following[(*history,label)]=matrix@evolved@matrix.T
            output=following
        return {history:float(np.trace(rho).real) for history,rho in output.items()}
    actual= laws(h);exchange=laws(a*lam)
    assert len(actual)==len(exchange)==27
    assert abs(sum(actual.values())-1)<1e-13
    assert abs(sum(exchange.values())-1)<1e-13
    assert all(actual[word]==exchange[word]==0 for word in actual if 'reserved_zero' in word)
    tv=sum(abs(actual[word]-exchange[word]) for word in actual)/2
    assert tv>1e-7
    budget=e.finite_history_bound(4,6,24,2,3,(1,2,1),(labels,)*3,times,floor=F(1,2))
    assert tv<=float(budget['complete_history_total_variation_upper'])
    assert budget['full_declared_history_count']==27


def test_complete_menu_refuses_missing_zero_and_nonzero_outcomes():
    _,_,labels,menus=_complete_local_menus()
    with pytest.raises(ValueError):e.complete_rational_instrument(labels,menus[0][:-1])
    with pytest.raises(ValueError):e.complete_rational_instrument(labels,menus[0][1:])
    changed=menus[0][:-1]+(('reserved_zero',menus[0][0][1]),)
    with pytest.raises(ValueError):e.complete_rational_instrument(labels,changed)
    with pytest.raises(ValueError):e.complete_rational_instrument(labels,((labels[0],((True,),)),))


def test_fixed_finite_complete_menu_bound_converges():
    labels=(('first','second','reserved_zero'),)*3
    values=[]
    for level in (12,24,100):
        n=4**level
        row=e.finite_history_bound(240*n,630*n,5400*n,80*n,6,(1,2,1),labels,(F(1),)*3)
        values.append(row['complete_history_total_variation_upper'])
    assert values[0]>values[1]>values[2]
    assert values[2]<F(1,10**20)
    with pytest.raises(ValueError):
        e.finite_history_bound(240,630,5400,80,6,(1,2,1),labels,(F(1),)*2)
    with pytest.raises(ValueError):
        e.finite_history_bound(240,630,5400,80,6,(1,), (('first','first'),),(F(1),))
