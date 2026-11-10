"""Exact coin, full-CAR parity, original probabilities and reference controls."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product
import json
import math

import numpy as np
import pytest
from scipy.linalg import expm

from . import model as m
from model.particle_local_wilson_density import density


def _annihilators(count):
    ops=[]
    for j in range(count):
        op=np.zeros((2**count,2**count))
        for bits in range(2**count):
            if bits>>j&1:op[bits^(1<<j),bits]=(-1)**((bits&((1<<j)-1)).bit_count())
        ops.append(op)
    return ops


def _polynomial_state(H):
    E,V=np.linalg.eigh(H);x=(max(E)-E)/16
    q=sum(x**j/math.factorial(j) for j in range(6))**8
    p=q*q;p/=p.sum()
    return E,V,p,(V*p)@V.conj().T


def test_machin_identity_and_outward_algebraic_coin():
    tangent4=2*F(5,12)/(1-F(5,12)**2)
    assert (tangent4-F(1,239))/(1+tangent4/F(239))==1
    lo,hi=m.pi_interval()
    import mpmath as mp
    with mp.workdps(80):
        assert mp.mpf(lo.numerator)/lo.denominator<mp.pi<mp.mpf(hi.numerator)/hi.denominator
        for n in (0,1,2,4):
            R=2**n
            p=mp.mpf(R)/2*sum((1 if u==R else 2)*(1 if v==R else 2)*
                (1-mp.cos(mp.pi*u*v/R**4)) for u in range(1,R+1) for v in range(1,R+1))
            lower,upper=m.coin_interval(n)
            assert mp.mpf(lower.numerator)/lower.denominator<=p<=mp.mpf(upper.numerator)/upper.denominator
    assert m.coin_interval(0)==(1,1)
    with pytest.raises(ValueError):m.coin_interval(-1)


def test_exact_periodic_profile_probability_and_open_bond_work_coefficients():
    for n in (0,1,2):
        R=2**n;L=2*R*R;a=1/R
        c=np.arange(L);t=np.maximum(1-np.minimum(c,L-c)/R,0)
        f=np.einsum('i,j,k->ijk',t,t,t)
        coefficients=[np.sum(1-np.cos(np.pi*a*(np.roll(f,-1,axis=i)-f))) for i in range(3)]
        assert max(coefficients)-min(coefficients)<1e-13
        lo,hi=m.coin_interval(n)
        assert float(lo)-1e-13<=sum(coefficients)/12<=float(hi)+1e-13
        # Every open bond sees the same average coefficient, even at a boundary.
        for x in ((0,0,0),(L-2,L-2,L-2)):
            y=(x[0]+1,x[1],x[2]);total=0.
            for X in product(range(L),repeat=3):
                fx=f[tuple((x[i]-X[i])%L for i in range(3))]
                fy=f[tuple((y[i]-X[i])%L for i in range(3))]
                total+=1-np.cos(np.pi*a*(fy-fx))
            assert total==pytest.approx(coefficients[0],abs=2e-12)


def test_full_car_parity_conjugation_and_nonstationary_work_identity():
    c=_annihilators(4);nums=[x.T@x for x in c];I=np.eye(16)
    P=[np.diag([(-1)**sum((b>>j)&1 for j in sites) for b in range(16)])
       for sites in ((0,1),(2,3))]
    T=c[0].T@c[2]+c[2].T@c[0]+.7j*(c[1].T@c[3]-c[3].T@c[1])
    onsite=2*nums[0]-3*nums[1]+nums[2]-nums[3];H=onsite+T
    for x in range(2):
        assert np.array_equal(P[x]@P[x],I)
        assert np.allclose(P[x]@H@P[x]-H,-2*T)
        for j in range(4):
            assert np.array_equal(P[x]@c[j]@P[x],(-1)**(j//2==x)*c[j])
            assert np.array_equal(P[x]@nums[j],nums[j]@P[x])
    psi=np.arange(1,17)+1j*np.arange(16,0,-1);psi=psi/np.linalg.norm(psi)
    rho=np.outer(psi,psi.conj())
    assert abs(np.trace(rho@(H@T-T@H)))>0
    W=sum(np.trace(rho@(U@H@U-H)).real for U in P)/2
    assert W==pytest.approx(-4*np.trace(rho@T).real/2)


def test_literal_dirac_cube_control_matches_both_work_forms():
    # Algebra control only: flat one-particle links are not the primary selected Gauss state.
    cells=list(product(range(2),repeat=3));index={x:i for i,x in enumerate(cells)}
    gamma=[np.array(g.tolist(),dtype=complex) for g in density.gamma_data()];beta=gamma[3]
    H=np.kron(np.eye(8),3*beta);hop=np.zeros_like(H)
    for x in cells:
        for axis in range(3):
            if x[axis]:continue
            y=list(x);y[axis]=1;y=tuple(y);s=index[x];t=index[y]
            block=(-beta+beta@gamma[axis])/2
            hop[4*s:4*s+4,4*t:4*t+4]=block
            hop[4*t:4*t+4,4*s:4*s+4]=block.conj().T
    H+=hop;_,_,_,rho=_polynomial_state(H)
    phase=.71;coefficient=2*(1-math.cos(phase));coin=coefficient/4
    Wphase=Wcontrol=0.
    for x in range(8):
        local=np.zeros(8);local[x]=1
        P=np.kron(np.diag(1-2*local),np.eye(4))
        U=np.kron(np.diag(np.exp(1j*phase*local)),np.eye(4))
        Wphase+=np.trace(rho@(U.conj().T@H@U+U@H@U.conj().T-2*H)).real/16
        Wcontrol+=coin*np.trace(rho@(P@H@P-H)).real/8
    assert Wphase==pytest.approx(Wcontrol,abs=2e-15)
    assert Wcontrol==pytest.approx(-coefficient*np.trace(rho@hop).real/8,abs=2e-15)
    assert Wcontrol>0


def test_complete_reference_output_and_coherent_pointer_bound():
    rng=np.random.default_rng(47);psi=rng.normal(size=8)+1j*rng.normal(size=8);psi/=np.linalg.norm(psi)
    rho=np.outer(psi,psi.conj());P=np.diag([1,-1,-1,1]);U=np.kron(P,np.eye(2))
    for p in (.003,.4,1.):
        delta=p*(U@rho@U.conj().T-rho)
        half=np.sum(abs(np.linalg.eigvalsh(delta)))/2
        assert half<=p+1e-14
        control=np.concatenate((math.sqrt(1-p)*psi,math.sqrt(p)*U@psi))
        idle=np.concatenate((math.sqrt(1-p)*psi,math.sqrt(p)*psi))
        trace=math.sqrt(max(0,1-abs(np.vdot(control,idle))**2))
        assert trace<=min(1,2*math.sqrt(p))+1e-14
        # Tracing the original purifier is unnecessary in either bound.
        assert np.linalg.matrix_rank(rho.reshape(4,2,4,2).trace(axis1=1,axis2=3),tol=1e-12)==2


def test_original_transfer_mixture_and_two_sided_moment_inequality():
    H=np.array([[0,1,0],[1,3,.6],[0,.6,7]],float);P=np.diag([1,-1,1])
    E,V,weights,rho=_polynomial_state(H)
    transition=abs(V.T@P@V)**2;law=weights[:,None]*transition.T
    omega=E[None,:]-E[:,None]
    assert law.sum()==pytest.approx(1)
    for p in (.01,.5):
        W=p*np.sum(law*omega);M2=p*np.sum(law*omega**2)
        assert W>0 and M2>=W*W/p-1e-14
        assert M2<=p*np.linalg.norm(P@H@P-H,2)**2+1e-14
        for K in (1.,5.,15.):
            tail=p*np.sum(law*omega*(omega>K))
            assert tail>=W-K*p-1e-14


def test_general_profile_work_coin_is_eventually_admitted():
    for shape,amplitude in (((1,1,1),1),((2,1,1),2),((2,2,2),F(3,2))):
        values=[m.coin_recipe(n,shape,amplitude) for n in (3,4,5,6)]
        assert values[-1]['probability_admitted_by_bound']
        assert all(values[j+1]['probability_upper']<values[j]['probability_upper'] for j in range(3))
        for x in values:
            assert x['same_profile_work_exactly_matched_when_probability_at_most_one']
            assert x['probability_upper']/x['spacing']<100


def test_actual_members_clock_ids_and_original_probability_ancestry():
    for n in (0,1,2,8):
        member=m.contract.member(n);bound=m.bind_member(member);r=bound['result']
        assert json.loads(json.dumps(bound))==bound
        m.verify_bound(member,bound)
        assert r['source_state_reference_H_and_clock_unchanged']
        assert not r['apparatus_identical_to_critical_phase_menu']
        assert not r['critical_phase_collapse_proved']
        assert not r['original_critical_current_covariance_evaluated']
        assert r['physical_mass_table_GeV']==[]
        lo,hi=map(F,r['original_second_moment_interval']);a=F(member['source']['spacing'])
        assert 0<lo<=hi<=182250/a
        assert F(r['complete_reported_work_lower'])>0
        assert r['retained_origin_coin_outcomes']==2*member['source']['side']**3
        assert F(r['original_nonelastic_probability_upper'])<=F(r['uniform_coin_upper'])
        bad=deepcopy(bound);bad['result']['critical_phase_collapse_proved']=True
        with pytest.raises(ValueError):m.verify_bound(member,bad)
    altered=deepcopy(member);altered['source']['spacing']='1'
    with pytest.raises(ValueError):m.bind_member(altered)


def test_probability_envelope_and_uniform_moment_constants():
    for n in range(1,15):
        a=F(1,2**n);g=m.critical.profile_geometry(n)
        assert a*g['three_axis_gradient_square_sum']==F(2,3)*(2+a*a)**2
        assert F(5,12)*a*a*g['three_axis_gradient_square_sum']<=F(45,32)*a<1
    assert 3600*6**2*F(45,32)==182250
    assert F(32,45)<1
