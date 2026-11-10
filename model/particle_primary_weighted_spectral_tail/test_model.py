"""Independent stationary moments, source incidence and exact bound controls."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product
import json
import math

import numpy as np
import pytest
from scipy.linalg import expm

from . import model as m


def _comm(A,B):return A@B-B@A


def _flat_dirac_control():
    """One-particle flat-link algebra control, not the primary physical state."""
    gamma=[np.array(x.tolist(),dtype=complex) for x in m.work.density.gamma_data()]
    beta=gamma[3];hop=(-beta+beta@gamma[0])/2
    H=np.kron(np.eye(3),3*beta)
    for x in range(2):
        H[4*x:4*x+4,4*x+4:4*x+8]=hop
        H[4*x+4:4*x+8,4*x:4*x+4]=hop.conj().T
    E,V=np.linalg.eigh(H);inverse_temperature=2.;scale=32;order=5
    y=inverse_temperature*(max(E)-E)/(2*scale)
    Q=sum(y**k/math.factorial(k) for k in range(order+1))**scale
    p=Q*Q; p/=sum(p);rho=(V*p)@V.conj().T
    F0=np.kron(np.diag([.15,-.3,.25]),np.eye(4))
    units=[expm(1j*F0),expm(-1j*F0)]
    T=sum(abs((V.conj().T@U@V).T)**2 for U in units)/2
    paths=p[:,None]*T;omega=E[None,:]-E[:,None]
    return H,E,V,rho,p,units,paths,omega


def test_exact_stationary_third_moment_and_sign_cubic_identity():
    H,E,V,rho,p,units,P,omega=_flat_dirac_control()
    assert np.linalg.norm(_comm(H,rho),2)<2e-15
    D=[U.conj().T@H@U-H for U in units];C=(D[0]+D[1])/2;S=(D[0]-D[1])/2
    m2=float(np.sum(P*omega**2));m3=float(np.sum(P*omega**3))
    expected2=np.trace(rho@(C@C+S@S)).real
    expected3=np.trace(rho@(.5*(_comm(C,_comm(H,C))+_comm(S,_comm(H,S)))+
                                  C@C@C+3*S@C@S+_comm(S,_comm(S,C)))).real
    assert m2==pytest.approx(expected2,abs=2e-14)
    assert m3==pytest.approx(expected3,abs=2e-13)
    assert np.linalg.norm((D[0]@D[0]@D[0]+D[1]@D[1]@D[1])/2-
            C@C@C-3*S@C@S-_comm(S,_comm(S,C)),2)<1e-14
    c=np.linalg.norm(C,2)
    A=(np.linalg.norm(_comm(C,_comm(H,C)),2)+np.linalg.norm(_comm(S,_comm(H,S)),2))/2+np.linalg.norm(_comm(S,_comm(S,C)),2)
    assert m3<=A+3*c*m2
    # Exact selected source stationarity is essential to the spectral identity.
    bad=units[0]@rho@units[0].conj().T
    cubic_sum=sum(np.trace(bad@(U.conj().T@np.linalg.matrix_power(H,3)@U-
            3*H@U.conj().T@H@H@U+3*H@H@U.conj().T@H@U-H@H@H)).real for U in units)/2
    stationary_expression=np.trace(bad@(.5*(_comm(C,_comm(H,C))+_comm(S,_comm(H,S)))+
            C@C@C+3*S@C@S+_comm(S,_comm(S,C)))).real
    assert abs(cubic_sum-stationary_expression)>1e-5


def test_actual_owner_incidence_proves_nested_commutator_coefficient():
    for side in (2,3):
        source=dict(m.contract.member(0)['source']);source['side']=side
        cells=list(product(range(side),repeat=3));terms={}
        for cell in cells:
            for t in m.stress._cell_terms(source,cell):terms[t['id']]=t
        edges=[set(t['owner_support']) for t in terms.values() if t['kind']=='bond']
        assert len(edges)==3*side*side*(side-1)
        assert max(sum(x in e for e in edges) for x in cells)<=6
        assert max(len(t['owner_support']) for t in terms.values())<=3
        kappa=max(sum((t['unallocated_norm_upper'] for t in terms.values() if x in t['owner_support']),F()) for x in cells)
        K=np.zeros((len(edges),len(edges)))
        for t in terms.values():
            support=set(t['owner_support']);weight=float(t['unallocated_norm_upper'])
            for c,edge in enumerate(edges):
                if edge.isdisjoint(support):continue
                union=support|edge
                for b,outer in enumerate(edges):
                    if not outer.isdisjoint(union):K[b,c]+=weight
        assert max(K.sum(axis=0))<=60*float(kappa)
        assert max(K.sum(axis=1))<=60*float(kappa)
        rng=np.random.default_rng(side);weights=rng.uniform(0,2,len(edges))
        assert 4*(weights@K@weights)<=240*float(kappa)*(weights@weights)
        analytic=270+F(9,7)+F(1827,2)
        assert kappa<=analytic


def test_local_cosine_current_bounds_use_exact_discrete_profile():
    for n in (0,1,2):
        R=2**n;side=2*R*R;a=1/R;coordinates=np.indices((side,side,side));profile=np.ones((side,side,side))
        for axis in range(3):profile*=np.maximum(1-np.minimum(coordinates[axis],side-coordinates[axis])/R,0)
        d2=0.;squares=0.;cosine_sum=0.
        for axis in range(3):
            gradient=np.roll(profile,-1,axis=axis)-profile;phase=np.pi*a*gradient
            assert np.max(abs(phase))<=np.pi*a*a+1e-15
            cosine=30/a*(1-np.cos(phase));sine=30/a*np.sin(phase)
            assert np.max(cosine)<=150*a**3
            d2+=np.sum(gradient**2);squares+=np.sum(cosine**2+sine**2);cosine_sum+=sum(cosine.ravel())
        member=m.contract.member(n);row=m.local_current_caps(member)
        # Origin averaging replaces periodic coefficient-list bonds by the
        # exact fraction of actual open bonds in each orientation.
        fraction=member['source']['interior_hopping_bonds']/(3*side**3)
        assert squares*fraction<=float(row['origin_averaged_cosine_plus_sine_bond_norm_squares_upper'])+1e-11
        if n:assert cosine_sum<=float(row['cosine_operator_norm_upper'])
        else:assert cosine_sum/2==pytest.approx(float(row['cosine_operator_norm_upper']))
        assert a*d2<=6


def test_selected_negative_moments_and_quadratic_closure():
    H,E,V,rho,p,units,P,omega=_flat_dirac_control()
    slope=m.work.modular.polynomial_slope(2,str(min(E)),str(max(E)),5,32)
    b=float(slope['modular_slope_lower']);neg=np.maximum(-omega,0)
    n1=np.sum(P*neg);n2=np.sum(P*neg**2);n3=np.sum(P*neg**3)
    assert n1<=1/(3*b) and n2<=2/b**2 and n3<=6/b**3
    W=np.sum(P*omega);M2=np.sum(P*omega**2);M3=np.sum(P*omega**3)
    D=[U.conj().T@H@U-H for U in units];C=(D[0]+D[1])/2;S=(D[0]-D[1])/2
    A=(np.linalg.norm(_comm(C,_comm(H,C)),2)+np.linalg.norm(_comm(S,_comm(H,S)),2))/2+np.linalg.norm(_comm(S,_comm(S,C)),2)
    B=3*np.linalg.norm(C,2);P1=W+1/(3*b)
    assert max(M2-2/b**2,0)**2<=P1*(A+B*M2+6/b**3)
    cap=m.quadratic_moment_cap(*map(str,(P1,A,B,2/b**2,6/b**3)))
    assert float(cap)>=M2 and M3<=A+B*M2


def test_reference_complete_signed_menu_preserves_original_normalization():
    H,E,V,rho,p,units,P,omega=_flat_dirac_control();root=(V*np.sqrt(p))@V.conj().T
    reference=np.zeros_like(rho);completeness=np.zeros_like(rho)
    for U in units:
        K=U/np.sqrt(2);A=K@root
        reference+=A.T@A.conj();completeness+=K.conj().T@K
    assert np.linalg.norm(reference-rho.T,2)<2e-15
    assert np.linalg.norm(completeness-np.eye(len(E)),2)<1e-15
    assert P.sum()==pytest.approx(1.)
    W=float(np.sum(P*omega));M2=float(np.sum(P*omega**2));band=m.work.local.positive_window(str(W),str(M2))
    left,right=map(float,band['original_transfer_interval'])
    assert P[(omega>=left)&(omega<=right)].sum()>=float(band['original_weight_lower'])


def test_rational_sqrt_and_quadratic_solver_round_outward():
    for x in (F(),F(1),F(2),F(17,19),F(10**90+7,3)):
        for bits in (0,8,32):
            q=m.sqrt_upper(x,bits);unit=F(1,2**bits)
            assert q*q>=x
            if q:assert (q-unit)**2<x
    with pytest.raises(ValueError):m.sqrt_upper(-1)
    with pytest.raises(ValueError):m.quadratic_moment_cap(0,1,1,1,1)
    with pytest.raises(ValueError):m.quadratic_moment_cap(1,-1,1,1,1)
    # Independent exact quadratics include the branch M2<=N2.
    for P,A,B,N2,N3 in ((1,0,0,0,0),(3,7,4,2,1),(F(1,9),9,0,4,2)):
        cap=m.quadratic_moment_cap(P,A,B,N2,N3)
        assert cap>=N2 and (cap-N2)**2>=P*(A+B*cap+N3)


def test_actual_source_caps_improve_the_large_capacity_weighted_bound():
    improvements=[]
    for n in (0,1,2,8):
        member=m.contract.member(n);bound=m.bind_member(member);m.verify_bound(member,bound);r=bound['result']
        a=F(member['source']['spacing']);ell=member['source']['per_link_word_degree']
        assert F(r['selected_second_moment_upper'])<=4000000/a
        assert F(r['selected_second_moment_upper'])<=F(r['inherited_operator_second_moment_upper'])
        improvements.append(float(F(r['second_moment_cap_improvement_factor'])))
        ideal=r['ideal_original_positive_window'];reported=r['finite_reported_positive_window']
        assert F(reported['reported_transfer_interval'][0])>0
        assert F(reported['original_probability_lower'])>F(ideal['original_weight_lower'])/2
        assert not reported['conditioned_or_postselected_probability']
        assert not r['uniform_selected_second_moment_proved']
        assert not r['uniform_selected_energy_weighted_tail_proved']
        assert not r['fixed_positive_energy_probability_limit_proved']
        assert r['physical_mass_table_GeV']==[]
        assert json.loads(json.dumps(bound,sort_keys=True,allow_nan=False))==bound
    assert improvements[1]>2 and improvements[2]>15 and improvements[3]>1000000
    bad=deepcopy(bound);bad['result']['selected_second_moment_upper']='1'
    with pytest.raises(ValueError):m.verify_bound(member,bad)
    damaged=deepcopy(member);damaged['source']['electric_coefficient']='1'
    with pytest.raises(ValueError):m.bind_member(damaged)


def test_global_growth_envelope_constants_do_not_hide_a_uniform_bound():
    assert 1+1519*444<700000
    assert 444*(14100000000+1520)<3000000**2
    previous=F()
    for n in (0,1,2,4,8):
        member=m.contract.member(n);caps=m.local_current_caps(member);a=F(member['source']['spacing']);ell=n+1
        assert caps['signed_third_constant_upper']<=19440000000*ell**2/a**2
        assert caps['selected_electric_third_constant_upper']<=14100000000/a**2
        assert caps['signed_third_M2_coefficient_upper']<=1519
        bound=4000000/a;assert bound>previous;previous=bound


def test_native_electric_derivatives_and_selected_energy_density():
    control=m.native_derivative_control();data=m.contract.characteristics.native_generators()
    assert sum(control['native_Casimir_diagonal'])==12
    for i,T in enumerate(data['generators']):
        matrix=np.array(T.tolist(),dtype=complex)/np.sqrt(float(data['gram'][i,i]))
        assert np.trace(matrix@matrix).real==pytest.approx(1.)
        assert np.linalg.norm(matrix,2)**2<=1+1e-15
    for n in (0,1,2,8):
        member=m.contract.member(n);s=member['source'];a=F(s['spacing']);V=s['side']**3
        row=m.local_current_caps(member)
        trial=F(-90*V)/a;mean_upper=trial+F(member['state']['target_ground_excess'])
        matter_lower=trial;hopping_lower=-F(30*s['interior_hopping_bonds'])/a
        magnetic_lower=-F(s['magnetic_coefficient'])*s['interior_plaquettes']/a
        electric_upper=mean_upper-matter_lower-hopping_lower-magnetic_lower
        assert electric_upper/V==row['selected_electric_energy_upper_per_cell_from_actual_box']
        assert electric_upper<91*V/a
        assert 2*a*electric_upper/105<F(182*V,105)
    assert F(12*182,35)<64


def test_electric_form_identity_survives_product_projection_without_false_commutation():
    rng=np.random.default_rng(91);L=np.diag([-3.,-2.,-.5,.2,1.,4.])
    raw=rng.normal(size=(6,6))+1j*rng.normal(size=(6,6));Q=raw+raw.conj().T
    P=np.diag([1.,1.,1.,1.,0.,0.]);cutL=(P@L@P)[:4,:4];cutQ=(P@Q@P)[:4,:4]
    assert np.linalg.norm(_comm(P,L),2)==0
    assert np.linalg.norm(_comm(P,Q),2)>1
    assert np.linalg.norm(_comm(cutL,cutQ)-(P@_comm(L,Q)@P)[:4,:4],2)<1e-15
    assert np.linalg.norm(_comm(cutL,cutQ),2)<=np.linalg.norm(_comm(L,Q),2)
    Z=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4));rho=Z@Z.conj().T;rho/=np.trace(rho)
    J=1j*_comm(cutL,cutQ);D=_comm(cutQ,_comm(cutL,cutQ))
    left=_comm(cutQ,_comm(cutL@cutL,cutQ));right=2*J@J+cutL@D+D@cutL
    assert np.linalg.norm(left-right,2)<1e-12
    expectation=np.trace(rho@left).real
    upper=2*np.linalg.norm(J,2)**2+2*np.linalg.norm(D,2)*np.sqrt(np.trace(rho@cutL@cutL).real)
    assert abs(expectation)<=upper
