"""Independent source-matrix, symmetry and complete-score checks."""
from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest

from model.full_fock_field.field import FullFock
from model.source_family_compatibility import composite, momentum_response as m


def matrix(poly,sites):
    space=FullFock(2*sites);identity=np.eye(1<<(2*sites),dtype=complex)
    gammas=[]
    for j in range(2*sites):
        a=space.annihilator(j).toarray()
        gammas.extend((a+a.conj().T,-1j*(a-a.conj().T)))
    result=np.zeros_like(identity)
    for mask,(re,im)in poly.items():
        value=identity
        for j,gamma in enumerate(gammas):
            if mask>>j&1:value=value@gamma
        result+=(float(re)+1j*float(im))*value
    return result


@pytest.mark.parametrize('sites',[2,3])
def test_exact_majorana_codec_against_independent_full_source_matrices(sites):
    data=composite.source((sites,1,1),F(3,2),F(4))
    for i in range(2*sites):
        np.testing.assert_array_equal(matrix(m.annihilator(i),sites),data['a'][i])
        np.testing.assert_array_equal(matrix(m.dressed(i),sites),data['b'][i])
    exact=m.hopping(tuple((x,y,F(3,4)*sign)for x,y,sign in data['edges']))
    T=data['H'].copy()
    np.fill_diagonal(T,0)
    np.testing.assert_array_equal(matrix(exact,sites),T/data['scale'])
    for x,y,_ in data['edges']:
        for spin in(0,1):
            i=2*x+spin;b=data['b'][i]
            edge=m.hopping(((x,y,F(1)),))
            t=matrix(edge,sites)
            np.testing.assert_array_equal(matrix(m.bond_derivative(x,y,spin),sites),b@t-t@b)


def test_disjoint_odd_anticommutator_is_zero_with_full_car_signs():
    left=m.bond_derivative(0,1);right=m.dagger(m.bond_derivative(2,3))
    assert m.anticommutator(left,right)=={}
    # A two-site odd generator anticommutes with an independent onsite odd probe.
    a=matrix(m.bond_derivative(0,1),3);b=matrix(m.annihilator(4),3)
    np.testing.assert_array_equal(a@b+b@a,np.zeros_like(a))


def test_actual_source_cubic_and_magnetic_lifts_cancel_mixed_axis_coefficients():
    row=m.source_symmetry_audit()
    assert row['proper_rotations']==24
    assert row['exact_link_translation_checks']==1125
    assert row['exact_link_rotation_checks']==row['exact_rotation_translation_normalizations']==9000
    assert len(row['mixed_axis_cancellations'])==12
    for witness in row['mixed_axis_cancellations']:
        d,p,s=witness['displacement'],witness['permutation'],witness['axis_signs']
        assert m.rotate(d,p,s)==d and m.rotation_phase(d,p)==-1
        assert m.permutation_sign(p)*np.prod(s)==1


def test_omitting_the_source_restoring_gauge_is_rejected(monkeypatch):
    monkeypatch.setattr(m,'rotation_phase',lambda x,p:1)
    with pytest.raises(ArithmeticError,match='proper cubic lift'):
        m.source_symmetry_audit()


def test_original_weyl_translation_has_exact_hidden_magnetic_and_spin_action():
    from model.field_symmetry.symmetry import onsite_frame, SIGMA
    for x in product(range(1,4),repeat=3):
        for j in range(3):
            y=list(x);y[j]-=1
            actual=onsite_frame(x)@onsite_frame(tuple(y)).conj().T
            expected=m.magnetic_phase(x,j)*1j*np.asarray(SIGMA[j])
            np.testing.assert_array_equal(actual,expected)


def test_spin_averaged_local_derivatives_equal_actual_two_orbital_correlators():
    inv=m.invariant_observables()  # exact seven-site polynomial identities checked here
    for name in('q1','q2'):
        sites=2 if name=='q1'else 3
        values=np.linalg.eigvalsh(matrix(inv[name],sites))
        assert values[0]>-.750000000000001 and values[-1]<.250000000000001
        assert set(np.round(values,8))=={-.75,0.,.25}
    D=[]
    for spin in(0,1):
        D.append(matrix(m.bond_derivative(0,1,spin),2))
    gram=sum(x@x.conj().T+x.conj().T@x for x in D)/2
    np.testing.assert_array_equal(gram,3*np.eye(16)-8*matrix(inv['q1'],2))


def test_axial_coefficient_replays_the_actual_full_three_site_ground():
    row=composite.finite_case(3)
    data=composite.source((3,1,1),row['c'],row['U'])
    numerator,denominator,_=composite.ground(data,row['ground_energy'])
    rho=numerator/denominator
    q2=matrix(m.invariant_observables()['q2'],3)
    expected=float(row['c']**2)*np.trace(rho@q2).real
    C=np.array(row['C'],float);R=np.array(row['residual_Gram'],float)
    G=R+C@C
    assert abs(G[0,4]-expected)<2e-14
    # Finite boundary breaks translation invariance; the middle diagonal differs.
    assert abs(G[0,0]-G[2,2])>1


def test_complete_local_score_plans_reconstruct_actual_observables():
    inv=m.invariant_observables()
    counts={'q1':6,'q2':6,'qstar':216}
    norms={'q1':F(3,4),'q2':F(3,4),'qstar':F(9,2)}
    for name,poly in inv.items():
        plan=m.score_plan(poly)
        assert len(plan['terms'])==counts[name] and plan['score_absolute_bound']==norms[name]
        assert sum(plan['setting_probabilities'])==1
        rebuilt={}
        for term in plan['terms']:
            inds=term['majorana_indices'];degree=len(inds);phase=degree*(degree-1)//2%4
            unit=(m._ONE,m._I,(-F(1),F(0)),(F(0),-F(1)))[phase]
            mask=sum(1<<i for i in inds)
            word={mask:unit}
            assert m.multiply(word,word)=={0:m._ONE} and m.dagger(word)==word
            rebuilt=m.add(rebuilt,m.scale(word,term['coefficient']))
        assert rebuilt==poly
    zero=m.score_plan({})
    assert zero['score_absolute_bound']==0 and zero['zero_observable_identity_record']


def test_readout_polynomials_are_detached_on_repeated_calls():
    first=m.local_observables();second=m.local_observables()
    first['g2_over_c_squared'].clear()
    assert second['g2_over_c_squared'] and m.local_observables()['g2_over_c_squared']


def test_finite_magnetic_bloch_matrix_and_range_two_symbol_match_exact_source():
    length=6;points=tuple(product(range(length),repeat=3));index={x:i for i,x in enumerate(points)}
    H=np.zeros((len(points),)*2);G=np.eye(len(points))*7
    for x in points:
        for axis in range(3):
            y=list(x);y[axis]=(y[axis]+1)%length
            H[index[x],index[tuple(y)]]+=m.signed_link(x,axis)/2
            H[index[tuple(y)],index[x]]+=m.signed_link(x,axis)/2
            for sign in(-1,1):
                y=list(x);y[axis]=(y[axis]+2*sign)%length
                G[index[x],index[tuple(y)]]+=.25
    parities=tuple(product(range(2),repeat=3))
    for k in((0.,0.,0.),(np.pi/3,0.,-np.pi/3),(-np.pi/3,np.pi/3,np.pi/3)):
        V=np.zeros((len(points),8),complex)
        for i,x in enumerate(points):
            s=tuple(t%2 for t in x)
            V[i,parities.index(s)]=np.exp(1j*np.dot(k,x))/np.sqrt(27)
        symbol=m.symbols(k,-1.,7.,.25)
        np.testing.assert_allclose(V.conj().T@V,np.eye(8),atol=2e-15)
        np.testing.assert_allclose(H@V,V@symbol['h'],atol=8e-16)
        np.testing.assert_allclose(G@V,V@symbol['Y_Gram'],atol=3e-15)


def test_full_momentum_positivity_uses_both_exact_endpoints():
    good=m.endpoint_bounds(F(-1),F(7),F(1,4))
    assert good['full_momentum_positivity']
    assert good['residual_endpoint_S0']==F(11,2)
    assert good['residual_endpoint_S3']==F(11,2)
    bad=m.endpoint_bounds(F(-2),F(7),F(0))
    assert not bad['full_momentum_positivity'] and bad['residual_minimum']==-5
    other=m.endpoint_bounds(F(0),F(1),F(1))
    assert not other['full_momentum_positivity'] and other['residual_endpoint_S0']==-5


def test_uniform_symbol_error_and_first_four_actual_composite_moments():
    k=(.13,-.41,.72);alpha=F(3,2)
    a=m.symbols(k,-1.2,10.,-.3);b=m.symbols(k,-1.19,10.02,-.301)
    budget=m.uniform_error(F(1,100),F(1,50),F(1,1000))
    assert np.linalg.norm(a['residual_Gram']-b['residual_Gram'],2)<=float(budget['residual_operator_error'])
    h,C,G=a['h'],a['C'],a['Y_Gram'];I=np.eye(8)
    block=np.block([[h,float(alpha)*I],[float(alpha)*I,C]])
    exact4=h@h@h@h+float(alpha**2)*(3-2.4)*(h@h)+float(alpha**4)*I+float(alpha**2)*G
    np.testing.assert_allclose(np.linalg.matrix_power(block,4)[:8,:8]+float(alpha**2)*a['residual_Gram'],exact4,atol=3e-14)


def test_tracial_control_is_distinguished_from_the_actual_ground():
    row=m.exact_algebra_audit()
    assert row['tracial_control_g0_over_c_squared']==F(9,2)
    assert row['tracial_control_g2_over_c_squared']==0
    assert not row['tracial_control_is_a_ground_state']
