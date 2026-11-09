"""Independent boundary, covariance, source and complete-detector checks."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest
import sympy as s

from . import model as m


def test_actual_open_flat_source_and_dilation():
    control=m.flat_control()
    assert control['flat_source_equality'] and control['dilation_identity']
    assert control['all_native_rows_retained_in_operator_proof']


def test_open_face_completion_and_physical_extremal_norm():
    for side in (1,2,3):
        shifts=m.open_shifts(side);data=m.shift_data(shifts);I=s.eye(side**3)
        missing=s.zeros(side**3)
        for i,t in enumerate(shifts):
            assert s.trace(data['holes'][i])==side**2
            assert data['L'][i][i]==2*I-t-t.H
            missing+=data['holes'][i]
        # Thirty positive and thirty negative beta rows per site.
        assert 30*s.trace(missing)/2==45*side**2
        for C in (s.eye(3),s.Matrix([[2,1,0],[1,2,0],[0,0,1]])):
            L=sum((C[i,j]*data['L'][i][j] for i,j in product(range(3),repeat=2)),s.zeros(side**3))
            assert min(np.linalg.eigvalsh(np.asarray(L,dtype=complex)))>=-1e-12


def test_local_gauge_and_internal_frame_covariance():
    shifts=m.open_shifts(2,(s.I,-1,1));V=s.Matrix([[1,s.Rational(1,3),0],[0,1,s.Rational(1,5)],[0,0,1]])
    h=m.matter(V,shifts)
    gauge=s.diag(1,s.I,-1,-s.I,s.I,1,-s.I,-1)
    lifted=s.kronecker_product(s.eye(4),gauge)
    changed=tuple(gauge*t*gauge.H for t in shifts)
    assert m.matter(V,changed)==m.geometry.clean(lifted*h*lifted.H)
    gamma=m.density.gamma_data();beta=gamma[3]
    alpha=[s.I*beta*g for g in gamma[:3]]
    rotation=s.Matrix([[0,-1,0],[1,0,0],[0,0,1]])
    spin=s.eye(4)-alpha[0]*alpha[1]
    U=s.kronecker_product(spin,s.eye(8))
    assert m.matter(rotation*V,shifts)==m.geometry.clean(U*h*U.H/2)


def test_six_strains_by_independent_coefficient_derivatives():
    t=s.symbols('t',real=True);shifts=m.open_shifts(2,(s.I,1,-1))
    data=m.shift_data(shifts);gamma=m.density.gamma_data();beta=gamma[3]
    alpha=[s.I*beta*g for g in gamma[:3]]
    for _,S in m.geometry.strain_basis():
        V=s.eye(3)+t*S;A=V.inv().T;C=(V.T*V).inv()
        B=(V.T*V)/V.det()
        assert B.diff(t).subs(t,0)==2*S-s.trace(S)*s.eye(3)
        direct=s.zeros(32)
        for i,j in product(range(3),repeat=2):
            direct+=s.diff(A[i,j],t).subs(t,0)*s.kronecker_product(alpha[i],data['D'][j])
            direct+=s.diff(C[i,j],t).subs(t,0)*s.kronecker_product(beta,data['L'][i][j])/2
        assert m.strain(S,shifts)==m.geometry.clean(direct)
        assert m.strain(S,shifts).H==m.strain(S,shifts)


def test_native_magnetic_missing_planes_and_gauge_covariance():
    loops=m.geometry.link_data(m.geometry.native_links())['loops']
    gauge=m.geometry.native_links()[0]
    for flags in product((False,True),repeat=3):
        selected=tuple(loop if present else None for loop,present in zip(loops,flags))
        Q=m.magnetic_gram(selected)
        assert Q==Q.T
        for i,loop in enumerate(selected):
            expected=0 if loop is None else -s.re(s.trace(loop))/15
            assert s.simplify(Q[i,i]-expected)==0
        changed=tuple(None if loop is None else gauge*loop*gauge.H for loop in selected)
        assert m.magnetic_gram(changed)==Q
        gram=Q+s.diag(*(int(x) for x in flags))
        assert min(np.linalg.eigvalsh(np.asarray(gram,dtype=float)))>=-1e-12


def test_complete_mean_detector_keeps_entangled_reference_and_disturbance():
    K=np.array([[.2,.3j],[-.3j,-.7]],complex);R=2.
    effects=[(np.eye(2)+sign*K/R)/2 for sign in (-1,1)]
    roots=[]
    for E in effects:
        val,vec=np.linalg.eigh(E);roots.append((vec*np.sqrt(val))@vec.conj().T)
    psi=np.array([np.sqrt(.7),0,0,np.sqrt(.3)],complex)
    rho_joint=np.outer(psi,psi.conj());rho=np.diag([.7,.3])
    blocks=[np.kron(B,np.eye(2))@rho_joint@np.kron(B.conj().T,np.eye(2)) for B in roots]
    assert np.allclose(sum(B.conj().T@B for B in roots),np.eye(2))
    assert abs(sum(np.trace(x) for x in blocks)-1)<1e-14
    assert abs(sum(sign*R*np.trace(x) for sign,x in zip((-1,1),blocks))-np.trace(rho@K))<1e-14
    complete=sum(blocks)
    reference=np.trace(complete.reshape(2,2,2,2),axis1=0,axis2=2)
    assert np.allclose(reference,rho)
    assert np.linalg.norm(complete-rho_joint)>.01


def test_full_primary_strain_bounds_and_source_forgery_rejection():
    member=m.contract.member(0);receipt=m.bind_member(member);m.verify_bound(member,receipt)
    result=receipt['result'];assert result['flat_primary_action_identical']
    assert len(result['all_six_spatial_strains'])==6
    assert F(result['all_six_spatial_strains'][0]['full_strain_norm_upper'])==F(71340,7)
    assert result['uncompleted_flat_boundary_error']['full_Fock_and_physical_norm_exact']=='180'
    zero=m.strain_bounds(member,s.zeros(3))
    assert zero['complete_record_values']==[0,0]
    assert zero['effects']=='I/2 for both zero-score records'
    assert not result['proper_metric_recovered'] and not result['symmetric_Lorentz_stress_derived']
    assert result['physical_particle_mass_table_GeV']==[]
    damaged=deepcopy(receipt);damaged['result']['proper_metric_recovered']=True
    with pytest.raises(ValueError):m.verify_bound(member,damaged)
    damaged=deepcopy(member);damaged['source']['electric_coefficient']='1'
    with pytest.raises(ValueError):m.bind_member(damaged)


def test_invalid_coframe_shift_and_strain_rejected():
    for shifts in ((s.eye(2),)*2,(2*s.eye(2),)*3):
        with pytest.raises(ValueError):m.shift_data(shifts)
    for frame in (s.diag(-1,1,1),s.zeros(3),s.eye(2)):
        with pytest.raises(ValueError):m.matter(frame,m.open_shifts(2))
    with pytest.raises(ValueError):m.strain(s.Matrix([[0,1,0],[0,0,0],[0,0,0]]),m.open_shifts(2))


def test_flat_action_does_not_uniquely_fix_spatial_response():
    shifts=m.open_shifts(1);h=m.matter(s.eye(3),shifts)
    A=s.diag(1,-1,1,-1);S=s.diag(1,0,0);t=s.symbols('t',real=True)
    V=s.eye(3)+t*S;extra=(s.trace(V.T*V)-3)*A
    assert extra.subs(t,0)==s.zeros(4)
    assert extra.diff(t).subs(t,0)==2*A
    assert h==(h+extra).subs(t,0)
