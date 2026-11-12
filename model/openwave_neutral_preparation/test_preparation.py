"""Native local states, complete filters and independent polar controls."""
import numpy as np
import scipy.linalg as la
import sympy as s
from model.coherent_defect_wave.chiral_joint.ground_composite_response.model import source,target
from model.openwave_neutral_stationary.isolation import data
from .local import line
from .filtering import parameters,certificate
from .records import certificate as records
from .replay import certificate as replay


def test_native_gram_householder_prepares_each_actual_line():
    for kind in ('neutral','ground_up','ground_down'):
        row=line(kind);block=target(2,0) if kind=='neutral' else source(1 if kind=='ground_up' else -1)
        g=np.array(block['gram'],complex);h=np.array(block['H'],complex)
        energies,vectors=la.eigh(g@h,g)
        v=vectors[:,0];j=row['seed_basis_column'];b=np.eye(len(g),dtype=complex)[:,j]/np.sqrt(g[j,j])
        overlap=np.vdot(v,g@b);v=v*(overlap/abs(overlap))
        w=b-v;u=np.eye(len(g))-2*np.outer(w,w.conj()@g)/np.vdot(w,g@w)
        assert np.linalg.norm(u.conj().T@g@u-g)<1e-9
        assert np.linalg.norm(u@b-v)<1e-12
        assert np.linalg.norm(h@(u@b)-energies[0]*(u@b))<1e-11
        p=abs(np.vdot(b,g@v))**2
        e=s.Symbol('e');root=sum(data()['en' if kind=='neutral' else 'e0'])/2
        numerator=np.array([complex(s.sympify(z).subs(e,root).evalf(65)) for z in row['line_numerator_coefficients_in_e']])
        denominator=complex(s.sympify(row['projector_denominator']).subs(e,root).evalf(65))
        explicit=numerator/(denominator*np.sqrt(g[j,j]*p))
        assert np.linalg.norm(explicit-v)<1e-9
        assert .01<p<.99


def test_rational_filter_bounds_and_resource_ancestry():
    p=parameters();c=certificate()
    assert p['pass_loss']**p['flattening_power']<=p['epsilon']
    assert p['flattening_power']*p['stop']/(1-p['flattening_power']*p['stop'])<=p['epsilon']
    assert p['success_lower']>s.Rational('0.99999999999639')
    assert p['trace_bound']<s.Rational('2.000001e-16')
    assert p['degree']==5091 and p['operator_evaluations']==30546
    assert not c['Stinespring_gate_count_or_runtime_derived']
    assert not c['stationarity_claim_about_arbitrary_joint_reference']


def test_independent_polar_instrument_and_omitted_correction_control():
    # Three unequal principal angles prevent a plain filtered product
    # mixture from being uniform on the spectral band.
    angles=np.array([.0007,.002,.004]);r=3;n=6;eps=1e-9
    p=np.diag([1.]*r+[0.]*r)
    columns=np.vstack((np.diag(np.cos(angles)),np.diag(np.sin(angles))))
    e=columns@columns.T;m=e+eps*(np.eye(n)-e)
    a=(p@m@m@p)[:r,:r];q0=(1-max(np.sin(angles))**2)*(1-eps)**2
    b=(3*np.eye(r)-a)/2;k=np.sqrt(q0)*m[:,:r]@b
    acceptance=k.T@k;vals=np.linalg.eigvalsh(acceptance)
    tau=(1-q0)**2/q0
    assert vals[0]>=q0*(1-tau)**2-1e-15 and vals[-1]<=q0+1e-15
    kfull=np.zeros((n,n));kfull[:,:r]=k
    failure=la.sqrtm(np.eye(n)-kfull.T@kfull)
    assert np.linalg.norm(kfull.T@kfull+failure.T@failure-np.eye(n))<1e-14
    rho=k@k.T;rho/=np.trace(rho)
    bound=4*tau/(1-tau)**2+2*eps/q0
    assert np.linalg.norm(rho-e/r,ord='nuc')<=bound
    plain=m@p@m;plain/=np.trace(plain)
    assert np.linalg.norm(plain-e/r,ord='nuc')>1000*bound
    # Completeness includes an arbitrary entangled untouched reference.
    psi=np.eye(n).reshape(-1)/np.sqrt(n)
    assert abs(np.linalg.norm(np.kron(kfull,np.eye(n))@psi)**2+
               np.linalg.norm(np.kron(failure,np.eye(n))@psi)**2-1)<1e-13


def test_reduced_stationarity_does_not_imply_joint_reference_stationarity():
    psi=np.array([1,0,0,1],complex)/np.sqrt(2)
    evolved=np.kron(np.diag([1,np.exp(-.7j)]),np.eye(2))@psi
    assert np.linalg.norm(np.outer(psi,psi.conj())-np.outer(evolved,evolved.conj()))>.4
    for state in (psi,evolved):
        matrix=state.reshape(2,2)
        assert np.allclose(matrix@matrix.conj().T,np.eye(2)/2)


def test_complete_recipe_momentum_energy_position_records():
    c=records()
    assert c['complete_records_per_recipe']==401
    assert c['complete_original_record_count']==410624
    assert c['archive_not_discarded_to_prove_channel_completeness']
    assert not c['position_and_momentum_joint_sharp_values_claimed']


def test_all_native_local_filter_modes_and_failure_records():
    c=replay();rows=c['all_modes']
    assert sum(row['multiplicity'] for row in rows)==16
    assert sum(row['multiplicity'] for row in rows if row['selected_neutral_line'])==1
    assert all('complete_failure_probability_decimal' in row for row in rows)
    assert not c['full_many_cell_H_materialized']
