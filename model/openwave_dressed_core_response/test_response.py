from fractions import Fraction as F
import numpy as np
from scipy.linalg import eigh
from .response import (exact_polynomial, action, norm_square, exact_residual,
                       ground_member, response_records, probe_records, CHARGE)
from model.openwave_dirac_maxwell_refinement.operators import polynomial_member


def test_exact_integer_polynomial_and_both_spin_blocks():
    cap=3;m=F(1);r=F(1,16)
    o=polynomial_member(cap,CHARGE,m,regulator=r)
    basis,rows,den=exact_polynomial(cap,CHARGE,m,r);n=len(basis)
    h=np.zeros((n,n),complex)
    for i,row in enumerate(rows):
        for j,(a,b) in row.items():h[i,j]=complex(float(F(a,den)),float(F(b,den)))
    assert np.max(abs(h-o['H_minus_mass'][:n,:n]))<1e-14
    assert np.array_equal(o['H_minus_mass'][n:,n:],o['H_minus_mass'][:n,:n].conj())
    assert np.count_nonzero(o['magnetic_commutator'])>0
    assert np.array_equal(o['electron_charge_relative_dipole_current']*3,o['relative_dipole_current'])


def test_padded_residual_rejects_boundary_eigenvector_shortcut():
    cap=3;basis,rows,den=exact_polynomial(cap,CHARGE,F(1),F(1,16))
    v=[(0,0)]*len(basis);v[basis.index((0,cap))]=(1,0)
    # A unit vector at the cutoff has nonzero actual electric/magnetic residual.
    assert exact_residual(rows,den,v,1,F(0))>1
    assert norm_square(action(rows,v),den)>1


def test_certified_response_radius_and_static_contact():
    d=ground_member('electron',F(1));r=response_records(d);p=probe_records(d)
    assert abs(r['static_Ward_residual_diagnostic'])<1e-12
    assert F(p['radius_x_cells_interval_exact'][1])<F(1,8)
    assert F(p['largest_first_Brillouin_zone_suppression_interval_exact'][1])<F(1,16)
    assert not p['order_one_core_onset_resolved']
    for row in p['records']:
        total=sum(sum(sum(spin)for spin in pointer)for pointer in row['complete_joint_position_flux_probabilities'])
        assert abs(total-.5)<2e-14
        assert abs(row['final_mean_source_energy_with_failure_electron_units']-.125-p['original_preparation_mean_work_electron_units']-row['original_input_probe_mean_work_electron_units'])<2e-13


def test_complete_probe_matches_literal_kraus_energy_balance():
    o=polynomial_member(2,CHARGE,F(1));n=len(o['basis'])//2;h=o['H_minus_mass'][:n,:n]+np.eye(n)
    _,u=eigh(h);v=u[:,0]
    # Vacuum input, full charged block and complete preparation complement.
    H=np.zeros((n+1,n+1),complex);H[0,0]=.125;H[1:,1:]=h
    Ks=np.zeros_like(H);Ks[1:,0]=v/np.sqrt(2)
    Kf=np.eye(n+1);Kf[0,0]=1/np.sqrt(2)
    assert np.max(abs(Ks.conj().T@Ks+Kf.T@Kf-np.eye(n+1)))<1e-14
    rho=np.zeros_like(H);rho[0,0]=1
    prepared=Ks@rho@Ks.conj().T+Kf@rho@Kf.T
    x=np.array([o['plaquette_coordinates'][z[1]][0]for z in o['basis'][:n]])
    phase=np.exp(-1j*(1-np.cos(.7*x)))
    L=[]
    for factor in (1,-1j):
        for sign in (-1,1):
            k=np.zeros_like(H);k[1:,1:]=np.diag((1+sign*factor*phase)/np.sqrt(8));L.append(k)
    failure=np.zeros_like(H);failure[0,0]=1;L.append(failure)
    assert np.max(abs(sum(k.conj().T@k for k in L)-np.eye(n+1)))<1e-14
    phase_after=sum(k@prepared@k.conj().T for k in L)
    after=np.diag(np.diag(phase_after))
    prepwork=np.trace(H@(prepared-rho)).real;phasework=np.trace(H@(phase_after-prepared)).real
    positionwork=np.trace(H@(after-phase_after)).real
    assert positionwork>1e-5
    assert abs(prepwork+phasework+positionwork-np.trace(H@(after-rho)).real)<1e-14
    # Omitting the failure's energy violates the complete balance by exactly b/2.
    assert abs(np.trace(H@after).real-.5*np.trace(h@after[1:,1:]*2).real-.0625)<1e-13
