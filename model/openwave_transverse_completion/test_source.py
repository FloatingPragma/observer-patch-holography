from fractions import Fraction as F
import numpy as np
from scipy.linalg import eigh,expm
from .source import matrices,gauss_checks,tail_bound,energy_certificate,exact_cap_enclosure,finite_position_certificate,exact_phase,HOP,MAGNETIC

def test_all_native_charges_and_complete_positive_compressions():
    assert gauss_checks()['all_hop_columns_tested']==504
    for q in (-3,-2,-1,0,1,2,3):
        d=matrices(5,q)
        assert np.linalg.eigvalsh(d['K'].astype(float)).min()>-1e-12
        assert np.linalg.eigvalsh(d['V'].astype(float)).min()>-1e-12
        assert d['A'].shape==(44,44)

def test_exact_tail_controls_independent_larger_carrier_records():
    small=matrices(6,-3);large=matrices(15,-3);t=F(3,2)
    def state(d):
        e,U=eigh(d['A'].astype(float)/32);v=np.zeros(len(e));v[d['index'][0,0]]=1
        return U@(np.exp(-1j*e*float(t))*(U.T@v))
    u=state(small);v=state(large);embedded=np.zeros_like(v)
    for i,z in enumerate(small['basis']):embedded[large['index'][z]]=u[i]
    assert np.linalg.norm(embedded-v)<float(tail_bound(6,-3,t))
    assert abs(u[small['index'][1,0]])>1e-3

def test_finite_numerical_error_is_separate_from_the_exact_tail():
    from fractions import Fraction as F
    import mpmath as mp
    mp.mp.dps=90
    for x in (F(0),F(1234567,1024),F(-7,3)):
        re,im,error=exact_phase(x)
        convert=lambda a:mp.mpf(a.numerator)/a.denominator
        z=mp.exp(-1j*convert(x))
        assert abs(z-(convert(re)+1j*convert(im)))<convert(error)
    d=matrices(2,0);initial=np.zeros(len(d['A']));index=d['index'][0,0];initial[index]=1
    budgets=[]
    for driver in ('ev','evd','evr','evx'):
        cert,e,U=exact_cap_enclosure(d['A'],eigh(d['A'].astype(float)/32,driver=driver))
        psi=U@(np.exp(-2j*e)*(U.T@initial))
        numerical=finite_position_certificate(d['A'],e,U,psi,index,2,cert)
        budgets.append(numerical)
        error=float(F(numerical['conditional_vector_error_upper_exact']))
        assert np.linalg.norm(psi-expm(-2j*d['A']/32)@initial)<error<.001
        assert numerical['complete_position_flux_table_numerically_enclosed'] is True
    assert all(x==budgets[0]for x in budgets)
    import pytest
    wrong=psi.copy();wrong[0]+=.001
    with pytest.raises(AssertionError):finite_position_certificate(d['A'],e,U,wrong,index,2,cert)

def test_actual_infinite_energy_not_unshifted_anchor():
    c=energy_certificate(12);low,high=map(F,c['charged_rest_minus_full_field_vacuum_electron_units_interval_exact'])
    assert 1<low<high<F(6,5)
    d=matrices(3,-3);assert not np.array_equal(d['A'],np.diag(np.diag(d['A'])))
    # Suppressing hopping preserves position exactly; suppressing both off-diagonal
    # mechanisms preserves the complete initially declared flux state.
    no_hop=d['A']-d['K'];x=np.repeat(np.arange(4),7)
    assert np.count_nonzero(no_hop*(x[:,None]-x[None,:]))==0

def test_complete_preparation_and_charged_absorption_energy_account():
    # Vacuum, matter, and retained charged-detector sectors. Every flux value
    # and every success/failure outcome of this complete compression remains.
    from scipy.linalg import block_diag
    cap=3;mass=F(5,3);d=matrices(cap,-3);M=len(d['A']);N=2*cap+1
    fluxes=np.arange(-cap,cap+1)
    h0=(np.diag(4*fluxes**2+4)+np.diag(np.full(N-1,-2),1)+np.diag(np.full(N-1,-2),-1))/32
    h1=d['A']/32+float(mass)*np.eye(M)
    detector=(d['A']-d['K'])/32
    H=block_diag(h0,h1,detector);initial=np.zeros(N+2*M);initial[cap]=1
    success=np.zeros_like(H);success[N+d['index'][0,0],cap]=1/np.sqrt(2)
    failure=np.eye(len(H));failure[cap,cap]=1/np.sqrt(2)
    assert np.linalg.norm(success.T@success+failure.T@failure-np.eye(len(H)))<4e-16
    u=success@initial;f=failure@initial
    energy=lambda v:float(np.vdot(v,H@v).real)
    before=energy(initial);prepared=energy(u)+energy(f)
    assert abs(before-float(MAGNETIC))<1e-15
    work=prepared-before
    assert abs(work-float((mass+2*HOP)/2))<1e-15
    U=expm(-1.25j*H);u=U@u;f=U@f
    assert abs(energy(f)-float(MAGNETIC)/2)<1e-14
    p=abs(u[N:N+M])**2
    final=energy(f)+float(np.dot(p,np.diag(detector)))
    readout=final-energy(u)-energy(f)
    electric=float(np.dot(p,d['electric'])/32)
    assert abs(final-before-electric)<2e-14
    assert abs(work+readout-electric)<2e-14
    conditional=np.sqrt(2)*u[N:N+M]
    K=float(np.vdot(conditional,d['K']@conditional).real/32)
    V=float(np.vdot(conditional,d['V']@conditional).real/32)
    assert abs(readout-.5*(-float(mass)-K+float(MAGNETIC)-V))<2e-14
