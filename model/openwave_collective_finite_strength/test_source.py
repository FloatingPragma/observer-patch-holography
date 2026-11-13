from fractions import Fraction as F
from math import comb
import numpy as np
import pytest
from scipy.linalg import eigh_tridiagonal
from .bounds import moments, complete_gap
from .source import complete_matrix, uniform_certificate, quench_certificate, coarse_record_certificate


def test_all_spin_blocks_and_full_bit_spectrum():
    for n in range(1,7):
        H,q=complete_matrix(n,.7,1.3,.4)
        actual=np.linalg.eigvalsh(H); reconstructed=[]
        for j2 in range(n,-1,-2):
            m=np.arange(-j2,j2+1,2)/2
            diagonal=n/2+.7*np.exp(.25*m/n)+1.3*np.exp(.1*m/n)
            off=-.5*np.sqrt(np.arange(1,j2+1)*np.arange(j2,0,-1))
            values=eigh_tridiagonal(diagonal,off,eigvals_only=True) if j2 else diagonal
            k=(n-j2)//2
            multiplicity=comb(n,k)-(comb(n,k-1) if k else 0)
            reconstructed.extend(np.repeat(values,multiplicity))
        assert len(reconstructed)==2**n
        assert np.max(abs(np.sort(reconstructed)-actual))<2e-12


def test_rational_gap_encloses_complete_matrix():
    for n in (1,3,5):
        cert=complete_gap(n,[(F(1),(F(2),F(2))),(F(1,3),(F(1),F(1)))])
        values=np.linalg.eigvalsh(complete_matrix(n,2.,1.,1/3)[0])
        assert float(cert['gap'][0])-1e-12<=values[1]-values[0]<=float(cert['gap'][1])+1e-12
        assert cert['gap'][0]>0


def test_product_moments_against_all_outcomes():
    for n in range(1,8):
        for m in (F(-3,5),F(0),F(1,3)):
            p=[F(comb(n,k))*((1+m)/2)**k*((1-m)/2)**(n-k) for k in range(n+1)]
            q=[F(2*k-n,8*n) for k in range(n+1)]
            assert tuple(sum(a*b**power for a,b in zip(p,q,strict=True)) for power in (2,4))==moments(n,F(1,8),m*m)


def test_finite_energy_mutual_bound_and_time_records():
    # A moderate carrier makes the bound testable without using an enormous
    # physical calibration; the exact full Hamiltonian supplies the control.
    n=8; es=.2; ep=.3; f=.7
    cert=uniform_certificate(n,F(1),F(1,5))
    ground=lambda a,b:np.linalg.eigvalsh(complete_matrix(n,a,b,f)[0])[0]
    mutual=ground(es,ep)-ground(es,0)-ground(0,ep)+ground(0,0)
    leading=-2*(1/8)**2*f*es*ep/n
    assert abs(mutual/leading-1)<=float(F(cert['mutual_energy_relative_error_upper_exact']))+1e-10
    H,q=complete_matrix(n,es,ep,f); e,V=np.linalg.eigh(H); initial=np.ones(2**n)/np.sqrt(2**n)
    B=es+f*ep; h=B/(8*n); omega=np.sqrt(1+4*h*h)
    for t in (1.,3.):
        state=V@(np.exp(-1j*t*e)*(V.T@initial))
        m=-2*h/omega**2*(1-np.cos(omega*t))
        product=(np.cosh(f/(8*n))+m*np.sinh(f/(8*n)))**n
        actual=float((abs(state)**2)@np.exp(f*q))
        bound=quench_certificate(n,F(1),F(t))
        normalizer=np.cosh(f/(8*n))**n
        assert abs(actual-product)/normalizer<=f*float(F(bound['normalized_lapse_error_per_profile_upper_exact']))+2e-12


def test_exact_huge_carrier_bounds_are_not_rounded_to_zero():
    n=6*10**42
    cert=uniform_certificate(n,F(7),F(4,5))
    assert 0<F(cert['mutual_energy_relative_error_upper_exact'])<F(1,10**18)
    assert F(cert['full_geometry_sector_gap_lower_exact'])>F(999,1000)
    sampled=coarse_record_certificate(n,F(7),F(4,5))
    assert 0<F(sampled['relative_error_upper_exact'])<F(1,10**11)


def test_invalid_energy_domain_is_not_accepted():
    with pytest.raises(ValueError): uniform_certificate(1,F(10**6),F(1))
    with pytest.raises(ValueError): complete_gap(0,[])
