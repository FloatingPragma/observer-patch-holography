from fractions import Fraction as F
from pathlib import Path
import importlib.util
import math
import numpy as np
import pytest
from scipy.linalg import expm

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('ground_refinement_under_test',HERE/'refinement.py')
r=importlib.util.module_from_spec(spec);spec.loader.exec_module(r)


@pytest.fixture(scope='module')
def certificate():return r.scientific()


def test_uniform_ground_full_gap_and_actual_records(certificate):
    b={k:F(v) for k,v in certificate['uniform_ground_bounds_exact'].items()}
    assert b['strict_Schur_margin_lower']>F(7,10**6)
    assert b['uniform_internal_code_gap_lower_over_m']>F(1328,1000)
    assert b['uniform_full_physical_gap_lower']>F(25,10**11)
    assert b['uniform_complete_filter_tree_success_lower']>F(9993,10000)
    assert b['uniform_physical_preparation_and_T111_joint_probability_lower']>F(45,10**8)
    assert b['ideal_ground_two_field_record_error_at_K1024_upper']<F(31,10**10)
    assert b['fixed_pointer_complete_prep_and_two_field_tests_error_upper']<F(1,10**8)


def test_actual_vacuum_moments_all_caps_and_boundary_countercontrol():
    # One collective oscillator is exact for each fixed edge under total-cap
    # unitary mode rotations; this test does not substitute product cutoffs.
    import sympy as s
    for K in range(3,9):
        a=s.zeros(K+1)
        for n in range(1,K+1):a[n-1,n]=s.sqrt(n)
        q=a+a.T
        assert [(q**j)[0,0] for j in (2,4,6)]==[1,3,15]
        comm=a*a.T-a.T*a
        assert comm[K,K]==-K
        assert comm!=s.eye(K+1)
    a=s.zeros(4)
    for n in range(1,4):a[n-1,n]=s.sqrt(n)
    assert ((a+a.T)**8)[0,0]!=105


def test_fixed_cayley_ground_chord_without_excited_norm_bound():
    B=F(7);gap=F(1,100)
    # Exact squared chord comparison, including raw energies far outside B.
    for ground in (-B,-B/2,F(0),B/2,B):
        for d in (gap,2*gap,F(1),2*B,10**6*B,10**30*B):
            energy=ground+d
            chord_squared=4*B*B*d*d/((B*B+ground*ground)*(B*B+energy*energy))
            assert chord_squared>=gap*gap/(4*B*B)


def test_range_reduced_outward_exponential():
    import mpmath as mp
    mp.mp.dps=110
    for x in (F(0),F(1,1000),F(1),F(16),F(650)):
        lo,hi=r.exp_interval(x)
        v=mp.exp(mp.mpf(x.numerator)/x.denominator)
        assert mp.mpf(lo.numerator)/lo.denominator<=v
        assert v<=mp.mpf(hi.numerator)/hi.denominator
        assert (hi-lo)/max(F(1),lo)<F(1,10**75)


def test_low_sector_literal_exponential_bound_and_no_compression_identity():
    radius=.025
    for K in (7,8,11):
        L=K//2;full=45
        q=np.diag(np.sqrt(np.arange(1,full)),1)
        q=q+q.T
        infinity=expm(1j*radius*q)
        finite=expm(1j*radius*q[:K+1,:K+1])
        embedded=np.zeros((full,L+1),complex);embedded[:K+1]=finite[:,:L+1]
        error=np.linalg.norm(embedded-infinity[:,:L+1],2)
        d=K+1-L;rho=6*radius*math.sqrt(K+1)/d
        bound=2*rho**d/(1-rho)
        assert error<=bound+1e-14
        assert bound<(K+1)**-2
    assert np.linalg.norm(finite-infinity[:K+1,:K+1])>1e-5


def test_complete_filter_tree_ground_fidelity_with_high_energy():
    # Exact three-energy spectral law and all M outcomes. High-energy alias
    # states do not invalidate the ground chord, and all outcomes sum to 1.
    B=7.;M=257;energies=np.array([-3.,-2.,1e8]);p=np.array([.995,.003,.002])
    C=(1-1j*energies/B)/(1+1j*energies/B)
    W=C[:,None]*C.conj()[None,:]
    outcomes=np.fft.fft(W[...,None]**np.arange(M),axis=-1)/M
    weights=p[:,None,None]*p[None,:,None]*np.abs(outcomes)**2
    assert abs(weights.sum()-1)<1e-12
    zero=weights[:,:,0];conditional_ground=zero[0].sum()/zero.sum()
    epsilon=4*B/(M*1.)
    odds=(1-p[0])/p[0]
    lower=1/(1+odds**2+2*epsilon**2*odds)
    assert conditional_ground>=lower-1e-13
    assert zero.sum()>=p[0]**2-1e-13


def test_guards_and_resource_ancestry(certificate):
    c=certificate
    assert c['fixed_complete_energy_clock']['pointer_M']==315179904752653412
    assert c['fixed_complete_energy_clock']['same_pointer_alphabet_for_all_K_and_mathematical_comparison']
    assert not c['fixed_complete_energy_clock']['imaginary_heat_time_is_physical_clock']
    assert c['constructive_complete_preparation_history']['K']==2**96
    assert c['constructive_ideal_ground_history']['K']==1024
    assert not c['boundaries']['uniform_spatial_or_relativistic_cone']


def test_raw_ground_bound_rounds_outward(certificate):
    parent,source=r.admit()
    mu=F(parent['bounds_exact']['free_vacuum_shift_over_m_upper'])
    m=F(source['clock_and_capacity']['source_selected_excitation_interval'][1])
    raw=32*(480+F(360,2**34))+m*(3+r.bounds.root(6)+mu)
    recorded=F(certificate['uniform_ground_bounds_exact']['uniform_raw_ground_absolute_bound'])
    assert recorded>=raw
    assert recorded-raw<=F(1,10**70)
