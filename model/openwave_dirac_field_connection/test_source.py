from fractions import Fraction as F
import numpy as np
from .source import block,spectral_intervals,records,run,HOP,ELECTRIC


def test_connection_identity_by_independent_differentiation():
    import sympy as s
    theta=s.symbols('theta',real=True);c=F(4,5);ss=F(3,5);q=3
    u=s.Matrix([[s.Rational(c),-s.Rational(ss)*s.exp(-s.I*q*theta)],[s.Rational(ss)*s.exp(s.I*q*theta),s.Rational(c)]])
    a=s.simplify(-s.I*u.conjugate().T*u.diff(theta))
    expected=q*s.Matrix([[s.Rational(ss**2),s.Rational(c*ss)*s.exp(-s.I*q*theta)],[s.Rational(c*ss)*s.exp(s.I*q*theta),-s.Rational(ss**2)]])
    assert s.simplify(a-expected)==s.zeros(2)
    assert s.simplify(a*a-q*q*s.Rational(ss**2)*s.eye(2))==s.zeros(2)


def test_all_block_eigenvalues_have_rational_enclosures():
    import mpmath as mp
    with mp.workdps(90):
        m=lambda x:mp.mpf(F(x).numerator)/F(x).denominator
        for mass in (F(1),F(207),F(3477)):
            for n in range(-5,6):
                for q in (-3,3):
                    lo,hi=spectral_intervals(n,q,mass)
                    mm,aa,t=m(mass),m(ELECTRIC),m(HOP)
                    eps=mp.sqrt(mm*mm+t*t);cc=aa*(n*n+(n+q)**2)/2;dd=mm+aa*(n*n-(n+q)**2)/2
                    rad=mp.sqrt(dd*dd+t*t)
                    assert m(lo[0])<=eps+cc-rad<=m(lo[1])
                    assert m(hi[0])<=eps+cc+rad<=m(hi[1])


def test_missing_connection_changes_dynamics_and_spectrum():
    for n in (-1,0,1):
        bare,u,h,naive=block(n,3,F(1))
        np.testing.assert_allclose(np.linalg.eigvalsh(bare),np.linalg.eigvalsh(h),atol=1e-14)
        assert np.linalg.norm(h-naive)>1e-4
    r=records(0,3,F(1),.5)
    assert r['complete_band_records'][1]['probability']>1e-5
    assert abs(r['terminal_readout_work'])>1e-5
    assert abs(r['after_readout_energy']-r['initial_source_energy']-r['terminal_readout_work'])<1e-12


def test_full_reference_and_scope():
    r=run();assert len(r['cases'])==6
    assert r['all_flux_blocks_retained_analytically']
    assert not r['source']['native_Fock_source_equivalence_proved']
    assert not r['scope']['new_physical_mass_prediction']
