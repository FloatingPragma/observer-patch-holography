"""Independent complete-channel and finite-observation controls."""
from fractions import Fraction as F

import numpy as np
import pytest
from scipy.linalg import expm
from scipy.sparse import csr_matrix

from . import recorded_response as r


def control():
    rng=np.random.default_rng(53)
    def herm(n):
        a=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
        return (a+a.conj().T)/3
    a,_=np.linalg.qr(rng.normal(size=(4,4))+1j*rng.normal(size=(4,4)))
    b,_=np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))
    return herm(3),herm(4),a[:,:2]@b[:,:2].conj().T,a[:,1:3]@b[:,1:3].conj().T


@pytest.mark.parametrize('sign',[-1,1])
@pytest.mark.parametrize('quadrature',['X','Y'])
def test_complete_kraus_channel_against_full_exponentials(sign,quadrature):
    h0,h1,b,c=control();h=np.zeros((7,7),complex);h[:3,:3]=h0;h[3:,3:]=h1
    def w(b,s):
        g=np.zeros((7,7),complex);g[3:,:3]=b;g[:3,3:]=b.conj().T
        return expm(-.5j*s*np.pi*g)
    u=expm(-.37j*h);left=w(b,1)@u;right=u@w(c,sign)
    phase=1 if quadrature=='X' else -1j
    expected=((left+phase*right)/2,(left-phase*right)/2)
    actual=r.ramsey_kraus_actions(h0,h1,b,c,np.eye(7),.37,sign,quadrature)
    for a,e in zip(actual,expected):np.testing.assert_allclose(a,e,atol=2e-14)
    np.testing.assert_allclose(sum(a.conj().T@a for a in actual),np.eye(7),atol=2e-14)
    # Columns retain an untouched entangled reference, including coherences.
    rho=np.eye(7)/7
    np.testing.assert_allclose(sum(a@rho@a.conj().T for a in actual).trace(),1,atol=2e-14)


def test_record_probabilities_match_postmeasurement_norms():
    h0,h1,b,c=control();times=[0.,.37,1.2]
    law=r.record_law(h0,h1,[b,c],times)
    state=np.r_[law['ground'],np.zeros(4)]
    for ti,t in enumerate(times):
        for si,sign in enumerate((1,-1)):
            for qi,q in enumerate(('X','Y')):
                out=r.ramsey_kraus_actions(h0,h1,b,c,state,t,sign,q)
                np.testing.assert_allclose([np.vdot(x,x).real for x in out],
                                           law['probabilities'][ti,si,0,1,qi],atol=2e-14)


def test_invalid_physical_inputs_and_forged_spectrum_rejected():
    h0,h1,b,c=control()
    with pytest.raises(ValueError,match='Hermitian'):
        r.record_law(h0+1j*np.eye(3),h1,[b],[0.])
    with pytest.raises(ValueError,match='partial isometry'):
        r.ramsey_kraus_actions(h0,h1,.5*b,c,np.eye(7),.2)
    e,u=np.linalg.eigh(h1)
    with pytest.raises(ValueError,match='eigensystem'):
        r.record_law(h0,h1,[b],[0.],spectral1=(e+.1,u))
    with pytest.raises(ValueError):r.record_law(h0,h1,[b],[float('nan')])
    with pytest.raises(ValueError):r.ramsey_kraus_actions(h0,h1,b,c,np.eye(7),.2,True)
    with pytest.raises(ValueError):r.record_law(h0,h1,[b],np.array([1+2j]))
    with pytest.raises(ValueError):r.spectral_matrix(h0,np.eye(3),float('nan'))


def test_duplicate_sparse_entries_cannot_inflate_hermiticity_tolerance():
    h=csr_matrix(([1e12,-1e12,1.],[0,0,1],[0,3,3]),shape=(2,2))
    with pytest.raises(ValueError,match='Hermitian'):r.hermitian(h)
    assert len(h.data)==3  # Validation leaves the caller's representation intact.


def test_shot_counts_and_complete_packed_records():
    law=r.record_law(np.array([[0.]]),np.array([[.2]]),[np.ones((1,1))],[0.,.5])
    for bad in (True,0,-1,1.5):
        with pytest.raises(ValueError):r.sample_diagonal_records(law,bad)
    sample=r.sample_diagonal_records(law,37)
    bits=np.unpackbits(sample['packed_binary_records'],axis=-1,count=37)
    np.testing.assert_array_equal(bits.sum(axis=-1),sample['plus_counts'])


def test_exact_finite_shot_information_bound():
    times=tuple(F(i,10) for i in range(129))
    bound=r.spectral_shift_tv_squared(times,4096,24,F(3,10**8))
    assert bound==F(3,10**8)**2*4096*24*F(128*129*257,600)/4
    assert bound<F(1,2500)**2
    # This channel has a nonzero line and one retained coherent detector.
    # Its exact overlap loss is sin²(delta*t/2), below delta²*t²/4.
    for delta,t in ((.1,.2),(.4,.8),(1.,2.)):
        loss=np.sin(delta*t/2)**2
        assert loss <= delta**2*t**2/4+1e-15

