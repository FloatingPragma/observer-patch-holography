from pathlib import Path
from fractions import Fraction as F
import sys
import numpy as np
import sympy as s
from scipy.linalg import expm

HERE=Path(__file__).resolve().parent
from . import dynamics as d


def test_analytic_heat_and_dyadic_constants():
    assert F(d.constants()['analytic_coefficient_upper'])<8
    r=s.Symbol('r',integer=True,nonnegative=True)
    for R in (0,1,4,16):
        assert s.summation((r+2)**2/2**r,(r,R+1,s.oo))==s.Rational(R*R+8*R+18,2**R)
    assert d.hopping_tail(128)<d.hopping_tail(64)<d.hopping_tail(32)
    tail=sum(F(1,s.factorial(j)) for j in range(25,70))
    assert 2*tail<d.dyson_tail(1,24)


def test_literal_one_form_heat_and_coordinate_covariance_tail():
    # Actual mixed Dirichlet/Neumann product eigenbasis; no inferred
    # periodic continuum density of states is used.
    for N in (2,3,4,8,16):
        k=np.arange(N);lam=4*np.sin(np.pi*k/(2*N))**2
        vertex=np.sqrt(2/N)*np.cos(np.pi*(np.arange(N)[:,None]+.5)*k[None,:]/N)
        vertex[:,0]=1/np.sqrt(N)
        edge=np.sqrt(2/N)*np.sin(np.pi*(np.arange(N-1)[:,None]+1)*k[None,1:]/N)
        assert np.linalg.norm(edge.T@edge-np.eye(N-1))<1e-13
        x=(N-2)//2;y=(N-1)//2
        weights=edge[x,:,None,None]**2*vertex[y,None,:,None]**2*vertex[y,None,None,:]**2
        energy=lam[1:,None,None]+lam[None,:,None]+lam[None,None,:]
        for t in (.125,1.,2.,10.,100.):
            assert float(np.sum(weights*np.exp(-2*t*energy)))<=8*t**(-1.5)
            # The actual transverse projection only decreases this
            # positive one-form tail. erfc is the exact heat integral.
            from scipy.special import erfc
            covariance_tail=float(np.sum(weights*erfc(np.sqrt(t*energy))/(2*np.sqrt(energy))))
            assert covariance_tail<=12/t


def car(count,mode):
    a=np.zeros((2**count,2**count),complex)
    for word in range(2**count):
        if word>>mode&1:
            a[word^(1<<mode),word]=(-1)**((word&((1<<mode)-1)).bit_count())
    return a


def test_full_Fock_Coulomb_phase_and_mixed_insertions():
    # Two sites, two spin modes each; full 16 states, including nonzero
    # total-charge sectors. G is the actual two-vertex pseudoinverse.
    annihilators=[car(4,j) for j in range(4)]
    numbers=[a.conj().T@a for a in annihilators]
    I=np.eye(16);rho=[numbers[0]+numbers[1]-I,numbers[2]+numbers[3]-I]
    G=np.array([[1,-1],[-1,1]],float)/4;g=.17;t=.31
    Hc=sum(g*g*G[x,y]*rho[x]@rho[y]/2 for x in range(2) for y in range(2))
    hopping=annihilators[2].conj().T@annihilators[0]
    phase=g*g*((rho[1]-rho[0])/2+I/2)
    actual=expm(1j*t*Hc)@hopping@expm(-1j*t*Hc)
    assert np.linalg.norm(actual-hopping@expm(1j*t*phase))<1e-13
    assert np.linalg.norm(actual-hopping@expm(1j*t*(phase-g*g*I/2)))>1e-4
    # A centered two-site Slater state, with each occupied orbital evenly
    # spread between sites, tests the insertion-state bound directly.
    vacuum=np.zeros(16,complex);vacuum[0]=1
    omega=(annihilators[0]+annihilators[2]).conj().T@(annihilators[1]+annihilators[3]).conj().T@vacuum/2
    assert abs(np.linalg.norm(omega)-1)<1e-14
    a=np.array([.23,-.37]);X=a[0]*rho[0]+a[1]*rho[1]
    assert np.linalg.norm(X@omega)<=2*np.linalg.norm(a)
    for word in (annihilators[0],hopping,hopping@annihilators[1]):
        k={id(annihilators[0]):1,id(hopping):2}.get(id(word),3)
        assert np.linalg.norm(X@word@omega)<=(2+k)*np.linalg.norm(a)
