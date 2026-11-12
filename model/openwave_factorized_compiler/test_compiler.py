from fractions import Fraction as F
import numpy as np
import sympy as s
from scipy.linalg import expm
from .compiler import block,step,trace_bound,rounds_for_error,complete_control
from model.openwave_kinetic_phase_cycle.compiler import annihilator

def test_exact_centered_block_norm_and_phase_mean():
 c=[annihilator(4,i)for i in range(4)];u=s.Matrix([s.Rational(3,5),s.Rational(4,5)]);w=u*u.T
 V=sum((w[i,j]*c[i].H*c[2+j].H for i in range(2)for j in range(2)),s.zeros(16));V+=V.H
 assert V**3==V
 N=sum((a.H*a for a in c),s.zeros(16));assert set((N-2*s.eye(16)).diagonal())=={-2,-1,0,1,2}

def test_global_pulse_factorization_is_same_operator():
 A,_=block(.2,.3);B,_=block(.7,-.8);h=np.kron(A[2],np.eye(16))+np.kron(np.eye(16),B[2])
 assert np.linalg.norm(expm(-.031j*h)-np.kron(expm(-.031j*A[2]),expm(-.031j*B[2])),2)<1e-12

def test_volume_linear_bound_and_integer_threshold():
 n=rounds_for_error(F(1,10**8),sites=64)
 assert trace_bound(n,sites=64)<=F(1,10**8)<trace_bound(n-1,sites=64)
 assert trace_bound(n,sites=64)==64*trace_bound(n,sites=1)

def test_complete_reference_and_invalid_coupling_control():
 r=complete_control();assert sum(map(sum,r['all_original_ancillary_records']))>1-1e-12
 assert r['coupled_block_mutation_factorization_defect']>1e-4

def test_actual_nonaxial_momentum_factorization():
 from .compiler import generic_momentum_factorization
 assert max(r['maximum_intertwiner_error']for r in generic_momentum_factorization())<2e-13

def test_state_energy_bound_with_arbitrary_reference():
 from .compiler import state_energy_error
 rng=np.random.default_rng(714)
 for _ in range(25):
  U,_=np.linalg.qr(rng.normal(size=(5,5))+1j*rng.normal(size=(5,5)));e=np.arange(5.)
  A=(U*e)@U.conj().T;v=rng.normal(size=15)+1j*rng.normal(size=15);v/=np.linalg.norm(v)
  G=rng.normal(size=(15,15))+1j*rng.normal(size=(15,15));G=(G+G.conj().T)/2;T=expm(-.01j*G)
  actual=T@v;full=np.kron(A,np.eye(3));mean=float(np.vdot(v,full@v).real);difference=abs(float(np.vdot(actual,full@actual).real)-mean)
  delta=float(np.linalg.norm(T-np.eye(15),2));bound=float(state_energy_error(F(str(delta+1e-14)),4,F(str(mean+1e-14))))
  assert difference<=bound
