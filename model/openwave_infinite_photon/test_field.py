import numpy as np
import sympy as s
from .field import analytic_trial,padded_operator,certify
from model.openwave_quantum_rotational.enclosure import exact_sparse

def test_exact_isolated_infinite_spin_band_trial():
 d=analytic_trial();assert s.Rational(d['vacuum_spin_trial_shift_exact'])<s.Rational(1,16)
 assert d['exact_dressed_Casimir']=='3/4'

def test_sparse_real_embedding_matches_exact_source_and_retains_boundary():
 M,R,number,basis,maxrow=padded_operator(2);B,_,_,_=exact_sparse(2)
 n=len(number);exact=s.zeros(n)
 for ij,z in B.items():exact[ij]=z
 a=np.asarray(exact.evalf(18),complex);real=np.block([[a.real,-a.imag],[a.imag,a.real]])
 assert np.max(np.abs(M.toarray()-real)-R.toarray())<=0
 assert maxrow<256
 # Finite B squared omits a positive boundary contribution, so its square root is not a compression.
 big,_,_,large,_=padded_operator(3);idx=[i*len(large)+large.index(p)for i in range(20)for p in basis];idx+= [i+20*len(large)for i in idx]
 difference=(big@big)[np.ix_(idx,idx)].toarray()-(M@M).toarray()
 assert np.linalg.norm(difference)>1e-5
 assert np.linalg.eigvalsh(difference).min()>-1e-12

def test_cap_three_certifies_infinite_energy_and_tiny_tail():
 d=certify(3);lo,hi=d['infinite_photon_ground_energy_interval']
 assert lo<hi and hi-lo<3e-12
 assert d['infinite_ground_photon_tail_above_trial_cap_upper']<1e-14
 assert d['padded_exact_moment_cap']==17
