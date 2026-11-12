from fractions import Fraction as F
import numpy as np
from .clocks import native_clock_factorization,joint_instrument,outputs,weak_probe_certificate,exact_npt_witness

def test_native_factorization_and_complete_joint_instrument():
 d=native_clock_factorization();assert d['exact_quadratic_reconstruction']
 assert joint_instrument()['complete_native_outcomes']==32
 rng=np.random.default_rng(71);v=rng.normal(size=(4,7))+1j*rng.normal(size=(4,7));v*=2/np.linalg.norm(v)
 p,s,z=outputs(v);assert abs(p.sum()-1)<1e-14 and abs(s-1)<1e-14
 assert len(z)==2

def test_exact_entanglement_witness_and_probe_bound():
 # A Bell pair with a trivial environment checks the partial transpose
 # convention and exact-fraction witness independently of time evolution.
 v=np.zeros((4,1),complex);v[0,0]=v[3,0]=np.sqrt(2)
 r=exact_npt_witness(v,F(1,10**12));assert r['source_expectation_interval'][1]<-.49
 d=weak_probe_certificate();assert d['all_environment_states_and_passive_quantum_references']
 assert F(d['finite_probe_bounds'][-1]['two_clocks_same_environment_normalized_phase_disagreement_upper'])<F(35,10**8)
