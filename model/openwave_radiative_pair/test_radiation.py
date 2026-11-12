import numpy as np
from scipy import sparse as sp
from scipy.sparse.linalg import expm_multiply
from .radiation import apply,lift,construction,run_probe
from .exact import certificate

def test_exact_current_coefficient_and_positive_full_H_bound():
 c=certificate();assert c['original_two_photon_radiation_probability_lower_float']>4e-38
 assert c['exact_opposite_momentum_two_photon_coefficients'][0]['H_squared_amplitude_divided_by_eta_squared']=='0'

def test_independent_full_Fock_lift_and_charge_sector():
 from itertools import combinations
 rng=np.random.default_rng(15);x=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4));x+=x.conj().T
 c=[]
 for j in range(4):
  a=np.zeros((16,16))
  for mask in range(16):
   if mask>>j&1:a[mask^(1<<j),mask]=(-1)**((mask&((1<<j)-1)).bit_count())
  c.append(a)
 H=sum(x[i,j]*c[i].T@c[j]for i in range(4)for j in range(4))
 masks=tuple(sum(1<<i for i in comb)for comb in combinations(range(4),2))
 assert np.max(abs(H[np.ix_(masks,masks)]-lift(x,masks).toarray()))<1e-14


def test_native_momentum_and_zero_coupling_controls():
 c=construction();H=c['H'];q=H.tocoo()
 assert np.max(abs((H-H.conj().T).data),initial=0)==0
 assert np.all(c['total_momentum'][q.row]==c['total_momentum'][q.col])
 r=run_probe(1,0,(2,))[0];assert r['two_photon_no_fermion_probability']==0


def test_actual_second_order_matches_exact_formula():
 c=construction();v=c['H']@(c['H']@c['initial']);out=v.reshape(924,16)[c['ground_index']]
 expected=np.zeros(16,complex);expected[c['photon_basis'].index((1,0,0,1))]=2j*np.sqrt(2/5)*.05**2
 expected[c['photon_basis'].index((0,1,1,0))]=-2j*np.sqrt(2/5)*.05**2
 assert np.max(abs(out-expected))<2e-15


def test_complete_records_energy_and_vacuum_control():
 pair=run_probe(1,.05,(5,))[0];vac=run_probe(1,.05,(5,),initial_kind='vacuum')[0]
 assert pair['two_photon_no_fermion_probability']>1e-4
 assert abs(pair['full_energy_before_readout']-2)<1e-9
 assert abs(vac['full_energy_before_readout'])<1e-9
 assert abs(sum(map(sum,pair['complete_source_ground_or_complement_and_photon_probabilities']))-1)<1e-9
 assert abs(pair['free_source_energy']+pair['free_photon_energy']+pair['interaction_energy']-2)<1e-9
