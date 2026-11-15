import numpy as np
import sympy as s
from .spectrum import native_tensors,exact_symbol,countercontrols,pair_records
from model.lepton_mass_compiler import compiler as native

def test_native_full_tensor_removes_only_charged_zeros():
 d=native_tensors();assert d['P'].rank()==14 and d['Q']*(s.eye(15)-d['P'])==s.zeros(15)
 assert exact_symbol()['exact_three_axis_Clifford_mass_anticommutators']==3

def test_full_spin_pair_Fock_charge_floor():
 c=[native._annihilator(4,i)for i in range(4)];N=[a.conj().T@a for a in c]
 for charge in(-1,2/3,-1/3):
  V=c[0].conj().T@c[3].conj().T-c[1].conj().T@c[2].conj().T;H=V+V.conj().T
  Q=charge*(N[0]+N[1]-N[2]-N[3]);assert np.linalg.norm(H@Q-Q@H)<1e-13
  ground=np.linalg.eigvalsh(H)[0];absQ=np.diag(abs(np.diag(Q)))
  assert np.linalg.eigvalsh(H-ground*np.eye(16)-absQ)[0]>-1e-12

def test_missing_quark_channel_fails_global_charged_floor():
 r=countercontrols();assert r['charged_lepton_only_charged_zero_internal_count']==12
 assert all(len(x['charged_zero_internal_indices'])==6 for x in r['removed_channel_controls'])

def test_complete_source_creation_energy_charge_records():
 r=pair_records();assert len(r['pair_orientation_records'])==21
 for row in r['pair_orientation_records']:
  assert abs(sum(row['all_four_occupation_probabilities'])+row['failure_probability']-1)<1e-13
