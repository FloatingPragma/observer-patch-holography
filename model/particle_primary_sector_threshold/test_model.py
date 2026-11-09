"""Physical gauge kernels, conserved-sector bounds and complete memory control."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest
import sympy as s

from . import model as m


def test_actual_native_chart_and_gauss_singlets_both_signs():
 assert [(q,len(indices)) for q,indices in m.native_blocks()]==[(-4,3),(-3,2),(1,6),(2,3),(6,1)]
 for _,delta,_ in m.SECTORS:
  for sign in (1,-1):
   d=tuple(sign*x for x in delta);w=m.local_singlet(d)
   assert w['Gauss_generators_checked']==12
   assert w['CAR_parity']==sum(delta)%2
   assert sum(q*x for (q,_),x in zip(m.native_blocks(),d))==0
   assert w['defect_count']==sum(map(abs,d))
   assert m._norm(w['state'])==w['norm_squared']
 with pytest.raises(ValueError):m.local_singlet((0,0,0,0,1))
 # Zero central charge alone is insufficient: a lone color basis mask
 # in the nine-dimensional candidate space is moved by a color current.
 w=m.local_singlet((1,0,0,2,0));mask=next(iter(w['state']))
 moved=[]
 for g in m.native.native_currents():
  terms=[(15*r+i,15*r+j,v) for r in range(4) for (i,j),v in g.todok().items() if v]
  moved.append(m._act({mask:s.S.One},terms))
 assert any(moved)


def test_all_occupation_defect_bound_and_parity():
 for d,_,_ in m.refinement.native_control():
  for particles,holes in product(range(16*d+1),repeat=2):
   delta=particles-holes;D=particles+holes
   assert D>=abs(delta)
   assert (D-delta)%2==0


def test_exact_first_action_uses_all_actual_axes_and_native_rows():
 for _,delta,_ in m.SECTORS:
  plus=m.hopping_packets(delta);minus=m.hopping_packets(tuple(-x for x in delta))
  assert plus==minus
  assert sum(r['unchanged_defect_packet_norm_squared']+r['two_added_defect_packet_norm_squared'] for r in plus)==180
  for r,n in zip(plus,delta):
   assert r['unchanged_defect_packet_norm_squared']==F(3*abs(n),4)
   assert len(r['affected_axis_norms'])==3
   for axis in r['affected_axis_norms']:
    assert axis[0]+axis[2]==r['dimension']


def test_Ritz_value_against_independent_normalized_dense_compression():
 for _,delta,_ in m.SECTORS:
  trial=m.sector_trial(delta);packets=trial['hopping_packets'];k=len(packets);g=float(trial['loop_norm_squared'])
  H=np.zeros((k+6,k+6));v=np.ones(k+6)
  for i,row in enumerate(packets,1):
   H[0,i]=H[i,0]=np.sqrt(float(row['norm_squared']));H[i,i]=float(row['relative_H0_energy'])
   v[i]=-H[0,i]*float(trial['hopping_coefficients'][i-1])
  H[0,k+1:]=np.sqrt(g);H[k+1:,0]=np.sqrt(g)
  H[k+1:,k+1:]=np.array(trial['loop_packet_matrix'],float)
  v[k+1:]=-np.sqrt(g)*np.array(trial['loop_coefficients'],float)
  measured=v@H@v/(v@v)
  assert abs(measured-float(trial['trial_Rayleigh_shift']))<1e-12
  assert abs(v@v-float(trial['trial_norm_squared']))<1e-13
  assert np.linalg.eigvalsh(H)[0]<=measured+1e-12
  assert trial['packet_dimension']==13
  assert not trial['full_residual_or_isolated_line_computed']


def test_paired_memory_preserves_nonclassical_reference_and_needs_parity_correction():
 # S,A,B,R in that order. All |j>_S|j>_R terms have even combined parity.
 shape=(3,3,3,3);initial=np.zeros(shape,complex)
 amplitudes=np.array([1,2,1j])/np.sqrt(6)
 for j,c in enumerate(amplitudes):initial[j,0,0,j]=c
 actual=np.zeros(shape,complex);uncorrected=np.zeros_like(actual)
 parity=(0,1,1)
 for x,a,b,r in product(range(3),repeat=4):
  value=initial[x,a,b,r]
  if (a,b)==(0,0):aa,bb=1,2
  elif (a,b)==(1,2):aa,bb=0,0
  else:aa,bb=a,b
  phase=(-1)**(parity[x]*parity[aa])
  uncorrected[aa,x,bb,r]+=phase*value
  actual[aa,x,bb,r]+=phase*(-1)**parity[x]*value
 expected=np.zeros(shape,complex)
 for j,c in enumerate(amplitudes):expected[1,j,2,j]=c
 assert np.linalg.norm(actual-expected)<1e-15
 assert np.linalg.norm(uncorrected-expected)>1
 assert abs(actual[1,0,2,0]*actual[1,1,2,1].conjugate())>0
 assert m.paired_memory_control()['arbitrary_reference_basis_identity_checked']


def test_Markov_window_does_not_require_a_stationary_prepared_state():
 rng=np.random.default_rng(731)
 H=np.diag([-3.,-2.5,1.,8.]);lower=-3.
 for _ in range(40):
  v=rng.normal(size=4)+1j*rng.normal(size=4);v/=np.linalg.norm(v)
  rho=np.outer(v,v.conj());mean=np.trace(rho@H).real;high=2*mean-lower
  probability=sum(abs(v[i])**2 for i,E in enumerate(np.diag(H)) if E<=high)
  assert probability>=.5-1e-14
  assert np.linalg.norm(rho@H-H@rho)>1e-3


@pytest.fixture(scope='module')
def bound():
 member=m.contract.member(0)
 return member,m.bind_member(member)


def test_full_sector_support_has_a_positive_lower_edge_and_complete_readout(bound):
 member,value=bound;rows=value['result']['sector_results']
 for row,minimum in zip(rows,[F(3,2),F(1,2),F(7,2)]):
  plus,minus=row['signed_sectors']
  assert plus['trial']['sector_ground_energy_upper']==minus['trial']['sector_ground_energy_upper']
  for signed in row['signed_sectors']:
   lo,hi=map(F,signed['sector_bottom_excitation_interval'])
   assert lo>minimum and hi>lo
   p=signed['prepared_spectroscopy']
   assert F(p['reported_excitation_window'][0])>0
   assert F(p['reported_window_original_weight_lower'])>=F(255,512)
   assert p['additional_memory_qubit_upper']==1392
   assert p['original_state_and_reference_correlations_retained_in_memory']
   assert p['compensating_opposite_occupation_sector_retained_in_second_memory']
   assert not p['autonomous_controller_or_efficient_compilation_proved']
  assert F(row['selected_original_probability_of_either_sign_upper'])<F(1,16)/minimum
 assert value['result']['physical_mass_table_GeV']==[]
 assert not value['result']['charge_conjugation_spectral_equality_assumed']
 m.verify_bound(member,value)


def test_source_and_rehashed_threshold_forgery_rejected(bound):
 member,value=bound;bad=deepcopy(value)
 bad['result']['sector_results'][0]['signed_sectors'][0]['physical_sector_energy_lower']='0'
 bad['result_sha256']=m.contract.digest(bad['result'])
 with pytest.raises(ValueError):m.verify_bound(member,bad)
 with pytest.raises(ValueError):m.bind_member(m.contract.member(1))
