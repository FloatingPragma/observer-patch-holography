"""Literal CAR transport, cutoff admission and full-source error controls."""
from copy import deepcopy
from fractions import Fraction as F

import numpy as np
import pytest
import sympy as s
from scipy.linalg import expm

from . import model as m


def test_source_chart_isomorphism_and_actual_Wilson_spin_blocks():
 assert m.native_chart()==(0,1,4,2,5,7,10,3,6,8,11,9,12,13,14)
 for h in m.spin_blocks():
  assert h[:2,:2]==-s.eye(2)/2
  assert h[2:,2:]==s.eye(2)/2
 R=m.nontrivial_link();p=m.native_chart()
 inverse=[p.index(i) for i in range(15)]
 wrong=R.extract(inverse,inverse)
 w=m.sector.local_singlet((0,2,0,0,1))
 assert m.cubic_between_vectors(w['state'],w['state'],0,wrong)!=-s.Rational(3,4)


def test_literal_full_CAR_Haar_transport_all_axes_and_both_signs():
 for _,delta,_ in m.sector.SECTORS[:2]:
  for sign in (1,-1):
   actual=tuple(sign*x for x in delta);control=m.moment_control(actual)
   assert control['neighbor_moments_through_order_three']==[0,0,0,F(-3,4)]
   assert control['CAR_parity']==1
   for axis in control['full_Wilson_axis_controls']:
    assert axis['identity_link_cubic_moment']==axis['nontrivial_native_link_cubic_moment']==F(-3,4)
 with pytest.raises(ValueError):m.moment_control((0,0,0,3,-1))


def test_independent_dense_CAR_triple_and_hole_signs():
 # Three moving modes and one occupied spectator at each site. The
 # spectator makes the hole anchor even, as in the full thirty-mode anchor.
 n=8;dimension=1<<n;creation=[]
 for j in range(n):
  C=np.zeros((dimension,dimension))
  for mask in range(dimension):
   if not mask>>j&1:C[mask|1<<j,mask]=(-1)**((mask&((1<<j)-1)).bit_count())
  creation.append(C)
 transfer=sum((creation[4+j]@creation[j].T for j in range(3)),np.zeros((dimension,dimension)))
 for holes in (False,True):
  H=(.5 if holes else -.5)*(transfer+transfer.T)
  if holes:initial=(1<<3)|sum(1<<j for j in range(4,8));target=sum(1<<j for j in range(4))|(1<<7)
  else:initial=7;target=7<<4
  v=np.eye(dimension)[:,initial];w=np.eye(dimension)[:,target]
  assert np.dot(w,v)==0
  assert np.dot(w,H@v)==0
  assert np.dot(w,H@H@v)==0
  assert np.dot(w,H@H@H@v)==-.75
  R2=np.linalg.norm(H@H@v)*np.linalg.norm(H@H@w)
  t=.1
  amplitude=w@expm(-1j*t*H)@v
  assert abs(amplitude+1j*t**3/8)<=R2*t**4/24
  assert abs(amplitude)**2>=t**6/256


def test_exact_cutoff_third_moment_without_invariance():
 H=s.zeros(6)
 for i,j in ((0,2),(1,3),(2,3),(2,4),(3,4),(4,5)):H[i,j]=H[j,i]=1
 P=s.diag(1,1,1,1,0,0);C=P*H*P
 for j in (0,1):assert P*H[:,j]==H[:,j]
 assert (C**3)[1,0]==(H**3)[1,0]==1
 assert (C**4)[1,0]!=(H**4)[1,0]
 assert P*H!=H*P


def test_decoupled_resolvent_denominators_use_actual_native_Casimirs():
 expected=(
  -1/(4*F(189,8)**2)-1/(2*F(189,8)*F(63,4)),
  -1/(4*F(147,4)**2)-1/(2*F(147,4)*F(42)),
 )
 for (_,delta,_),value in zip(m.sector.SECTORS[:2],expected):
  for sign in (1,-1):
   d=tuple(sign*x for x in delta);w=m.sector.local_singlet(d)
   for row in m.moment_control(d)['full_Wilson_axis_controls']:
    assert row['third_order_decoupled_resolvent_transport_coefficient']==value<0
   assert m.cubic_between_vectors(w['state'],w['state'],0,m.nontrivial_link(),electric_denominators=True)/w['norm_squared']==s.Rational(value.numerator,value.denominator)


def test_exact_full_source_moments_and_controlled_error():
 for _,delta,_ in m.sector.SECTORS[:2]:
  result=m.moment_control(delta)
  assert result['bare_source_mean']==-711
  assert result['centered_variance']==F(132301,735)
  assert result['centered_diagonal_third_moment']==F(369869271,42875)
  t=result['positive_transfer_time_upper'];r=result['amplitude_Taylor_remainder_upper_coefficient']
  assert t>0 and r*t**4==t**3/16
  assert result['at_time_upper_transition_probability_lower']==t**6/256
  assert 4.94e-8<float(t)<4.95e-8
  assert not result['prepared_vector_is_original_selected_density']
  assert not result['invariant_eight_state_compression_assumed']


def test_complete_cube_moment_matrix_is_not_compressed_Hamiltonian():
 control=m.cube_control();A=s.Matrix(control['nearest_neighbor_adjacency'])
 W=s.Matrix(control['orthogonal_Walsh_numerator']);assert W.T*W==8*s.eye(8)
 d3=s.Rational(369869271,42875)
 actual=W.T*(d3*s.eye(8)-s.Rational(3,4)*A)*W/8
 assert actual.is_diagonal()
 for row in control['sector_rows']:
  values=[F(actual[i,i]) for i in range(8)]
  assert [r['exact_centered_third_energy_moment'] for r in row['Walsh_rows']]==values
  assert len(set(values))==4
  assert max(values)-min(values)==F(9,2)
 assert not control['Walsh_label_is_conserved_physical_momentum']


def test_all_six_physical_spin_components_are_transported():
 result=m.multiplet_transport_control()
 assert len(result['signed_multiplets'])==4
 for row in result['signed_multiplets']:
  assert row['source_spin_Casimir_multiplicities']=={'3/4':2,'15/4':4}
  for M in row['three_exact_transport_matrices']:
   assert s.Matrix(M)==-s.eye(6)*s.Rational(3,4)
 assert not result['spin_chart_is_Lorentz_little_group']


def test_actual_member_binding_rejects_forged_scope_and_other_capacity():
 member=m.contract.member(0);bound=m.bind_member(member);m.verify_bound(member,bound)
 assert bound['result']['physical_mass_table_GeV']==[]
 assert not bound['result']['isolated_dispersion_band_proved']
 bad=deepcopy(bound);bad['result']['isolated_dispersion_band_proved']=True
 bad['result_sha256']=m.contract.digest(bad['result'])
 with pytest.raises(ValueError):m.verify_bound(member,bad)
 changed=deepcopy(member);changed['source']['electric_coefficient']='106'
 with pytest.raises(ValueError):m.bind_member(changed)
 with pytest.raises(ValueError):m.bind_member(m.contract.member(1))
