"""Actual Haar/CAR moments, Stieltjes controls and source-bound record scopes."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest
import sympy as s

from . import model as m


def test_all_native_Haar_walks_all_cells_and_both_spin_rows():
 assert m.slater_hopping_cumulant('anchor')['five_Haar_trace_terms']==(72,36,48,24,15)
 for kind,cell,spin in product(('particle','hole'),range(8),range(2)):
  x=m.slater_hopping_cumulant(kind,cell,spin)
  assert x['retained_backtracking_closed_walks']==120
  assert x['vanishing_plaquette_closed_walks']==48
  assert x['one_native_dimension_variance']==12
  assert x['fourth_cumulant_per_native_dimension']==-18
 with pytest.raises(ValueError):m.slater_hopping_cumulant('particle',8)


def test_independent_tuple_exponent_CAR_second_action():
 sites,where,_,_=m.cube_data();terms=[];edges=0
 for x in sites:
  for axis in range(3):
   if x[axis]:continue
   y=list(x);y[axis]=1;y=tuple(y);a,b=where[x],where[y]
   h=m.transport.spin_blocks()[axis]
   for i,j in product(range(4),repeat=2):
    if h[i,j]:
     terms.append((4*a+i,4*b+j,edges,1,h[i,j]))
     terms.append((4*b+j,4*a+i,edges,-1,s.conjugate(h[i,j])))
   edges+=1
 def act(vector):
  result={}
  for (mask,exponents),value in vector.items():
   for i,j,edge,direction,c in terms:
    if not mask>>j&1:continue
    sign=(-1)**((mask&((1<<j)-1)).bit_count());target=mask^(1<<j)
    if target>>i&1:continue
    sign*=(-1)**((target&((1<<i)-1)).bit_count());target|=1<<i
    power=list(exponents);power[edge]+=direction;key=target,tuple(power)
    result[key]=result.get(key,0)+sign*c*value
  return {key:s.expand(value) for key,value in result.items() if value!=0}
 omega=sum(1<<(4*x+i) for x in range(8) for i in (2,3))
 for mask in (omega,omega|1,omega^(1<<2)):
  first=act({(mask,(0,)*12):s.S.One});second=act(first)
  assert s.simplify(sum(s.conjugate(x)*x for x in first.values()))==12
  assert s.simplify(sum(s.conjugate(x)*x for x in second.values()))==414
 assert all(row['fourth_moment']==414 for row in m.direct_character_CAR()['rows'])


def test_exact_mixed_fourth_moment_includes_all_native_cross_terms():
 x=m.fourth_moment()
 assert x['hopping_fourth_moment']==96930==3*180**2-18*15
 assert x['magnetic_fourth_moment']==F(41,6002500)
 assert x['twice_H0_magnetic_cubic_cross_term']==F(-188,18375)
 assert x['actual_cutoff_centered_fourth_moment_upper']==F(9738478674383,18007500)
 assert not x['actual_cutoff_fourth_moment_equality_claimed']
 assert x['transport_time_bound_improvement_factor']>56
 assert 2.77366e-6<float(x['improved_prepared_transport_time_upper'])<2.77367e-6
 assert x['normalized_complement_coupling_energy_variance_upper']>0
 # Dropping the between-native squared-variance terms gives a different
 # moment; representation blocks are not independently Haar sampled.
 wrong=sum(3*(12*d)**2-18*d for d,_,_ in m.transport.sector.refinement.native_control())
 assert wrong!=x['hopping_fourth_moment']


def test_quadratic_majorant_and_two_point_Radau_saturation():
 L=s.Rational(-2);z=s.Rational(-5);mu=s.Rational(3);tau=s.Rational(7);d=mu-L;c=mu+tau/d
 q=s.Symbol('q');f=1/(q-z)
 p=s.cancel(f+(q-L)*(q-c)**2/((L-z)*(c-z)**2*(q-z)))
 assert s.Poly(p,q).degree()==2 and s.Poly(p,q).LC()>0
 for value in (L,L+s.Rational(1,3),mu,c,c+2,c+100):assert s.simplify(p.subs(q,value)-f.subs(q,value))>=0
 weight=tau/(tau+d*d)
 exact=weight/(L-z)+(1-weight)/(c-z)
 asserted=((L-z)*d+tau)/((L-z)*(d*(mu-z)+tau))
 assert s.simplify(exact-asserted)==0


def test_full_complement_enclosure_against_exact_nontrivial_matrices():
 # The Q block has three different energies and all couplings nonzero.
 # H>=L follows from the positive Schur complement, then its full inverse
 # independently checks both self-energy and source-resolvent bounds.
 L=F(-2);q=[F(1),F(3),F(8)];u=[s.Rational(1),s.Rational(2),s.Rational(3)]
 weights=[F(x*x) for x in u];v=sum(weights);mu=sum(w*x for w,x in zip(weights,q))/v
 variance=sum(w*(x-mu)**2 for w,x in zip(weights,q))/v
 b=L+sum(w/(x-L) for w,x in zip(weights,q))+F(1,3)
 H=s.diag(s.Rational(b.numerator,b.denominator),*map(s.Rational,q))
 for j,value in enumerate(u,1):H[0,j]=H[j,0]=value
 for distance in (F(1,100),F(1,2),F(2),F(20)):
  z=L-distance;S=sum(w/(x-z) for w,x in zip(weights,q))
  for cap in (variance,variance+F(2)):
   bounds=m.scalar_self_energy_enclosure(L,z,b,v,mu,cap)
   assert bounds['lower']<=S<=bounds['upper']<b-z
   exact=F((H-s.Rational(z.numerator,z.denominator)*s.eye(4)).inv()[0,0])
   assert 1/(b-z-bounds['lower'])<=exact<=1/(b-z-bounds['upper'])
 with pytest.raises(ValueError):m.scalar_self_energy_enclosure(L,L,b,v,mu,variance)


def test_low_window_polynomial_controls_full_measure_without_invariance():
 L=F(-18055527,25000);x=m.low_window(L);v=m.fourth_moment()['centered_variance'];d3=m.fourth_moment()['centered_third_moment']
 assert x['polynomial_mean']==v*v/3-3*v-d3==F(13236625211,8103375)
 b=F(-711);r=x['polynomial_root'];d=b-L
 for k in (L-b,(L-b-3)/2,F(-3),F(0),r,r+1,F(10000)):
  p=(-3-k)*(k-r)**2
  if k<=-3:assert p<=x['polynomial_maximum_on_low_interval']
  else:assert p<=0
 assert x['spectral_weight_lower']>F(1,26)
 assert x['reported_original_weight_lower']>F(1,27)
 assert not x['original_selected_state_weight_claimed'] and not x['isolated_line_claimed']


def test_complete_retained_allocation_method_ceiling_is_not_actual_no_gap():
 result=m.retained_method_ceiling()
 assert len(result['cases'])==2
 assert all(row['comparison_floor_strict_upper']<F(-7154,10) for row in result['cases'])
 assert result['full_native_allocation_domain_covered']
 assert not result['proved_actual_complement_spectrum_upper']


def test_original_density_dominance_uses_complete_positive_event_and_marginal_only():
 # A nontrivial retained purification has an off-diagonal source density.
 # The source operator comparison holds, while its tempting extension
 # to a specified source-reference witness has a negative principal minor.
 U=s.Matrix([[3,-4],[4,3]])/5
 reference0=s.Matrix([1,0]);reference1=s.Matrix([0,1])
 psi=(2*s.kronecker_product(U[:,0],reference0)+s.kronecker_product(U[:,1],reference1))/s.sqrt(5)
 joint=psi*psi.H;rho=U*s.diag(s.Rational(4,5),s.Rational(1,5))*U.H
 c=s.Rational(1,5);A=s.Matrix([[0,0],[1,0]]);phi=s.Matrix([0,1])
 E=s.Matrix([[s.Rational(1,2),s.Rational(1,8)],[s.Rational(1,8),s.Rational(1,4)]])
 assert all(x>=0 for x in E.eigenvals()) and all(x>=0 for x in (s.eye(2)-E).eigenvals())
 difference=A*rho*A.H-c*phi*phi.H
 assert all(x>=0 for x in difference.eigenvals())
 out=s.kronecker_product(A,s.eye(2))*joint*s.kronecker_product(A.H,s.eye(2))
 probability=s.trace(s.kronecker_product(E,s.eye(2))*out)
 assert probability==s.trace(E*A*rho*A.H)>c*(phi.H*E*phi)[0]
 invalid=out-c*s.kronecker_product(phi*phi.H,reference0*reference0.H)
 assert invalid.extract((2,3),(2,3)).det()<0
 fail=s.diag(0,1)
 assert s.trace(out+s.kronecker_product(fail,s.eye(2))*joint*s.kronecker_product(fail,s.eye(2)))==1
 # The clock allowance belongs inside the prepared event, then c scales
 # the whole implemented effect; an unscaled subtraction can erase it.
 tiny=s.Rational(1,1000);error=s.Rational(1,512)
 prepared=s.Rational(m.low_window(F(-18055527,25000))['spectral_weight_lower'])
 assert tiny*(prepared-error)>tiny/27>0
 assert tiny*prepared-error<0


def test_actual_member_clock_resolvent_and_mutation_controls():
 member=m.contract.member(0);bound=m.bind_member(member);m.verify_bound(member,bound)
 result=bound['result'];window=result['low_energy_window'];clock=result['complete_low_window_clock_plan']
 assert result['physical_mass_table_GeV']==[]
 assert not result['isolated_primary_odd_band_proved']
 assert not result['uniform_dimensionless_gap_proved']
 assert F(clock['energy_resolution'])==F(1,64)
 assert F(clock['clock_plan']['complete_clock_tail_upper'])+F(clock['one_complete_clock_error_upper'])<=F(1,512)
 assert F(window['reported_excitation_window'][0])>0
 event=result['original_selected_finite_final_event']
 assert result['finite_original_selected_final_energy_event_proved']
 assert event['full_source_CAR_modes']==480 and event['complete_uniform_address_spin_count']==32
 assert event['actual_word_on_filled_beta_anchor_has_unit_norm']
 assert abs(event['actual_word_on_filled_beta_anchor_phase'])==1
 assert event['physical_density_floor_factored']==m.regional.common_contract_payload(member,1)['physical_density_floor_factored']
 assert F(event['simple_fixed_local_original_weight_strict_lower_factored']['multiplier_of_physical_density_floor'])==F(1,27)
 assert F(event['simple_retained_uniform_menu_row_weight_strict_lower_factored']['multiplier_of_physical_density_floor'])==F(1,864)
 assert F(event['retained_uniform_menu_row_weight_lower_factored']['multiplier_of_physical_density_floor'])==F(window['reported_original_weight_lower'])/32
 assert not event['initial_energy_record_applied'] and not event['energy_transfer_measure_claimed']
 assert not event['source_reference_joint_density_domination_claimed']
 assert not event['division_by_success_probability'] and not event['clock_error_subtracted_after_density_floor_scaling']
 for row in result['resolvent_rows']:
  lo,hi=map(F,row['complete_effect_success_interval']);assert 0<lo<=hi<=1
 bad=deepcopy(bound);bad['result']['original_selected_odd_probe_low_energy_residue_proved']=True
 bad['result_sha256']=m.contract.digest(bad['result'])
 with pytest.raises(ValueError):m.verify_bound(member,bad)
 bad=deepcopy(bound);bad['result']['original_selected_finite_final_event']['initial_energy_record_applied']=True
 bad['result_sha256']=m.contract.digest(bad['result'])
 with pytest.raises(ValueError):m.verify_bound(member,bad)
 with pytest.raises(ValueError):m.bind_member(m.contract.member(1))
 other=deepcopy(member);other['source']['electric_coefficient']='106'
 with pytest.raises(ValueError):m.bind_member(other)
