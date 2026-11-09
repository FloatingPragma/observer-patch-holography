"""Exact interacting odd-sector resolvent bounds with full complement control."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product

import sympy as s

from model.particle_common_contract import contract, provenance
from model.particle_primary_composite_transport import model as transport
from model.particle_primary_residual_refinement import model as parent
from model.particle_primary_odd_residual import model as car
from model.particle_spatial_modular_response import regional
from model.particle_primary_selected_odd_response import model as selected_odd

SOURCE_INPUTS=(
 'model/particle_primary_odd_resolvent/__init__.py',
 'model/particle_primary_odd_resolvent/__main__.py',
 'model/particle_primary_odd_resolvent/model.py',
 'model/particle_primary_odd_resolvent/test_model.py',
 'model/particle_primary_odd_resolvent/PROOF.md',
 'model/particle_primary_odd_resolvent/README.md',
 'model/particle_primary_composite_transport/PROOF.md',
 'model/particle_primary_sector_threshold/PROOF.md',
 'model/particle_primary_branch_resolution/PROOF.md',
 'model/particle_spatial_modular_response/PROOF.md',
 *regional.SOURCE_INPUTS,*selected_odd.SOURCE_INPUTS,
 *transport.SOURCE_INPUTS,*parent.SOURCE_INPUTS,*car.SOURCE_INPUTS,
)
DELTA=(0,2,0,0,1)


@lru_cache(maxsize=1)
def cube_data():
 sites=tuple(product(range(2),repeat=3));where={x:i for i,x in enumerate(sites)}
 neighbors={i:tuple(j for j,y in enumerate(sites) if sum(a!=b for a,b in zip(x,y))==1) for i,x in enumerate(sites)}
 blocks={}
 for i,x in enumerate(sites):
  for j in neighbors[i]:
   y=sites[j];axis=next(k for k in range(3) if x[k]!=y[k]);h=transport.spin_blocks()[axis]
   blocks[i,j]=2*(h if x[axis]==0 else h.H)
 return sites,where,neighbors,blocks


@lru_cache(maxsize=None)
def slater_hopping_cumulant(kind='anchor',origin=0,spin_row=0):
 """Exact gauge-Haar fourth cumulant per native dimension on all32 spins."""
 if kind not in ('anchor','particle','hole'):raise ValueError('anchor, particle or hole required')
 if origin not in range(8) or spin_row not in (0,1):raise ValueError('an actual cube cell and beta spin row are required')
 _,_,neighbors,blocks=cube_data()
 P=[]
 for i in range(8):
  diagonal=[0,0,1,1]
  if i==origin and kind=='particle':diagonal[spin_row]=1
  if i==origin and kind=='hole':diagonal[2+spin_row]=0
  P.append(s.diag(*diagonal))
 terms=[s.S.Zero]*5;retained=discarded=0
 variance=s.S.Zero
 for i in range(8):
  for j in neighbors[i]:
   variance+=s.trace(P[i]*blocks[i,j]*(s.eye(4)-P[j])*blocks[j,i])/4
   for k in neighbors[j]:
    for l in neighbors[k]:
     if i not in neighbors[l]:continue
     stack=[]
     for a,b in ((i,j),(j,k),(k,l),(l,i)):
      if stack and stack[-1]==(b,a):stack.pop()
      else:stack.append((a,b))
     if stack:discarded+=1;continue
     retained+=1;a,b,c,d=blocks[i,j],blocks[j,k],blocks[k,l],blocks[l,i]
     words=(P[i]*a*b*c*d,P[i]*a*P[j]*b*c*d,P[i]*a*b*P[k]*c*d,
            P[i]*a*P[j]*b*P[k]*c*d,P[i]*a*P[j]*b*P[k]*c*P[l]*d)
     terms=[old+s.trace(word)/16 for old,word in zip(terms,words)]
 terms=tuple(F(s.simplify(x)) for x in terms);variance=F(s.simplify(variance))
 cumulant=sum((c*x for c,x in zip((1,-4,-3,12,-6),terms)),F())
 if variance!=12 or cumulant!=-18:raise ArithmeticError('actual Slater/Haar hopping cumulant changed')
 return {'occupation_kind':kind,'origin':origin,'beta_spin_row':spin_row,'one_native_dimension_variance':variance,
  'five_Haar_trace_terms':terms,'fourth_cumulant_per_native_dimension':cumulant,
  'retained_backtracking_closed_walks':retained,'vanishing_plaquette_closed_walks':discarded,
  'nontrivial_native_representations_have_zero_single_entry_Haar_mean':True}


@lru_cache(maxsize=1)
def direct_character_CAR():
 """Independent literal32-mode CAR action with all12 Haar exponents."""
 sites,where,_,_=cube_data();terms=[];edge=0;phases={s.S.One:0,s.I:1,-s.S.One:2,-s.I:3}
 for x in sites:
  for axis in range(3):
   if x[axis]:continue
   y=list(x);y[axis]=1;y=tuple(y);a,b=where[x],where[y]
   h=2*transport.spin_blocks()[axis]
   for i,j in product(range(4),repeat=2):
    if h[i,j]:
     phase=phases[h[i,j]]
     terms.extend(((4*a+i,4*b+j,5**edge,phase),(4*b+j,4*a+i,-5**edge,(-phase)%4)))
   edge+=1
 omega=sum(1<<(4*x+i) for x in range(8) for i in (2,3));rows=[]
 for label,mask in (('anchor',omega),('particle',omega|1),('hole',omega^(1<<2))):
  phi={(mask,0):(1,0)};first=car._act(phi,terms);second=car._act(first,terms)
  a,b=car._haar_inner(first,first);c,d=car._haar_inner(second,second)
  if b or d or F(a,4)!=12 or F(c,16)!=414:raise ArithmeticError('literal full-CAR fourth moment differs')
  rows.append({'occupation_kind':label,'first_action_terms':len(first),'second_action_terms':len(second),
               'variance':F(a,4),'fourth_moment':F(c,16)})
 return {'CAR_modes':32,'independent_link_exponents':edge,'rows':rows,
  'all_occupation_outputs_retained':True,'scalar_character_is_an_independent_control_not_a_replacement_source':True}


@lru_cache(maxsize=1)
def fourth_moment():
 slater_hopping_cumulant('anchor')
 for kind,origin,spin_row in product(('particle','hole'),range(8),range(2)):
  slater_hopping_cumulant(kind,origin,spin_row)
 tr=transport.moment_control(DELTA);trial=transport.sector.sector_trial(DELTA)
 vT=F(180);T4=3*vT*vT-18*15
 magnetic=parent.mixed_parent_moment(2,[0]*5,[1]*5)
 B2=magnetic['magnetic_variance'];B4=magnetic['exact_magnetic_mixed_fourth']
 H0square=tr['centered_H0_first_action_norm_squared'];g=trial['loop_norm_squared']
 C=transport.sector.branch.variational.character_control()['cubic_loop_singlet_matrix']
 energies=[210*c for _,c,_ in transport.sector.refinement.native_control()]
 cross=2*g*sum((-energies[i]*F(C[i][j],420) for i,j in product(range(5),repeat=2)),F())
 V4=T4+6*vT*B2+B4
 cap=H0square+V4+cross
 v=tr['centered_variance'];d3=tr['centered_diagonal_third_moment']
 spread=cap/v-v-(d3/v)**2
 if spread<=0 or cap>=tr['centered_second_action_norm_upper']**2:
  raise ArithmeticError('the exact parent fourth moment did not improve the full-space bound')
 return {'hopping_variance':vT,'hopping_fourth_moment':T4,'magnetic_variance':B2,
  'magnetic_fourth_moment':B4,'exact_parent_perturbation_fourth_moment':V4,
  'H0_first_action_norm_squared':H0square,'twice_H0_magnetic_cubic_cross_term':cross,
  'exact_parent_centered_fourth_moment':cap,'actual_cutoff_centered_fourth_moment_upper':cap,
  'actual_cutoff_fourth_moment_equality_claimed':False,
  'centered_variance':v,'centered_third_moment':d3,
  'normalized_complement_coupling_energy_mean':F(-711)+d3/v,
  'normalized_complement_coupling_energy_variance_upper':spread,
  'improved_prepared_transport_time_upper':F(3,2)/cap,
  'transport_time_bound_improvement_factor':tr['centered_second_action_norm_upper']**2/cap,
  'inherited_global_J_second_action_norm_squared_upper':tr['centered_second_action_norm_upper']**2}


def scalar_self_energy_enclosure(L,z,b,v,mu,var):
 """Positive full-complement Stieltjes bounds from admitted exact moments."""
 L,z,b,v,mu,var=map(F,(L,z,b,v,mu,var))
 if not z<L<b or mu<=L or v<=0 or var<0:raise ValueError('positive complement, coupling and moment domain required')
 d=mu-L;e=L-z;anchor=b-L
 lower=v/(mu-z)
 radau=v*(e*d+var)/(e*(d*(e+d)+var))
 positive_floor=v*anchor/(v+e*anchor)
 upper=min(radau,positive_floor,v/e)
 if not 0<lower<=upper<b-z:raise ArithmeticError('actual self-energy interval failed positivity')
 return {'lower':lower,'upper':upper,'Gauss_Radau_upper':radau,'full_operator_floor_upper':positive_floor}


def resolvent_enclosure(sector_floor,energy):
 """Actual rank-one odd-sector Feshbach function for z below its floor."""
 L=F(sector_floor);z=F(energy);b=F(-711)
 m=fourth_moment();v=m['centered_variance'];mu=m['normalized_complement_coupling_energy_mean'];var=m['normalized_complement_coupling_energy_variance_upper']
 bounds=scalar_self_energy_enclosure(L,z,b,v,mu,var)
 lower,upper=bounds['lower'],bounds['upper'];e=L-z
 Glo=1/(b-z-lower);Ghi=1/(b-z-upper)
 return {'real_energy':z,'complement_spectrum_lower':L,'coupling_norm_squared':v,
  'complement_coupling_energy_mean':mu,'complement_coupling_energy_variance_upper':var,
  'self_energy_lower':lower,'self_energy_upper':upper,
  'Gauss_Radau_self_energy_upper':bounds['Gauss_Radau_upper'],
  'full_operator_floor_self_energy_upper':bounds['full_operator_floor_upper'],
  'source_resolvent_expectation_interval':[Glo,Ghi],
  'complete_effect_success_interval':[e*Glo,e*Ghi],
  'complete_effect':'(L-z) P_Delta (H_Delta-z)^(-1) P_Delta, with full complementary outcome',
  'uses_primary_interaction_coefficient':1,'all_orthogonal_physical_sector_directions_retained':True}


def low_window(sector_floor):
 L=F(sector_floor);b=F(-711);a=F(3);v=fourth_moment()['centered_variance'];d3=fourth_moment()['centered_third_moment'];r=v/a;d=b-L
 numerator=v*v/a-a*v-d3;maximum=(d-a)*(d+r)**2
 if not d>a or numerator<=0:raise ArithmeticError('low-energy polynomial witness is not positive')
 weight=numerator/maximum;eta=F(1,64);loss=F(1,512)
 if weight<=F(1,26) or weight-loss<=F(1,27):raise ArithmeticError('source-derived low-band probability failed its public rational floor')
 return {'source_energy_window':[L,b-a],'spectral_weight_lower':weight,
  'simple_spectral_weight_lower':F(1,26),'polynomial_root':r,'polynomial_mean':numerator,
  'polynomial_maximum_on_low_interval':maximum,'clock_resolution':eta,'complete_clock_error_upper':loss,
  'reported_source_energy_window':[L-eta,b-a+eta],
  'reported_original_weight_lower':weight-loss,'simple_reported_weight_lower':F(1,27),
  'prepared_state':'bare weak gauge singlet at one source corner',
  'original_selected_state_weight_claimed':False,'isolated_line_claimed':False}


def retained_method_ceiling():
 dims=transport.sector.refinement.native_control();rows=[]
 for label,shift,penalty in (('r_min<=1/2',F(7),F(12)),('r_min>=1/2',F(12),F(15))):
  loss=12*sum((F(d)/(F(105,2)*C-shift) for d,C,_ in dims),F())+F(6,7)
  ceiling=-720-loss+penalty
  if ceiling>=F(-7154,10):raise ArithmeticError('retained-penalty method ceiling lost separation')
  rows.append({'parameter_case':label,'necessary_scalar_loss_strict_lower':loss,
               'bare_complement_penalty_upper':penalty,'comparison_floor_strict_upper':ceiling})
 return {'cases':rows,'full_native_allocation_domain_covered':True,
  'full_bare_minimal_sector_complement':True,'extra_neutral_pair_is_physical_and_orthogonal_to_minimal_defects':True,
  'proved_actual_complement_spectrum_upper':False}


def original_final_event(member,window):
 """Actual original-density final-clock event, without an initial energy readout."""
 transport.sector.branch._member_zero(member)
 V=member['source']['side']**3
 payload=regional.common_contract_payload(member,1)
 floor=payload['physical_density_floor_factored']
 if floor['base']!=2 or floor['operator']!='P_G' or type(floor['negative_exponent']) is not int or floor['negative_exponent']<=0:
  raise ArithmeticError('authenticated physical source-density floor changed')
 selected=selected_odd.exact_controls()
 blocks=dict(transport.sector.native_blocks());word=tuple(blocks[-3])+tuple(blocks[6])
 omega=sum(1<<(60*x+i) for x in range(V) for i in range(30,60))
 target,phase=transport.spin.apply_operators(omega,tuple((i,True) for i in word))
 local=transport.sector.local_singlet(DELTA)
 if len(local['state'])!=1 or local['norm_squared']!=1 or abs(phase)!=1:
  raise ArithmeticError('actual creation word does not supply a unit Slater witness')
 expected=next(iter(local['state']))|(omega&~((1<<60)-1))
 if target!=expected or selected['spin_choices_per_cell']!=4:
  raise ArithmeticError('actual original-probe row differs from the prepared witness')
 prepared=F(window['reported_original_weight_lower'])
 if prepared<=F(1,27):raise ArithmeticError('complete prepared clock event lacks its positive margin')
 return {'physical_density_floor_factored':floor,'physical_floor_payload_sha256':contract.digest(payload),
  'fixed_local_original_weight_lower_factored':{'multiplier_of_physical_density_floor':prepared},
  'simple_fixed_local_original_weight_strict_lower_factored':{'multiplier_of_physical_density_floor':F(1,27)},
  'retained_uniform_menu_row_weight_lower_factored':{'multiplier_of_physical_density_floor':prepared/(4*V)},
  'simple_retained_uniform_menu_row_weight_strict_lower_factored':{'multiplier_of_physical_density_floor':F(1,108*V)},
  'complete_uniform_address_spin_count':4*V,'full_source_CAR_modes':60*V,
  'actual_local_creation_word_modes':word,'actual_word_on_filled_beta_anchor_phase':phase,
  'actual_word_on_filled_beta_anchor_has_unit_norm':True,
  'reported_source_energy_window':window['reported_source_energy_window'],
  'reported_excitation_window':window['reported_excitation_window'],
  'final_clock_is_the_same_complete_implemented_positive_event':True,
  'density_domination_is_on_source_marginal_only':True,
  'source_reference_joint_density_domination_claimed':False,
  'original_source_input_replaced_by_prepared_witness':False,
  'initial_energy_record_applied':False,'energy_transfer_measure_claimed':False,
  'division_by_success_probability':False,'clock_error_subtracted_after_density_floor_scaling':False,
  'compensating_memory_and_full_failures_retained':True,
  'uniform_or_practically_resolvable_original_weight_claimed':False}


def bind_member(member):
 transport.sector.branch._member_zero(member)
 inherited=transport.sector.bind_member(member);transport.sector.verify_bound(member,inherited)
 row=inherited['result']['sector_results'][0];authenticated_floor=F(row['signed_sectors'][0]['physical_sector_energy_lower'])
 L=F(-18055527,25000)
 if L>authenticated_floor:raise ArithmeticError('the short rational sector floor is not outward')
 mean=F(row['signed_sectors'][0]['trial']['sector_ground_energy_upper'])
 obstruction=retained_method_ceiling()
 if any(F(x['comparison_floor_strict_upper'])>=mean for x in obstruction['cases']):
  raise ArithmeticError('the allocation obstruction does not cover the authenticated trial comparison')
 clock=transport.sector.spectroscopy.clock_plan(member,'1/64','1/1024','1/1024')
 loss=F(clock['clock_plan']['complete_clock_tail_upper'])+F(clock['one_complete_clock_error_upper'])
 window=low_window(L)
 if loss>window['complete_clock_error_upper'] or F(clock['energy_resolution'])!=window['clock_resolution']:
  raise ArithmeticError('complete clock budgets do not support the reported low window')
 U0=F(inherited['result']['global_ground_upper']);L0=F(inherited['result']['global_ground_lower'])
 window['reported_excitation_window']=[window['reported_source_energy_window'][0]-U0,
                                       window['reported_source_energy_window'][1]-L0]
 if window['reported_excitation_window'][0]<=0:raise ArithmeticError('reported window is not above the full ground')
 result=contract.serial({'member_id':member['member_id'],**member['ids'],
  'native_occupation_increment':DELTA,'opposite_sign_has_same_source_energy_law':True,
  'opposite_sign_comparison_uses_literal_conjugate_preparation':True,
  'sector_bound_receipt_sha256':contract.digest(inherited),'physical_sector_floor':L,
  'authenticated_physical_sector_floor':authenticated_floor,
  'Slater_Haar_controls':[slater_hopping_cumulant(k) for k in ('anchor','particle','hole')],
  'all32_Slater_position_and_spin_controls':[
   slater_hopping_cumulant(k,x,r) for k,x,r in product(('particle','hole'),range(8),range(2))],
  'independent_character_CAR_control':direct_character_CAR(),'fourth_moment':fourth_moment(),
  'resolvent_rows':[resolvent_enclosure(L,L-e) for e in (F(1,2),F(1),F(2),F(4),F(8),F(16))],
  'low_energy_window':window,'retained_allocation_method_obstruction':obstruction,
  'original_selected_finite_final_event':original_final_event(member,window),
  'complete_low_window_clock_plan':clock,
  'primary_contract_replaced':False,'prepared_vector_experiment_explicit':True,
  'source_marginal_changes_by_explicit_preparation':True,
  'original_selected_odd_probe_low_energy_residue_proved':False,
  'original_selected_odd_probe_low_energy_residue_flag_scope':'uniform physical particle residue, not finite final-event positivity',
  'finite_original_selected_final_energy_event_proved':True,
  'isolated_primary_odd_band_proved':False,'uniform_dimensionless_gap_proved':False,
  'physical_mass_table_GeV':[]})
 return {'schema':'oph.primary_odd_resolvent.v1','result':result,'result_sha256':contract.digest(result)}


def verify_bound(member,value):
 if value!=bind_member(member):raise ValueError('source odd resolvent or full-complement scope changed')


def certificate():
 paths=provenance.definition_inputs(SOURCE_INPUTS);pins=provenance.definition_hashes(paths)
 result={'schema':'oph.primary_odd_resolvent_certificate.v1','selected_member_bound':bind_member(contract.member(0)),
         'source_sha256':pins}
 if provenance.definition_hashes(paths)!=pins:raise ValueError('source changed during odd resolvent replay')
 return result
