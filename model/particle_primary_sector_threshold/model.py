"""Full physical thresholds in conserved native-occupation sectors."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product
from math import gcd

import sympy as s

from model.particle_common_contract import contract, provenance
from model.particle_primary_branch_resolution import model as branch
from model.particle_primary_refinement_response import model as refinement
from model.particle_primary_occupation_bound import model as occupation
from model.particle_primary_capacity import capacity
from model.particle_primary_joint_spectroscopy import spectroscopy
from model.native_quantum_gauss import algebra as native
from model.spin_packet_readout.bounds import sqrt_upper

SOURCE_INPUTS=(
 'model/particle_primary_sector_threshold/__init__.py',
 'model/particle_primary_sector_threshold/__main__.py',
 'model/particle_primary_sector_threshold/model.py',
 'model/particle_primary_sector_threshold/test_model.py',
 'model/particle_primary_sector_threshold/PROOF.md',
 'model/particle_primary_sector_threshold/README.md',
 'model/particle_primary_branch_resolution/PROOF.md',
 'model/particle_primary_occupation_bound/PROOF.md',
 'model/particle_primary_refinement_response/PROOF.md',
 'model/particle_primary_spectral_resolution/PROOF.md',
 'model/particle_primary_capacity/PROOF.md',
 'model/particle_primary_joint_spectroscopy/PROOF.md',
 'model/native_quantum_gauss/algebra.py',
 'model/particle_common_contract/PROOF.md',
)
# Rational choices optimize a source inequality, never a measured mass.
SECTORS=(
 ('weak_pair_charge_six',(0,2,0,0,1),(F(),F(73,100),F(),F(),F(72,100))),
 ('color_triple',(1,0,0,2,0),(F(52,100),F(),F(),F(78,100),F())),
 ('color_determinant_charge_six_hole',(0,0,0,3,-1),(F(),F(),F(),F(81,100),F(66,100))),
)


@lru_cache(maxsize=1)
def native_blocks():
 blocks=tuple(sorted((int(-2*s.Rational(center)),indices) for center,indices in native.native_blocks()))
 if tuple((q,len(indices)) for q,indices in blocks)!=tuple((q,d) for d,_,q in refinement.native_control()):raise ArithmeticError('native current chart and Q6 order disagree')
 return blocks


def _bilinear(mask,i,j):
 return refinement._bilinear(mask,i,j)


def _act(vector,terms):
 out={}
 for mask,value in vector.items():
  for i,j,c in terms:
   moved=_bilinear(mask,i,j)
   if moved:out[moved[0]]=out.get(moved[0],0)+value*c*moved[1]
 return {k:s.simplify(v) for k,v in out.items() if v!=0}


def _norm(vector):
 return s.simplify(sum(s.conjugate(v)*v for v in vector.values()))


def _defects(mask,sites=1):
 return sum(sum(bool(mask&(1<<(60*x+j))) for j in range(30))+
            sum(not bool(mask&(1<<(60*x+j))) for j in range(30,60)) for x in range(sites))


@lru_cache(maxsize=None)
def local_singlet(delta):
 """Solve the actual local Gauss kernel on explicit particle/hole masks."""
 delta=tuple(delta);blocks=native_blocks();generators=native.native_currents()
 if len(delta)!=5:raise ValueError('five native occupation increments required')
 if sum(q*n for (q,_),n in zip(blocks,delta))!=0:raise ValueError('total native charge must vanish')
 omega=sum(1<<j for j in range(30,60));choices=[]
 for (_,indices),number in zip(blocks,delta):
  choices.append(tuple(combinations(indices,abs(number))))
 masks=[]
 for selected in product(*choices):
  mask=omega
  for labels,number in zip(selected,delta):
   for i in labels:
    if number>0:mask|=1<<i
    elif number<0:mask^=1<<(30+i)
  masks.append(mask)
 masks=tuple(sorted(set(masks)));where={mask:i for i,mask in enumerate(masks)}
 matrices=[]
 for g in generators:
  terms=[(15*row+i,15*row+j,value) for row in range(4) for (i,j),value in g.todok().items() if value]
  A=s.zeros(len(masks))
  for column,mask in enumerate(masks):
   for target,value in _act({mask:s.S.One},terms).items():
    if target not in where:raise ArithmeticError('the fixed-occupation local basis is not gauge invariant')
    A[where[target],column]=value
  matrices.append(A)
 kernel=s.Matrix.vstack(*matrices).nullspace()
 if len(kernel)!=1:raise ArithmeticError('the declared local singlet is not one dimensional')
 v=kernel[0];den=s.ilcm(*[x.q for x in v]) if len(v)>1 else v[0].q
 v=[s.expand(x*den) for x in v]
 if any(not x.is_Integer for x in v):raise ArithmeticError('exact integer singlet coordinates expected')
 common=0
 for x in v:common=gcd(common,int(x))
 v=[int(x)//common for x in v]
 if next(x for x in v if x)<0:v=[-x for x in v]
 state={mask:s.Integer(x) for mask,x in zip(masks,v) if x}
 for A in matrices:
  if A*s.Matrix(v)!=s.zeros(len(masks),1):raise ArithmeticError('actual local Gauss generator does not annihilate the witness')
 norm=_norm(state);D=sum(map(abs,delta))
 if any(_defects(mask)!=D for mask in state):raise ArithmeticError('witness has incorrect beta defects')
 for (_,indices),n in zip(blocks,delta):
  if any(sum(bool(mask&(1<<(15*r+i))) for r in range(4) for i in indices)!=2*len(indices)+n for mask in state):
   raise ArithmeticError('witness native occupation mismatch')
 return {'delta':delta,'state':state,'norm_squared':int(norm),'local_basis_dimension':len(masks),
         'Gauss_generators_checked':len(generators),'defect_count':D,'CAR_parity':D%2}


@lru_cache(maxsize=None)
def hopping_packets(delta):
 """Actual Haar/CAR first-action norms on a localized physical singlet."""
 witness=local_singlet(tuple(delta));omega=sum(1<<j for j in range(30,60))
 vector={mask|(omega<<60):value for mask,value in witness['state'].items()}
 norm=witness['norm_squared'];D=witness['defect_count'];pairs=branch.response.beta_pairs()
 chart=s.Matrix.hstack(pairs[0][0],pairs[1][0],pairs[0][1],pairs[1][1])
 gamma=branch.response.density.gamma_data();beta=gamma[3]
 rows=[]
 for (charge,indices),(d,C,_) in zip(native_blocks(),refinement.native_control()):
  totals={0:s.S.Zero,2:s.S.Zero};axis_rows=[]
  for axis in range(3):
   hop=branch.apparatus.clean(chart.H*((-beta+beta*gamma[axis])/2)*chart)
   parts={0:s.S.Zero,2:s.S.Zero}
   for i,j in product(indices,repeat=2):
    forward=[(15*a+i,60+15*b+j,hop[a,b]) for a,b in product(range(4),repeat=2) if hop[a,b]]
    reverse=[(y,x,s.conjugate(c)) for x,y,c in forward]
    for terms in (forward,reverse):
     output=_act(vector,terms)
     for mask,value in output.items():
      shift=_defects(mask,2)-D
      if shift not in parts:raise ArithmeticError('unaccounted first-action defect sector')
      parts[shift]+=s.conjugate(value)*value/(d*norm)
   parts={k:F(s.simplify(v)) for k,v in parts.items()}
   axis_rows.append(parts)
   for k,v in parts.items():totals[k]+=v
  n0=F(totals[0]);n2=F(9*d)+F(totals[2]) # nine bonds do not meet the chosen corner
  expected=F(3*abs(delta[len(rows)]),4)
  if n0!=expected or n2!=12*d-expected:raise ArithmeticError(f'exact localized Haar/CAR norms changed: delta={delta},charge={charge},n0={n0},n2={n2},expected={expected},axes={axis_rows}')
  rows.append({'Q6':charge,'dimension':d,'Casimir':C,'affected_axis_norms':axis_rows,
               'unchanged_defect_packet_norm_squared':n0,'two_added_defect_packet_norm_squared':n2,
               'electric_energy':F(105,2)*C})
 return rows


@lru_cache(maxsize=None)
def local_allocation(retained):
 retained=tuple(map(F,retained));native_rows=refinement.native_control();onsite=tuple(1-r for r in retained)
 if len(retained)!=5 or any(not 0<=r<1 for r in retained):raise ValueError('five fractions in [0,1) required')
 delta=4*sum((d*(sqrt_upper(m*m+F(1,4))-m) for (d,_,_),m in zip(native_rows,onsite)),F())
 if delta>=F(63,4):raise ArithmeticError('nonconstant link block is not positive')
 rows=[]
 for (d,C,q),m,r in zip(native_rows,onsite,retained):
  w=1/(F(105,2)*C-delta)
  proof=refinement.exact_q6_loss(m,w) if d==1 else occupation.occupation_loss(d,m,w)
  rows.append({'Q6':q,'dimension':d,'onsite':m,'retained_fraction':r,'inverse_electric_margin':w,**proof})
 loss=12*sum((F(row['occupation_loss_upper']) for row in rows),F())+F(6,7)
 return {'retained_fractions':retained,'pointwise_fermion_shift_upper':delta,
         'nonconstant_link_margin_lower':F(63,4)-delta,'native_rows':rows,
         'full_source_scalar_loss_upper':loss,'operator_floor_constant':-720-loss}


@lru_cache(maxsize=None)
def sector_trial(delta):
 """Exact finite Ritz expectation; trial coefficients do not alter the action."""
 delta=tuple(delta);witness=local_singlet(delta);rows=hopping_packets(delta)
 packets=[]
 for row in rows:
  for shift,key in ((0,'unchanged_defect_packet_norm_squared'),(2,'two_added_defect_packet_norm_squared')):
   n=row[key]
   if n:packets.append({'Q6':row['Q6'],'added_defects':shift,'norm_squared':n,'relative_H0_energy':row['electric_energy']+3*shift})
 control=branch.variational.character_control();C=control['cubic_loop_singlet_matrix']
 native_rows=refinement.native_control();g=F(6,22050)
 A=s.diag(*[s.Rational(210*c.numerator,c.denominator) for _,c,_ in native_rows])-s.Matrix(C)/420
 u=(A+4*s.eye(5)).inv()*s.ones(5,1);u=tuple(F(x) for x in u)
 t=tuple(1/(row['relative_H0_energy']+4) for row in packets)
 square=sum((row['norm_squared']*x*x for row,x in zip(packets,t)),F())+g*sum((x*x for x in u),F())
 q=sum((row['norm_squared']*x for row,x in zip(packets,t)),F())+g*sum(u)
 numerator=-2*q+sum((row['norm_squared']*row['relative_H0_energy']*x*x for row,x in zip(packets,t)),F())
 numerator+=g*sum((u[i]*F(A[i,j])*u[j] for i,j in product(range(5),repeat=2)),F())
 shift=numerator/(1+square);bare=-720+3*witness['defect_count']
 return {'native_occupation_increment':delta,'local_singlet':{k:v for k,v in witness.items() if k!='state'},
  'local_singlet_integer_coordinates':[[str(mask),int(value)] for mask,value in witness['state'].items()],
  'hopping_norm_control':rows,'hopping_packets':packets,'hopping_coefficients':t,
  'loop_norm_squared':g,'loop_packet_matrix':[[F(x) for x in row] for row in A.tolist()],
  'loop_coefficients':u,'packet_dimension':1+len(packets)+5,'trial_norm_squared':1+square,
  'bare_trial_source_energy':bare,'trial_Rayleigh_shift':shift,'sector_ground_energy_upper':bare+shift,
  'full_residual_or_isolated_line_computed':False,'coefficients_use_measured_mass':False}


def prepared_record(member,trial,lower,ground_lower,ground_upper):
 mean=trial['sector_ground_energy_upper'];eps=F(member['state']['target_ground_excess'])
 plan=spectroscopy.clock_plan(member,'1/64','1/1024','1/1024');eta=F(plan['energy_resolution'])
 loss=F(plan['clock_plan']['complete_clock_tail_upper'])+F(plan['one_complete_clock_error_upper'])
 high=2*mean-lower
 if not lower-eta>ground_upper or mean<lower:raise ArithmeticError('the prepared sector window is not above the full ground')
 cap=capacity.member_capacity(member)
 return {'intervention':'paired-sector Householder on two idle physical memories; graded exchange; parity correction; one energy clock',
  'source_post_intervention_state':'normalized source sector Ritz trial',
  'source_marginal_changes':True,'original_state_and_reference_correlations_retained_in_memory':True,
  'compensating_opposite_occupation_sector_retained_in_second_memory':True,
  'total_source_and_memory_native_occupations_preserved':True,'combined_CAR_parity_preserved':True,
  'additional_memory_qubit_upper':2*cap['allocated_source_qubit_count_upper'],
  'source_work_mean_interval':[mean-ground_upper-eps,mean-ground_lower],
  'exact_spectral_window':[lower,high],'original_window_weight_lower':F(1,2),
  'reported_source_energy_window':[lower-eta,high+eta],
  'reported_excitation_window':[lower-eta-ground_upper,high+eta-ground_lower],
  'reported_window_original_weight_lower':F(1,2)-loss,
  'complete_clock_plan':plan,'complete_energy_outcomes':plan['clock_plan']['clock_size'],
  'all_energy_outcomes_and_quantum_references_retained':True,
  'autonomous_controller_or_efficient_compilation_proved':False}


def bind_member(member):
 branch._member_zero(member)
 ground=branch.occupation.bind_member(member);branch.occupation.verify_bound(member,ground)
 lo=F(ground['result']['zero_retention_physical_ground_lower']);up=F(ground['result']['physical_ground_upper'])
 eps=F(member['state']['target_ground_excess']);rows=[]
 for label,delta,retained in SECTORS:
  allocation=local_allocation(retained)
  lower=allocation['operator_floor_constant']+3*sum((r*abs(n) for r,n in zip(retained,delta)),F())
  gap=lower-up
  if gap<=0:raise ArithmeticError('the full conserved-sector threshold did not separate from the ground')
  signs=[]
  for sign in (1,-1):
   actual=tuple(sign*n for n in delta);trial=sector_trial(actual)
   signs.append({'sign':sign,'delta':actual,'trial':trial,
    'physical_sector_energy_lower':lower,
    'sector_bottom_excitation_interval':[gap,trial['sector_ground_energy_upper']-lo],
    'prepared_spectroscopy':prepared_record(member,trial,lower,lo,up)})
  rows.append({'label':label,'allocation':allocation,'signed_sectors':signs,
               'selected_original_probability_of_either_sign_upper':min(F(1),eps/gap),
               'selected_sector_probability_strictly_positive_from_faithfulness':True})
 result=contract.serial({'member_id':member['member_id'],**member['ids'],
  'native_charge_order':[-4,-3,1,2,6],'anchor_native_occupations':[48,32,96,48,16],
  'global_ground_lower':lo,'global_ground_upper':up,'ground_receipt_sha256':contract.digest(ground),
  'sector_results':rows,'source_state_reference_clock_unchanged_for_bounds':True,
  'source_preparation_intervention_explicit':True,'charge_conjugation_spectral_equality_assumed':False,
  'isolated_particle_branch_or_relativistic_dispersion_proved':False,'physical_mass_table_GeV':[]})
 return {'schema':'oph.primary_sector_threshold.v1','result':result,'result_sha256':contract.digest(result)}


def verify_bound(member,value):
 if value!=bind_member(member):raise ValueError('source sector threshold or retained-state result changed')



@lru_cache(maxsize=1)
def paired_memory_control():
 """Exact graded exchange with a retained opposite odd sector and reference."""
 parity=(0,1,1);charge=(0,1,-1)
 def gate(state,correct=True):
  x,a,b,r=state
  if (a,b)==(0,0):a,b=1,2
  elif (a,b)==(1,2):a,b=0,0
  sign=(-1)**(parity[x]*parity[a]);x,a=a,x
  if correct:sign*=(-1)**parity[a]
  return (x,a,b,r),sign
 # The full gate is bijective, even, and preserves combined number.
 outputs=[]
 for state in product(range(3),repeat=4):
  target,sign=gate(state);outputs.append(target)
  if sum(charge[x] for x in state[:3])!=sum(charge[x] for x in target[:3]):raise ArithmeticError('paired memory changes total occupation')
  if sum(parity[x] for x in state[:3])%2!=sum(parity[x] for x in target[:3])%2:raise ArithmeticError('paired memory changes total parity')
 if len(set(outputs))!=81:raise ArithmeticError('memory control is not unitary')
 for x,r in product(range(3),repeat=2):
  if gate((x,0,0,r))!=((1,x,2,r),1):raise ArithmeticError('arbitrary source/reference coherence was not retained')
 uncorrected=[gate((x,0,0,x),False)[1] for x in range(3)]
 if uncorrected!=[1,-1,-1]:raise ArithmeticError('odd-sector exchange phase was omitted')
 return {'complete_control_dimension':81,'source_and_reference_dimension':3,
  'paired_even_Householder_exchange_checked':True,'full_unitary_checked':True,
  'total_occupation_and_parity_preserved':True,'arbitrary_reference_basis_identity_checked':True,
  'uncorrected_source_parity_phases':uncorrected,'parity_correction_required_for_odd_target':True,
  'original_source_information_discarded':False}


def certificate():
 paths=provenance.definition_inputs((*SOURCE_INPUTS,*branch.SOURCE_INPUTS,*refinement.SOURCE_INPUTS,*capacity.SOURCE_INPUTS))
 pins=provenance.definition_hashes(paths)
 result={'schema':'oph.primary_sector_threshold_certificate.v1','selected_member_bound':bind_member(contract.member(0)),
         'paired_memory_control':paired_memory_control(),'source_sha256':pins}
 if provenance.definition_hashes(paths)!=pins:raise ValueError('source inputs changed during threshold replay')
 return result
