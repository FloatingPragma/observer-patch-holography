"""Finite corrected source evolution, complete quadratures and interval inference."""
from pathlib import Path
import json
from fractions import Fraction as F
from functools import lru_cache
from math import pi,sqrt,factorial
import numpy as np
from model.common_source import repair_continuation as rc,parameter_selection as sel,record_code as code
from model.observer.source_scalar import compiler
from model.native_effect_synthesis import synthesis as syn
from model.characteristic_constraints.bounds import trig_interval
from model.openwave_native_process.source import raw_local_branches
from model.openwave_coherent_native_architecture.source import construction_descriptor
ROOT=Path(__file__).resolve().parents[2]

def strings(x):return [str(v)for v in x]
def pair_vectors(spectator,axis='X'):
 if spectator not in range(3) or axis not in ('X','Y'):raise ValueError('three spectator sectors and X/Y required')
 a,b=2*spectator,2*spectator+1;u=np.ones(6)/sqrt(6);v=np.zeros((36,2),complex);v[6*a:6*a+6,0]=u;v[6*b:6*b+6,1]=u
 phase=1 if axis=='X'else 1j
 return [(v[:,0]+sign*phase*v[:,1])/sqrt(2)for sign in (1,-1)]

def effects(spectator,axis):
 vectors=pair_vectors(spectator,axis);P=[np.outer(x,x.conj())for x in vectors];return P+[np.eye(36)-sum(P)]

@lru_cache(maxsize=32)
def powers(law,time,steps):
 geo=sel.code_geometry();energy=sel.diagonal_cost(2)['values'];values=np.array([[float(energy[k])for k in row]for row in geo['supports']]);out={}
 for a in range(6):
  for b in range(a,6):
   phase=np.exp(-1j*float(time/steps)*(values[a,:,None]-values[b,None,:])).ravel();out[a,b]=np.linalg.matrix_power(rc._pair_transfer(law,a,b)*phase,steps)
 return out

def evolve(rho,law,time,steps):
 out=np.zeros_like(rho);maps=powers(law,time,steps)
 for a in range(6):
  for b in range(a,6):
   block=(maps[a,b]@rho[6*a:6*a+6,6*b:6*b+6].ravel()).reshape(6,6);out[6*a:6*a+6,6*b:6*b+6]=block
   if a!=b:out[6*b:6*b+6,6*a:6*a+6]=block.conj().T
 return out

def exact_readout_control():
 # A supplied coherent chart maps the two code quadrature rays to literal
 # native states (0,2) and (1,2). E-law signed parity distinguishes them.
 states=(5,6);records=[];weights=[F(0),F(0)]
 for slot,branch in enumerate(raw_local_branches('E')):
  row={'slot':slot,'record':list(branch.record),'input_probabilities_exact':[],'native_images':[]}
  mapping=dict(branch.pairs)
  for index,mask in enumerate(states):
   local=mask&3;weight=branch.weight if local in mapping else F(0)
   if weight:assert branch.record[-1]==(1 if index==0 else -1)
   weights[index]+=weight;row['input_probabilities_exact'].append(str(weight));row['native_images'].append((mask&~3)|mapping[local]if local in mapping else None)
  records.append(row)
 assert weights==[F(1),F(1)]
 # Reversible record-selected source routing followed by inverse chart
 # restores the measured quadrature eigenstate. Its record is retained.
 return {'all_native_E_raw_slots':len(raw_local_branches('E')),'fine_records':records,'probability_per_input_exact':strings(weights),'outside_code_outcome_retained':True,'record_selected_coherent_routing_supplied':True,'measurement_reset_is_not_free':True,'ideal_instrument':'P_plus rho P_plus, P_minus rho P_minus, Q rho Q; each native fine branch retained','full512_complement':'Q=I512-P_plus-P_minus; no projection or renormalization of Q'}

def compiled_quadrature_control():
 outputs=[]
 for axis in ('X','Y'):
  word=([compiler.pulse('Z',pi/4,[0])]if axis=='Y'else[])+[compiler.pulse('H',pi/2,[0])]
  primitive=[syn.certify_primitive(x['kind'],x['angle'])for x in word]
  U=compiler.execute(word,1);phase=1 if axis=='X'else 1j
  rays=np.array([[1,1],[phase,-phase]],complex)/sqrt(2);p=abs(U@rays)**2
  assert np.max(abs(p-np.eye(2)))<2e-11
  outputs.append({'axis':axis,'word':word,'complete_binary_input_probabilities_diagnostic':p.tolist(),'primitive_interval_checks':len(primitive),'native_pulses':len(word),'proper_load_export_services':2*len(word),'exact_phase_target_angle_error_upper':'1/2251799813685248','processor_complement_outcomes_retained':4*len(word)})
 return {'controls':outputs,'exact_source_response_identity_checked':True,'stored_primitive_interval_gate':'2e-11 per primitive; inherited complete gain/error composition','arbitrary_phase_clock_or_physical_precision_inferred':False,'code_load_and_export_are_declared_architecture':True}

def interval_joint(gap,time,axis,epsilon,twirled=False):
 c=trig_interval(gap*time,True);v=c if axis=='X'else trig_interval(gap*time)
 # rho01 evolves exp(+i gap*t), so Y expectation is -sin(gap*t).
 if axis=='Y':v=(-v[1],-v[0])
 if twirled:v=(F(0),F(0))
 out=[]
 for i,sign in enumerate((1,-1)):
  first=tuple((1+sign*z)/2 for z in v);first=tuple(sorted(first))
  row=[]
  for j,second in enumerate((1,-1)):
   cond=tuple(sorted((F(1,2)+sign*second*z/2 for z in c)));corners=[x*y for x in first for y in cond];row.append([max(F(0),min(corners)-epsilon),min(F(1),max(corners)+epsilon)])
  row.append([F(0),min(F(1),epsilon)]);out.append(row)
 out.append([[F(0),min(F(1),epsilon)]for _ in range(3)])
 return out

@lru_cache(maxsize=128)
def infer_gap(sine_interval,time):
 # Inverse sees record bounds and the declared phase window only.
 lo,hi=-F(3,2),F(3,2)
 def inverse(y):
  if y<=-1:return lo
  if y>=1:return hi
  a,b=lo,hi
  for _ in range(48):
   m=(a+b)/2
   center=sum(((-1)**(k//2)*m**k/F(factorial(k))for k in range(1,26,2)),F(0));error=abs(m)**26/F(factorial(26));l,u=center-error,center+error
   if u<y:a=m
   elif l>y:b=m
   else:return a,b
  return a,b
 a=inverse(sine_interval[0]);b=inverse(sine_interval[1]);a=(a,a)if isinstance(a,F)else a;b=(b,b)if isinstance(b,F)else b
 return a[0]/time,b[1]/time

def charge_recovery_join():
 d=json.loads((ROOT/'model/openwave_native_charge_readout/results.json').read_text())
 j=d['matching_clock_recovery_join'];w=j['nonzero_charge_witness'];services=d['complete_transport_compiler']['compiled_services']
 assert w['law']=='C' and w['charge_character']==1 and w['logical_labels']==[0,2] and w['coherence_multiplier_real_exact']=='85/96'
 assert all(x['source_memory_input_columns']==4096 and x['actual_emitted_forward_inverse_identity_on_all_admitted_columns']for x in services)
 return {'common_hardware':construction_descriptor(),'uncompensated_counterexample':w,'exact_compensated_identity':j['compensated_identity'],'compensated_policy':j['compensated_policy'],'arbitrary_reference_clock_channel_transfers_under_compensation':True,'compiled_feedback_services':len(services),'maximum_added_native_pulses_per_event':2*max(x['native_pulses']for x in services),'maximum_added_proper_code_services_per_event':2*max(x['proper_code_service_events']for x in services),'all_source_memory_columns':4096,'dirty_scratch_and_processor_failures_retained':True,'coherent_charge_dependent_frequency_derived':False,'active_charge_signal_timing':j['active_charge_record_timing'],'stored_compiled_error_budget_additional':True}

def run():
 h=tuple(F(x)for x in rc.logical_hamiltonian().diagonal());gaps=[h[2*j+1]-h[2*j]for j in range(3)];assert gaps==[F(1,256),-F(1,1024),-F(1,3072)]
 readout=exact_readout_control();compiled=compiled_quadrature_control();rows=[]
 energy=np.diag([float(sel.diagonal_cost(2)['values'][k])for support in sel.code_geometry()['supports']for k in support])
 for law in ('C','E'):
  for time,steps in [(F(64),4096),(F(64),65536),(F(256),65536)]:
   err=F(rc.error_bound(law,time,steps)['trace_norm_upper_exact']);joint_error=min(F(2),2*err)
   for spectator,gap in enumerate(gaps):
    for axis in ('X','Y'):
     P=effects(spectator,axis);psi=pair_vectors(spectator,'X')[0];rho=np.outer(psi,psi.conj());first=evolve(rho,law,time,steps);joint=np.zeros((3,3))
     for i,p in enumerate(P):
      second=evolve(p@first@p,law,time,steps)
      for j,q in enumerate(P):joint[i,j]=np.trace(q@second).real
     after1=sum(p@first@p for p in P);second_unread=evolve(after1,law,time,steps);after2=sum(p@second_unread@p for p in P)
     source_work=[float(np.trace(energy@(x-y)).real)for x,y in [(first,rho),(after1,first),(second_unread,after1),(after2,second_unread)]]
     assert abs(sum(source_work)-np.trace(energy@(after2-rho)).real)<2e-12
     box=interval_joint(gap,time,axis,joint_error)
     assert abs(joint.sum()-1)<2e-10 and joint.min()>-2e-12
     assert all(float(box[i][j][0])-2e-11<=joint[i,j]<=float(box[i][j][1])+2e-11 for i in range(3)for j in range(3))
     sine=trig_interval(gap*time);inferred=infer_gap((max(F(-1),sine[0]-err),min(F(1),sine[1]+err)),time)
     assert inferred[0]<=gap<=inferred[1]
     rows.append({'law':law,'spectator':spectator,'axis':axis,'time_per_block_exact':str(time),'steps_per_block':steps,'signed_gap_exact':str(gap),'two_time_probabilities_diagnostic':joint.tolist(),'complete_joint_probability_intervals_exact':[[strings(z)for z in line]for line in box],'two_time_L1_bound_exact':str(joint_error),'single_block_trace_norm_bound_exact':str(err),'gap_interval_from_sine_record_exact':strings(inferred),'phase_branch': 'abs(gap*time)<=1<pi/2; no winding ambiguity','inference_input':'predicted complete-record interval plus finite-source bound; no independent experimental sample','source_energy_before_readout_diagnostic':float(np.trace(energy@first).real),'complete_nonselective_source_work_diagnostic':dict(zip(['correction_block1','readout1','correction_block2','readout2'],source_work))})
 # Diagonal logical preparation remains unchanged in the ideal effective
 # evolution; the same finite source bound controls its measured departure.
 controls=[]
 for j,gap in enumerate(gaps):
  P=effects(j,'X');a,b=pair_vectors(j,'X');rho=(np.outer(a,a.conj())+np.outer(b,b.conj()))/2
  out=evolve(rho,'E',F(256),65536);contrast=np.trace((P[0]-P[1])@out).real
  assert abs(contrast)<2e-11
  controls.append({'spectator':j,'twirled_initial_X_contrast_diagnostic':float(contrast),'source_preparation_invariant_under_logical_phase':True,'self_starting_phase_observed':False})
 return {'schema':'oph.complete_native_coherent_clock.v1','verdict':'PASS_THREE_SIGNED_GAPS_AND_COMPLETE_REPEATED_NATIVE_QUADRATURES','shared_charge_clock_program':charge_recovery_join(),'signed_gaps_exact':strings(gaps),'same_source_logical_energies_exact':strings(h),'readout':readout,'compiled_quadrature':compiled,'repeated_records':rows,'twirled_controls':controls,'member':{'hardware':'common C/E source dilation plus declared coherent M6/load/export services','program':'two-hop H2 phase then weighted one-proposal native C/E repair, committed recovery, timeout as idle','logical_layout':'a=2*s+n; one binary degree of freedom with three spectator values','program_equal_to_uniform_one_load_spatial_walk':False},'resource_contract':{'source_dimension':512,'executed_invariant_rows':36,'two_blocks_without_repreparation':True,'measurement_backaction_retained':True,'all_classical_history_expanded':False,'phase_and_continuing_schedule_supplied':True,'error_scope':'Exact declared source/reference marginal; actual complete stored primitives require their additional gain/error budget. No physical hardware accuracy inferred.','record_storage':'Fresh selector/coin/raw record per interval; no overwritten archive. Source conditioning and apparatus reset work are declared resources.','physical_work_energy_scale_supplied':False},'scope':{'three_spectator_values_are_three_lepton_species':False,'physical_GeV_or_seconds_derived':False,'self_starting_or_sustained_clock':False,'literal_OpenWave_row3_complete':False,'new_Nature_validation':False}}
