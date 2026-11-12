"""Complete clock/arrival records and event calibration from one source."""
from fractions import Fraction as F
import math
import sympy as sy
import numpy as np
from mpmath import iv
from model.openwave_geometry_refinement.source import compression,initial,upper,scalar,ends
from model.openwave_geometry_refinement.propagation import advance,step_error,accumulated,probability_error
PHASES=np.array([1,1j,-1,-1j])

def phase_bounds(z,error):
 """An outward argument enclosure of a certified complex disk."""
 iv.dps=80;error=F(error);square=F(float(z.real))**2+F(float(z.imag))**2
 normlower=F(math.isqrt(square.numerator*square.denominator),square.denominator)
 assert normlower>error
 angle=iv.atan2(scalar(F(float(z.imag))),scalar(F(float(z.real))));lo,hi=ends(angle);radius=error/(normlower-error)
 return lo-radius,hi+radius

def floats(bounds):return [float(np.nextafter(float(bounds[0]),-np.inf)),float(np.nextafter(float(bounds[1]),np.inf))]
def cycles(bounds):
 plo,phi=ends(2*iv.pi);lo,hi=bounds
 values=[lo/plo,lo/phi,hi/plo,hi/phi];return min(values),max(values)
def overlap_error(d0,d1):return upper(d0+d1+d0*d1+F(1,2**35))

def phase_records(v):
 # The setting coin chooses X orY with probability1/2. Within that
 # setting the native0/1 clock is measured sharply; other native levels
 # have their own retained zero-probability outcomes on this preparation.
 values=np.array([abs(v[0]+np.conj(phi)*v[1])**2/8 for phi in PHASES]);p=(np.vdot(v[0],v[0]).real+np.vdot(v[1],v[1]).real)/2
 assert abs(values.sum()-p)<1e-10
 return values,float(p),complex(np.vdot(v[0],v[1]))

def execute(background):
 wave=[compression(background,j,4,True)for j in(0,1)];vacuum=[compression(background,j,4,False)for j in(0,1)]
 w=[initial()for _ in range(2)];v=[initial(wave=False)for _ in range(2)];ordering_residual=0.;mask=np.tile(np.repeat(np.arange(7)>=5,7),25)
 ew=[step_error(x)for x in wave];ev=[step_error(x)for x in vacuum];hits=[];references=[];Z=0j;probabilities=[];raw_ref=[]
 for event in range(1,17):
  before=[advance(source,state)for source,state in zip(wave,w)];v=[advance(source,state)for source,state in zip(vacuum,v)];hit=[np.where(mask,x,0)for x in before];w=[np.where(mask,0,x)for x in before]
  p,success,z=phase_records(hit);pv,pvac,zvac=phase_records(v);assert abs(pvac-1)<1e-10
  before_clock,_,_=phase_records(before);ordering_residual=max(ordering_residual,float(np.max(abs(p-before_clock*mask[None,:]))));assert ordering_residual<1e-14
  d=max(accumulated(x,event)for x in ew);de=overlap_error(*(accumulated(x,event)for x in ew));vr=overlap_error(*(accumulated(x,event)for x in ev));Z+=z;probabilities.append(success);raw_ref.append(zvac)
  before_energy=sum(np.vdot(x,a['H']@x).real for a,x in zip(wave,before))/2
  nohit_energy=sum(np.vdot(x,a['H']@x).real for a,x in zip(wave,w))/2
  mixed_diagonal=(np.diag(wave[0]['H'])+np.diag(wave[1]['H']))/2;hit_energy=float(p.sum(axis=0)@mixed_diagonal)
  hits.append({'event':event,'all_setting_sign_geometry_wave_site_probabilities':(p.reshape(4,25,49)/2).tolist(),
   'unconditional_hit_probability':success/2,'clock_overlap_real_imag':[z.real,z.imag],'conditional_clock_visibility':abs(z)/success,
   'clock_phase_interval_radians':floats(phase_bounds(z,de)),'clock_overlap_error_upper':str(de),
   'native_clock_readout_work_conditioned_on_initial_wave':hit_energy+nohit_energy-before_energy,
   'work_numeric_error_upper':float(128*probability_error(d)),
   'other_native_background_or_clock_complement_outcomes_each_probability':'0'})
  references.append({'event':event,'all_setting_sign_geometry_number_probabilities':pv.tolist(),
   'vacuum_clock_overlap_real_imag':[zvac.real,zvac.imag],'clock_visibility':abs(zvac),
   'clock_phase_interval_radians':floats(phase_bounds(zvac,vr)),'clock_overlap_error_upper':str(vr),
   'independent_vacuum_preparation_for_this_calibration_event':True})
 timeout,timeout_probability,z_timeout=phase_records(w);vacuum_final,_,_=phase_records(v);assert abs(sum(probabilities)+timeout_probability-1)<1e-10
 histories_error=[accumulated(x,16)for x in ew];Zerror=overlap_error(*histories_error)
 # A separately prepared vacuum clock calibrates each event. Its measured
 # phase has a unique increasing lift, certified by adjacent interval gaps.
 lifted=np.unwrap(-np.angle(raw_ref));reference_cycles=[];last=(F(0),F(0))
 for phase,record,z in zip(lifted,references,raw_ref):
  lo,hi=phase_bounds(z,F(record['clock_overlap_error_upper']));k=int(round((phase+np.angle(z))/(2*np.pi)));plo,phi=ends(2*iv.pi)
  shiftlo,shifthi=(k*plo,k*phi)if k>=0 else(k*phi,k*plo);bounds=(-hi+shiftlo,-lo+shifthi)
  assert bounds[0]>last[1] and bounds[1]-last[0]<ends(iv.pi)[0]
  last=bounds;cycle=cycles(bounds);record['unwrapped_clock_cycle_interval']=floats(cycle);reference_cycles.append(cycle)
 return {'background_native_level':background,'clock_native_levels':[0,1],
  'all_arrival_clock_records':hits,'all_timeout_setting_sign_geometry_wave_site_probabilities':(timeout.reshape(4,25,49)/2).tolist(),
  'all_vacuum_timeout_setting_sign_geometry_probabilities':(vacuum_final/2).tolist(),
  'vacuum_calibration_controls':references,'conditional_wave_hit_probability':sum(probabilities),'conditional_wave_timeout_probability':timeout_probability,
  'total_arrival_clock_overlap_real_imag':[Z.real,Z.imag],'total_arrival_clock_overlap_error_upper':str(Zerror),
  'total_arrival_clock_visibility':abs(Z)/sum(probabilities),'complete_history_probability_error_upper':str(probability_error(max(histories_error))),
  'source_dimension':80000,'setting_choice_is_retained':True,'same_event_readout_reordering_numeric_residual':ordering_residual,'all_native_and_field_complements_retained':True,
  'wave_and_vacuum_are_distinct_declared_preparations_of_the_same_source':True},reference_cycles,probabilities

def instrument_certificate():
 identity=sy.eye(8);effects=[]
 for phi in(1,sy.I):
  for sign in(1,-1):
   v=sy.zeros(8,1);v[0]=1/sy.sqrt(2);v[1]=sign*phi/sy.sqrt(2);effects.append(v*v.conjugate().T/2)
  for j in range(2,8):
   v=identity[:,j];effects.append(v*v.T/2)
 assert sum(effects,sy.zeros(8))==identity
 field=sy.diag(1,0);geometry=sy.diag(0,1)
 for effect in effects:
  C=sy.kronecker_product(effect,sy.eye(2),sy.eye(2));P=sy.kronecker_product(identity,field,geometry);assert C*P==P*C
 return {'native_clock_levels':8,'setting_and_native_outcomes':len(effects),'exact_effect_sum':'I_8',
  'exact_same_event_field_geometry_clock_effect_commutators':'0',
  'event_zero_clock_probabilities_in_active_order':['1/2','1/4','0','1/4'],
  'phase_origin':'The prepared native equal-phase state supplies event zero; every calibration uses its relative phase.',
  'clock_phase_reference_is_declared_apparatus':True,
  'unitary_dilation':'The complete finite Kraus word has V*V=I on the entire source. An orthonormal extension implements V with finite setting, event and observer-record registers.'}

def run():
 first,calibration,p0=execute(0);second,_,p1=execute(1);records=[first,second];means=[];errors=[]
 for r,p in zip(records,(p0,p1)):
  law=F(r['complete_history_probability_error_upper']);success=F(float(sum(p)));pguard=F(1,2**35)
  # The clock-coordinate weights are measured vacuum Ramsey phases from
  # the background-zero preparation, shared by both arrival experiments.
  numerator_lo=sum(F(float(x))*bounds[0]for x,bounds in zip(p,calibration));numerator_hi=sum(F(float(x))*bounds[1]for x,bounds in zip(p,calibration));maxcycle=calibration[-1][1]
  interval=((numerator_lo-maxcycle*law-pguard)/(success+law+pguard),(numerator_hi+maxcycle*law+pguard)/(success-law-pguard));r['mean_arrival_in_reference_clock_cycles_interval']=floats(interval);means.append(interval)
  z=complex(*r['total_arrival_clock_overlap_real_imag']);e=F(r['total_arrival_clock_overlap_error_upper']);bounds=phase_bounds(z,e);r['aggregate_arrival_clock_phase_interval_radians']=floats(bounds);errors.append(e)
 delta=(means[1][0]-means[0][1],means[1][1]-means[0][0]);assert delta[0]>0
 z0=complex(*first['total_arrival_clock_overlap_real_imag']);z1=complex(*second['total_arrival_clock_overlap_real_imag']);z=z1*z0.conjugate();product_error=upper(errors[1]*(F(abs(z0))+errors[0])+errors[0]*F(abs(z1))+F(1,2**35));phase=phase_bounds(z,product_error);phasecycles=cycles((-phase[1],-phase[0]))
 # A finite repetition theorem uses unconditional arrival/sign frequencies.
 # Z=4[(q_X+ -q_X-) +i(q_Y+ -q_Y-)], so eight frequency
 # errors bounded by eta imply a coherence disk of radius12eta.
 sample_count=2**40;eta=F(1,200000);se=[x+12*eta for x in errors]
 sampling_product_error=upper(se[1]*(F(abs(z0))+se[0])+se[0]*F(abs(z1))+F(1,2**35));sample_phase=phase_bounds(z,sampling_product_error);sample_cycles=cycles((-sample_phase[1],-sample_phase[0]));assert sample_cycles[0]>0
 failure=ends(16*iv.exp(-2*sample_count*scalar(eta*eta)))[1];assert failure<F(1,10**20)
 total=[]
 for r in records:
  probs=[x['unconditional_hit_probability']for x in r['all_arrival_clock_records']];timeout=float(np.array(r['all_timeout_setting_sign_geometry_wave_site_probabilities']).sum());vac=float(np.array(r['all_vacuum_timeout_setting_sign_geometry_probabilities']).sum());total.append(sum(probs)+timeout+vac)
 assert max(abs(x-1)for x in total)<1e-10
 return {'schema':'oph.native_clock_wave_event_join.v1','verdict':'PASS_OPERATIONAL_NATIVE_CLOCK_ARRIVAL_JOIN',
  'same_master_source':'model/openwave_geometry_refinement/source.py','complete_experiments':records,'exact_complete_clock_instrument':instrument_certificate(),
  'delay_in_calibrated_native_clock_cycles_interval':floats(delta),'actual_aggregate_arrival_clock_phase_response_cycles_interval':floats(phasecycles),
  'complete_probability_totals':total,
  'finite_sampling_guarantee':{'independent_prepared_trials_per_background':sample_count,'each_unconditional_frequency_tolerance':str(eta),'simultaneous_failure_probability_upper':str(upper(failure,160)),
   'positive_arrival_phase_response_cycles_interval_with_sampling':floats(sample_cycles),'large_repetition_experiment_executed':False,
   'rule':'Two independent ensembles, four arrival/sign frequencies per ensemble; Hoeffding union bound16 exp(-2N eta^2). Complete failure and timeout records count in every trial.'},
  'clock_observable':'The retained setting coin chooses X orY. Differences of their sign probabilities reconstruct the clock coherence; no simultaneous sharpX/Y outcome is assigned.',
  'calibration_rule':'Unique increasing vacuum-clock phase lift from independently prepared reference runs; the same measured event-to-cycle dictionary weights both arrival laws.',
  'event_partial_order':{'evolution':'sixteen copies of the declared finite source unitary, with retained no-hit branches',
   'at_hit':'field hit/site, geometry number and native clock-setting/sign effects act on distinct tensor factors and commute',
   'record_copying':'orthogonal setting/sign/site/event records may be copied to two observers; no unknown quantum state is copied',
   'all_topological_orders_of_each_measurement_layer_have_identical_complete_records':True,
   'evolution_and_measurement_events_are_not_freely_reordered':True},
  'phase_statistic_is_a_circular_coherence_statistic_not_the_mean_arrival_time':True,
  'physical_seconds_Einstein_proper_time_or_fundamental_global_clock_assumed':False,
  'scope':'Operational native-clock cycles and joint clock/wave observations in one finite autonomous source. The experiment word specifies causal event order. No physical proper-time calibration or relativistic continuum identification is supplied.'}
