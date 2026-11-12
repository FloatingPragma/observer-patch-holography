"""Joint Ramsey clocks, wave arrivals and certified common-source responses."""
from fractions import Fraction as F
from itertools import product
import numpy as np
import sympy as sy
from mpmath import iv
from model.openwave_geometry_refinement.source import compression,initial,upper,ends,native
from model.openwave_geometry_refinement.propagation import advance,step_error,accumulated,probability_error
from model.openwave_relational_clock.clock import phase_bounds,floats,cycles,overlap_error
# a,c are the smallest/largest positive Jordan-Wigner normal frequencies.
# The middle-frequency spectator is unoccupied in the prepared clock state.
LEVELS=(0,3,1,5)
X=np.array([[1,1],[1,-1]],complex)/np.sqrt(2)
Y=np.array([[1,-1j],[1,1j]],complex)/np.sqrt(2)
TRANSFORMS=tuple(np.kron(a,b)for a,b in product((X,Y),repeat=2))

def native_clock_factorization():
 from model.openwave_native_force.force import source_matrices
 H=source_matrices()[2];I=sy.eye(2);sx=sy.Matrix([[0,1],[1,0]]);syy=sy.Matrix([[0,-sy.I],[sy.I,0]]);sz=sy.diag(1,-1)
 majorana=[]
 for j in range(3):
  for end in(sx,syy):majorana.append(sy.kronecker_product(*([sz]*j+[end]+[I]*(2-j))))
 assert all(majorana[a]*majorana[b]+majorana[b]*majorana[a]==2*int(a==b)*sy.eye(8)for a,b in product(range(6),repeat=2))
 A=sy.zeros(6)
 for j in range(3):A[2*j,2*j+1]=1;A[2*j+1,2*j]=-1
 for j in range(2):A[2*j+1,2*j+2]=sy.Rational(2,7);A[2*j+2,2*j+1]=-sy.Rational(2,7)
 reconstructed=3*sy.eye(8)+sy.I*sum((A[a,b]*majorana[a]*majorana[b]for a,b in product(range(6),repeat=2)),sy.zeros(8))/4
 assert H==reconstructed
 x=sy.Symbol('x');polynomial=sy.factor((sy.I*A).charpoly(x).as_expr());levels=native()[2]
 # The exact quadratic identity gives additive normal occupations; the
 # disjoint native spectral boxes fix their sorted energy labels.
 a,b,c=levels[1],levels[2],levels[3]
 expected=[(F(0),F(0)),(a.lo,a.hi),(b.lo,b.hi),(c.lo,c.hi),
  (a.lo+b.lo,a.hi+b.hi),(a.lo+c.lo,a.hi+c.hi),(b.lo+c.lo,b.hi+c.hi),(a.lo+b.lo+c.lo,a.hi+b.hi+c.hi)]
 assert all(expected[j][1]<expected[j+1][0]for j in range(7))
 assert all(max(lo,e.lo)<=min(hi,e.hi)for(lo,hi),e in zip(expected,levels))
 return {'exact_native_matrix':'D/2-(X X I+I X X)/7','exact_Majorana_CAR':True,'exact_quadratic_reconstruction':True,
  'positive_normal_frequency_characteristic_polynomial':str(polynomial),'native_energy_basis_labels_for_abits_cbits_spectator0':list(LEVELS),
  'spectator1_labels':[2,6,4,7],'bare_clock_gap_intervals':[list(map(str,(a.lo,a.hi))),list(map(str,(c.lo,c.hi)))],
  'meaning':'The exact finite Jordan-Wigner transform makes the8-level source three additive normal-occupation bits. Both prepared clocks occupy the same native source cell.'}

def joint_instrument():
 effects=[]
 for a,b in product((sy.Matrix([[1,1],[1,-1]])/sy.sqrt(2),sy.Matrix([[1,-sy.I],[1,sy.I]])/sy.sqrt(2)),repeat=2):
  U=sy.kronecker_product(a,b)
  for spectator,labels in enumerate((LEVELS,(2,6,4,7))):
   for k in range(4):
    v=sy.zeros(8,1)
    for j,label in enumerate(labels):v[label]=sy.conjugate(U[k,j])
    effects.append(v*v.conjugate().T/4)
 assert sum(effects,sy.zeros(8))==sy.eye(8)
 return {'native_dimension':8,'retained_independent_setting_coin_values':4,'complete_native_outcomes':32,'exact_effect_sum':'I8',
  'independent_initial_clocks':'|+>_a tensor |+>_c; spectator b in0','terminal_field_geometry_and_clock_effects_commute':True,
  'clock_and_geometry_evolution_not_freely_reordered':True,'unknown_quantum_states_not_copied':True}

def outputs(v):
 # v has four equally weighted source branches; the initial amplitude is1/2.
 # Each of four setting pairs has probability1/4, so amplitudes are Uv/4.
 p=np.stack([abs(U@v)**2/16 for U in TRANSFORMS]);success=float(sum(np.vdot(x,x).real for x in v)/4)
 assert abs(p.sum()-success)<1e-10
 shape=v.reshape(2,2,-1);z=[sum(np.vdot(shape[0,b],shape[1,b])for b in range(2))/2,
                         sum(np.vdot(shape[a,0],shape[a,1])for a in range(2))/2]
 return p,success,z

def exact_npt_witness(v,error):
 # An exact dyadic witness for the reduced two-clock state. All entries of
 # the accepted numerical state vector are binary rationals.
 rho=v@v.conj().T/4;pt=rho.reshape(2,2,2,2).transpose(0,3,2,1).reshape(4,4);e,U=np.linalg.eigh(pt)
 wr=np.rint(U[:,0].real*2**24).astype(int);wi=np.rint(U[:,0].imag*2**24).astype(int)
 w=[(F(int(a)),F(int(b)))for a,b in zip(wr,wi)];norm=sum(a*a+b*b for a,b in w)
 def mul(x,y):return x[0]*y[0]-x[1]*y[1],x[0]*y[1]+x[1]*y[0]
 def conj(x):return x[0],-x[1]
 entries=[[(F(float(z.real)),F(float(z.imag)))for z in row]for row in v]
 value=F(0)
 for i,j in product(range(4),repeat=2):
  a,b=divmod(i,2);c,d=divmod(j,2);k,l=2*a+d,2*c+b
  weight=mul(conj(w[i]),w[j]);inner=(F(0),F(0))
  for x,y in zip(entries[k],entries[l]):
   term=mul(x,conj(y));inner=(inner[0]+term[0],inner[1]+term[1])
  value+=mul(weight,inner)[0]/(4*norm)
 law=probability_error(error);bounds=(value-law,value+law)
 assert bounds[1]<0
 return {'integer_witness_real':wr.tolist(),'integer_witness_imag':wi.tolist(),'exact_dyadic_partial_transpose_expectation':str(value),
  'source_expectation_interval':floats(bounds),'native_clock_entanglement_certified':True,
  'physical_entanglement_or_equivalence_principle_test_claimed':False}

def execute(background):
 wave=[compression(background,j,4,True)for j in LEVELS];vacuum=[compression(background,j,4,False)for j in LEVELS]
 # Columns are independent free-wave and stopped-arrival experiments.
 w=[np.column_stack((initial(),initial()))for _ in LEVELS];v=[initial(wave=False)for _ in LEVELS]
 ew=[step_error(x)for x in wave];ev=[step_error(x)for x in vacuum];mask=np.tile(np.repeat(np.arange(7)>=5,7),25)
 hits=[];controls=[];single=[];allz=[];arrivals=[];raw=[[],[]]
 for event in range(1,17):
  w=[advance(h,x)for h,x in zip(wave,w)];v=[advance(h,x)for h,x in zip(vacuum,v)]
  hit=np.stack([np.where(mask,x[:,1],0)for x in w]);p,success,z=outputs(hit);arrivals.append(success);allz.append(z)
  for x in w:x[mask,1]=0
  dv=max(accumulated(x,event)for x in ev);dw=max(accumulated(x,event)for x in ew)
  vp,_,vz=outputs(np.stack(v));ce=overlap_error(dv,dv)
  for i in range(2):raw[i].append(vz[i])
  controls.append({'event':event,'all_settings_signs_geometry_number_probabilities':vp.tolist(),
   'clock_coherences_real_imag':[[z.real,z.imag]for z in vz],
   'clock_phase_intervals':[floats(phase_bounds(z,ce))for z in vz],
   'clock_overlap_error_upper':str(ce),'independent_same_source_vacuum_preparation':True})
  hits.append({'event':event,'all_settings_clock_sign_pair_probabilities':(p.sum(axis=2)/2).tolist(),
   'unconditional_receiver_probability':success/2,'clock_coherences_real_imag':[[z.real,z.imag]for z in z],
   'field_and_geometry_postmeasurement_quantum_factors_retained':True,'native_spectator1_outcomes_each_probability':'0',
   'complete_history_state_error_upper':str(dw)})
  if event in(8,16):
   free=np.stack([x[:,0]for x in w]);wp,_,wz=outputs(free)
   before=sum(np.vdot(x,h['H']@x).real for x,h in zip(free,wave))/4
   diagonal=sum(np.diag(h['H'])for h in wave)/4
   after=float(wp.sum(axis=(0,1))@diagonal)
   single.append({'event':event,'all_settings_clock_signs_geometry_wave_site_probabilities':(wp.reshape(4,4,25,49)/2).tolist(),
    'all_vacuum_settings_clock_signs_geometry_probabilities':(vp/2).tolist(),
    'wave_clock_coherences_real_imag':[[z.real,z.imag]for z in wz],
    'wave_clock_phase_intervals':[floats(phase_bounds(z,overlap_error(dw,dw)))for z in wz],
    'unconditional_wave_readout_mean_work':(after-before)/2,'work_numeric_error_upper':float(128*probability_error(dw)),
    'vacuum_joint_clock_NPT_witness':exact_npt_witness(np.stack(v),dv)})
 timeout,tp,tz=outputs(np.stack([x[:,1]for x in w]));vp,_,_=outputs(np.stack(v));assert abs(sum(arrivals)+tp-1)<1e-10
 # Lift each measured clock phase from the same product preparation.
 calibrated=[[],[]]
 for clock in range(2):
  lifted=np.unwrap(-np.angle(raw[clock]));last=(F(0),F(0))
  for phase,z,row in zip(lifted,raw[clock],controls):
   lo,hi=phase_bounds(z,F(row['clock_overlap_error_upper']));k=int(round((phase+np.angle(z))/(2*np.pi)));plo,phi=ends(2*iv.pi)
   shift=(k*plo,k*phi)if k>=0 else(k*phi,k*plo);bounds=(-hi+shift[0],-lo+shift[1]);assert bounds[0]>last[1]and bounds[1]-last[0]<ends(iv.pi)[0]
   last=bounds;calibrated[clock].append(cycles(bounds))
  for row,bounds in zip(controls,calibrated[clock]):row.setdefault('unwrapped_clock_cycle_intervals',[]).append(floats(bounds))
 total=sum(r['unconditional_receiver_probability']for r in hits)+float(timeout.sum()/2)+float(vp.sum()/2);assert abs(total-1)<1e-10
 return {'background_native_level':background,'source_dimension':80000,'right_clock_preparation_levels':list(LEVELS),
  'vacuum_clock_controls':controls,'complete_arrival_clock_records':hits,'complete_single_time_joint_records':single,
  'complete_timeout_setting_clock_sign_probabilities':(timeout.sum(axis=2)/2).tolist(),
  'vacuum_timeout_setting_clock_sign_geometry_probabilities':(vp/2).tolist(),
  'conditional_wave_hit_probability':sum(arrivals),'complete_original_probability_total':total,
  'final_history_law_error_upper':str(probability_error(max(accumulated(e,16)for e in ew))),
  'all64_native_and_field_number_complements_retained':True},calibrated,arrivals

def weak_probe_certificate():
 # H_lambda=H_env+lambda B(epsilon_a n_a+epsilon_c n_c), with
 # B=exp(Q_right) on the identical geometry/wave carrier.
 B=F(121,100);a_min=F(4,5);frequency=F(61,50);pair_sum=F(13,4);T=F(8)
 assert ends(iv.exp(iv.mpf('0.19')))[1]<B
 g=native()[2];assert g[1].lo>a_min and g[3].hi<frequency and g[1].hi+2*g[3].hi<pair_sum
 remainder=B*B*T*T*(pair_sum*pair_sum+frequency*frequency)/2;assert remainder<600
 rows=[]
 for lam in(F(1,1000),F(1,2**20),F(1,2**32)):
  R=600*lam*lam;single=1500*lam/(1-R)
  rows.append({'probe_multiplier':str(lam),'complete_clock_coherence_disk_upper':str(R),
   'two_clocks_same_environment_normalized_phase_disagreement_upper':str(single),
   'two_backgrounds_normalized_response_disagreement_upper':str(2*single)})
 return {'source_family':'H_env+lambda exp(Q_right)(epsilon_a n_a+epsilon_c n_c); lambda1 is the fixed master on the spectator-zero preparation',
  'operator_B_upper':str(B),'event_source_parameter_upper':'8','probe_domain':'0<lambda<=1/1000',
  'all_environment_states_and_passive_quantum_references':True,
  'initial_clocks_independent_of_environment':'|+>_a tensor |+>_c, spectator b=0',
  'exact_common_first_derivative':'-arg(z_j)/(lambda epsilon_j) tends to integral_0^t Tr[rho exp(iH_env s) B exp(-iH_env s)] ds for both clocks',
  'rational_second_order_constant_upper':'600','finite_probe_bounds':rows,
  'scaled_probe_family_executed_numerically':False,'physical_proper_time_identification':False}

def run():
 factor=native_clock_factorization();instrument=joint_instrument();a,cal,p0=execute(0);b,_,p1=execute(1);responses=[]
 gaps=[tuple(map(F,x))for x in factor['bare_clock_gap_intervals']]
 for clock in range(2):
  x,y=a['vacuum_clock_controls'][-1],b['vacuum_clock_controls'][-1]
  z0=complex(*x['clock_coherences_real_imag'][clock]);z1=complex(*y['clock_coherences_real_imag'][clock]);e0=F(x['clock_overlap_error_upper']);e1=F(y['clock_overlap_error_upper'])
  er=upper(e1*(F(abs(z0))+e0)+e0*F(abs(z1))+F(1,2**35));phase=phase_bounds(z1*z0.conjugate(),er);gap=gaps[clock]
  response=(phase[0]/gap[1],phase[1]/gap[0]);assert response[0]>0;responses.append(response)
 difference=(responses[0][0]-responses[1][1],responses[0][1]-responses[1][0]);assert difference[0]>0
 ratio=(responses[0][0]/responses[1][1]-1,responses[0][1]/responses[1][0]-1)
 means=[]
 for clock in range(2):
  intervals=[]
  for record,prob in((a,p0),(b,p1)):
   law=F(record['final_history_law_error_upper']);success=F(float(sum(prob)));guard=F(1,2**35)
   lo=sum(F(float(p))*iv[0]for p,iv in zip(prob,cal[clock]));hi=sum(F(float(p))*iv[1]for p,iv in zip(prob,cal[clock]));largest=cal[clock][-1][1]
   intervals.append(((lo-largest*law-guard)/(success+law+guard),(hi+largest*law+guard)/(success-law-guard)))
  means.append(floats((intervals[1][0]-intervals[0][1],intervals[1][1]-intervals[0][0])))
 normalized=[]
 for intervals,gap in zip(means,gaps):normalized.append((F(intervals[0])/gap[1],F(intervals[1])/gap[0]))
 common=(max(r[0]for r in normalized),min(r[1]for r in normalized));assert 0<common[0]<common[1]
 return {'schema':'oph.two_native_clocks_same_quantum_source.v1','verdict':'PASS_COMPLETE_JOINT_CLOCKS_AND_FINITE_NONUNIVERSAL_RESPONSE',
  'native_clock_factorization':factor,'complete_joint_instrument':instrument,'complete_experiments':[a,b],
  'event16_bare_gap_normalized_vacuum_phase_responses':[floats(r)for r in responses],
  'event16_normalized_phase_response_difference_interval':floats(difference),
  'event16_relative_response_difference_ppm_interval':floats((ratio[0]*10**6,ratio[1]*10**6)),
  'same_arrival_law_delays_in_each_measured_clock_cycles_intervals':means,
  'bare_gap_normalized_arrival_delay_intervals':[floats(r)for r in normalized],
  'intersection_of_normalized_arrival_delay_intervals':floats(common),
  'weak_probe_common_clock_theorem':weak_probe_certificate(),
  'same_fixed_master_geometry_wave_and_native_source':True,'clock_backreaction_and_correlations_retained':True,
  'exact_finite_strength_clock_universality':False,'physical_clock_G_or_equivalence_principle_violation_identified':False,
  'scope':'Two independently prepared native normal-occupation clocks are observed jointly on one fixed quantum geometry and wave source. Their complete operational records include geometry-mediated clock correlations. The finite response difference is a property of this specified source; no laboratory time, physical G or equivalence-principle violation is inferred.'}
