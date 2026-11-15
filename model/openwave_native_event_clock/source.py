"""Three exact event-decay modes reconstructed from actual C/E currents."""
from fractions import Fraction as F
from collections import defaultdict
from hashlib import sha256
import json
from .field import Q
from model.openwave_native_process import source as process
from model.native_program import interpreter as native

SEAMS=tuple(native.local.SEAMS)
BASIS=tuple(sum(1<<j for j,z in enumerate(x)if z)for x in native.local.basis())
INDEX={x:i for i,x in enumerate(BASIS)}

def graph_vectors():
 adjacency=[[int((min(a,b),max(a,b))in SEAMS)for b in range(12)]for a in range(12)]
 distance=[None]*12;distance[0]=0
 for level in range(3):
  for a in range(12):
   if distance[a]==level:
    for b in range(12):
     if adjacency[a][b]and distance[b]is None:distance[b]=level+1
 assert sorted(distance)==[0]+[1]*5+[2]*5+[3]
 vectors=[]
 for name,adj,values in [('icosahedral_plus',Q(0,1),(Q(F(1,4)),Q(0,F(1,20)),Q(0,F(-1,20)),Q(F(-1,4)))),('icosahedral_minus',Q(0,-1),(Q(F(1,4)),Q(0,F(-1,20)),Q(0,F(1,20)),Q(F(-1,4)))),('icosahedral_minus_one',Q(-1),(Q(F(5,12)),Q(F(-1,12)),Q(F(-1,12)),Q(F(5,12))))]:
  v=tuple(values[d]for d in distance)
  assert sum(v)==Q() and all(sum(adjacency[a][b]*v[b]for b in range(12))==adj*v[a]for a in range(12))
  vectors.append((name,Q(5)-adj,Q(1)-(Q(5)-adj)/60,v))
 return adjacency,vectors

def population_and_record_check(vectors):
 T=[defaultdict(F)for _ in BASIS];count=0;rawcount=0
 summaries=[]
 for x in BASIS:
  i=INDEX[x]
  for edge in SEAMS:
   for branch in process.raw_local_branches():
    # Original uniform selector has ONE global zero hold, not one per seam.
    if branch.record[0]=='hold':continue
    rawcount+=1;y=process.apply_branch(x,edge,branch)
    if y is None or not branch.weight:continue
    weight=branch.weight/len(SEAMS);T[i][INDEX[y]]+=weight;count+=1
    a,b=edge
    assert ((x>>a)&1)-((y>>a)&1)==branch.current
    for _,_,_,v in vectors:
     delta=sum(v[j]*(((y>>j)&1)-((x>>j)&1))for j in range(12))
     assert delta==(v[b]-v[a])*branch.current
  assert sum(T[i].values())==1
 for i in range(66):
  for j,w in T[i].items():assert T[j][i]==w
 serial=[[[j,str(w)]for j,w in sorted(row.items())]for row in T]
 return T,{'basis_columns':66,'raw_nonhold_slots_checked':rawcount,'single_global_zero_hold_per_column':True,'nonzero_transition_contributions':count,'population_kernel_sha256':sha256(json.dumps(serial,separators=(',',':')).encode()).hexdigest(),'full_record_current_equals_source_mode_increment':True}

def forward(T,p):
 out=[F(0)]*66
 for i,row in enumerate(T):
  for j,w in row.items():out[j]+=p[i]*w
 return out

def backward(T,f):return[sum(w*f[j]for j,w in row.items())for row in T]

def run():
 A,modes=graph_vectors();T,control=population_and_record_check(modes)
 initial_mask=(1<<SEAMS[0][0])|(1<<SEAMS[0][1]);start=INDEX[initial_mask]
 rows=[]
 for name,mu,lam,v in modes:
  f=[sum(v[j]for j in range(12)if x>>j&1)for x in BASIS]
  assert all(z==lam*y for z,y in zip(backward(T,f),f))
  p=[F(int(i==start))for i in range(66)];means=[]; distributions=[]
  for ordinal in range(4):
   mean=sum(p[i]*w*(f[j]-f[i])for i,row in enumerate(T)for j,w in row.items())
   assert mean==(lam-Q(1))*(lam**ordinal)*f[start]
   means.append(mean);p=forward(T,p);distributions.append([str(x)for x in p])
  assert all(means[i+1]/means[i]==lam for i in range(3))
  V=sum(z*z for z in f)/66;expected=Q(F(5,33))*sum(z*z for z in v);assert V==expected
  inc2=sum(w*(f[j]-f[i])**2 for i,row in enumerate(T)for j,w in row.items())/66
  assert inc2==2*V*(Q(1)-lam)
  B=backward(T,f);B=[x-y for x,y in zip(B,f)];lag=[]
  for n in range(1,4):
   cov=sum(w*(f[j]-f[i])*B[j] for i,row in enumerate(T)for j,w in row.items())/66
   assert cov==-V*(Q(1)-lam)**2*(lam**(n-1));lag.append(cov);B=backward(T,B)
  rows.append({'id':name,'laplacian_eigenvalue_exact':mu.pair(),'event_multiplier_exact':lam.pair(),'native_vertex_mode_exact':[z.pair()for z in v],
   'four_actual_record_means_exact':[z.pair()for z in means],'recovered_event_multiplier_exact':(means[1]/means[0]).pair(),
   'stationary_mode_variance_exact':V.pair(),'stationary_current_variance_exact':inc2.pair(),'record_lag_covariances_exact':[z.pair()for z in lag],
   'four_complete_population_distributions_exact':distributions})
 # The original whole-history evaluator's finite host stop is reported.
 command={'op':'repair','selector':{'kind':'uniform'},'law':'mixture','repeat':4}
 program=process.admitted_local_program(commands=[command],capacity=2)
 original=native.execute(program,{(start,0,start,0):F(1)},max_steps=1)
 assert original.status=='RESOURCE_STOP'and original.attempted_steps==1
 exact_one={}
 for block in original.branches:
  for(i,a,j,b),weight in block.density.items():
   assert i==j and a==b==0;exact_one[i]=exact_one.get(i,F(0))+weight
 assert exact_one==dict(T[start])
 # Bound12 mean estimates from independent four-event preparations. The
 # raw current mode increment lies in[-1,1]; this is a conditional sample plan.
 from model.openwave_collective_finite_strength.bounds import exp_bounds
 repetitions=2**24;epsilon=F(1,1000);exponent=F(repetitions)*epsilon**2/2
 failure=24*exp_bounds(-exponent/16)[1]**16;assert failure<F(1,100)
 return{'schema':'oph.native_current_event_clock.v1','verdict':'PASS_NATIVE_CURRENT_EVENT_MODES_AND_COMPLETE_MOMENT_CONTRACTION',
  'quadratic_field_encoding':'[a,b] denotes a+b*sqrt(5), with exact rational a,b','native_source':{'shell':'H66','ports':12,'seams':list(map(list,SEAMS)),'uniform_seam_probability_exact':'1/30','C_E_weights_exact':['1/2','1/2'],'endpoint_coin_weights_exact':['1/2','1/2'],'repair_probability_exact':'1','source_parameters_fitted_to_particle_masses':False},
  'initial_occupation_mask':initial_mask,'mode_records':rows,'complete_basis_current_control':control,
  'event_clock':{'admitted_program':program,'two_recording_blocks':[[0,1],[2,3]],'native_events':4,'raw_slots_per_event':721,'complete_raw_record_words':721**4,
   'full_history_kernel':'Product of four existing complete native instruments with record-addressed outputs; all zero raw slots remain. Exact moment contraction traces only the records not requested by each moment, with no success normalization.',
   'exact_four_event_moment_contraction_executed':True,'expanded_full_history_tree_executed':False,'original_prefix_status':original.status,'original_prefix_attempts':original.attempted_steps,'original_prefix_raw_blocks':len(original.branches),'original_one_event_distribution_matches':True},
  'sampling_plan':{'independent_preparations':repetitions,'native_events_per_preparation':4,'relative_assessment':'absolute record-mean error, not relative decay-frequency error','absolute_mean_error_exact':str(epsilon),'simultaneous_mean_count':12,'simultaneous_failure_upper_exact':str(failure),'independent_initialization_and_record_resolution_supplied':True},
  'scope':{'actual_native_current_records':True,'population_and_record_moments_exact':True,'three_decay_modes_are_particle_families':False,'two_oscillator_cycles_executed':False,'physical_mass_phase_or_GeV_assignment':False,'absolute_physical_time_derived':False,'new_nature_match':False}}
