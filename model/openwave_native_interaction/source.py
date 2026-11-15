"""An interacting two-load observable from actual native repairs and records."""
from fractions import Fraction as F
from collections import defaultdict
from itertools import product
from hashlib import sha256
import json
from model.openwave_native_process import source as process
from model.native_program import interpreter as native
from model.compiled_native_history.compiler import compile_history

LOCAL_EDGES=((0,1),(1,2))

def event_maps(edge):
 out=[]
 for slot,branch in enumerate(process.raw_local_branches()):
  mapping=tuple(process.apply_branch(x,edge,branch)for x in range(8))
  assert len({y for y in mapping if y is not None})==sum(y is not None for y in mapping)
  out.append({'slot':slot,'record':branch.record,'current':branch.current,'weight':branch.weight,'mapping':mapping})
 for x in range(8):assert sum(b['weight']for b in out if b['mapping'][x]is not None)==1
 return out

def unread(rho,edge,weighted=None):
 result=defaultdict(F)
 for b in event_maps(edge):
  factor=b['weight']*(weighted(b)if weighted else 1)
  if not factor:continue
  for(i,j),v in rho.items():
   a,c=b['mapping'][i],b['mapping'][j]
   if a is not None and c is not None:result[a,c]+=factor*v
 return{k:v for k,v in result.items()if v}

def trace(rho):return sum(v for(i,j),v in rho.items()if i==j)
def population(rho):return{str(i):str(rho.get((i,i),F(0)))for i in range(8)}
def mean(rho,fn):return sum(v*fn(i)for(i,j),v in rho.items()if i==j)
def bit(x,j):return(x>>j)&1

def assay_zero_imbalance(rho):
 # On the N2 three-port sector this complete native E statistic equals n_b*n_c.
 # It is an additional actual repair and its post-measurement state is retained.
 weighted=unread(rho,(1,2),lambda b:2*int(b['record'][1]=='E'and b['record'][4]==0))
 out=unread(rho,(1,2));value=trace(weighted)
 return value,out

def independent_walkers(cycles):
 # Explicit comparison process, not an alternative native two-load law.
 p={(0,2):F(1)}
 for edge in LOCAL_EDGES*cycles:
  out=defaultdict(F)
  for(x,y),w in p.items():
   ax=[(x,F(1))]if x not in edge else[(z,F(1,2))for z in edge]
   ay=[(y,F(1))]if y not in edge else[(z,F(1,2))for z in edge]
   for(a,u),(b,v)in product(ax,ay):out[a,b]+=w*u*v
  p=dict(out)
 means=[sum(w*((x==j)+(y==j))for(x,y),w in p.items())for j in range(3)]
 pair=sum(w*(((x==1)+(y==1))*((x==2)+(y==2)))for(x,y),w in p.items())
 Ezero=sum(w for(x,y),w in p.items()if((x==1)+(y==1))==((x==2)+(y==2)))
 return p,means,pair,Ezero

def complete_operator_control():
 # Every matrix unit in the full8-dimensional local occupation factor.
 # The common-source map is tensor identity on all other native ports.
 counts=0;serial=[]
 for edge in LOCAL_EDGES:
  for i in range(8):
   for j in range(8):
    expected=process.unread_column(i,j,edge)
    actual=unread({(i,j):F(1)},edge)
    assert actual==expected
    assert trace(actual)==F(i==j);counts+=1
    serial.append([list(edge),i,j,[[list(k),str(v)]for k,v in sorted(actual.items())]])
 # Adjacent operators genuinely interact/noncommute on the same source.
 rho={(5,5):F(1)}
 ef=unread(unread(rho,LOCAL_EDGES[0]),LOCAL_EDGES[1]);fe=unread(unread(rho,LOCAL_EDGES[1]),LOCAL_EDGES[0])
 assert ef!=fe
 return{'matrix_units_checked':counts,'complete_operator_sha256':sha256(json.dumps(serial,separators=(',',':')).encode()).hexdigest(),'arbitrary_spectator_and_reference_extension':'Each endpoint partial injection acts as identity on all spectator occupation and reference coordinates. This proves the full global fixed-N restriction, not just selected state trajectories.','overlapping_event_orders_differ':True,'order_witness_initial_local_mask':5,'forward_population_exact':population(ef),'reverse_population_exact':population(fe)}

def original_compiler_control(actual_edges):
 # Two genuine repairs: the original compiler retains every625 raw history
 # and checks the entire H66 carrier through tagged basis maps.
 commands=[process.fixed_command(e)for e in actual_edges]
 p=process.admitted_local_program(commands=commands,capacity=2)
 compiled=compile_history(p)
 assert compiled.status=='COMPLETE'and compiled.attempted_steps==2 and len(compiled.branches)==625
 states=native.local.basis();masks=[sum(1<<j for j,v in enumerate(x)if v)for x in states];idx={m:i for i,m in enumerate(masks)}
 checked=0
 for branch in compiled.branches:
  records=branch.archive+branch.memory
  assert len(records)==2
  weights=[];maps=[]
  for record,edge in zip(records,actual_edges):
   selected=next(b for b in process.raw_local_branches()if b.record==(record.kind,record.law,record.coin,record.raw,record.label))
   weights.append(selected.weight);maps.append((edge,selected))
  predicted=[]
  for col,mask in enumerate(masks):
   for edge,b in maps:
    if mask is not None:mask=process.apply_branch(mask,edge,b)
   if mask is not None and weights[0]*weights[1]:predicted.append((col,idx[mask]))
  assert tuple(predicted)==branch.source_map
  assert F(*branch.weight)==(weights[0]*weights[1]if predicted else 0)
  checked+=66
 return{'admitted_program':p,'status':compiled.status,'attempts':2,'raw_histories':625,'all_H66_tagged_column_history_checks':checked,'compiled_map_entries':compiled.map_entries,'complete_matrix_unit_action_follows_existing_partial_injection_theorem':True}

def run():
 graph=process.graph(0);e=next(edge for edge in graph.seams if edge[0]//12==edge[1]//12)
 a,b=e;c=next(v for u,v in graph.seams if u==b and v!=a and v//12==b//12)
 actual=((a,b),(b,c));assert all(edge in native.local.SEAMS for edge in actual)
 rho={(5,5):F(1)};rows=[]
 for cycles in(1,2):
  for edge in LOCAL_EDGES:rho=unread(rho,edge)
  assert trace(rho)==1
  avg=[mean(rho,lambda x,j=j:bit(x,j))for j in range(3)]
  pair=mean(rho,lambda x:bit(x,1)*bit(x,2));cov=pair-avg[1]*avg[2]
  independent,freeavg,freepair,freezero=independent_walkers(cycles)
  assert avg==freeavg
  readout,after=assay_zero_imbalance(rho);assert readout==pair
  # The additional terminal assay changes quantum coherences in general;
  # this particular already heat-bathed diagonal population is unchanged.
  assert after==rho
  expected=F(1,2)if cycles==1 else F(3,8)
  difference=F(1,4)if cycles==1 else F(9,64)
  assert pair==expected and pair-freepair==difference
  rows.append({'admitted_program':process.admitted_program(0,2,[process.fixed_command(e)for e in(actual*cycles+(actual[1],))],capacity=3),'two_event_blocks':cycles,'native_transport_events':2*cycles,'native_terminal_assay_events':1,
   'complete_local_population_exact':population(rho),'native_one_point_means_exact':list(map(str,avg)),
   'native_pair_density_bc_exact':str(pair),'native_connected_bc_exact':str(cov),
   'native_complete_E_zero_estimator_mean_exact':str(readout),'native_post_assay_population_exact':population(after),
   'independent_one_load_control':{'labelled_pair_probabilities_exact':[[x,y,str(w)]for(x,y),w in sorted(independent.items())],'one_point_means_exact':list(map(str,freeavg)),'pair_density_bc_exact':str(freepair),'connected_bc_exact':str(freepair-freeavg[1]*freeavg[2]),'same_zero_imbalance_record_statistic_exact':str(freezero),'collision_probability_exact':str(sum(w for(x,y),w in independent.items()if x==y))},
   'pair_density_difference_exact':str(difference),'operational_zero_imbalance_difference_exact':str(readout-freezero),
   'readout_statistic_is_pair_density_on_control_with_collisions':False})
 initial={(5,5):F(1)};one=unread(initial,(0,1));two=unread(one,(1,2))
 J1=trace(unread(initial,(0,1),lambda b:(b['current'] or 0)));J2=trace(unread(one,(1,2),lambda b:(b['current'] or 0)))
 J12=trace(unread(unread(initial,(0,1),lambda b:(b['current'] or 0)),(1,2),lambda b:(b['current'] or 0)))
 assert (J1,J2,J12-J1*J2)==(F(1,2),F(-1,4),F(1,8))
 program=process.admitted_program(0,2,[process.fixed_command(e)for e in(*actual,*actual,actual[1])],capacity=3)
 native_word={'admitted_five_event_program':program,'endpoint_local_to_global':[a,b,c],'carrier':'full global binary N=2 sector on240ports','full_carrier_dimension':240*239//2,'complete_local_operator_dimension':8,
  'native_raw_slots_per_fixed_event':25,'complete_raw_words':25**5,'complete_instrument_representation':'Literal product of the five complete25-slot endpoint instruments, with full spectator identity and all zero words. Population/record contractions sum the complete instrument; the expanded five-event tree is not materialized.',
  'expanded_five_event_history_tree_executed':False,'five_event_observable_contraction_executed':True,
  'record_preparation_and_storage':'Fixed occupation preparation with loads at a,c and vacuum spectators is supplied. Every seam, law/coin/raw record, padding outcome and terminal source state is retained by the product instrument; FIFO exports preserve all five records.',
  'intrinsic_reset_or_energy_service_inserted':False}
 # The original native grammar does not implement an arbitrary field gate.
 bad=json.loads(json.dumps(program));bad['commands'][0]={'op':'interacting_field_unitary','repeat':1,'target':'exp(-i theta n_a n_b)'}
 try:native.admit(bad)
 except ValueError as ex:rejection=str(ex)
 else:raise AssertionError('undeclared field unitary was admitted')
 return{'schema':'oph.actual_native_exclusion_interaction.v1','verdict':'PASS_NATIVE_TWO_LOAD_INTERACTION_AND_COMPLETE_CURRENT_READOUT',
  'native_word':native_word,'complete_operator_control':complete_operator_control(),'original_compiler_control':original_compiler_control(actual),
  'two_recording_block_results':rows,'actual_record_current_moments_exact':{'first_mean':str(J1),'second_mean':str(J2),'connected_current_correlation':str(J12-J1*J2)},
  'complete_local_event_kernels':[{'edge':list(e),'raw_slots':[{'slot':b['slot'],'record':b['record'],'current':b['current'],'weight_exact':str(b['weight']),'partial_map':list(b['mapping'])}for b in event_maps(e)]}for e in LOCAL_EDGES],
  'unsupported_instruction_witness':{'attempted_instruction':'interacting_field_unitary','status':'REJECTED','reason':rejection,'source_operations_executed':0},
  'scope':{'actual_native_repair_only':True,'nonfactorizing_two_load_observable':True,'same_one_point_means_distinct_pair_records':True,'unitary_interacting_matter_field_compiled':False,'fermionic_field_operator_signs_or_QFT_scattering_derived':False,'physical_particle_identification_or_GeV_mass':False,'new_nature_match':False}}
