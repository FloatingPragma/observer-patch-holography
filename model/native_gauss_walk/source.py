"""Charged coherent walk of the source, tag and integer electric links.

Exact amplitudes use integer Q(zeta) coordinates with one common denominator.
The stored native event matrices are independently applied to the same full
sparse carrier; their padding and normalization defects are retained.
"""
from dataclasses import dataclass
from collections import defaultdict,deque
from fractions import Fraction as F
from math import gcd,lcm,sqrt
from functools import lru_cache
import base64,hashlib,json,zlib
import numpy as np
from model.native_electric_gauss import source as electric
from model.common_source.record_code import MATCHINGS
from model.native_observer_program import source as program

SEAMS=electric.SEAMS
EVENTS=60
PADDED=64
LINK_LIMIT=4
NORM=720
ZERO=(0,0)
ONE=(1,0)


def add(a,b):return(a[0]+b[0],a[1]+b[1])
def scale(c,a):return(c*a[0],c*a[1])
def norm(a):return a[0]*a[0]+a[0]*a[1]+a[1]*a[1]
def rotate(a,p):
    x,y=a;p%=6
    return((x,y),(-y,x+y),(-x-y,x),(-x,-y),(y,-x-y),(x+y,-x))[p]

def add_arc(out,key,value):
    result=add(out.get(key,ZERO),value)
    if result==ZERO:out.pop(key,None)
    else:out[key]=result


# A configuration is (pair,tagged_port,flux tuple); q and active character
# are explicit wave-level labels and never inferred from a particle name.
@dataclass
class Wave:
    amplitudes:dict
    denominator:int=1
    charge:int=1
    character:int=1
    dressed:bool=True

    def __post_init__(self):
        if type(self.denominator)is not int or self.denominator<=0 or type(self.charge)is not int or type(self.character)is not int or not 0<=self.character<6 or type(self.dressed)is not bool:
            raise ValueError('exact carrier labels and a positive denominator required')
        for(c,e),a in self.amplitudes.items():
            if type(e)is not int or not 0<=e<PADDED or len(a)!=2 or any(type(x)is not int for x in a):
                raise ValueError('complete padded event and exact cyclotomic amplitude required')
            electric.validate(state(c,self.charge),LINK_LIMIT)

    def probability(self):return F(sum(norm(v)for v in self.amplitudes.values()),NORM*self.denominator**2)

    def reduce(self):
        factor=self.denominator
        for a in self.amplitudes.values():factor=gcd(factor,gcd(abs(a[0]),abs(a[1])))
        if factor>1:
            self.denominator//=factor
            self.amplitudes={k:(a[0]//factor,a[1]//factor)for k,a in self.amplitudes.items()}
        return self


def config(s):return(s.pair,s.tagged_port,s.flux)
def state(c,q):return electric.ElectricState(c[0],c[1],q,0,c[2])


@lru_cache(None)
def shortest_path(pair):
    origin,target=pair
    adjacent={v:[]for v in range(12)}
    for a,b in SEAMS:adjacent[a].append(b);adjacent[b].append(a)
    parents={origin:None};queue=deque([origin])
    while queue:
        a=queue.popleft()
        for b in sorted(adjacent[a]):
            if b not in parents:parents[b]=a;queue.append(b)
    path=[target]
    while path[-1]!=origin:path.append(parents[path[-1]])
    return tuple(reversed(path))


def canonical_config(pair,q=1):return config(electric.pair_string(shortest_path(pair),q))


def prepare(q=1,character=None):
    if type(q)is not int or q not in(0,1):raise ValueError('primary q1 or neutral q0 control required')
    character=q%6 if character is None else character
    if type(character)is not int or not 0<=character<6:raise ValueError('one declared memory character required')
    pairs=tuple(MATCHINGS[0]+MATCHINGS[1]);assert len(set(pairs))==12
    return Wave({(canonical_config(pair,q),e):ONE for pair in pairs for e in range(EVENTS)},1,q,character)


@lru_cache(maxsize=200000)
def event_column(c,event,q,character,dressed=True):
    if not 0<=event<PADDED:raise ValueError('complete64 event register required')
    if event>=EVENTS:return c,0,0,False
    seam=SEAMS[event//2];s=state(c,q)
    hop=electric.apply_native_current(s,seam,event%2,LINK_LIMIT if dressed else None)
    target=hop.state
    if not dressed:target=electric.ElectricState(target.pair,target.tagged_port,q,target.memory,s.flux)
    return config(target),(character*target.memory)%6,hop.realized_signed_current,hop.overflow


def reflect(w):
    groups=defaultdict(dict)
    for(c,e),a in w.amplitudes.items():groups[c][e]=a
    out={}
    for c,row in groups.items():
        total=ZERO
        for e,a in row.items():
            if e<EVENTS:total=add(total,a)
        for e in range(EVENTS):add_arc(out,(c,e),add(scale(2,total),scale(-EVENTS,row.get(e,ZERO))))
        for e in range(EVENTS,PADDED):
            if e in row:add_arc(out,(c,e),scale(-EVENTS,row[e]))
    return Wave(out,w.denominator*EVENTS,w.charge,w.character,w.dressed).reduce()


def swap(w):
    out={}
    for(c,e),a in w.amplitudes.items():
        target,p,_,_=event_column(c,e,w.charge,w.character,w.dressed)
        add_arc(out,(target,e),rotate(a,p))
    return Wave(out,w.denominator,w.charge,w.character,w.dressed)


def phase(w,inverse=False):
    sign=1 if inverse else-1
    return Wave({(c,e):rotate(a,sign*sum(v*v for v in c[2]))for(c,e),a in w.amplitudes.items()},w.denominator,w.charge,w.character,w.dressed)


def step(w,inverse=False,electric_phase=False):
    if inverse:
        if electric_phase:w=phase(w,True)
        return reflect(swap(w))
    result=swap(reflect(w))
    return phase(result)if electric_phase else result


def equivalent(a,b):
    if(a.charge,a.character,a.dressed)!=(b.charge,b.character,b.dressed):return False
    keys=set(a.amplitudes)|set(b.amplitudes)
    return all(scale(b.denominator,a.amplitudes.get(k,ZERO))==scale(a.denominator,b.amplitudes.get(k,ZERO))for k in keys)


def decode_canonical_dressing(w):
    """D* plus its orthogonal complement, never a source/link reset."""
    admitted=defaultdict(lambda:ZERO);outside=0
    for(c,e),a in w.amplitudes.items():
        if c==canonical_config(c[0],w.charge):admitted[(c[0],e)]=add(admitted[(c[0],e)],a)
        else:outside+=norm(a)
    denominator=NORM*w.denominator**2
    return {'canonical_dressing_probability':str(F(sum(norm(a)for a in admitted.values()),denominator)),
            'complete_outside_dressing_probability':str(F(outside,denominator)),
            'reset_performed':False}


@lru_cache(None)
def source_energies():
    from model.common_source import parameter_selection as selection
    costs=selection.diagonal_cost(2)['values'];states=selection.code_geometry()['states']
    result={}
    for loads,cost in zip(states,costs,strict=True):
        pair=tuple(i for i,n in enumerate(loads)if n)
        if len(pair)==2:result[pair]=F(int(cost.p),int(cost.q))
    assert set(result)==set(electric.PAIRS)
    return result


def complete_records(w):
    denominator=NORM*w.denominator**2
    pair={p:0 for p in electric.PAIRS};event=[0]*PADDED;current=defaultdict(int)
    gauss=defaultdict(int);field=0;source_energy=0;overflow=0;idle=0;rows=[]
    for(c,e),a in sorted(w.amplitudes.items()):
        weight=norm(a);pair[c[0]]+=weight;event[e]+=weight
        s=state(c,w.charge);g=electric.gauss(s);gauss[g]+=weight
        field+=F(weight)*electric.electric_energy(s);source_energy+=weight*source_energies()[c[0]]
        # S is an involution. Its reverse column has the opposite realized
        # forward current. This is a terminal spectral readout of the last
        # recoded event, not a secretly committed intermediate record.
        _,_,backward,failed=event_column(c,e,w.charge,w.character,w.dressed)
        current[-backward]+=weight;overflow+=weight*failed
        idle+=weight*int(e>=EVENTS or e%2==0)
        rows.append([list(c[0]),c[1],list(c[2]),e,list(a),weight])
    raw=json.dumps(rows,separators=(',',':')).encode()
    assert sum(pair.values())==sum(event)==sum(current.values())==sum(gauss.values())
    return {'norm_squared':str(w.probability()),'nonzero_joint_arcs':len(rows),
      'nonzero_source_tag_flux_configurations':len({c for c,e in w.amplitudes}),
      'joint_probability_denominator':denominator,
      'complete_joint_table_encoding':'zlib+base64 JSON rows[pair,tag,flux30,event,integer_Qzeta_amplitude,norm_numerator]; omitted columns have exact zero',
      'complete_joint_table_zlib_base64':base64.b64encode(zlib.compress(raw,9)).decode(),
      'complete_joint_table_sha256':hashlib.sha256(raw).hexdigest(),
      'source_pair_probabilities':[{'pair':list(p),'probability':str(F(v,denominator))}for p,v in sorted(pair.items())],
      'all64_event_probabilities':[str(F(v,denominator))for v in event],
      'signed_current_probabilities':{str(k):str(F(v,denominator))for k,v in sorted(current.items())},
      'gauss_sector_probabilities':[{'G':list(k),'probability':str(F(v,denominator))}for k,v in sorted(gauss.items())],
      'electric_energy_mean':str(field/denominator),'source_H2_energy_mean':str(source_energy/denominator),
      'source_plus_electric_energy_mean':str((source_energy+field)/denominator),'overflow_probability':str(F(overflow,denominator)),
      'idle_or_padding_event_probability':str(F(idle,denominator)),
      'event_padding_probability':str(F(sum(event[EVENTS:]),denominator)),
      'active_memory_character':w.character,'memory_coordinate_distribution':'each full joint arc has six equal coordinate weights in its retained character; no character projection or normalization',
      'canonical_dressing_readout':decode_canonical_dressing(w)}


def unpack_records(record):
    raw=zlib.decompress(base64.b64decode(record['complete_joint_table_zlib_base64']))
    assert hashlib.sha256(raw).hexdigest()==record['complete_joint_table_sha256']
    return json.loads(raw)


def exact_protocol(q=1,character=None,dressed=True,electric_phase=False):
    initial=prepare(q,character);initial.dressed=dressed;states=[initial]
    for _ in range(2):
        moved=step(states[-1],electric_phase=electric_phase)
        assert moved.probability()==1
        states.append(moved)
    records=[complete_records(x)for x in states]
    inverses=[]
    for n in(1,2):
        value=states[n]
        for _ in range(n):value=step(value,inverse=True,electric_phase=electric_phase)
        assert equivalent(value,initial)
        inverses.append({'steps':n,'full_state_exactly_restored':True,'decoded':decode_canonical_dressing(value),
                         'source_link_tag_event_reset':False})
    if dressed:
        assert all(r['gauss_sector_probabilities']==[{'G':[0]*12,'probability':'1'}]for r in records)
    return {'charge_number':q,'active_character':initial.character,'flux_dressing':dressed,
      'electric_phase':electric_phase,'phase_duration':'16*pi/3 source-time units'if electric_phase else None,
      'records':records,'inverse_controls':inverses,
      'electric_work_per_step':[str(F(records[n+1]['electric_energy_mean'])-F(records[n]['electric_energy_mean']))for n in(0,1)],
      'source_plus_electric_work_per_step':[str(F(records[n+1]['source_plus_electric_energy_mean'])-F(records[n]['source_plus_electric_energy_mean']))for n in(0,1)],
      'work_scope':'common declared H2 plus electric energy; tracking and walk supplied by apparatus, not an autonomous energy-conservation claim'},states


def complex_value(a):return complex(a[0]+a[1]/2,sqrt(3)*a[1]/2)


def stored_prepare(q=1):
    U,_=program.compiled_event_matrices();out={}
    for pair in MATCHINGS[0]+MATCHINGS[1]:out[canonical_config(pair,q)]=U[:,0].copy()/sqrt(12)
    return out


def stored_reflect(value,inverse=False):
    R=program.compiled_inverse_reflection() if inverse else program.compiled_event_matrices()[1]
    keys=list(value);matrix=np.stack([value[k]for k in keys],axis=1)
    result=R@matrix
    return{k:result[:,i].copy()for i,k in enumerate(keys)}


def stored_swap(value,q=1,character=1):
    out={}
    for c,row in value.items():
        for e,a in enumerate(row):
            if a==0:continue
            target,p,_,_=event_column(c,e,q,character,True)
            if target not in out:out[target]=np.zeros(PADDED,dtype=complex)
            out[target][e]+=complex_value(rotate(ONE,p))*a
    return out


def stored_phase(value,inverse=False):
    sign=1 if inverse else-1
    return{c:complex_value(rotate(ONE,sign*sum(x*x for x in c[2])))*v for c,v in value.items()}


def stored_step(value,inverse=False,electric_phase=False):
    if inverse:
        if electric_phase:value=stored_phase(value,True)
        return stored_reflect(stored_swap(value),True)
    value=stored_swap(stored_reflect(value))
    return stored_phase(value)if electric_phase else value


def stored_error(value,exact):
    error=0.;norm_squared=0.;padding=0.;gauss_failure=0.;energy=0.;outside_dressing=0.
    for c in set(value)|{x for x,e in exact.amplitudes}:
        actual=value.get(c,np.zeros(PADDED,complex));target=np.zeros(PADDED,complex)
        for e in range(PADDED):target[e]=complex_value(exact.amplitudes.get((c,e),ZERO))/(sqrt(NORM)*exact.denominator)
        error+=float(np.vdot(actual-target,actual-target).real)
        weight=float(np.vdot(actual,actual).real);norm_squared+=weight;padding+=float(np.vdot(actual[EVENTS:],actual[EVENTS:]).real)
        gauss_failure+=weight*bool(any(electric.gauss(state(c,1))))
        energy+=weight*float(electric.electric_energy(state(c,1)))
        if c!=canonical_config(c[0],1):outside_dressing+=weight
    return {'vector_error_diagnostic':sqrt(error),'unnormalized_norm_squared':norm_squared,
      'all_event_padding_mass_diagnostic':padding,'Gauss_failure_mass':gauss_failure,
      'electric_energy_diagnostic':energy,'configuration_count':len(value),
      'outside_canonical_dressing_mass_diagnostic':outside_dressing,
      'normalization_or_amplitude_pruning_performed':False}


def stored_protocol(exact_states,electric_phase=False):
    value=stored_prepare();rows=[stored_error(value,exact_states[0])]
    for n in(1,2):value=stored_step(value,electric_phase=electric_phase);rows.append(stored_error(value,exact_states[n]))
    for _ in range(2):value=stored_step(value,inverse=True,electric_phase=electric_phase)
    inverse=stored_error(value,exact_states[0])
    from pathlib import Path
    parent=json.loads(Path(program.__file__).with_name('results.json').read_text())
    contract=next(r for r in parent['event_reflection_compilation']if r['valid_events']==60)
    prep=F(contract['preparation_bound']['full_word_operator_error_upper'])
    ref=F(contract['reflection_bound']['full_word_operator_error_upper'])
    # These are worst-case full-reference operator bounds from the executed
    # source word, separate from the observed floating vector residuals.
    operator_bounds=[(1+prep)*(1+ref)**n-1 for n in range(3)]
    for row,bound in zip(rows,operator_bounds,strict=True):
        assert row['vector_error_diagnostic']<float(bound)
        assert row['Gauss_failure_mass']==0
        row['inherited_operator_error_upper']=str(bound)
    bound_inverse=(1+prep)*(1+ref)**4-1
    assert inverse['vector_error_diagnostic']<float(bound_inverse)
    d0=F(contract['preparation_bound']['complete_halfdiamond_upper']);dr=F(contract['reflection_bound']['complete_halfdiamond_upper'])
    g0=F(contract['preparation_bound']['complete_map_gain_upper']);gr=F(contract['reflection_bound']['complete_map_gain_upper'])
    f0=F(contract['preparation_bound']['all_first_processor_failure_mass_upper']);fr=F(contract['reflection_bound']['all_first_processor_failure_mass_upper'])
    composed=[{'event_reflections':n,'complete_halfdiamond_upper':str(d0*gr**n+dr*sum((gr**j for j in range(n)),F(0))),
               'all_first_processor_failure_mass_upper':str(f0+g0*fr*sum((gr**j for j in range(n)),F(0)))}for n in(0,1,2,3,4)]
    return {'stored_forward_records':rows,'stored_inverse_record':inverse,'complete_word_composition':composed,
      'inverse_operator_error_upper':str(bound_inverse),
      'source_word_bound':{'preparation':contract['preparation_bound'],'reflection':contract['reflection_bound']},
      'charged_link_event':'exact declared integer/active-memory permutation; no inferred physical hardware precision',
      'source_and_character_preparation':'same declared matching0+1 and character1 services as the common program',
      'complete_processor_failures':'inherited native-word bound retained separately; success component is never renormalized'}


def run():
    primary,states=exact_protocol();phased,phase_states=exact_protocol(electric_phase=True)
    neutral,_=exact_protocol(q=0);undressed,_=exact_protocol(dressed=False)
    assert undressed['records'][1]['gauss_sector_probabilities']!=[{'G':[0]*12,'probability':'1'}]
    assert primary['records'][2]['complete_joint_table_sha256']!=phased['records'][2]['complete_joint_table_sha256']
    p=primary['records'][2];v=phased['records'][2]
    observable_difference=sum(abs(F(a['probability'])-F(b['probability']))for a,b in zip(p['source_pair_probabilities'],v['source_pair_probabilities'],strict=True))/2
    event_difference=sum(abs(F(a)-F(b))for a,b in zip(p['all64_event_probabilities'],v['all64_event_probabilities'],strict=True))/2
    assert observable_difference==0 and event_difference==0
    third=step(states[2]);third_phased=step(phase_states[2],electric_phase=True)
    third_record=complete_records(third);third_phase_record=complete_records(third_phased)
    third_tv=sum(abs(F(a['probability'])-F(b['probability']))for a,b in zip(third_record['source_pair_probabilities'],third_phase_record['source_pair_probabilities'],strict=True))/2
    assert third_tv>0 and third.probability()==1 and third_phased.probability()==1
    assert third_record['gauss_sector_probabilities']==third_phase_record['gauss_sector_probabilities']==[{'G':[0]*12,'probability':'1'}]
    stored=stored_protocol(states);stored_phase_result=stored_protocol(phase_states,True)
    delta3=F(next(x for x in stored['complete_word_composition']if x['event_reflections']==3)['complete_halfdiamond_upper'])
    native_lower=third_tv-2*delta3
    assert native_lower>0
    return {'schema':'oph.joint_native_charged_integer_field_walk.v1',
      'verdict':'PASS: joint native charged source and integer field evolution with complete records',
      'common_program_sha256':program.descriptor_sha256(),'integer_extension':electric.extension_descriptor(),
      'preparation':{'matching_labels':[0,1],'source_pairs':[list(p)for p in MATCHINGS[0]+MATCHINGS[1]],
       'canonical_paths':[list(shortest_path(p))for p in MATCHINGS[0]+MATCHINGS[1]],
       'uniform_events':EVENTS,'padded_events':PADDED,'link_limit':LINK_LIMIT,
       'tagged_port':'smaller input endpoint, then transported with its load','charge_number':1,'active_character':1,
       'canonical_dressing':'declared lexicographic shortest native strings; no empirical selection'},
      'primary_joint_walk':primary,'electric_phase_joint_walk':phased,
      'observed_electric_phase_response':{'two_step_source_pair_total_variation':str(observable_difference),'two_step_event_total_variation':str(event_difference),
       'three_step_source_pair_total_variation':str(third_tv),'three_step_response_lower_with_native_event_word_bound':str(native_lower),
       'bound_scope':'executed source-word event/preparation gates, exact declared signed-link permutation and electric phase services; no physical hardware precision inferred',
       'three_step_primary':third_record,'three_step_with_phase':third_phase_record,
       'implementation':'exact full-state third-step interference probe; stored native matrices executed through2steps andtheirfullinverse'},
      'neutral_control':neutral,'undressed_negative_control':undressed,
      'stored_native_primary':stored,'stored_native_with_electric_phase':stored_phase_result,
      'scope':{'same_evolved_source_tag_links_event_memory':True,'all_quantum_references':'exact event involutions and reflection tensor identity; complete source-word bounds',
       'freely_evolving_local_clock_proved':False,'physical_photon_or_electron_identified':False,
       'undressed_source_spectrum_transferred':False,'reset_in_inverse_test':False,'OpenWave_awards':0}}
