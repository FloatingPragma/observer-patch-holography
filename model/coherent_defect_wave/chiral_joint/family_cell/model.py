"""Exact full-family local action and its complete vacuum-connected graph."""
from collections import defaultdict,deque
from functools import lru_cache
from fractions import Fraction
from math import factorial,isqrt
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.sparse import coo_matrix
from ..periodic_family.model import yukawa_words


@lru_cache(maxsize=1)
def terms():
    by_scalar=defaultdict(list)
    for channel,row,col,word,coefficient in yukawa_words():
        assert coefficient.is_Integer
        by_scalar[col].append((channel,tuple(word),row,int(coefficient)))
    return dict(by_scalar)


def column(key,channels=('up','down','charged_lepton')):
    mask,scalar=key
    if type(mask) is not int or not 0<=mask<2**30 or scalar not in range(5):
        raise ValueError('complete F30 times scalar5 basis required')
    out=defaultdict(int)
    if scalar:out[key]=1
    for channel,word,row,value in terms().get(scalar,()):
        if channel not in channels:continue
        m=mask
        for mode,creation in reversed(word):
            if bool((m>>mode)&1)==creation:break
            if (m&((1<<mode)-1)).bit_count()%2:value=-value
            m^=1<<mode
        else:out[m,row]+=value
    return {key:v for key,v in out.items() if v}


def apply(state):
    out=defaultdict(int)
    for key,value in state.items():
        for target,weight in column(key).items():out[target]+=value*weight
    return {key:v for key,v in out.items() if v}


@lru_cache(maxsize=1)
def graph():
    seen={(0,0)};front=seen.copy();distances={(0,0):0};depth=0
    while front:
        depth+=1
        more={target for key in front for target in column(key)}-seen
        for key in more:distances[key]=depth
        seen.update(more);front=more
    keys=tuple(sorted(seen));index={key:j for j,key in enumerate(keys)}
    rows=[];cols=[];values=[]
    for j,key in enumerate(keys):
        for target,value in column(key).items():
            assert target in index
            rows.append(index[target]);cols.append(j);values.append(value)
    h=coo_matrix((np.array(values,dtype=np.int64),(rows,cols)),shape=(len(keys),len(keys))).tocsr()
    h.sort_indices();assert (h-h.T).nnz==0
    phase={index[(0,0)]:1};todo=deque(phase)
    while todo:
        j=todo.popleft()
        for offset in range(h.indptr[j],h.indptr[j+1]):
            i=int(h.indices[offset]);value=int(h.data[offset])
            if i==j:continue
            desired=-phase[j]*(1 if value>0 else -1)
            if i in phase:assert phase[i]==desired,'frustrated graph cannot use this certificate'
            else:phase[i]=desired;todo.append(i)
    assert len(phase)==len(keys)
    phases=np.array([phase[j] for j in range(len(keys))],dtype=np.int64)
    payload=json.dumps({'keys':keys,'rows':h.indptr.tolist(),'cols':h.indices.tolist(),'values':h.data.tolist()},separators=(',',':')).encode()
    digest=hashlib.sha256(payload).hexdigest()
    return keys,h,phases,max(distances.values()),digest


def exact_ground():
    from .quotient import spectrum
    result=spectrum()
    keys,h,phase,depth,digest=graph()
    return {'full_local_dimension':5*2**30,'preparation_sector_dimension':len(keys),
        'Hamiltonian_nonzero_entries':h.nnz,'vacuum_graph_radius':depth,
        'source_graph_sha256':digest,'full_source_action_leaves_graph_invariant':True,
        'switchable_stoquastic_and_connected':True,
        'ground_in_declared_preparation_sector_is_unique':True,
        'ground_energy_interval':result['ground_energy_interval'],
        'ground_preparation_probability_lower_bound':result['ground_filter_success_interval'][0],
        'ground_preparation_probability_interval':result['ground_filter_success_interval'],
        'ground_preparation_failure_interval':result['ground_filter_failure_interval'],
        'ground_preparation_failure_retained':True,
        'ground_is_full_source_Gauss_invariant':True,
        'global_ground_across_other_conserved_sectors_proved':False,
        'source_energy_is_particle_rest_mass':False}


def preparation_records():
    """Exact finite-time family-production records, with a full-source remainder."""
    powers=[{(0,0):1}]
    for _ in range(5):powers.append(apply(powers[-1]))
    t=Fraction(1,1000);real=defaultdict(Fraction);imag=defaultdict(Fraction)
    for n in range(5):
        target=real if n%2==0 else imag
        coefficient=(-1 if n%4 in (1,2) else 1)*t**n/factorial(n)
        for key,value in powers[n].items():target[key]+=coefficient*value
    fourth_remainder_norm2=sum(v*v for v in powers[5].values())
    scale=2**64;n=isqrt(fourth_remainder_norm2*scale*scale)
    if n*n<fourth_remainder_norm2*scale*scale:n+=1
    remainder=Fraction(n,scale)*t**5/factorial(5)
    # |<p,Ep>-<u,Eu>| <= delta*(2+delta), ||u||=1.
    error=remainder*(2+remainder)
    blocks=((0,1,4),(2,3,5,6,7,8),(9,),(10,11),(12,13,14))
    menus={'up':(1,1,0,0,0),'down':(0,1,0,0,1),'charged_lepton':(0,0,1,1,0)}
    support=set(real)|set(imag)
    def counts(mask):
        return tuple(sum((mask>>(15*spin+i))&1 for spin in (0,1) for i in block) for block in blocks)
    # Complete irreducible-block number sectors commute with every gauge
    # current; individual color-coordinate projectors would not suffice.
    groups={name:{key for key in support if key[1]!=0 and counts(key[0])==menu}
            for name,menu in menus.items()}
    assert sum(map(len,groups.values()))==len(set().union(*groups.values()))
    values={}
    for channel,keys in groups.items():
        center=sum(real[k]**2+imag[k]**2 for k in keys)
        values[channel]={'full_detector_rank':32 if channel=='charged_lepton' else 288,
            'native_irreducible_occupation_counts':list(menus[channel]),'probability_interval':
            [str(max(Fraction(0),center-error)),str(min(Fraction(1),center+error))]}
        assert Fraction(values[channel]['probability_interval'][0])>0
    total=sum(sum(real[k]**2+imag[k]**2 for k in keys) for keys in groups.values())
    values['all_other_outputs']={'probability_interval':
        [str(max(Fraction(0),1-total-3*error)),str(min(Fraction(1),1-total+3*error))]}
    return {'input':'literal full-family fermion/scalar vacuum, all Yukawas and scalar penalty one',
        'source_time':str(t),'Taylor_degree':4,'full_source_H5_input_squared_norm':fourth_remainder_norm2,
        'state_error_upper_bound':str(remainder),'complete_grouped_records':values,
        'outcomes_are_gauge_number_records_not_observed_particle_labels':True,
        'energy_or_ground_reference_substituted':False}
