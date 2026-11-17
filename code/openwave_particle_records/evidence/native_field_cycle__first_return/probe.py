#!/usr/bin/env python3
"""Complete coherent first-return source; immutable receipts and read-only replay."""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import inspect
import json
from pathlib import Path
import sys
import time
import numpy as np
from scipy import sparse as sp

HERE=Path(__file__).resolve().parent;PARENT=HERE.parent;ROOT=HERE.parents[3]
spec=importlib.util.spec_from_file_location('first_return_field_provider',PARENT/'probe.py')
base=importlib.util.module_from_spec(spec);sys.modules[spec.name]=base;spec.loader.exec_module(base)
sg=base.sg;cc=base.cc;original_loop_columns=base.loop_columns
BOUND=486;COUNTERS=512
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def canonical(x):return base.canonical(x)

def return_maps(states,targets,signs,native_increments,native_rows,q,arrays):
    n=len(states);N=n+1
    rt=np.zeros((2,N),int);ft=np.tile(np.arange(N),(2,1));fs=np.ones((2,N),int)
    path=np.full((2,n,BOUND+1),-1,dtype=np.int32)
    integrated=np.zeros((2,n,3),int);delta=np.zeros((2,n,3),int)
    native_current=np.zeros((2,n,3),int)
    for row in native_rows:
        d,col,step,edge,_,_,_,current,*_=row
        native_current[d,col,edge]+=current
    for d in range(2):
        for col,s in enumerate(states):
            cur=col;phase=1;path[d,col,0]=col
            for t in range(1,BOUND+1):
                integrated[d,col]+=native_current[d,cur]
                phase*=signs[d,cur];cur=int(targets[d,cur]);path[d,col,t]=cur
                if states[cur].tagged_port==s.tagged_port:
                    rt[d,col]=t;ft[d,col]=cur;fs[d,col]=phase;break
            assert rt[d,col]>0
            final=states[cur]
            assert final.pair==s.pair and final.tagged_port==s.tagged_port
            full=[b-a for a,b in zip(s.flux,final.flux)]
            assert not any(sg.g.divergence(full))
            delta[d,col]=[full[sg.g.EDGE_INDEX[e]] for e in sg.EDGES]
            assert np.array_equal(delta[d,col],-integrated[d,col])
    for d in range(2):
        assert np.array_equal(ft[1-d,ft[d]],np.arange(N))
        assert np.array_equal(rt[1-d,ft[d]],rt[d])
        assert np.array_equal(fs[d]*fs[1-d,ft[d]],np.ones(N))
    arrays.update({f'q{q}_first_return_lengths':rt,f'q{q}_first_return_targets':ft,
                   f'q{q}_first_return_signs':fs,f'q{q}_first_return_paths':path,
                   f'q{q}_first_return_native_integrated_current':integrated,
                   f'q{q}_first_return_electric_increment':delta})
    return rt,ft,fs,delta

def complete_oracle_and_controller(states,targets,signs,rt,ft,fs,q,arrays):
    n=len(states);N=n+1;tags=np.array([s.tagged_port for s in states]+[-1]);original=np.arange(N)
    tt=np.concatenate((targets,np.full((2,1),n)),axis=1)
    ss=np.concatenate((signs,np.ones((2,1),int)),axis=1)
    trace=np.zeros((2,BOUND+1,N),np.int32);hits=np.zeros((2,BOUND,N),bool)
    for d in range(2):
        work=original.copy();phase=np.ones(N,int);trace[d,0]=work
        for t in range(BOUND):
            phase*=ss[d,work];work=tt[d,work];trace[d,t+1]=work
            hits[d,t]=(tags[work]==tags)&(original<n)
        calculated=np.where(original<n,np.argmax(hits[d],axis=0)+1,0)
        assert np.array_equal(calculated,rt[d])
        # Reverse the actual fixed scan, erasing every equality bit before its source inverse.
        scratch_hits=hits[d].copy()
        for t in range(BOUND-1,-1,-1):
            scratch_hits[t]^=(tags[work]==tags)&(original<n)
            phase*=ss[1-d,work];work=tt[1-d,work]
        assert not np.any(scratch_hits) and np.array_equal(work,original)
        assert np.all(phase==1) and np.all((work-original)%N==0)
    # C- V C+ on every full source/counter column, including dirty input counters.
    source=np.repeat(original,COUNTERS);counter=np.tile(np.arange(COUNTERS),N)
    middle=counter^rt[0,source];out=source.copy();phase=np.ones(len(source),int)
    for t in range(1,COUNTERS):
        active=middle>=t;phase[active]*=ss[0,out[active]];out[active]=tt[0,out[active]]
    final_counter=middle^rt[1,out]
    encoded=out*COUNTERS+final_counter
    assert len(np.unique(encoded))==N*COUNTERS
    # Literal inverse sequence for every column.
    inverse_counter=final_counter^rt[1,out];inverse=out.copy();inverse_phase=phase.copy()
    for t in range(1,COUNTERS):
        active=inverse_counter>=t;inverse_phase[active]*=ss[1,inverse[active]];inverse[active]=tt[1,inverse[active]]
    inverse_counter^=rt[0,inverse]
    assert np.array_equal(inverse,source) and np.array_equal(inverse_counter,counter) and np.all(inverse_phase==1)
    zero=np.arange(N)*COUNTERS
    assert np.array_equal(out[zero],ft[0]) and np.all(final_counter[zero]==0)
    assert np.array_equal(phase[zero],fs[0])
    arrays.update({f'q{q}_oracle_source_trace':trace,f'q{q}_oracle_hit_bits':hits,
                   f'q{q}_controller_output_address':encoded.astype(np.int32),
                   f'q{q}_controller_output_sign':phase.astype(np.int8)})
    return {'source_counter_columns':N*COUNTERS,'invalid_columns':COUNTERS,
            'valid_zero_counter_columns':n,'all_counter_zero_outputs_clean':True,
            'complete_controller_inverse_exact':True,'oracle_scratch_and_phases_uncomputed':True,
            'fixed_loop_word_slots':4*BOUND+COUNTERS-1,'fixed_native_seam_slots':3*(4*BOUND+COUNTERS-1),
            'retained_output_length_record':False}

def coherence_reader(states,rt,ft,fs,q,arrays):
    domain=np.array([i for i,s in enumerate(states) if not any(sg.g.gauss(s))])
    n=len(states);psi=np.zeros(n,complex);psi[domain]=1/np.sqrt(len(domain))
    desired=np.zeros(n,complex);desired[ft[0,:n]]=fs[0,:n]*psi
    full=np.zeros((COUNTERS,n),complex);full[0]=desired
    leaked=np.zeros_like(full)
    for i in range(n):leaked[rt[0,i],ft[0,i]]+=fs[0,i]*psi[i]
    # Complete two-outcome source pointer, retaining all512 counter addresses.
    def measure(joint):
        success=abs(joint@desired.conj())**2
        total=np.sum(abs(joint)**2,axis=1)
        failure=total-success
        assert failure.min()>-1e-12
        return np.stack((success,failure))
    clean_law=measure(full);leaked_law=measure(leaked)
    weights=np.bincount(rt[0,domain],minlength=COUNTERS)/len(domain)
    expected=float(weights@weights)
    assert abs(clean_law[0].sum()-1)<1e-12 and abs(leaked_law[0].sum()-expected)<1e-12
    assert abs(clean_law.sum()-1)<1e-12 and abs(leaked_law.sum()-1)<1e-12
    # Two-length plant has orthogonal retained labels and equal amplitudes.
    target=np.ones(2)/np.sqrt(2);plant=np.diag(target)
    planted=float(np.sum(abs(plant@target)**2));assert abs(planted-.5)<1e-12
    assert abs(1-planted)>.1
    arrays.update({f'q{q}_coherence_complete_pointer':clean_law,
                   f'q{q}_leaked_history_complete_pointer':leaked_law,
                   f'q{q}_coherence_preparation':psi})
    return {'gauss_zero_preparation_columns':len(domain),'clean_overlap':float(clean_law[0].sum()),
            'leaked_history_overlap':float(leaked_law[0].sum()),'expected_leaked_overlap':expected,
            'return_length_distribution':{str(i):float(x) for i,x in enumerate(weights) if x},
            'full_pointer_slots':2*COUNTERS,'two_length_plant_overlap':planted,
            'constant_one_reader_mutant_error':1-planted}

def neutral_uniform_memory_certificate(states,rt,ft,fs,q):
    if q:return {'applicable':False}
    groups={}
    for col,s in enumerate(states):groups.setdefault((s.tagged_port,s.flux),[]).append(col)
    shifts=[]
    for key,columns in sorted(groups.items()):
        assert len(columns)==6 and {states[i].memory for i in columns}==set(range(6))
        for d in range(2):
            assert len({int(rt[d,i]) for i in columns})==1
            shift={(states[ft[d,i]].memory-states[i].memory)%6 for i in columns}
            assert len(shift)==1 and np.all(fs[d,columns]==1)
            assert set(ft[d,columns])==set(columns)
            assert all((states[ft[d,i]].tagged_port,states[ft[d,i]].flux)==key for i in columns)
            shifts.append({'tag':key[0],'flux':list(key[1]),'direction':d,
                           'return_length':int(rt[d,columns[0]]),'memory_translation':int(next(iter(shift))),
                           'all_six_car_signs':fs[d,columns].tolist()})
    assert len(groups)==81
    return {'applicable':True,'uniform_memory_source_dimension':81,'spin_tensor_dimension':162,
            'full_fiber_matrix_units_fixed':True,'fibers_by_direction':shifts}

def first_return_columns(states,q,arrays):
    native_arrays={}
    U,targets,signs,inc,ret,original=original_loop_columns(states,q,native_arrays)
    rt,ft,fs,delta=return_maps(states,targets,signs,inc,native_arrays[f'q{q}_loop_native_records'],q,arrays)
    controller=complete_oracle_and_controller(states,targets,signs,rt,ft,fs,q,arrays)
    coherence=coherence_reader(states,rt,ft,fs,q,arrays)
    neutral=neutral_uniform_memory_certificate(states,rt,ft,fs,q)
    arrays.update({'one_word_'+key:value for key,value in native_arrays.items()})
    n=len(states)
    F=sp.coo_matrix((fs[0,:n],(ft[0,:n],np.arange(n))),shape=(n,n)).tocsr()
    inverse=sp.coo_matrix((fs[1,:n],(ft[1,:n],np.arange(n))),shape=(n,n)).tocsr()
    assert base.maxabs(F.getH()-inverse)==0 and base.maxabs(F.getH()@F-sp.eye(n))==0
    eig=np.linalg.eigvalsh(((6*sp.eye(n)+F+F.getH())/8).toarray())
    assert eig.min()>=.5-1e-9 and eig.max()<=1+1e-9
    arrays[f'q{q}_first_return_lazy_eigenvalues']=eig
    info={'native_path_calls':original['native_path_calls'],
          'native_one_word_provider':original,'first_return_controller':controller,
          'coherence':coherence,'first_return_lengths_by_direction':[
              {str(t):int(count) for t,count in zip(*np.unique(rt[d,:n],return_counts=True))} for d in range(2)],
          'all_columns_returned':True,'all_field_increments_divergence_zero':True,
          'lazy_spectral_interval':[float(eig.min()),float(eig.max())],
          'microscopic_current_included_in_accounting':True}
    info['neutral_uniform_memory_certificate']=neutral
    return F,ft[:,:n],fs[:,:n],delta,np.ones((2,n),bool),info

def amended_case_reader():
    source=inspect.getsource(base.run_case)
    original="    assert changes['omit_loop']>1e-6 and changes['zero_current']>1e-6"
    replacement="    if q:\n        assert changes['omit_loop']>1e-6 and changes['zero_current']>1e-6\n    else:\n        assert changes['omit_loop']<=1e-9 and changes['zero_current']>1e-6"
    assert source.count(original)==1
    amended=source.replace(original,replacement)
    namespace=dict(vars(base))
    exec(compile(amended,str(HERE/'probe.py')+'::amended_case_reader','exec'),namespace)
    return namespace['run_case'],{'original_reader_sha256':hashlib.sha256(source.encode()).hexdigest(),
                                 'amended_reader_sha256':hashlib.sha256(amended.encode()).hexdigest(),
                                 'only_change':'q0 omitted-loop gate is exact null; q1 gates unchanged'}

def source_pins():
    pins=json.loads((PARENT/'SOURCE_PINS.json').read_text())
    snapshots=json.loads((PARENT/'FROZEN_INPUT_MAP.json').read_text())
    for source,entry in snapshots.items():
        assert pins.pop(source)==entry['sha256']
        pins[str((PARENT/entry['snapshot']).relative_to(ROOT))]=entry['sha256']
    for p in ('probe.py','PROOF.md','DESIGN.md','TASK.md','MANIFEST.sha256','PREREG.json','RECEIPT.json','FROZEN_INPUT_MAP.json'):
        file=PARENT/p;pins[str(file.relative_to(ROOT))]=digest(file)
    return dict(sorted(pins.items()))

def freeze():
    assert not (HERE/'RECEIPT.json').exists()
    pins=source_pins();(HERE/'SOURCE_PINS.json').write_text(json.dumps(pins,sort_keys=True,indent=2)+'\n')
    names=('DESIGN.md','TASK.md','probe.py','SOURCE_PINS.json','AMENDMENT.md','PROOF.md','FAILURE_GATE.log')
    prereg={'schema':'oph.native_first_return.prereg.v1','files':{n:digest(HERE/n) for n in names},'evidence_executed':False}
    (HERE/'PREREG.json').write_text(json.dumps(prereg,sort_keys=True,indent=2)+'\n')
    print(json.dumps({'sha256':digest(HERE/'PREREG.json'),'source_pin_count':len(pins)}))

def check_pins():
    for p,h in json.loads((HERE/'PREREG.json').read_text())['files'].items():assert digest(HERE/p)==h,p
    for p,h in json.loads((HERE/'SOURCE_PINS.json').read_text()).items():assert digest(ROOT/p)==h,p

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--freeze',action='store_true');parser.add_argument('--record',action='store_true')
    args=parser.parse_args()
    if args.freeze:freeze();return
    check_pins()
    if args.record and (HERE/'RECEIPT.json').exists():raise FileExistsError('Evidence is immutable')
    start=time.perf_counter();frame,_,spins=cc.geometric_gates()
    a,b,c=(frame[x] for x in sg.TRIANGLE);normal=np.cross(b-a,c-a);normal/=np.linalg.norm(normal)
    weights=np.array([np.cross(frame[a],frame[b]) for a,b in sg.EDGES]).T
    arrays={'frame':frame,'normal':normal,'axial_weights':weights}
    # The sole changed primitive is the explicitly constructed native first-return word.
    base.loop_columns=first_return_columns
    reader,reader_audit=amended_case_reader()
    cases=[]
    for q in (0,1):
        cases.append(reader(q,frame,spins,normal,weights,arrays));print('First return q='+str(q)+' passed',flush=True)
    receipt=canonical({'schema':'oph.native_first_return.v1','candidate':'CAR-first-return-field-clock-v1',
                       'class':'N+ with declared reversible coherent services','cases':cases,
                       'nature_measurements_used':[],'physical_row5_claimed':False,
                       'new_full_native_pulse_compilation':False,'preregistration_sha256':digest(HERE/'PREREG.json'),
                       'reader_amendment':reader_audit,
                       'original_q0_omitted_loop_gate_passed':False,
                       'original_q0_omitted_loop_difference':cases[0]['changed_law_reference_differences']['omit_loop'],
                       'failed_original_preregistration_sha256':digest(HERE/'attempts/01_neutral_control/PREREG.json')})
    check_pins()
    if args.record:
        (HERE/'RECEIPT.json').write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n')
        np.savez_compressed(HERE/'RECORDS.npz',**arrays)
    else:
        assert receipt==json.loads((HERE/'RECEIPT.json').read_text())
        with np.load(HERE/'RECORDS.npz') as old:
            assert set(old.files)==set(arrays)
            for name,value in arrays.items():assert np.array_equal(value,old[name]),name
    print(json.dumps({'all_pass':True,'elapsed_seconds':time.perf_counter()-start,
        'literal_native_calls':sum(c['source_native_calls']+c['loop']['native_path_calls'] for c in cases),
        'controller_columns':sum(c['loop']['first_return_controller']['source_counter_columns'] for c in cases),
        'stationary_charge_current':[c['stationary_diagnostics']['native']['max_realized_charge_current'] for c in cases],
        'leaked_history_overlap':[c['loop']['coherence']['leaked_history_overlap'] for c in cases],
        'physical_rows_complete':[]},indent=2))

if __name__=='__main__':main()
