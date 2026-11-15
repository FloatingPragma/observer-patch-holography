#!/usr/bin/env python3
"""Frozen literal signed source, coherent compact preparation and native readout."""
from __future__ import annotations
from collections import Counter, defaultdict
from dataclasses import asdict
from fractions import Fraction as F
import gzip
import hashlib
import json
import math
from pathlib import Path
import random
import sys
import time

sys.set_int_max_str_digits(0)
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
RER=ROOT/'reverse-engineering-reality'
sys.path.insert(0,str(RER))
CALLS=0

def canonical(value):
    return json.dumps(value,sort_keys=True,indent=2,ensure_ascii=False)+'\n'

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def check_pins():
    for file in ('PREREG.txt','SOURCE_PINS.json'):
        data=json.loads((HERE/file).read_text())
        for path,sha in (data['files'] if file=='PREREG.txt' else data).items():
            assert digest(ROOT/path)==sha,path
    return digest(HERE/'PREREG.txt')

def imports():
    global np,native,source,compiler,pointer,configuration
    import numpy as np
    from model.observer import native
    original=native.native
    def counted(*args):
        global CALLS
        CALLS+=1
        return original(*args)
    native.native=counted
    from model.federated_observer import source
    from model.observer.source_scalar import compiler
    from model.native_interface import pointer
    from model.native_clock import configuration

def loads(p,r,V):
    q=[0]*V
    if p>=0:q[p]+=1
    if r>=0:q[r]-=1
    return tuple(q)

def coordinates(q):
    pos=[i for i,x in enumerate(q) if x==1]
    neg=[i for i,x in enumerate(q) if x==-1]
    assert all(x in (-1,0,1) for x in q) and len(pos)<=1 and len(neg)<=1
    return (pos[0] if pos else -1,neg[0] if neg else -1)

def kernel(evidence):
    graph=source.federation(0)
    V,M=graph.ports,len(graph.seams)
    assert (V,M,max(map(len,graph.neighbours)))==(240,630,6)
    seams=graph.seams
    f={4:1,5:-1,6:-1,7:1}
    g={16:1,17:-1,18:-1,19:1}
    assert not(set(f)&set(g))
    assert not any((a in f and b in g) or (b in f and a in g) for a,b in seams)
    L=np.zeros((V,V),dtype=int)
    for a,b in seams:
        L[a,a]+=1;L[b,b]+=1;L[a,b]-=1;L[b,a]-=1
    for vec in (f,g):
        x=np.array([vec.get(i,0) for i in range(V)])
        assert np.array_equal(L@x,6*x)
    endpoints=[]
    for x in (-1,0,1):
        for y in (-1,0,1):
            for c in (0,1):
                out,record=native.native((x,y),(0,1),c)
                j=(x-y-(2*c-1)*(abs(x-y)%2))//2
                assert out==(x-j,y+j) and record==(c,j)
                endpoints.append([x,y,c,j,list(out)])
    vectors={
        'positive':{(p,-1):a for p,a in f.items()},
        'negative':{(-1,r):a for r,a in g.items()},
        'pair':{(p,r):a*b for p,a in f.items() for r,b in g.items()},
    }
    output={arm:defaultdict(int) for arm in vectors}
    touched={arm:set() for arm in vectors}
    records=[]
    continuity=0
    for arm,vector in vectors.items():
        for (p,r),amp in vector.items():
            q=loads(p,r,V)
            for ei,(a,b) in enumerate(seams):
                for c in (0,1):
                    out,(_,j)=native.native(q,(a,b),c)
                    pp,rr=coordinates(out)
                    expected=list(q);expected[a]-=j;expected[b]+=j
                    assert tuple(expected)==out
                    for region in (set(f),set(g)):
                        change=sum(out[i]-q[i] for i in region)
                        assert change==-j*(int(a in region)-int(b in region))
                        continuity+=1
                    records.append([arm,p,r,ei,c,j,pp,rr])
                    output[arm][pp,rr]+=amp
                    touched[arm].add((pp,rr))
        energy=12 if arm=='pair' else 6
        wanted={key:(2*M-energy)*val for key,val in vector.items()}
        actual={key:val for key,val in output[arm].items() if val}
        assert actual==wanted,(arm,actual,wanted)
        assert output[arm].get((-1,-1),0)==0
    contacts=[]
    for ei,(a,b) in enumerate(seams):
        for p,r in ((a,b),(b,a)):
            q=loads(p,r,V)
            for c in (0,1):
                out,(_,j)=native.native(q,(a,b),c)
                assert not any(out) and j==(1 if p==a else -1)
                contacts.append(['contact',p,r,ei,c,j,-1,-1])
    allrecords=records+contacts
    with gzip.open(evidence/'branches.jsonl.gz','wt') as fobj:
        for row in allrecords:fobj.write(json.dumps(row,separators=(',',':'))+'\n')
    (evidence/'graph.json').write_text(canonical({'ports':V,'seams':seams,'f':f,'g':g,
        'endpoint_table':endpoints,'basis':'ordered signed positions; -1 means absent'}))
    # A two-cell constant contrast has source Rayleigh quotient below6.
    trial=np.zeros(V,dtype=int);trial[:12]=1;trial[12:24]=-1
    rayleigh=F(int(trial@L@trial),int(trial@trial))
    assert 0<rayleigh<6 and int(trial.sum())==0
    return {'vertices':V,'seams':M,'signed_endpoint_cases':18,
        'source_branches':len(records),'contact_branches':len(contacts),
        'regional_continuity_checks':continuity,
        'full_output_coordinates':{a:len(v) for a,v in touched.items()},
        'zero_exterior_amplitude_coordinates':{a:len(v-set(vectors[a])) for a,v in touched.items()},
        'compact_eigenvalues':{'positive':6,'negative':6,'pair':12},
        'pair_transfer_fraction':str(F(2*M-12,2*M)),
        'single_transfer_fraction':str(F(2*M-6,2*M)),
        'positive_pair_floor':str(1-F(24,2*M)),
        'lower_nonconstant_trial_rayleigh':str(rayleigh),
        'extended_constant_single_energy':0,'vacuum_output_amplitude':0}

def tensor_gate(one,bit,bits=8):
    result=np.array([[1.]],complex)
    for j in reversed(range(bits)):
        result=np.kron(result,one if j==bit else np.eye(2))
    return result

def prep_word(bit):
    return [compiler.pulse('X',math.pi/2,[bit]),
            compiler.pulse('Y',math.pi/4,[0]),compiler.pulse('Y',math.pi/4,[1])]

def actual_effects(law):
    matrix=pointer.local_matrix()
    counts=pointer.native_effects(law)
    result=[]
    for b in (0,1):
        E=np.zeros((2,2),complex)
        for z in range(4):
            block=matrix[np.ix_([z,z+4],[pointer.READY_CODE,pointer.READY_CODE+4])]
            E+=(counts[b][z]/2)*(block.conj().T@block)
        result.append(E)
    return np.array(result)

def distribution(rho,effects):
    probabilities=[]
    for label in range(256):
        effect=np.array([[1.]],complex)
        for bit in reversed(range(8)):
            effect=np.kron(effect,effects[(label>>bit)&1])
        value=np.einsum('ij,ji->',rho,effect)
        assert abs(value.imag)<1e-12
        probabilities.append(value.real)
    return np.array(probabilities)

def source_reader(evidence):
    words={name:prep_word(bit) for name,bit in (('f',2),('g',4))}
    unitaries={name:compiler.execute(word,8) for name,word in words.items()}
    inverse={name:compiler.execute(compiler.inverse(word),8) for name,word in words.items()}
    X=np.array([[0,1j],[1j,0]])
    Y=np.array([[1,1],[-1,1]])/math.sqrt(2)
    ideal={name:tensor_gate(Y,1)@tensor_gate(Y,0)@tensor_gate(X,bit)
           for name,bit in (('f',2),('g',4))}
    errors={name:float(np.linalg.norm(unitaries[name]-ideal[name],'fro')) for name in words}
    inverse_errors={name:float(np.linalg.norm(inverse[name]-ideal[name].conj().T,'fro')) for name in words}
    defects={name:float(np.linalg.norm(unitaries[name].conj().T@unitaries[name]-np.eye(256),'fro')) for name in words}
    assert max((*errors.values(),*inverse_errors.values(),*defects.values()))<=1e-11
    certificate=pointer.certificate()
    effects={law:actual_effects(law) for law in ('C','E')}
    pointer_error=0.
    blank_reset_error=0.
    witnesses=[np.array([1,0]),np.array([0,1]),np.array([1,1])/math.sqrt(2),np.array([1,1j])/math.sqrt(2)]
    for law in ('C','E'):
        for v in witnesses:
            rho=np.outer(v,v.conj()).astype(complex)
            out=pointer.read_bit(rho,law,bit=0,work_qubits=1)
            wanted=[np.trace(E@rho) for E in effects[law]]
            pointer_error=max(pointer_error,max(abs(np.trace(a)-b) for a,b in zip(out,wanted)))
        # Complete one-bit measure/conditional-source-X blank preparation.
        local=pointer.local_matrix();counts=pointer.native_effects(law)
        actual_x=compiler.execute([compiler.pulse('X',math.pi/2,[0])],1)
        superoperator=np.zeros((2,2,2,2),complex)
        for bit in (0,1):
            for z in range(4):
                block=local[np.ix_([z,z+4],[pointer.READY_CODE,pointer.READY_CODE+4])]
                if bit:block=actual_x@block
                superoperator+=(counts[bit][z]/2)*np.einsum('ac,bd->abcd',block,block.conj())
        target=np.zeros((2,2,2,2),complex)
        target[0,0,0,0]=target[0,0,1,1]=1
        blank_reset_error=max(blank_reset_error,float(np.max(abs(superoperator-target))))
    assert pointer_error<=5e-11
    assert blank_reset_error<=5e-11
    vectors={name:U[:,0] for name,U in unitaries.items()}
    rhos={name:np.outer(v,v.conj()) for name,v in vectors.items()}
    dephased={name:np.diag(np.diag(rho)) for name,rho in rhos.items()}
    plant=vectors['f'].copy();plant[4]*=-1
    plant=np.outer(plant,plant.conj())
    laws={};arrays={};summaries={}
    for law in ('C','E'):
        ordinary={name:distribution(rho,effects[law]) for name,rho in rhos.items()}
        mixed={name:distribution(rho,effects[law]) for name,rho in dephased.items()}
        for name in words:assert np.max(abs(ordinary[name]-mixed[name]))<=5e-11
        resolved={}
        for mode,inputs in (('coherent',rhos),('dephased',dephased),('sign_plant',{'f':plant,'g':rhos['g']})):
            probs={name:distribution(inverse[name]@rho@inverse[name].conj().T,effects[law])
                   for name,rho in inputs.items()}
            joint=np.outer(probs['f'],probs['g'])
            arrays[f'{law}_energy_{mode}']=joint
            resolved[mode]=float(joint[0,0])
        for mode,wanted in (('coherent',1.),('dephased',1/16),('sign_plant',1/4)):
            assert abs(resolved[mode]-wanted)<=5e-10,(law,mode,resolved[mode])
        position=np.outer(ordinary['f'],ordinary['g'])
        dephased_position=np.outer(mixed['f'],mixed['g'])
        arrays[f'{law}_position_coherent']=position
        arrays[f'{law}_position_dephased']=dephased_position
        region_mass=float(position[np.ix_(range(4,8),range(16,20))].sum())
        assert abs(region_mass-1)<=5e-10
        assert max(abs(array.sum()-1) for key,array in arrays.items() if key.startswith(law))<=5e-10
        summaries[law]={'energy_ray_probability':resolved,'regional_charge_pair_mass':region_mass,
            'complete_outcomes_per_arm':65536,'position_dephasing_max_error':float(np.max(abs(position-dephased_position)))}
    np.savez_compressed(evidence/'raw_source_arrays.npz',
        U_f=unitaries['f'],U_g=unitaries['g'],inverse_f=inverse['f'],inverse_g=inverse['g'],
        ideal_f=ideal['f'],ideal_g=ideal['g'],pointer_matrix=pointer.local_matrix(),
        effects_C=effects['C'],effects_E=effects['E'],
        full_M6_X=compiler.native_pulse('X',math.pi/2),
        full_M6_Y=compiler.native_pulse('Y',math.pi/4),
        full_M6_inverse_X=compiler.native_pulse('X',-math.pi/2),
        full_M6_inverse_Y=compiler.native_pulse('Y',-math.pi/4))
    np.savez_compressed(evidence/'complete_outcome_laws.npz',**arrays)
    (evidence/'source_words.json').write_text(canonical({'preparation':words,
        'inverse':{name:compiler.inverse(word) for name,word in words.items()},
        'pointer_local_word':pointer._local_word(),'blank_reset_conditional_pulse':'X(pi/2) after observed1'}))
    return {'preparation_source_pulses':6,'preparation_full_columns_each':256,
        'preparation_frobenius_errors':errors,'inverse_frobenius_errors':inverse_errors,
        'unitarity_frobenius_defects':defects,'pointer_local_columns':8,
        'pointer_superoperator_entries':certificate['exact_instrument']['logical_superoperator_entries_checked'],
        'pointer_operator_error_upper':str(certificate['operator_error_upper']),
        'pointer_half_diamond_upper':str(certificate['half_diamond_upper']),
        'pointer_factorization_crosscheck_error':float(pointer_error),'laws':summaries,
        'blank_reset_full_superoperator_entry_error':blank_reset_error,
        'raw_arrays_normalized_or_projected':False,'proper_code_gluing_supplied':True}

def reset_service(evidence):
    target=pointer.ready_state()
    probabilities=configuration.success_probabilities(full=True,mode='reset',target=target)
    bound=configuration.budget(full=True,mode='reset',target=target)
    exact_timeout=1-min(probabilities)
    assert exact_timeout<=bound['failure_upper']
    records=[];totals=Counter();checks=0
    for law in ('C','E'):
        for index,initial in enumerate(native.states()):
            seed=int.from_bytes(hashlib.sha256(f'localized-rest-join/{law}/{index}'.encode()).digest(),'big')
            rng=random.Random(seed)
            coins=(rng.randrange(2) for _ in range(bound['native_attempts_max']))
            execution=configuration.execute(initial,coins,full=True,mode='reset',target=target,law=law)
            state=list(initial)
            for seam,coin,label in execution.records:
                a,b=seam;d=state[a]-state[b]
                current=(d-(2*coin-1)*(abs(d)%2))//2
                assert label==(current if law=='C' else d)
                state[a]-=current;state[b]+=current;checks+=1
            assert tuple(state)==execution.output
            outcome=asdict(execution.outcome)
            if execution.outcome.success:
                assert execution.output==target
                assert execution.outcome.initial==initial
            # Controller's own terminal schema is retained without recoding.
            records.append({'law':law,'index':index,'initial':initial,'execution':asdict(execution)})
            totals[execution.outcome.reason]+=1
    with gzip.open(evidence/'auxiliary_reset_tapes.json.gz','wt') as fobj:
        fobj.write(canonical(records))
    (evidence/'auxiliary_reset_probabilities.json').write_text(canonical([str(x) for x in probabilities]))
    delta=F(pointer.certificate()['half_diamond_upper'])
    return {'input_states_per_law':78,'laws':['C','E'],'tape_endpoint_checks':checks,
        'terminal_outcome_counts':dict(totals),'timeout_uniform_bound':str(bound['failure_upper']),
        'exact_timeout_max':str(exact_timeout),'single_reset_attempt_budget':bound['native_attempts_max'],
        'sixteen_bit_complete_error_bound':str(16*(bound['failure_upper']+delta)),
        'blank_preparation_plus_read_complete_error_bound':str(32*(bound['failure_upper']+delta)),
        'timeout_retained':True,'sample_success_used_as_exact_probability':False}

def log_interval(x,n=32):
    partial=sum((x**k/k for k in range(1,n+1)),F())
    return partial,partial+x**(n+1)/((n+1)*(1-x))

def energy():
    intervals={str(lam):[1260*v for v in log_interval(F(lam,1260))] for lam in (6,12)}
    E6,E12=(float(sum(intervals[str(lam)])/2) for lam in (6,12))
    return {'exact_log_intervals':{k:[str(v) for v in vs] for k,vs in intervals.items()},
        'E6':E6,'E12':E12,'E12_over_E6':E12/E6,'finite_clock_nonadditivity':E12-2*E6,
        'optional_anchor_keV':511,'anchored_pair_keV':511*E12/E6,
        'held_source_energy':0,'mode_selected_by_anchor':False,
        'bounded_energy_continuum_limit_claimed':False,
        'physical_mass_clock_identity_claimed':False}

def main():
    started=time.perf_counter();prereg=check_pins();imports()
    evidence=HERE/'evidence'
    evidence.mkdir(exist_ok=False)
    receipt={'schema':'oph.localized_signed_rest.v1','prereg_sha256':prereg,
        'candidate':'A-signed-compact','kernel':kernel(evidence)}
    receipt['reader']=source_reader(evidence)
    receipt['auxiliary_reset']=reset_service(evidence)
    receipt['energy']=energy()
    receipt['native_calls']=CALLS
    receipt['elapsed_seconds']=time.perf_counter()-started
    receipt['finite_component_pass']=True
    receipt['official_row_award']=False
    receipt['full_row_completion_claimed']=False
    receipt['artifact_sha256']={str(p.relative_to(HERE)):digest(p) for p in sorted(evidence.iterdir())}
    (HERE/'RECEIPT.json').write_text(canonical(receipt))
    print(canonical({k:receipt[k] for k in ('candidate','native_calls','elapsed_seconds','finite_component_pass')}))

if __name__=='__main__':main()
