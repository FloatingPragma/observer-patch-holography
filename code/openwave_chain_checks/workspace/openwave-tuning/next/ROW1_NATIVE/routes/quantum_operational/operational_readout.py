"""Coarse C10 regional-charge apparatus; explicit coherent code-transfer scope."""
from pathlib import Path
import sys,math
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[4]
sys.path.insert(0,str(ROOT/'reverse-engineering-reality'))
from model.observer.source_scalar import compiler
from model.reference_recovery import recovery
from model.native_program import interpreter


def phase(theta,bit):
    return [compiler.pulse('I',theta/2,[bit]),compiler.pulse('Z',-theta/2,[bit])]


def toffoli(a,b,t):
    # Standard decomposition; two source iH pulses require one -I correction.
    w=[compiler.pulse('H',math.pi/2,[t])]
    w+=compiler.cnot(b,t)+phase(-math.pi/4,t)
    w+=compiler.cnot(a,t)+phase(math.pi/4,t)
    w+=compiler.cnot(b,t)+phase(-math.pi/4,t)
    w+=compiler.cnot(a,t)+phase(math.pi/4,b)+phase(math.pi/4,t)
    w+=[compiler.pulse('H',math.pi/2,[t])]
    w+=compiler.cnot(a,b)+phase(math.pi/4,a)+phase(-math.pi/4,b)+compiler.cnot(a,b)
    return w+[compiler.pulse('I',math.pi,[t])]


def sum_word():
    # Three source bits; pointer low3/high4. Each controlled increment uses
    # carry from the OLD low bit, then flips low. Count0..3 never overflows.
    w=[]
    for control in range(3):w+=toffoli(control,3,4)+compiler.cnot(control,3)
    return w


def compiled_oracle(source_bits=3):
    if source_bits not in (2,3):raise ValueError("two or three decision bits required")
    local=1<<source_bits;work_bits=source_bits+2
    initial=np.zeros((4*local,local),complex);initial[:local]=np.eye(local)
    expected=np.zeros_like(initial)
    for bits in range(local):expected[bits+(bits.bit_count()<<source_bits),bits]=1
    word=(sum_word() if source_bits==3 else toffoli(0,1,3)+compiler.cnot(0,2)+compiler.cnot(1,2));admitted=interpreter.source_word(word,work_bits)
    run=recovery.execute_complete(word,work_bits,initial)
    error=float(np.linalg.norm(run.success-expected,2))
    if error>1e-10:raise ArithmeticError(('actual sum oracle',error))
    return run,expected,admitted,error


def decoded_occupations(records,boundary,modulus=6):
    records=np.asarray(records);boundary=np.asarray(boundary)
    if (type(modulus) is not int or modulus!=6 or records.ndim!=2 or records.shape[1]!=6
            or records.shape[0]==0 or boundary.shape!=(records.shape[0],)
            or records.dtype.kind not in 'iu' or boundary.dtype.kind not in 'iu'
            or np.any(records<0) or np.any(records>=6) or np.any(boundary<-3) or np.any(boundary>3)):
        raise ValueError('complete canonical integer Z6 records and integer B in[-3,3] required')
    edges=((0,3),(1,3),(1,4),(2,4),(2,5),(0,5));g=np.array([1,1,1,0,0,0])
    occupations=[]
    for f,B in zip(records,boundary):
        div=np.zeros(6,int)
        for value,(a,b) in zip(f,edges):div[a]+=value;div[b]-=value
        n=(g-div-B*np.eye(6,dtype=int)[0])%modulus
        if any(x not in (0,1) for x in n):raise ValueError('Gauss record outside Pauli interval')
        if int(sum(n))!=3-int(B):raise ValueError('integer boundary inconsistent with complete number sector')
        occupations.append(n)
    return np.array(occupations)


def coarse_branches(rho,values):
    # Boolean masks implement whole degenerate projectors, not fine-state readout.
    return {int(q):rho*np.outer(values==q,values==q) for q in np.unique(values)}


def lift_local(matrix,states,region=(0,1,3)):
    """Declared coherent Gauss code reindexing; untouched fiber is spectator.

    Canonical states contain six cycle words for every occupation. The index
    in each mask's listed fiber is an explicit spectator register. This is a
    finite isometry, not a claim that raw native C/E synthesizes the transfer.
    """
    bymask={};lookup={}
    for i,(mask,f) in enumerate(states):
        cycle=len(bymask.setdefault(mask,[]));bymask[mask].append(i);lookup[mask,cycle]=i
    if set(map(len,bymask.values()))!={6}:raise ValueError('six complete fibers required')
    out=np.zeros((len(states),len(states)),complex)
    for mask,ids in bymask.items():
        local=sum(((mask>>p)&1)<<z for z,p in enumerate(region))
        for output in range(1<<len(region)):
            target=mask
            for z,p in enumerate(region):target=(target&~(1<<p))|(((output>>z)&1)<<p)
            for cycle,col in enumerate(ids):out[lookup[target,cycle],col]=matrix[output,local]
    return out


if __name__=='__main__':
 r,e,a,error=compiled_oracle();print({'pulses':len(r.word),'error':error,'failure_effect_norm':float(np.linalg.norm(r.failure_gram,2)),'admitted':a['schema']})


def native_pointer_instrument(oracle):
    """Two actual C/E native-pointer readouts; no ideal terminal projector.

    Auxiliary reset is prepended separately as the pinned finite reset law.
    Export failures stop with their exact source-word prefix family.
    """
    import importlib.util
    path=HERE.parent/'auditor_bridge/native_pointer_bridge.py'
    spec=importlib.util.spec_from_file_location('operational_native_pointer_provider',path)
    provider=importlib.util.module_from_spec(spec);spec.loader.exec_module(provider)
    outcomes,receipt=provider.actual_pointer()
    from model.native_interface import pointer
    initial=np.eye(8,dtype=complex)[:,[pointer.READY_CODE,pointer.READY_CODE+4]]
    pointer_run=recovery.execute_complete(pointer._local_word(),3,initial)
    result={};errors={};work_bits=oracle.bits;source_bits=work_bits-2;local=1<<source_bits
    for law,records in outcomes.items():
        branches={q:[] for q in range(4)};fine=[]
        for meta1,ks1 in records:
            for meta2,ks2 in records:
                q=meta1['decoded_pointer_bit']+2*meta2['decoded_pointer_bit'];indices=[]
                for k1 in ks1:
                    for k2 in ks2:
                        state=compiler._apply(oracle.success,k1,[source_bits],work_bits)
                        state=compiler._apply(state,k2,[source_bits+1],work_bits)
                        # Trace the two output pointer bits, retaining one Kraus
                        # branch per hidden native output and output word.
                        for output in range(4):
                            k=state[local*output:local*(output+1)]
                            indices.append((q,len(branches[q])));branches[q].append(k)
                fine.append({'first':meta1,'second':meta2,'coarse_sum':q,'kraus_indices':indices})
        worst=0.
        for q,ks in branches.items():
            actual=sum(np.kron(k.conj(),k) for k in ks)
            p=np.diag([int(i.bit_count()==q) for i in range(local)])
            worst=max(worst,float(np.max(abs(actual-np.kron(p,p)))))
        success_effect=sum(k.conj().T@k for ks in branches.values() for k in ks)
        # Sequential complete effects: compute oracle abort, first read abort,
        # then second read abort after all first-read continuing outcomes.
        lift=lambda a,bit:interpreter._lift(a,[bit],work_bits)
        first=oracle.success.conj().T@lift(pointer_run.failure_gram,source_bits)@oracle.success
        second=np.zeros((local,local),complex)
        for _,ks in records:
            for k in ks:
                state=compiler._apply(oracle.success,k,[3],5)
                second+=state.conj().T@lift(pointer_run.failure_gram,source_bits+1)@state
        total=success_effect+oracle.failure_gram+first+second
        errors[law]={'coarse_instrument_matrix_unit_error':worst,'complete_effect_error':float(np.max(abs(total-np.eye(local)))),'source_word_abort_effect_norm':float(np.linalg.norm(oracle.failure_gram,2)),'first_native_pointer_abort_effect_norm':float(np.linalg.norm(first,2)),'second_native_pointer_abort_effect_norm':float(np.linalg.norm(second,2))}
        result[law]={'branches':branches,'fine_native_records':fine,'failure_effects':{'oracle':oracle.failure_gram,'pointer_first':first,'pointer_second':second}}
    return result,errors,pointer_run,receipt


def encoded_indices(states,region=(0,1,3)):
    others=tuple(x for x in range(6) if x not in region);counts={};ids=np.zeros((8,48),int)
    for i,(mask,f) in enumerate(states):
        cycle=counts.get(mask,0);counts[mask]=cycle+1
        local=sum(((mask>>p)&1)<<z for z,p in enumerate(region))
        rest=sum(((mask>>p)&1)<<z for z,p in enumerate(others));ids[local,6*rest+cycle]=i
    if set(ids.flat)!=set(range(len(states))):raise ValueError('coherent code isometry not bijective')
    return ids


def apply_local_channel(rho,kraus,indices):
    flat=indices.ravel();blocks=rho[np.ix_(flat,flat)].reshape(8,48,8,48)
    superoperator=sum(np.einsum('ab,cd->acbd',k,k.conj()) for k in kraus)
    out=np.einsum('acbd,brds->arcs',superoperator,blocks,optimize=True).reshape(384,384)
    result=np.zeros_like(rho,dtype=complex);result[np.ix_(flat,flat)]=out
    return result
