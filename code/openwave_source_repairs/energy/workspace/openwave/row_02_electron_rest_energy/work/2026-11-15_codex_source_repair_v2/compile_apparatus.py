#!/usr/bin/env python3
"""Native coarse region monitor and complete uniform reinjection recipe."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,math,sys,time
import numpy as np
HERE=Path(__file__).resolve().parent;OW=HERE.parents[2]
sys.path.insert(0,str(OW/'shared/compact_logarithmic_flow'))
import native_words as nw
import reference
from model.observer.native import native
from model.native_interface import pointer
BITS=13

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,x):p.write_text(json.dumps(x,sort_keys=True,indent=2,allow_nan=False)+'\n')
def comparator():
    ops=[]
    def gate(name,targets):ops.append((name,tuple(targets)))
    def controlled(controls,zeros):
        for z in zeros:gate('flip',(z,))
        scratch=list(range(9,13));chain=[(controls[0],controls[1],scratch[0])]
        for j in range(2,len(controls)-1):chain.append((scratch[j-2],controls[j],scratch[j-1]))
        for c in chain:gate('ccx',c)
        gate('ccx',(scratch[len(controls)-3],controls[-1],8))
        for c in reversed(chain):gate('ccx',c)
        for z in reversed(zeros):gate('flip',(z,))
    controlled([3,4,5,6,7],[3,4,5,6,7])
    controlled([2,3,4,5,6,7],[2,4,5,6,7])
    return ops

def main():
    started=time.monotonic();freeze=json.loads((HERE/'PREREG_APPARATUS.json').read_text())
    for name,digest in freeze['files'].items():assert sha(HERE/name)==digest
    out=HERE/'apparatus';out.mkdir(exist_ok=False)
    words,raw,failed,metadata=reference.routines();ops=comparator()
    labels=np.arange(512);truth=labels.copy()
    for name,t in ops:
        if name=='flip':truth^=1<<t[0]
        else:truth^=(((truth>>t[0])&1)&((truth>>t[1])&1))<<t[2]
    wanted=labels^(((labels&255)<12).astype(int)<<8)
    assert np.array_equal(truth,wanted)
    # Full scratch evolution of three source/reference witnesses, with no
    # normalization after raw primitive products.
    state=np.zeros((1<<BITS,3),complex)
    state[:256,0]=1/16
    state[255,1]=1
    state[[4,5,6,7],2]=[.5,-.5,-.5,.5]
    original=state.copy();fail=np.zeros((3,3),complex)
    for name,t in ops:
        ix=nw.indices(tuple(reversed(t)),BITS);block=state[ix]
        fail+=np.einsum('iac,ij,jad->cd',block.conj(),failed[name],block,optimize=True)
        state[ix]=(raw[name]@block.reshape(len(raw[name]),-1)).reshape(block.shape)
    ideal=np.zeros_like(state);ideal[wanted]=original[:512]
    err=float(np.linalg.norm(state-ideal,2));assert err<1e-8
    pointer_effects=reference.pointer_effects('C')
    flag_indices=nw.indices((8,),BITS);blocks=state[flag_indices]
    record_grams=np.asarray([np.einsum('iac,ij,jad->cd',blocks.conj(),E,blocks,optimize=True)
                             for E in pointer_effects])
    expected=np.asarray([[61/64,1,0],[3/64,0,1]])
    probabilities=np.real(np.diagonal(record_grams,axis1=1,axis2=2))
    assert np.max(abs(probabilities-expected))<1e-8
    dirty=state[512:].conj().T@state[512:]
    # Exact uniform12 selector, every four-bit native tape retained.
    tapes=[];calls=0
    for value in range(16):
        q=(1,0);records=[];bits=[]
        for j in range(4):
            coin=(value>>j)&1;outq,record=native(q,(0,1),coin);calls+=1
            records.append([list(q),coin,list(outq),list(record)]);bits.append(record[0]);q=outq
        read=sum(b<<j for j,b in enumerate(bits));assert read==value
        tapes.append({'word':value,'accept':value<12,'records':records})
    timeout=F(1,4)**16
    routes=[];route_fail=[];route_meta=[]
    for p in range(12):
        word=[gate for bit in range(8) if p>>bit&1 for gate in nw.c.remap(words['flip'],[bit])]
        blank=np.zeros(256);blank[0]=1
        routed,bad=nw.apply_word(word,blank,8)
        routes.append(routed[:,0]);route_fail.append(float(bad[0,0].real))
        route_meta.append({'p':p,'pulses':len(word),'sha256':nw.digest_word(word)})
    routes=np.asarray(routes).T;seed=routes@routes.conj().T/12
    target=np.diag([1/12]*12+[0.]*244)
    seed_error=float(np.linalg.norm(seed-target,2));assert seed_error<1e-9
    pulse_count=sum(len(words[name]) for name,t in ops)
    comp=nw.complete_bound(pulse_count)
    reset_timeout=F(35,2147483648)
    pointer_error=F(200000000001,20000000000000000000000)
    # Thirteen blank/comparator bits and eight work reset bits are charged
    # every epoch, even when an escape does not require reinjection.
    blank_cost=21*(reset_timeout+pointer_error)
    route_cost=max(r['pulses'] for r in route_meta)*nw.PULSE_GATE
    complete=blank_cost+comp['half_diamond']+pointer_error+route_cost+timeout
    assert 512*complete<F(1,2000)
    np.savez_compressed(out/'maps.npz',native_routine_ccx=raw['ccx'],native_routine_flip=raw['flip'],
                        monitor_states=state,monitor_processor_failure=fail,dirty_scratch=dirty,
                        record_grams=record_grams,uniform_routes=routes,route_failure=np.asarray(route_fail))
    report={'candidate':'A-native-escape-reinjection-v1','all_clean_truth_columns':512,
        'dirty_scratch_dimension':1<<BITS,'raw_quantum_witnesses':3,'raw_monitor_operator_witness_error':err,
        'complete_native_pointer_probabilities':probabilities.tolist(),'native_selector_calls':calls,
        'selector_draws_max':16,'selector_coin_calls_max':64,'selector_timeout':str(timeout),
        'selector_conditional_address_probability':'1/12','selector_tapes':tapes,
        'source_reinjection_operator_error':seed_error,'source_route_words':route_meta,
        'monitor_source_pulses':pulse_count,'routine_metadata':metadata,'reset_bits_charged_per_epoch':21,
        'reset_attempts_per_bit':2678,'reset_timeout_per_bit':str(reset_timeout),
        'pointer_complete_error':str(pointer_error),'complete_epoch_apparatus_halfdiamond_upper':str(complete),
        '512_epoch_apparatus_halfdiamond_upper':str(512*complete),
        'native_reset_bound_inherited_not_reexecuted':True,'full_record_histories_factorized':True,
        'apparatus_service_assumptions':['proper-code load/export','blank reset routine','coherent address routing',
                                       'analogue M6 control precision','retained controller memory'],
        'full_row_solved':False}
    write(out/'RECEIPT.json',{'scientific':report,'artifacts':{'maps.npz':sha(out/'maps.npz')},
          'prereg_sha256':sha(HERE/'PREREG_APPARATUS.json'),'elapsed_seconds':time.monotonic()-started})
    print(json.dumps({'apparatus_complete':True,'512_epoch_error':float(512*complete)}))

if __name__=='__main__':main()
