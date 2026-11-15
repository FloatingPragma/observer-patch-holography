"""Controlled native source words for a probe of ONE supplied walk W.

This module supplies an apparatus, never a Hamiltonian or mass parameter.
Every scalar phase is controlled, including source I and J phases.
"""
from __future__ import annotations
from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import expm

HERE=Path(__file__).resolve().parent
WORKSPACE=HERE.parents[2]
RER=WORKSPACE/'reverse-engineering-reality'
sys.path.insert(0,str(RER))
from model.observer import native
from model.observer.source_scalar import compiler
from model.observer.source_discrimination import circuit
from model.reference_recovery import recovery
from model.native_interface import pointer


def controlled_local_word(kind,angle):
    """Target bits first, high bit is the retained probe control.

    Exact target diag(I, exp(i theta P)) on the proper code. The literal
    full-six-state responses and all export complements are executed.
    """
    if kind=='J':
        return compiler.compile_diagonal([0.,0.,0.,0.,angle,angle,angle,0.])
    if kind=='H':
        return ([compiler.pulse('Y',math.pi/8,[0])]
                +controlled_local_word('Z',angle)
                +[compiler.pulse('Y',-math.pi/8,[0])])
    if kind not in ('I','X','Y','Z'):
        raise ValueError('a source compiler primitive is required')
    return (compiler.pauli_rotation(kind+'I',angle/2)
            +compiler.pauli_rotation(kind+'Z',-angle/2))


def control_word(word,control):
    if type(control)is not int or control<0:
        raise ValueError('one existing nonnegative probe address required')
    output=[]
    for event in compiler._validate(word):
        if control in event['targets']:
            raise ValueError('probe control aliases a source target')
        output+=compiler.remap(controlled_local_word(event['kind'],event['angle']),
                              [*event['targets'],control])
    return output


def target_generator(kind):
    X=np.array([[0,1],[1,0]],complex)
    Y=np.array([[0,-1j],[1j,0]],complex)
    Z=np.diag([1.,-1.])
    return {'I':np.eye(2),'X':X,'Y':Y,'Z':Z,'H':(X+Z)/np.sqrt(2),
            'J':np.diag([1.,1.,1.,0.])}[kind]


def actual_native_pointer():
    initial=np.eye(8,dtype=complex)[:,[pointer.READY_CODE,pointer.READY_CODE+4]]
    run=recovery.execute_complete(pointer._local_word(),3,initial)
    outcomes={}
    for law in ('C','E'):
        groups=defaultdict(lambda:defaultdict(lambda:np.zeros((2,2),complex)))
        metadata={}
        for code in range(4):
            source=circuit.decode(code);K=run.success[[code,code+4]]
            for coin in (0,1):
                output,(_,current)=native.native(source,pointer.SEAM,coin)
                label=current if law=='C'else source[0]-source[1]
                key=(coin,label);bit=int(2*current+2*coin-1<0)
                groups[key][output]+=K/np.sqrt(2)
                metadata[key]={'coin':coin,'label':label,'current':current,'bit':bit}
        outcomes[law]=[(metadata[k],list(v.values()))for k,v in sorted(groups.items())]
        for bit in (0,1):
            actual=sum(np.kron(k.conj(),k)for meta,ks in outcomes[law]if meta['bit']==bit for k in ks)
            P=np.diag([int(bit==0),int(bit==1)])
            assert np.max(abs(actual-np.kron(P,P)))<2e-12
    return outcomes,run


def certificate():
    rows=[]
    for kind in ('I','X','Y','Z','H','J'):
        generator=target_generator(kind);d=len(generator);bits=2 if d==2 else 3
        for angle in (0.,.37,-.23,math.pi/2):
            word=controlled_local_word(kind,angle)
            run=recovery.execute_complete(word,bits,np.eye(2*d))
            target=np.zeros((2*d,2*d),complex)
            target[:d,:d]=np.eye(d);target[d:,d:]=expm(1j*angle*generator)
            error=float(np.linalg.norm(run.success-target,2))
            gram=float(np.linalg.norm(run.success.conj().T@run.success+run.failure_gram-np.eye(2*d),2))
            assert max(error,gram)<3e-11
            rows.append({'kind':kind,'angle':angle,'pulses':len(word),'all_input_columns':2*d,
                         'target_operator_error':error,'complete_gram_error':gram,
                         'export_failure_trace':float(np.trace(run.failure_gram).real),'word':word})
    outcomes,pointer_run=actual_native_pointer()
    # A source scalar phase must remain a physical relative probe phase.
    angle=.37
    controlled=recovery.execute_complete(controlled_local_word('I',angle),2,np.eye(4)).success
    plus=np.array([1,0,1,0],complex)/np.sqrt(2)
    state=controlled@plus
    coherence=2*state[2]*state[0].conjugate()
    assert abs(coherence-np.exp(1j*angle))<2e-12
    assert abs(coherence-1)>.3
    paths=[Path(__file__),RER/'model/observer/source_scalar/compiler.py',
        RER/'model/reference_recovery/recovery.py',RER/'model/native_interface/pointer.py',
        RER/'model/observer/native.py']
    return {'schema':'oph.common_walk.native_phase_probe.v1',
        'checks':{'all_six_source_primitives_controlled':True,'all_local_input_columns_and_export_failures_retained':True,
                  'native_pointer_full_matrix_units':True,'scalar_phase_is_not_dropped':True},
        'controlled_gate_cases':rows,
        'native_pointer_records':{law:[meta for meta,ks in out]for law,out in outcomes.items()},
        'native_pointer_export_failure_trace':float(np.trace(pointer_run.failure_gram).real),
        'common_source_contract':'controlled words must come from the one frozen W provider; no H, mass, eigenbasis or empirical calibration is accepted here',
        'source_word_control_identity':'C(U_N ... U_1)=C(U_N)...C(U_1), including every scalar phase',
        'full_phase_probe_readout':'prepare probe+, controlled W^n, X/Y analysis, actual native pointer; every source and apparatus complement retained',
        'spectral_interpretation':'unit-circle phase modulo2pi per W step; no absolute rest mass or physical seconds inferred',
        'missing_before_full_walk_claim':['frozen common W word/provider','same common preparation/reference','actual apparatus reset and complete common W implementation bounds'],
        'source_sha256':{str(p.relative_to(WORKSPACE)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},
        'official_rows_changed':[]}


if __name__=='__main__':
    result=certificate()
    (HERE/'NATIVE_PHASE_PROBE.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'checks':result['checks'],'maximum_operator_error':max(r['target_operator_error']for r in result['controlled_gate_cases'])},indent=2))
