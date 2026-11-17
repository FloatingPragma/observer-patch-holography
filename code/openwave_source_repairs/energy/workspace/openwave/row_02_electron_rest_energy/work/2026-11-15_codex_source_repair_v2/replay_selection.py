#!/usr/bin/env python3
"""Raw source-word propagator for a complete escape/reinjection channel.

Native primitive products are executed once; channel composition then reuses
those stored maps. Represented pump histories are never counted as native calls.
The region and reset ideal comparison instruments have separately listed
implementation obligations. This script does not claim a full-row pass.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,math,sys,time,argparse
import numpy as np
HERE=Path(__file__).resolve().parent
OW=HERE.parents[2]
FLOW=OW/'shared/compact_logarithmic_flow'
sys.path.insert(0,str(FLOW))
import native_words as nw
import certify
import run_flow

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,x):p.write_text(json.dumps(x,sort_keys=True,indent=2,allow_nan=False)+'\n')
def check_freeze():
    frozen=json.loads((HERE/'PREREG_REPLAY.json').read_text())
    for name,digest in frozen['files'].items():assert sha(HERE/name)==digest,name
    run_flow.check_pins()

def gram_sqrt_upper(squared,bits=100):
    q=1<<bits;squared=F(squared)
    n=(squared.numerator*q*q+squared.denominator-1)//squared.denominator
    return F(math.isqrt(n)+1,q)

def main():
    started=time.monotonic();check_freeze()
    archive=FLOW/'attempts/stage2_002/evidence'
    origin=json.loads((HERE/'PREREG_REPLAY.json').read_text())
    for name,digest in origin['archive_files'].items():assert sha(archive/name)==digest
    index=json.loads((archive/'NATIVE_PRIMITIVES.json').read_text())
    matrices=np.load(archive/'native_M6_primitives.npz')['matrices']
    for key,matrix,bound in zip(index['keys'],matrices,index['exact_squared_error_bounds']):
        nw.GATES[tuple(key)]=matrix;nw.RAW_BOUNDS[tuple(key)]=F(bound)
    archive_count=len(nw.GATES)
    print(json.dumps({'retained_native_primitive_matrices_reused':archive_count}),flush=True)
    out=HERE/'execution';out.mkdir(exist_ok=False)
    proposal,source=certify.load();calls=run_flow.recheck_source(source,proposal)
    cert=certify.certify('native');rows=proposal['laws']['native']
    state=np.eye(256,dtype=complex)[:,:12]
    processor=np.zeros((12,12),complex);padding=np.zeros_like(processor)
    pulse_count=0;tree=F();phase=F();word_records=[]
    for j,vector in enumerate(rows['columns']):
        prep,terr,worst=nw.preparation(list(vector)+[0]*16)
        inv=nw.c.inverse(prep);angle,perr=certify.phase(rows['eigenvalues'][j],4)
        diag=nw.address_phase(angle,8)
        for word in (inv,diag,prep):
            state,bad=nw.apply_word(word,state,8);processor+=bad
        padding+=state[240:].conj().T@state[240:];state[240:]=0
        pulses=2*len(prep)+len(diag);pulse_count+=pulses;tree+=4*terr;phase+=perr
        word_records.append({'column':j,'pulses':pulses,'preparer_sha256':nw.digest_word(prep),
                             'phase_sha256':nw.digest_word(diag),'phase_error':str(perr)})
        if j%12==0:print(json.dumps({'native_selection_words':j+1,'of':240}),flush=True)
    # Exact integer dark basis, supplied only to the terminal diagnostic.
    D=np.zeros((240,2),dtype=np.int64)
    D[[4,5,6,7],0]=[1,-1,-1,1];D[[5,6,8,11],1]=[-1,-1,1,1]
    L=1260*np.eye(240,dtype=np.int64)-source['counts']
    assert np.array_equal(L@D,6*D)
    gram=D.T@D;assert np.array_equal(gram,[[4,2],[2,4]])
    Pd=D[:12]@np.array([[1/3,-1/6],[-1/6,1/3]])@D[:12].T
    np.savez_compressed(out/'RAW_FLOW_CHECKPOINT.npz',raw_flow_columns=state,
                        processor_failure_gram=processor,padding_failure_gram=padding)
    Q=state[:12];outside=state[12:240]
    outside_effect=outside.conj().T@outside
    sigma=np.eye(12)/12;rho=sigma.copy();history=[]
    cumulative_escape=0.;cumulative_processor=0.;cumulative_padding=0.
    for n in range(513):
        history.append([n,float(np.trace(rho).real),float(np.trace(Pd@rho).real),
                        cumulative_escape,cumulative_processor,cumulative_padding])
        if n==512:break
        escape=float(np.trace(outside_effect@rho).real)
        cumulative_escape+=escape
        cumulative_processor+=float(np.trace(processor@rho).real)
        cumulative_padding+=float(np.trace(padding@rho).real)
        rho=Q@rho@Q.conj().T+escape*sigma
    weight=float(np.trace(Pd@rho).real);mass=float(np.trace(rho).real)
    original_raw_gate=weight>1-1e-7 and abs(mass-1)<1e-7
    print(json.dumps({'raw_weight':weight,'raw_mass':mass,'raw_noncompact_mass':mass-weight,
                      'original_unconditional_raw_gate':original_raw_gate}),flush=True)
    assert mass-weight<1e-7 and abs(mass-1)<1e-6
    orth=D@np.linalg.inv(np.linalg.cholesky(gram)).T
    selected=state@orth[:12]
    interval=certify.energy_interval(209,210);energy=float(sum(interval)/2)
    expected=np.exp(-4j*energy)*orth
    residual=float(np.linalg.norm(selected[:240]-expected,2))
    escape=float(np.linalg.norm(selected[12:240],2)**2)
    assert residual<1e-7 and escape<1e-12
    single=orth.conj().T@selected[:240]
    # Two independent spatial processors share one local epoch. This is a
    # represented tensor product of the actually executed raw source map.
    pair=np.kron(single,single)
    pair_residual=float(np.linalg.norm(pair-np.exp(-8j*energy)*np.eye(4),2))
    assert pair_residual<1e-7
    # Recompute a sharper all-pulse error from every measured full6x6 bound.
    b=max(gram_sqrt_upper(s) for s in nw.RAW_BOUNDS.values())
    x=pulse_count*b;assert x<1
    e=(x+tree+phase)/(1-x)
    failure=pulse_count*b*b/(1-x)**2+e*e
    word_half=(e*(2+e)+failure)/2
    spectral=4*F(cert['energy_operator_upper'])+F(cert['phase_projector_product_upper'])
    per_epoch=word_half+(1/(1-x)**2+failure)*spectral
    assert 512*per_epoch<F(1,1000)
    # Padded source matrices and every processor complement retain their mass.
    complete_effect=state.conj().T@state+processor+padding
    effect_residual=float(np.linalg.norm(complete_effect-np.eye(12),2))
    assert effect_residual<1e-7
    hold_weight=float(np.trace(Pd@sigma).real);assert abs(hold_weight-1/6)<1e-14
    np.savez_compressed(out/'maps.npz',raw_flow_columns=state,processor_failure_gram=processor,
                        padding_failure_gram=padding,dark_integer_basis=D,terminal_density=rho,
                        epoch_records=np.asarray(history),selected_phase=single,pair_phase=pair)
    keys,matrices,bounds=nw.gate_archive()
    np.savez_compressed(out/'native_primitives.npz',matrices=matrices)
    write(out/'native_primitives.json',{'keys':keys,'squared_error_upper':bounds,'counts':nw.COUNTS})
    write(out/'word_records.json',word_records)
    scientific={'candidate':'A-native-escape-reinjection-v1','source_columns':240,'local_seed_rank':12,
        'selected_rank':2,'selected_laplacian_eigenvalue':6,'epochs':512,'source_flow_time':4,
        'selected_weight':weight,'retained_source_mass':mass,
        'original_unconditional_raw_gate_passed':original_raw_gate,
        'numerical_noncompact_mass':mass-weight,'retained_primitive_matrices_reused':archive_count,'cumulative_escape_records':cumulative_escape,
        'cumulative_processor_failure':cumulative_processor,'cumulative_padding_failure':cumulative_padding,
        'selected_hold_residual':residual,'selected_escape_probability_upper_numeric':escape,
        'same_local_epoch_pair_phase_residual':pair_residual,'source_energy_interval':list(map(str,interval)),
        'hold_control_selected_weight':hold_weight,'native_pulse_applications':nw.COUNTS['pulse_applications'],
        'native_matrix_executions':nw.COUNTS['native_matrix_executions'],'fresh_native_endpoint_calls':calls,
        'full_primitive_error_upper':str(b),'flow_complete_halfdiamond_upper':str(per_epoch),
        '512_flow_invocation_error_upper':str(512*per_epoch),'complete_effect_residual':effect_residual,
        'region_monitor_compilation_executed':False,'reinjection_compilation_executed':False,
        'physical_energy_identified':False,'full_row_solved':False,'external_G6_confirmed':False}
    artifacts={str(p.name):sha(p) for p in out.iterdir() if p.is_file()}
    receipt={'scientific':scientific,'artifacts':artifacts,'prereg_sha256':sha(HERE/'PREREG_REPLAY.json'),
             'elapsed_seconds':time.monotonic()-started}
    write(out/'RECEIPT.json',receipt)
    print(json.dumps({'done':True,'selected_weight':weight,'seconds':receipt['elapsed_seconds']}),flush=True)

if __name__=='__main__':main()
