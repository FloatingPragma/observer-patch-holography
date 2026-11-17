#!/usr/bin/env python3
"""Native controlled-flow Ramsey records for one and two local processors.

Two local source inputs execute the same complete controlled word. The joint
two-source instrument is their factorized composition with one shared control;
its65536-dimensional retained source branches and native control-pointer
records are materialized. No source-energy value is substituted for a record.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,math,sys,time
import numpy as np
HERE=Path(__file__).resolve().parent;OW=HERE.parents[2]
sys.path.insert(0,str(OW/'shared/compact_logarithmic_flow'))
import native_words as nw
import reference,certify,run_flow

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,x):p.write_text(json.dumps(x,sort_keys=True,indent=2,allow_nan=False)+'\n')
def controlled_phase(angle):
    word=[];coefficient=angle/512;assert F(coefficient)*512==F(angle)
    for mask in range(512):
        label=['Z' if mask>>j&1 else 'I' for j in range(9)]
        word.extend(nw.c.pauli_rotation(label,coefficient*(-1 if mask&256 else 1)))
    return word

def read(branches):
    # branches is a source-by-control array, retaining its absolute norm.
    rho=branches.conj().T@branches
    effects=reference.pointer_effects('C');records={};failures={}
    for name,word in [('X',[nw.c.pulse('Y',math.pi/4,[0])]),
        ('Y',[nw.c.pulse('Z',math.pi/4,[0]),nw.c.pulse('Y',math.pi/4,[0])])]:
        beam,bad=nw.apply_word(word,np.eye(2),1)
        effects_read=[beam.conj().T@E@beam for E in effects]
        records[name]=[float(np.einsum('ab,ba->',rho,E).real) for E in effects_read]
        failures[name]=float(np.einsum('ab,ba->',rho,bad).real)
    return records,failures

def phase(records):return complex(records['X'][0]-records['X'][1],records['Y'][0]-records['Y'][1])

def main():
    started=time.monotonic();freeze=json.loads((HERE/'PREREG_PHASE.json').read_text())
    for name,digest in freeze['files'].items():assert sha(HERE/name)==digest
    run_flow.check_pins();out=HERE/'phase';out.mkdir(exist_ok=False)
    proposal,source=certify.load();cert=certify.certify('native');rows=proposal['laws']['native']
    # Two orthonormal modes of the process-selected two-dimensional sector.
    # Their preparation here is a calibration of the retained reader, not a
    # claim that these vectors were the initial state of the selection pump.
    D=np.zeros((256,2));D[[4,5,6,7],0]=[.5,-.5,-.5,.5]
    D[[4,5,6,7,8,11],1]=np.array([-1,-1,-1,-1,2,2])/math.sqrt(12)
    state=np.zeros((512,2),complex);state[:256]=D
    preparation_beam=[nw.c.pulse('Y',-math.pi/4,[8])]
    state,beamfail=nw.apply_word(preparation_beam,state,9)
    actual_beam=nw.native_gate('Y',-math.pi/4)[:2,:2][:,0]
    failures=beamfail.copy();padding=np.zeros((2,2),complex);meta=[];tree=F();phase_error=F();N=0
    for j,vector in enumerate(rows['columns']):
        prep,terr,worst=nw.preparation(list(vector)+[0]*16);inv=nw.c.inverse(prep)
        angle,perr=certify.phase(rows['eigenvalues'][j],4);diag=controlled_phase(angle)
        for word in (inv,diag,prep):
            state,bad=nw.apply_word(word,state,9);failures+=bad
        bad=np.concatenate((state[240:256],state[496:512]));padding+=bad.conj().T@bad
        state[240:256]=0;state[496:512]=0
        n=2*len(prep)+len(diag);N+=n;tree+=4*terr;phase_error+=perr
        meta.append({'ray':j,'pulses':n,'preparation_sha256':nw.digest_word(prep),
                     'controlled_phase_sha256':nw.digest_word(diag),'phase_error':str(perr)})
        if j%12==0:print(json.dumps({'native_controlled_words':j+1,'of':240}),flush=True)
    singles=[];single_fail=[]
    for j in range(2):
        rec,bad=read(np.column_stack((state[:256,j],state[256:,j])));singles.append(rec);single_fail.append(bad)
    # Factorize the two executed conditional local source maps. The shared
    # beam is applied once; division below removes only that known linear
    #input coefficient, never a success probability or an output norm.
    branches=np.column_stack([np.kron(state[256*k:256*(k+1),0]/actual_beam[k],
                                     state[256*k:256*(k+1),1]/actual_beam[k])*actual_beam[k]
                              for k in (0,1)])
    pair,pair_fail=read(branches)
    z=[phase(r) for r in singles];zpair=phase(pair)
    additivity=abs(zpair-z[0]*z[1]);assert additivity<1e-7
    lo,hi=certify.energy_interval(209,210);energy=float((lo+hi)/2)
    expected=np.exp(-4j*energy);single_residual=max(abs(v-expected) for v in z)
    assert single_residual<1e-7
    # Same physical detector on a complete HOLD source continuation, keeping
    # the control/source coherence and the source input fixed.
    hold=np.column_stack((actual_beam[0]*D[:,0],actual_beam[1]*D[:,0]))
    hold_records,hold_fail=read(hold);hold_phase=phase(hold_records)
    assert abs(hold_phase-1)<1e-7 and abs(hold_phase-z[0])>.1
    # A deliberate relative phase is only a reader-sensitivity plant.
    plant=hold.copy();plant[:,1]*=1j;plant_records,plant_fail=read(plant)
    assert abs(phase(plant_records)-1j)<1e-7
    raw=max(F(s) for s in nw.RAW_BOUNDS.values());q=1<<100
    b=F(math.isqrt((raw.numerator*q*q+raw.denominator-1)//raw.denominator)+1,q)
    x=N*b;e=(x+tree+phase_error)/(1-x);f=N*b*b/(1-x)**2+e*e
    flow=(e*(2+e)+f)/2+(1/(1-x)**2+f)*(4*F(cert['energy_operator_upper'])+F(cert['phase_projector_product_upper']))
    pointer=F(200000000001,20000000000000000000000)
    joint_bound=2*flow+4*pointer+12*nw.PULSE_GATE
    assert joint_bound<F(1,10000)
    np.savez_compressed(out/'records.npz',controlled_source_columns=state,
        processor_failure_gram=failures,padding_failure_gram=padding,
        joint_source_control_branches=branches,actual_control_beam=actual_beam)
    write(out/'words.json',meta)
    report={'candidate':'A-native-escape-reinjection-v1','native_controlled_flow_pulses_per_local_word':N,
        'native_pulse_applications':nw.COUNTS['pulse_applications'],'native_matrix_executions':nw.COUNTS['native_matrix_executions'],
        'source_inputs_executed':2,'controlled_source_dimension':512,'joint_source_dimension':65536,
        'same_local_epoch':4,'single_native_pointer_records':singles,'pair_native_pointer_records':pair,
        'single_reader_failure':single_fail,'pair_reader_failure':pair_fail,
        'single_complex_phase':[[w.real,w.imag] for w in z],'pair_complex_phase':[zpair.real,zpair.imag],
        'measured_additivity_complex_residual':additivity,'source_phase_prediction_residual':single_residual,
        'hold_control_records':hold_records,'hold_reader_failure':hold_fail,
        'plant_records':plant_records,'plant_reader_failure':plant_fail,
        'complete_single_controlled_flow_error_upper':str(flow),'complete_joint_phase_error_upper':str(joint_bound),
        'joint_source_instrument_execution':'factorized composition of two executedlocal controlled maps; joint output branches materialized',
        'control_input_domain':'two selected energy-sector basis vectors; selection-to-sector error supplied separately',
        'event_clock':'one declared local controlled U4 macrostep per source; no shared exclusive event',
        'physical_proper_time_claimed':False,'full_row_solved':False}
    write(out/'RECEIPT.json',{'scientific':report,'artifacts':{p.name:sha(p) for p in out.iterdir() if p.is_file()},
        'prereg_sha256':sha(HERE/'PREREG_PHASE.json'),'elapsed_seconds':time.monotonic()-started})
    print(json.dumps({'phase_done':True,'additivity_residual':additivity,'seconds':time.monotonic()-started}))

if __name__=='__main__':main()
