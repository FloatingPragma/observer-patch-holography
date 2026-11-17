#!/usr/bin/env python3
"""Source-count-derived invariant-sector phase, with complete native reader."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,math,sys,time
import numpy as np
HERE=Path(__file__).resolve().parent;OW=HERE.parents[2]
sys.path.insert(0,str(OW/'shared/compact_logarithmic_flow'))
import native_words as nw
import reference,certify,run_flow
from run_phase import controlled_phase,read,phase

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,x):p.write_text(json.dumps(x,sort_keys=True,indent=2,allow_nan=False)+'\n')
def sector(cell,first,second):
    integer_vectors=[first,[2*b-a for a,b in zip(first,second)]]
    source=np.zeros((512,12),complex);source[cell:cell+12]=np.eye(12)
    source,failed=nw.apply_word([nw.c.pulse('Y',-math.pi/4,[8])],source,9)
    pulses=1;meta=[]
    for integer_vector in integer_vectors:
        full=[0]*256;full[cell:cell+12]=integer_vector
        prep,error,worst=nw.preparation(full);reflection=controlled_phase(math.pi)
        for word in [nw.c.inverse(prep),reflection,prep]:
            source,bad=nw.apply_word(word,source,9);failed+=bad;pulses+=len(word)
        meta.append({'integer_vector':integer_vector,'preparation_pulses':len(prep),
                     'preparation_sha256':nw.digest_word(prep),'controlled_reflection_sha256':nw.digest_word(reflection),
                     'preparation_error_upper':str(error)})
    source,bad=nw.apply_word([nw.c.pulse('Y',math.pi/4,[8])],source,9);failed+=bad;pulses+=1
    raw_first=np.asarray(first);raw_second=np.asarray(second)
    D=np.column_stack((raw_first,raw_second));P=D@np.array([[2,-1],[-1,2]])@D.T/6
    target=np.zeros_like(source);target[cell:cell+12]=np.eye(12)-P;target[cell+256:cell+268]=-P
    residual=float(np.linalg.norm(source-target,2));assert residual<1e-7
    pointer=reference.pointer_effects('C');blocks=np.stack((source[:256],source[256:]))
    effects=np.asarray([np.einsum('iac,ij,jad->cd',blocks.conj(),E,blocks,optimize=True) for E in pointer])
    assert np.linalg.norm(effects[1]-P,2)<1e-7
    assert np.linalg.norm(effects.sum(axis=0)+failed-np.eye(12),2)<1e-7
    return source,failed,effects,{'cell_offset':cell,'word_pulses':pulses,'source_word_metadata':meta,
        'full12column_operator_residual':residual,'selected_projector':P.tolist(),
        'complete_effect_error':float(np.linalg.norm(effects.sum(axis=0)+failed-np.eye(12),2))}

def clock_phase(energy):
    # On the source-certified invariant sector, controlledH is diagonal on
    # the retained control. Every one of these is an actual fullM6 pulse.
    theta=-4*float(energy);word=[]
    for k in range(32):word.extend([nw.c.pulse('I',theta/64,[0]),nw.c.pulse('Z',-theta/64,[0])])
    word=[nw.c.pulse('Y',-math.pi/4,[0])]+word
    out,failed=nw.apply_word(word,np.array([1.,0.]),1)
    records,read_failed=read(out.T)
    return records,read_failed,float(failed[0,0].real),nw.digest_word(word)

def main():
    started=time.monotonic();reg=json.loads((HERE/'PREREG_POISSON_ENERGY.json').read_text())
    for name,digest in reg['files'].items():assert sha(HERE/name)==digest
    clock=HERE.parent/'2026-11-15_codex_poisson_source_clock'
    for name,digest in reg['clock_files'].items():assert sha(clock/name)==digest
    counted=json.loads((clock/'RECEIPT.json').read_text());assert counted['all_gates_pass']
    singles=list(counted['native_single_generator_values']);pairs=list(counted['native_pair_generator_values'])
    assert len(set(singles))==1 and len(set(pairs))==1
    one=F(singles[0]);two=F(pairs[0]);assert one>0 and two==2*one
    run_flow.check_pins();out=HERE/'poisson_energy';out.mkdir(exist_ok=False)
    A=[0]*12;B=[0]*12;C=[0]*12;D=[0]*12
    for p,v in zip([4,5,6,7],[1,-1,-1,1]):A[p]=v
    for p,v in zip([5,6,8,11],[-1,-1,1,1]):B[p]=v
    for p,v in zip([0,3,5,6],[-1,-1,1,1]):C[p]=v
    for p,v in zip([0,3,4,7],[-1,-1,1,1]):D[p]=v
    array={};sector_reports=[]
    for name,cell,a,b in [('first',0,A,B),('second',12,C,D)]:
        raw,failed,effects,report=sector(cell,a,b);sector_reports.append(report)
        array[name+'_sector_columns']=raw;array[name+'_processor_failure']=failed;array[name+'_effects']=effects
    records={};failures={};wordhashes={}
    for name,energy in [('single_source',one),('pair_source',two),('hold_source',F(0))]:
        rec,rf,pf,wh=clock_phase(energy);records[name]=rec;failures[name]={'reader':rf,'processor':pf};wordhashes[name]=wh
    old_one=certify.energy_interval(1260-int(one),1260)
    old_two=certify.energy_interval(1260-int(two),1260)
    for name,interval in [('fixed_count_single',old_one),('fixed_count_pair',old_two)]:
        rec,rf,pf,wh=clock_phase(sum(interval)/2);records[name]=rec;failures[name]={'reader':rf,'processor':pf};wordhashes[name]=wh
    z1,z2=phase(records['single_source']),phase(records['pair_source'])
    residual=abs(z2-z1*z1);assert residual<1e-7
    oldres=abs(phase(records['fixed_count_pair'])-phase(records['fixed_count_single'])**2)
    assert oldres>.05 and abs(phase(records['hold_source'])-1)<1e-7 and abs(z1-1)>.1
    # Detectorplant is explicitly separate from the changed-clock comparison.
    plant=np.array([[1,1j]],complex)/math.sqrt(2);plant_records,plant_failure=read(plant)
    assert abs(phase(plant_records)-1j)<1e-7
    np.savez_compressed(out/'source_sector_maps.npz',**array)
    total_sector_pulses=max(r['word_pulses'] for r in sector_reports)
    bound=nw.complete_bound(total_sector_pulses,F(1,10**8))
    reset=F(35,2147483648);pointer=F(200000000001,20000000000000000000000)
    # Two sector readers, two native clock blocks, source/flag blank resets,
    # and every source pointer have separate retained failure symbols.
    complete=2*bound['half_diamond']+20*(reset+pointer)+150*nw.PULSE_GATE
    assert complete<F(1,10000)
    result={'candidate':'A-native-poisson-compact-energy-v1','source_clock_candidate':counted['candidate'],
        'clock_receipt_sha256':sha(clock/'RECEIPT.json'),'native_single_generator':str(one),'native_pair_generator':str(two),
        'native_pair_generator_imported_from_independent_complete_source_census':True,
        'source_block_compilation_domain':'D0 and D1; signedpairD0tensorD1, allcomplements terminalunresolved',
        'coarse_sector_instruments':sector_reports,'native_pointer_records':records,
        'retained_failures':failures,'native_block_word_hashes':wordhashes,
        'measured_additivity_residual':residual,'fixed_count_changed_clock_nonadditivity':oldres,
        'plant_records':plant_records,'plant_reader_failure':plant_failure,
        'complete_sector_and_phase_halfdiamond_upper':str(complete),
        'native_pulse_applications':nw.COUNTS['pulse_applications'],'native_matrix_executions':nw.COUNTS['native_matrix_executions'],
        '511keV_used_in_construction':False,'pair_at_single_511keV_anchor':'1022keV (exactideal sourceclock)',
        'full240_UL_executed':False,'full57361_UL_executed':False,'preparation_mode':'native logarithmic escape/reinjection pump',
        'full_row_solved':False,'physical_proper_time_claimed':False}
    write(out/'RECEIPT.json',{'scientific':result,'artifacts':{'source_sector_maps.npz':sha(out/'source_sector_maps.npz')},
        'prereg_sha256':sha(HERE/'PREREG_POISSON_ENERGY.json'),'elapsed_seconds':time.monotonic()-started})
    print(json.dumps({'poisson_energy_pass':True,'additivity':residual,'old_clock_mismatch':oldres,
                      'complete_error':float(complete),'seconds':time.monotonic()-started}))

if __name__=='__main__':main()
