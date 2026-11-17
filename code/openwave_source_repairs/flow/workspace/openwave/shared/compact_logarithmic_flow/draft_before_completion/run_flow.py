#!/usr/bin/env python3
"""Complete native-source logarithmic-flow experiment; no ideal propagator."""
from fractions import Fraction as F
import argparse,hashlib,json,math
from pathlib import Path
import time
import numpy as np
import certify
import native_words as nw
import reference
from model.native_interface import pointer
from model.observer.native import native

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
PROVIDER=HERE.parent/'localized_rest_join'

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def serial(value):
    if isinstance(value,F):return str(value)
    if isinstance(value,np.generic):return value.item()
    raise TypeError(type(value).__name__)
def write(p,x):p.write_text(json.dumps(x,sort_keys=True,indent=2,default=serial,allow_nan=False)+'\n')
def progress(x):print(json.dumps(x,sort_keys=True,default=serial),flush=True)

def check_pins():
    for name in ('SOURCE_PINS.json','PREREG_PROPOSAL.json','SOURCE_PINS_FLOW.json','PREREG_FLOW.json'):
        for path,digest in json.loads((HERE/name).read_text())['files'].items():
            assert sha(ROOT/path)==digest,path

def recheck_source(source,proposal):
    table={};calls=0
    for a in (0,1):
        for b in (0,1):
            for coin in (0,1):
                table[a,b,coin]=native((a,b),(0,1),coin);calls+=1
    V=len(source['counts']);labels=np.arange(V);C=np.zeros((V,V),dtype=np.int64);D=np.zeros_like(C)
    held=proposal['held_seam_index']
    assert held==next(i for i,(a,b) in enumerate(source['seams']) if (int(a) in {4,5,6,7})!=(int(b) in {4,5,6,7}))
    for ei,(a,b) in enumerate(source['seams']):
        a,b=int(a),int(b)
        for coin in (0,1):
            q=[];j=[]
            for p in labels:
                out,record=table[int(p==a),int(p==b),coin]
                q.append(int(p) if p not in (a,b) else (a if out[0] else b));j.append(record[1])
            q=np.asarray(q);j=np.asarray(j);ec=2*ei+coin
            qc=labels if ei==held else q;jc=np.zeros(V,dtype=int) if ei==held else j
            assert np.array_equal(q,source['destination'][:,ec]) and np.array_equal(j,source['current'][:,ec])
            assert np.array_equal(qc,source['changed_destination'][:,ec]) and np.array_equal(jc,source['changed_current'][:,ec])
            np.add.at(C,(q,labels),1);np.add.at(D,(qc,labels),1)
    assert np.array_equal(C,source['counts']) and np.array_equal(D,source['changed_counts'])
    return calls

def run_law(law,times,proposal,initial,certificate):
    rows=proposal['laws'][law];dimension=proposal['source_dimension'];padded=initial.shape[0];bits=padded.bit_length()-1
    states=np.concatenate([initial.copy() for _ in times],axis=1)
    failures=np.zeros((len(times),2,2),complex);padding=[];words=[];phase_errors=[F() for _ in times]
    total_pulses=0;tree_errors=F();diag_count=None
    for j,vector in enumerate(rows['columns']):
        integers=list(vector)+[0]*(padded-dimension)
        preparation,tree_error,worst=nw.preparation(integers)
        inverse=nw.c.inverse(preparation)
        states,bad=nw.apply_word(inverse,states,bits)
        for k in range(len(times)):failures[k]+=bad[2*k:2*k+2,2*k:2*k+2]
        phase_metadata=[]
        for k,t in enumerate(times):
            angle,error=certify.phase(rows['eigenvalues'][j],t)
            diagonal=nw.address_phase(angle,bits)
            partial,bad=nw.apply_word(diagonal,states[:,2*k:2*k+2],bits)
            states[:,2*k:2*k+2]=partial;failures[k]+=bad
            phase_errors[k]+=error
            phase_metadata.append({'time':t,'angle_hex':angle.hex(),'source_word_sha256':nw.digest_word(diagonal),
                                   'scalar_log_phase_error':str(error)})
            if diag_count is None:diag_count=len(diagonal)
            assert len(diagonal)==diag_count
        states,bad=nw.apply_word(preparation,states,bits)
        for k in range(len(times)):failures[k]+=bad[2*k:2*k+2,2*k:2*k+2]
        padding.append(states[dimension:].copy())
        states[dimension:]=0 # Explicit retained source-padding abort, never renormalization.
        pulses=2*len(preparation)+diag_count;assert pulses<=65536
        total_pulses+=pulses;tree_errors+=4*tree_error
        words.append({'ray':j,'preparation_pulses':len(preparation),'phase_word_pulses':pulses,
                      'preparation_word_sha256':nw.digest_word(preparation),'inverse_word_sha256':nw.digest_word(inverse),
                      'tree_operator_error':str(tree_error),'max_split_component_error':str(worst),'phases':phase_metadata})
        if j%12==0:progress({'law':law,'source_rank_one_words_executed':j+1,'of':dimension})
    reports=[]
    for k,t in enumerate(times):
        bound=nw.complete_bound(total_pulses,tree_errors+phase_errors[k])
        spectral=abs(t)*F(certificate['energy_operator_upper'])+F(certificate['phase_projector_product_upper'])
        complete=bound['half_diamond']+bound['gain']*spectral
        assert complete<F(1,5000)
        reports.append({'time':t,'native_pulses':total_pulses,'word_count':dimension,
                        'native_word_complete_halfdiamond_upper':str(bound['half_diamond']),
                        'spectral_product_operator_upper':str(spectral),'complete_flow_halfdiamond_upper':str(complete),
                        'stored_map_gain_upper':str(bound['gain']),'processor_failure_upper':str(bound['processor_failure']),
                        'target_tree_phase_error_upper':str(tree_errors+phase_errors[k])})
    return states,failures,np.asarray(padding),{'law':law,'experiments':reports,'words':words}

def read_complete(state,law):
    local=reference.pointer_effects(law);padded=state.shape[0];bits=padded.bit_length()-1
    outcomes=[]
    for label in range(padded):
        effect=np.ones((1,1),complex)
        for bit in reversed(range(bits)):effect=np.kron(effect,local[(label>>bit)&1])
        outcomes.append(state.conj().T@effect@state)
    return np.asarray(outcomes)

def main(record=False):
    started=time.monotonic();check_pins()
    assert record,'full producer is explicit --record; use headline.py for read-only review'
    if (HERE/'RECEIPT.json').exists():raise FileExistsError('RECEIPT.json')
    evidence=HERE/'evidence';evidence.mkdir(exist_ok=True)
    proposal,source=certify.load();fresh_calls=recheck_source(source,proposal)
    certificates={law:certify.certify(law) for law in ('native','held_seam')}
    write(evidence/'SPECTRAL_CERTIFICATES.json',certificates)
    pointer_certificate=pointer.certificate()
    progress({'stage':'source and exact full spectral certificates checked'})
    ref_report,ref_arrays=reference.reference_preparation(source,progress)
    write(evidence/'REFERENCE_PREPARATION.json',ref_report)
    np.savez_compressed(evidence/'reference_preparation.npz',**ref_arrays)
    progress({'stage':'complete source swap-parity reference preparation checked'})
    inherited=np.load(PROVIDER/'evidence'/'raw_source_arrays.npz')
    initial=np.column_stack((inherited['U_f'][:,0],inherited['U_g'][:,0]))
    outputs={};flow_reports={};all_arrays={}
    for law,times in [('native',[1,4]),('held_seam',[4])]:
        out,failure,padding,report=run_law(law,times,proposal,initial,certificates[law])
        write(evidence/(law+'_WORDS.json'),report)
        all_arrays[law+'_terminal_source']=out;all_arrays[law+'_processor_failure_reference_grams']=failure
        all_arrays[law+'_padding_abort_amplitudes']=padding
        outputs[law]=(times,out);flow_reports[law]=report['experiments']
    pointer_records={};diagnostics={};E6lo,E6hi=certify.energy_interval(209,210)
    E6=float((E6lo+E6hi)/2)
    for law,(times,out) in outputs.items():
        for k,t in enumerate(times):
            state=out[:,2*k:2*k+2];key=law+'_t'+str(t)
            for read_law in ('C','E'):
                position=read_complete(state,read_law)
                energy=read_complete(inherited['inverse_f']@state,read_law)
                all_arrays[key+'_position_'+read_law]=position
                all_arrays[key+'_energy_ray_'+read_law]=energy
                outside=float(sum(position[z,0,0].real for z in range(12,256)))
                click=float(energy[0,0,0].real)
                pointer_records[key+'_'+read_law]={'outside_original_cell':outside,'energy_ray_click':click,
                   'f_energy_unresolved':float(energy[1:,0,0].real.sum()),
                   'bell_energy_unresolved':float(np.trace(np.sum(energy[1:],axis=0)).real/2),
                   'f_position_total_mass':float(position[:,0,0].real.sum()),
                   'bell_total_record_mass':float(np.trace(np.sum(position,axis=0)).real/2)}
                if law=='native':assert outside<1e-7 and click>1-1e-7
                else:assert outside>1e-3
            if law=='native':
                phase=complex(math.cos(t*E6),-math.sin(t*E6))
                residual=float(np.linalg.norm(state-phase*initial,'fro')/math.sqrt(2))
                assert residual<1e-7
                diagnostics[key]={'joint_bell_vector_phase_residual':residual,'time':t}
    plant_word=[nw.c.pulse('I',-math.pi/2,[4]),nw.c.pulse('X',math.pi/2,[4])]
    plant,plant_failure=nw.apply_word(plant_word,initial,8)
    plant_law=read_complete(plant,'C');plant_escape=float(plant_law[12:,0,0].real.sum())
    assert plant_escape>1-1e-7
    all_arrays['planted_source_escape']=plant;all_arrays['planted_position_C']=plant_law
    all_arrays['planted_processor_failure_reference_gram']=plant_failure
    all_arrays['initial_source_columns']=initial
    np.savez_compressed(evidence/'flow_records.npz',**all_arrays)
    keys,matrices,bounds=nw.gate_archive()
    np.savez_compressed(evidence/'native_M6_primitives.npz',matrices=matrices)
    write(evidence/'NATIVE_PRIMITIVES.json',{'keys':keys,'exact_squared_error_bounds':bounds,
          'raw_operator_gate':str(nw.RAW_GATE),'pulse_product_gate':str(nw.PULSE_GATE),'counts':nw.COUNTS})
    science={'candidate':'A-signed-compact-log-flow-v1','original_H_retained':True,
             'source_dimension':240,'invalid_source_labels':16,'source_denominator':1260,
             'source_trace_T':239,'fresh_native_endpoint_cases':fresh_calls,
             'source_native_columns_from_registered_proposal':302400,
             'fresh_selector_native_calls':ref_report['selector_native_calls'],
             'spectral_certificates':certificates,'flow_certificates':flow_reports,
             'complete_pointer_records':pointer_records,'diagnostics':diagnostics,
             'planted_escape_probability':plant_escape,'native_processor_counts':nw.COUNTS,
             'reference_preparation_halfdiamond_upper':ref_report['complete_preparation_halfdiamond_upper'],
             'reference_plus_probability_given_selectors':'239/240','reference_minus_probability_given_selectors':'1/240',
             'pointer_complete_certificate':pointer_certificate,
             'inherited_blank_preparation_and_terminal_read_error_bound':'6/10000000',
             'ideal_all_real_time_compact_ray_persistence':True,
             'all_time_stored_floating_word_error_uniformly_bounded':False,
             'physical_clock_calibration_claimed':False,'row3_self_start_claimed':False,
             'strict_official_status':'untested; unchanged','external_review':'pending; no G6'}
    artifacts={str(p.relative_to(HERE)):sha(p) for p in sorted(evidence.iterdir()) if p.is_file()}
    receipt={'scientific':science,'artifacts':artifacts,'prereg_sha256':sha(HERE/'PREREG_FLOW.json'),
             'elapsed_seconds':time.monotonic()-started}
    write(HERE/'RECEIPT.json',receipt)
    progress({'stage':'finished','elapsed_seconds':receipt['elapsed_seconds'],'native_processor_counts':nw.COUNTS,
              'primary_records':pointer_records,'planted_escape_probability':plant_escape})

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--record',action='store_true');args=ap.parse_args();main(args.record)
