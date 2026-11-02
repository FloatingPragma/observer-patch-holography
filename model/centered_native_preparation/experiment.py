"""Complete full-row native execution, without postselection normalization."""
from functools import lru_cache
import json
from pathlib import Path
import time
import numpy as np

from ..ground_spectral_source import preparation as frames
from ..observer.source_scalar import compiler
from ..cube_eta_protection import cube
from . import target,tree,execute

DIRECTORY=Path(__file__).resolve().parent


def frame_word():
    for x,site in enumerate(cube.GEOMETRY.sites):
        yield from compiler.remap(frames.emit(frames.inverse_frame_plan(site)),[2*x,2*x+1])


def ideal_ground():return cube.operators().frame.matrix().conj().T@target.ideal_hidden()


def preparation(progress=None):
    blank=np.zeros(65536);blank[0]=1
    hidden=execute.complete(tree.word(),16,blank,progress)
    physical=execute.complete(frame_word(),16,hidden.success)
    abort=hidden.failure_probability+physical.failure_probability
    diagnostic=float(np.sqrt(np.linalg.norm(physical.success[:,0]-ideal_ground())**2+abort))
    result={'hidden_success':hidden.success,'prepared_ground':physical.success,
            'tree_failure_weights':hidden.failure_weights,'frame_failure_weights':physical.failure_weights,
            'tree_block_offsets':np.array([x[:2] for x in hidden.schedule],dtype=np.int64)}
    summary={'native_tree_blocks':len(hidden.schedule),'native_tree_pulses':len(hidden.failure_weights),
             'native_frame_pulses':len(physical.failure_weights),'all_rows':65536,
             'tree_word_sha256':hidden.word_sha256,'frame_word_sha256':physical.word_sha256,
             'tree_native_failure_slots':hidden.native_failure_slots,'frame_native_failure_slots':physical.native_failure_slots,
             'tree_first_failure_probability':hidden.failure_probability,
             'frame_first_failure_probability':physical.failure_probability,
             'complete_preparation_failure_probability':abort,
             'tree_normalization_residual':hidden.normalization_error,
             'frame_normalization_residual':physical.normalization_error,
             'complete_preparation_error_diagnostic':diagnostic,
             'rounded_negative_count':hidden.rounded_negative_count+physical.rounded_negative_count,
             'rounded_negative_mass':hidden.rounded_negative_mass+physical.rounded_negative_mass}
    return result,summary


def main():
    start=time.monotonic()
    def progress(block,pulses):
        if block%1024==0:print(json.dumps({'blocks':block,'pulses':pulses,'seconds':time.monotonic()-start}),flush=True)
    arrays,summary=preparation(progress)
    np.savez_compressed(DIRECTORY/'preparation.npz',**arrays)
    (DIRECTORY/'preparation.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary),flush=True)


if __name__=='__main__':main()


def charge_and_clock(prepared,progress=None):
    from ..native_spectral_detector import execution as native
    from ..native_spectral_detector import circuits as old
    from . import circuits
    state=np.zeros((131072,1),complex);state[:65536]=prepared
    charged=native.complete(circuits.charge_blocks(),17,state)
    branches=np.column_stack((charged.success[:65536,0],charged.success[65536:,0]))
    clock_input=np.zeros((4*65536,2),complex);clock_input[:65536]=branches
    clock=native.complete(circuits.clock_blocks(),18,clock_input,progress)
    probabilities=np.sum(abs(clock.success.reshape(4,65536,2))**2,axis=1)
    return {'charge_success':charged.success,'charge_failure_weights':charged.failure_weights,
            'clock_success':clock.success,'clock_failure_weights':clock.failure_weights,
            'pointer_probabilities':probabilities}, {
            'charge_failure_probability':charged.failure_probability,
            'clock_failure_probability':clock.failure_probability,
            'charge_normalization_residual':charged.normalization_error,
            'clock_normalization_residual':clock.normalization_error,
            'clock_pulses':len(clock.failure_weights),'clock_blocks':len(clock.blocks),
            'clock_full_rows':4*65536,'classical_charge_columns':2,
            'pointer_probabilities':probabilities.tolist(),
            'rounded_negative_count':charged.rounded_negative_weight_count+clock.rounded_negative_weight_count,
            'rounded_negative_mass':charged.rounded_negative_weight_total+clock.rounded_negative_weight_total}


def clock_main():
    start=time.monotonic()
    with np.load(DIRECTORY/'preparation.npz',allow_pickle=False) as archive:
        prepared=archive['prepared_ground']
    def progress(block,total):
        if block%1024==0:print(json.dumps({'clock_blocks':block,'total':total,'seconds':time.monotonic()-start}),flush=True)
    arrays,summary=charge_and_clock(prepared,progress)
    np.savez_compressed(DIRECTORY/'clock.npz',**arrays)
    (DIRECTORY/'clock.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary),flush=True)


def pointer_records(state):
    if state.shape!=(262144,2):raise ValueError('complete four-bin field and two charge columns required')
    probabilities=np.sum(abs(state.reshape(4,65536,2))**2,axis=1)
    return [{'charge_record':tag,'clock_record':[low,high],
             'probability':float(probabilities[low+2*high,tag]) if tag<2 and low<2 and high<2 else 0.}
            for tag in range(6) for low in range(6) for high in range(6)]


def ideal_clock():
    """Independent exact-matrix evolution of the certified dyadic target."""
    from scipy.sparse.linalg import expm_multiply
    from ..native_spectral_detector.experiment import ideal_charge
    fields=ideal_charge(ideal_ground()[:,None])
    h=cube.operators().H.matrix(cube.operators().scale)
    lo,hi=target.ground()[2]['energy_intervals'][0];energy=float((lo+hi)/2)
    evolved=np.zeros((4,65536,2),complex)
    for tag,number in enumerate((9,7)):
        masks=cube.weights(number);restricted=h[masks][:,masks]
        for j in range(4):
            evolved[j,masks,tag]=expm_multiply(-1j*j*restricted,fields[masks,tag])*np.exp(1j*j*(number/2+energy))/2
    fourier=np.exp(2j*np.pi*np.arange(4)[:,None]*np.arange(4)[None,:]/4)/2
    return np.tensordot(fourier,evolved,axes=(1,0)).reshape(262144,2)


def report_from(arrays,prep,clock):
    from fractions import Fraction as F
    from . import bounds,circuits
    output=arrays['clock_success'];pointer=pointer_records(output)
    stop=(prep['complete_preparation_failure_probability']+clock['charge_failure_probability']+clock['clock_failure_probability'])
    ideal=ideal_clock();ideal_pointer=pointer_records(ideal)
    deviation=float(np.sqrt(np.linalg.norm(output-ideal)**2+stop))
    variation=(sum(abs(a['probability']-b['probability']) for a,b in zip(pointer,ideal_pointer,strict=True))+stop)/2
    error=F(bounds.certificate()['complete_state_and_total_variation_error'])
    total=float(np.vdot(output,output).real)+stop
    if abs(total-1)>execute.ADMISSION or max(deviation,variation)>float(error)+execute.ADMISSION:
        raise ArithmeticError('complete native detector exceeds the source-derived bound')
    for row in bounds.spectral_probabilities():
        index=row['clock_index'];tag=row['charge_record'];actual=clock['pointer_probabilities'][index][tag]
        lo,hi=map(F,row['native_probability_interval'])
        if not float(lo)-execute.ADMISSION<=actual<=float(hi)+execute.ADMISSION:
            raise ArithmeticError('native probability escapes the certified spectral enclosure')
    coefficients=np.concatenate([c for _,_,c in tree.ground_design()]);raw=np.concatenate([a for _,a,_ in tree.ground_design()])
    arrays['literal_Gray_coefficients']=coefficients;arrays['checked_conditional_angles']=raw
    return {'preparation':prep,'clock':clock,'all_charge_clock_outcomes':pointer,
            'independent_matrix_evolution_outcomes':ideal_pointer,
            'independent_reference':'full original-H matrix evolution of the certified normalized dyadic ground approximation, with centered number phase',
            'complete_first_failure_probability':stop,'complete_probability':total,
            'complete_output_error_diagnostic':deviation,'total_variation_diagnostic':variation,
            'host_arithmetic_admission':execute.ADMISSION,
            'full_centered_ground_preparation_executed':True,
            'complete_M4_centered_phase_readout_executed':True,
            'source_clock_centering_explicit':True,
            'coherent_tag_reused_only_after_classical_readout':True,
            'all228_lines_separately_resolved':False,
            'spectral_or_number_projection_used_in_native_run':False,
            'native_success_postselected_or_renormalized':False,
            'autonomous_ground_selection_or_bare_repair_protection_proved':False,
            'physical_particle_species_or_masses_identified':False,
            'finite_shot_outcome_samples_generated':False}


def run(progress=None):
    arrays,prep=preparation(None if progress is None else lambda b,p:progress('preparation',b,p))
    clock_arrays,clock=charge_and_clock(arrays['prepared_ground'],None if progress is None else lambda b,p:progress('clock',b,p))
    arrays.update(clock_arrays)
    return report_from(arrays,prep,clock),arrays


def block_events(blocks):
    for instrument,targets in blocks:
        for kind,angle,local_targets in instrument.tape:
            yield compiler.pulse(kind,angle,[targets[t] for t in local_targets])


def programs():
    from . import circuits
    return {'tree_failure_weights':tree.word(),'frame_failure_weights':frame_word(),
            'charge_failure_weights':block_events(circuits.charge_blocks()),
            'clock_failure_weights':block_events(circuits.clock_blocks())}


def program_metadata():
    import hashlib
    output={}
    for name,events in programs().items():
        digest=hashlib.sha256();pulses=0;pairs=0
        for event in events:
            pulses+=1;pairs+=event['kind']=='J'
            digest.update((event['kind']+'|'+event['angle'].hex()+'|'+','.join(map(str,event['targets']))+';').encode())
        output[name]={'native_pulses':pulses,'word_sha256':digest.hexdigest(),
                      'native_first_failure_slots':4*pulses-2*pairs,'padding_slots':2*pairs}
    return output
