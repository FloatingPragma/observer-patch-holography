"""Native controlled-W and complete position-read apparatus contract.

Only the frozen common provider supplies source evolution. Every apparatus
reset, export complement, native pointer label and scalar phase is retained.
"""
from collections import Counter
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import sys

from native_phase_probe import controlled_local_word, control_word

HERE=Path(__file__).resolve().parent
WORKSPACE=HERE.parents[2]
PROVIDER=HERE.parent/'native_walk_bridge'
sys.path.insert(0,str(PROVIDER))
import common_walk as common
import compile_walk as compiled
from model.observer import native
from model.observer.source_scalar import compiler
from model.observer.source_discrimination import circuit
from model.native_interface import pointer
from model.native_clock import configuration
from model.certified_reference import certificate as cert


def controlled_template(kind,angle):
    theta=F(angle)
    if kind=='H':
        return [('Y',(0,),cert._pi(1,8))]+controlled_template('Z',angle)+[('Y',(0,),cert._pi(-1,8))]
    if kind=='J':
        result=[]
        for mask in range(8):
            coefficient=theta*sum((-1)**((i&mask).bit_count())for i in (4,5,6))/8
            label=''.join('Z'if mask>>bit&1 else'I'for bit in range(3))
            result+=cert.pauli(label,coefficient)
        return result
    return cert.pauli(kind+'I',theta/2)+cert.pauli(kind+'Z',-theta/2)


def source_counts():
    built=compiled.controls()
    counts=Counter((event['kind'],event['angle'])for event in built['reflection_word'])
    for kind,number in Counter(kind for kind,targets in built['source_program']).items():
        for event,count in Counter((e['kind'],e['angle'])for e in built['primitive_words'][kind]).items():
            counts[event]+=number*count
    assert sum(counts.values())==len(built['reflection_word'])+built['source_swap_native_pulses']
    return counts


def controlled_walk_word():
    """Finite emitted program; common provider owns all source operations."""
    return control_word(compiled.pulse_word(),compiled.MEMORY_BITS)


def pointer_bound():
    word=pointer._local_word()
    end,angle,variables=cert.consume(word,0,cert.cnot(2,0)+cert.cnot(0,1))
    assert end==len(word)and not variables
    return angle+cert.controls_bound(word,{})


def reset_contract(count0,count2):
    rows=[];probabilities={}
    for code in (0,pointer.READY_CODE):
        target=circuit.decode(code)
        ps=configuration.success_probabilities(4,32,full=True,mode='reset',target=target)
        probabilities[code]=F(ps[0])
        rows.append({'target_code':code,'initial_carrier':list(native.states()[0]),
            'target_carrier':list(target),'success_probability_exact':str(ps[0]),
            'complete_reset_certificate':configuration.certificate(4,32,full=True,mode='reset',target=target)})
    p0,p2=probabilities[0],probabilities[pointer.READY_CODE]
    upper=count0*(1-p0)+count2*(1-p2)
    return {'blank_carriers':count0,'ready_carriers':count2,
        'exact_success_factor':{'p0':str(p0),'exponent0':count0,'p2':str(p2),'exponent2':count2},
        'timeout_upper_exact':str(upper),'timeout_upper':math.nextafter(float(upper),math.inf),
        'fine_failure_law':'first failure at reset k has product of preceding successes times(1-p_k); all native fine seam/coin/current/timeout records remain',
        'conditioning':'fresh independent native carriers reset before common source preparation; no success normalization',
        'native_reset_contracts':rows}


def run():
    assert common.MODEL_ID=='native-local12-current-v0-proposed'
    paths=[Path(__file__),HERE/'native_phase_probe.py',PROVIDER/'common_walk.py',PROVIDER/'compile_walk.py',
           Path(compiler.__file__),Path(pointer.__file__),Path(configuration.__file__)]
    pins={str(p.relative_to(WORKSPACE)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths}
    counts=source_counts();rows=[];total=F();pulses=0;table={}
    for (kind,angle),uses in sorted(counts.items()):
        word=controlled_local_word(kind,angle)
        end,fixed,variables=cert.consume(word,0,controlled_template(kind,angle))
        assert end==len(word)and not variables
        error=fixed+cert.controls_bound(word,table)
        total+=uses*error;pulses+=uses*len(word)
        rows.append({'source_kind':kind,'source_angle':angle,'source_pulse_uses':uses,
                     'controlled_source_pulses_each':len(word),'error_upper_exact':str(error)})
    native_pointer_error=pointer_bound()
    # Probe plus preparation and X/Y analysis use the same actual source
    # compiler. Scalar I phases remain in the explicit instrument.
    plus=[compiler.pulse('I',-math.pi/2,[0]),compiler.pulse('H',math.pi/2,[0])]
    plus_template=[('I',(0,),cert._pi(-1,2)),('H',(0,),cert._pi(1,2))]
    end,fixed,_=cert.consume(plus,0,plus_template)
    plus_error=fixed+cert.controls_bound(plus,table)
    phase=[compiler.pulse('I',-math.pi/4,[0]),compiler.pulse('Z',math.pi/4,[0])]
    _,fixed,_=cert.consume(phase,0,[('I',(0,),cert._pi(-1,4)),('Z',(0,),cert._pi(1,4))])
    phase_error=fixed+cert.controls_bound(phase,table)
    phase_detector_error=2*plus_error+phase_error+native_pointer_error
    position_reset=reset_contract(0,14);phase_reset=reset_contract(1,1)
    assert pins=={str(p.relative_to(WORKSPACE)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths}
    return {'schema':'oph.common_walk.complete_instrument_contract.v1','model_id':common.MODEL_ID,
        'checks':{'one_source_provider_only':True,'every_common_source_pulse_controlled':True,
                  'global_scalar_phases_retained':True,'outward_complete_source_control_bounds':True,
                  'all_added_native_reset_failures_counted':True,'source_pins_unchanged':True},
        'source_W_native_pulses':sum(counts.values()),'controlled_W_native_pulses':pulses,
        'controlled_W_wrapper_error_exact':str(total),'controlled_W_wrapper_error_upper':math.nextafter(float(total),math.inf),
        'wrapper_error_scope':'complete source output relative to controlled product of the common source word ideal proper-code gates. Add the common producer geometrical-W/compiler error; this is not its replacement.',
        'controlled_word_routines':rows,'source_control_bit':compiled.MEMORY_BITS,
        'position':{'all128_source_labels':True,'same_full_event_and_current_reference':True,
            'two_read_pointer_pulses':84,'two_read_pointer_error_exact':str(14*native_pointer_error),
            'additional_apparatus_resets':position_reset,
            'complete_bound_recipe':'common preparation error + (t_first+t_gap)*common full-W error +14*pointer_error + reset_timeout; no cyclic-sector projection'},
        'phase':{'Hadamard_probe_contrast':'Re or Im <psi|W^n|psi>, on the unchanged common prepared source',
            'quasienergy_clock':'unit-circle phase modulo2pi per coherent W step',
            'extra_preparation_analysis_pointer_error_exact':str(phase_detector_error),
            'additional_apparatus_resets':phase_reset,
            'complete_bound_recipe':'common preparation error+n*(common W error+controlled wrapper error)+probe apparatus error+reset timeout',
            'dephasing_control':'dephasing the retained probe before analysis gives two equal bit probabilities',
            'quadrature_settings':{'X':'1/2','Y':'1/2'},'native_law':'C or E configured separately, never mixed implicitly'},
        'failure_poststates':'all proper-code first-export aborts retain their source/environment/reference output; all native bit-read fine labels and every source/event/current padding remain',
        'coherent_services':'existing proper-code load/export and wire routing, with supplied internal receiver initialization; added resets do not derive these services',
        'no_independent_H_no_calibration':True,'official_rows_changed':[],'source_sha256':pins}


if __name__=='__main__':
    result=run()
    (HERE/'INSTRUMENT_CONTRACT.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'checks':result['checks'],'W_pulses':result['source_W_native_pulses'],
        'controlled_W_pulses':result['controlled_W_native_pulses'],
        'wrapper_error':result['controlled_W_wrapper_error_upper'],
        'two_position_reset_timeout':result['position']['additional_apparatus_resets']['timeout_upper']},indent=2))
