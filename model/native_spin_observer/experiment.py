"""Actual native Choi and complete local observer-history executions."""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
import numpy as np
from scipy.linalg import expm
from ..certified_reference import certificate as cert
from ..centered_native_preparation import tree
from . import algebra,circuits,execute


def plans():
    return {'encoded_spin_choi':circuits.encode(2),
            'exchange_choi':circuits.exchange(F(1,3)),
            'C_repair_choi':circuits.repair('C'),
            'E_repair_choi':circuits.repair('E'),
            'C_crossing_reference':circuits.repair('C',3,(0,2)),
            'observed_Neel_history':circuits.xgate(2)+circuits.encode(2)+circuits.exchange(F(1))+circuits.repair('C')+circuits.probe('Z',2,0,9)}


def initial_and_target(name):
    p=algebra.inclusion(2)
    if name=='encoded_spin_choi':
        initial=np.zeros((16,4),complex)
        for q in range(4):initial[(q&1)+((q>>1)<<2),q]=.5
        return 4,initial,p/2
    if name=='exchange_choi':
        h=(sum(algebra.tensor(a*2) for a in 'XYZ')-np.eye(4))/4
        return 4,p/2,p@expm(-1j*h/3)/2
    if name in ('C_repair_choi','E_repair_choi'):
        initial=np.zeros((512,16),complex);initial[:16]=np.eye(16)/4
        return 9,initial,algebra.target_repair(name[0])
    if name=='C_crossing_reference':
        fields=np.zeros((64,4),complex)
        for j,mask in enumerate((1,5,16,20)):fields[mask,j]=.5
        initial=np.zeros((2048,4),complex);initial[:64]=fields
        return 11,initial,algebra.target_repair('C',3,(0,2),fields)
    if name=='observed_Neel_history':
        initial=np.zeros((1024,1),complex);initial[0,0]=1
        h=(sum(algebra.tensor(a*2) for a in 'XYZ')-np.eye(4))/4
        state=p@expm(-1j*h)[:,2];target=np.zeros_like(initial)
        for coin in range(2):
            for mask,a in enumerate(state):target[mask+(coin<<4)+((mask&1)<<9),0]=a/np.sqrt(2)
        return 10,initial,target
    raise ValueError('declared experiment required')


@lru_cache(maxsize=1)
def analytic_errors():
    output={}
    for name,template in plans().items():
        word=list(circuits.emit(template));_,angle,variables=cert.consume(word,0,template)
        if variables:raise ArithmeticError('unaccounted free angle')
        native=tree.port_bound(Counter((e['kind'],e['angle']) for e in word))
        output[name]={'angle_error':str(angle),'port_generator_error':str(native),
                      'complete_isometry_error':str(angle+native),'complete_trace_error':str(2*(angle+native)),
                      'pulses':len(word),'word_sha256':cert.word_digest(word)}
    return output


def pointer_probabilities(name,state):
    if name in ('C_repair_choi','E_repair_choi','C_crossing_reference'):
        sites=3 if name=='C_crossing_reference' else 2;dim=1<<(2*sites)
        return [{'coin':coin,'raw':raw,'probability':float(np.sum(abs(state[(raw*2+coin)*dim:(raw*2+coin+1)*dim])**2))}
                for coin in range(2) for raw in range(16)]
    if name=='observed_Neel_history':
        return [{'coin':coin,'raw':raw,'spin_pointer':pointer,
                 'probability':float(np.sum(abs(state[(pointer*32+raw*2+coin)*16:(pointer*32+raw*2+coin+1)*16])**2))}
                for coin in range(2) for raw in range(16) for pointer in range(2)]
    return []


def report_case(name,run,target):
    error=F(analytic_errors()[name]['complete_isometry_error'])
    distance=float(np.sqrt(np.linalg.norm(run.success-target)**2+run.failure_probability))
    if distance>float(error)+execute.ADMISSION:raise ArithmeticError('native source misses certified target')
    record={'work_bits':run.bits,'work_rows':len(run.success),'quantum_reference_dimension':run.success.shape[1],
            'native_pulses':len(run.failure_weights),'native_first_failure_slots':run.native_failure_slots,
            'word_sha256':run.word_sha256,'complete_output_error_diagnostic':distance,
            'complete_first_failure_probability':run.failure_probability,'normalization_residual':run.normalization_residual,
            'rounded_negative_count':run.rounded_negative_count,'rounded_negative_mass':run.rounded_negative_mass,
            'logical_pointer_records':pointer_probabilities(name,run.success),
            'target_pointer_records':pointer_probabilities(name,target),
            'reference_columns_dephased':False,'success_postselected_or_renormalized':False,
            'native_failure_weights_archive':name+'_failure_weights','success_archive':name+'_success'}
    pointer_bits=6 if name=='observed_Neel_history' else 5 if record['logical_pointer_records'] else 0
    record['physical_pointer_outcomes']=6**pointer_bits if pointer_bits else 0
    record['reserved_pointer_outcomes_with_zero_probability']=6**pointer_bits-2**pointer_bits if pointer_bits else 0
    return record


def run(progress=None):
    records={};arrays={}
    for name,template in plans().items():
        bits,initial,target=initial_and_target(name)
        actual=execute.complete(circuits.emit(template),bits,initial)
        records[name]=report_case(name,actual,target)
        arrays[name+'_success']=actual.success;arrays[name+'_failure_weights']=actual.failure_weights
        if progress is not None:progress(name,records[name])
    return records,arrays
