"""Complete native preparation and X-pointer readout of the finite cube ground."""
from collections import Counter
from fractions import Fraction as F
import math
from pathlib import Path
import numpy as np
from model.centered_native_preparation import tree
from model.certified_reference import certificate as cert,interval as iv
from model.native_spin_observer import circuits,execute
from . import source,certificate

HERE=Path(__file__).parent


def witness():
    with np.load(HERE/'witness.npz',allow_pickle=False) as data:
        return {k:data[k] for k in data.files}


def plan():
    proposal=witness();checked=certificate.verify(proposal)
    values=tuple(map(int,proposal['vectors'][:,0]));word=[];angle=F()
    for bit,raw,coeff in tree.design(values,8):
        a,b=tree.level_error(values,bit,raw,coeff,8);angle+=a+b
        for event in tree.level_word(bit,coeff,8):
            word.append(dict(event,targets=[2*j for j in event['targets']]))
    angle+=254*(2*iv.pi_angle_error(F(math.pi/2),F(1,2))+iv.pi_angle_error(F(-math.pi),F(-1)))
    suffix=[]
    for site in source.cube()[1]:
        suffix += [('I',(2*site,),cert._pi(-1,2)),('Z',(2*site,),cert._pi(1,2))]
    suffix+=circuits.encode(8)+circuits.probe('X',8,0,16)
    actual=list(circuits.emit(suffix));_,error,variables=cert.consume(actual,0,suffix)
    if variables:raise ArithmeticError('unbound native angle')
    word+=actual;angle+=error
    native=tree.port_bound(Counter((e['kind'],e['angle']) for e in word))
    total=angle+native+F(checked['normalized_ground_vector_error'])
    if len(word)!=1090 or total>=F(1,10**8):raise ArithmeticError('native ground preparation budget failed')
    initial=np.zeros((1<<17,1),complex);initial[0]=1
    v=np.array(values,dtype=float);v/=np.linalg.norm(v)
    odd=source.cube()[1]
    v*=np.array([(-1)**sum((mask>>j)&1 for j in odd) for mask in range(256)])
    # The pointer records bar-X eigenvalues: bit0=+1, bit1=-1.
    plus=(v+v[np.arange(256)^1])/2;minus=(v-v[np.arange(256)^1])/2
    target=np.zeros_like(initial)
    for mask in range(256):
        physical=sum(1<<(2*j+(0 if (mask>>j)&1 else 1)) for j in range(8))
        target[physical,0]=plus[mask];target[physical+(1<<16),0]=minus[mask]
    report={'pulses':len(word),'work_bits':17,'full_work_rows':1<<17,
      'word_sha256':cert.word_digest(word),'angle_error':str(angle),
      'native_generator_error':str(native),'complete_ground_isometry_error':str(total),
      'complete_ground_trace_error':str(2*total),'source_law_changed':False,
      'ground_postselected':False,'success_renormalized':False,
      'pointer_bit_zero_ideal_probability':float(np.linalg.norm(plus)**2)}
    return word,initial,target,checked,report


def run():
    word,initial,target,checked,report=plan()
    result=execute.complete(word,17,initial)
    distance=float(np.sqrt(np.linalg.norm(result.success-target)**2+result.failure_probability))
    if distance>float(F(report['complete_ground_isometry_error']))+execute.ADMISSION:
        raise ArithmeticError('complete native preparation disagrees with certified ground')
    report.update({'complete_output_distance_diagnostic':distance,
      'complete_first_failure_probability':result.failure_probability,
      'normalization_residual':result.normalization_residual,
      'native_first_failure_slots':result.native_failure_slots,
      'rounded_negative_count':result.rounded_negative_count,
      'rounded_negative_mass':result.rounded_negative_mass,
      'pointer_bit_zero_success_probability':float(np.linalg.norm(result.success[:1<<16])**2),
      'pointer_bit_one_success_probability':float(np.linalg.norm(result.success[1<<16:])**2)})
    return checked,report,{'success':result.success,'failure_weights':result.failure_weights}
