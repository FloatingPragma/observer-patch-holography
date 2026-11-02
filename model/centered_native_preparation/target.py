"""Exact expansion of the admitted dyadic centered-ground witness."""
from functools import lru_cache
from fractions import Fraction as F
import hashlib
import math
import numpy as np

from ..centered_particle_spectrum import packet as prior, certificate


@lru_cache(maxsize=1)
def ground():
    witness=prior.load_witness()
    result=certificate.eigensystem('ground',witness['ground_columns'],witness['ground_energies'],witness['ground_vectors'])
    block=result['block']
    values=block['full_basis'].astype(object)@result['vector_numerators'][:,0] if 'full_basis' in block else block['q'].astype(object)@result['orbit_vectors'][:,0]
    from ..centered_particle_spectrum.source import masks
    vector=[0]*65536
    for mask,value in zip(masks(4,4),values,strict=True):vector[int(mask)]=int(value)
    divisor=math.gcd(*vector)
    vector=tuple(v//divisor for v in vector)
    if next(v for v in vector if v)<0:vector=tuple(-v for v in vector)
    norm=sum(v*v for v in vector)
    if not norm:raise ArithmeticError('empty centered ground witness')
    return vector, norm, result


def ideal_hidden():
    values,norm,_=ground()
    return np.array(values,dtype=float)/math.sqrt(norm)


def summary():
    values,norm,result=ground()
    digest=hashlib.sha256(','.join(map(str,values)).encode()).hexdigest()
    return {'number':8,'hidden_up':4,'hidden_down':4,'complete_field_rows':65536,
            'nonzero_integer_amplitudes':sum(bool(v) for v in values),
            'integer_norm_squared':str(norm),'integer_amplitude_sha256':digest,
            'normalized_target_to_exact_ground_error':str(result['vector_errors'][0]),
            'ground_energy_interval':[str(x) for x in result['energy_intervals'][0]],
            'input_witness_sha256':hashlib.sha256(prior.WITNESS.read_bytes()).hexdigest(),
            'target_normalized_before_compilation':True,
            'native_success_output_renormalized':False}
