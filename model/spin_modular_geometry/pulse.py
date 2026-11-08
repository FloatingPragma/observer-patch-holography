"""Native current pulse on the complete selected eight-spin source."""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
import math

import numpy as np

from model.certified_reference import certificate as cert
from model.native_spin_observer import algebra, circuits
from model.native_ordered_ground import source, experiment
from . import native, regional


TABLE = {('X', 'Y'): ('Z', 1j), ('Y', 'Z'): ('X', 1j), ('Z', 'X'): ('Y', 1j),
         ('Y', 'X'): ('Z', -1j), ('Z', 'Y'): ('X', -1j), ('X', 'Z'): ('Y', -1j)}


def multiply(a, b):
    phase, result = 1, []
    for x, y in zip(a, b, strict=True):
        if x == 'I': result.append(y)
        elif y == 'I': result.append(x)
        elif x == y: result.append('I')
        else:
            value, factor = TABLE[x, y]
            result.append(value)
            phase *= factor
    return ''.join(result), phase


def current_terms():
    edges, odd = source.cube()
    h, e = defaultdict(F), defaultdict(F)
    for x, y in edges:
        for axis in 'XYZ':
            label = ['I']*8
            label[x] = label[y] = axis
            h[''.join(label)] += F(1, 4)
            if 0 in (x, y): e[''.join(label)] += F(1, 8)
    for x in range(8):
        label = ['I']*8
        label[x] = 'X'
        coefficient = F(1 if x in odd else -1, 8)
        h[''.join(label)] += coefficient
        if x == 0: e[''.join(label)] += coefficient
    result = defaultdict(F)
    for a, ca in h.items():
        for b, cb in e.items():
            label, phase = multiply(a, b)
            result[label] += -2*int(complex(phase).imag)*ca*cb
    result = {k: v for k, v in sorted(result.items()) if v}
    matrix = sum((float(v)*algebra.tensor(k) for k, v in result.items()),
                 np.zeros((256, 256), complex))
    if not np.array_equal(matrix, 1j*regional.state()['current_imaginary_integer']/64):
        raise ArithmeticError('source-current Pauli synthesis changed the physical generator')
    return result


@lru_cache(maxsize=1)
def run():
    parameter = F(1, 10000)
    terms = current_terms()
    ordered = list(terms.items())+list(reversed(list(terms.items())))
    template = [event for label, value in ordered
                for event in cert.pauli(''.join(axis+('X' if axis in 'XY' else 'I')
                                                for axis in label), -parameter*value/2)]
    data = regional.state()
    state = np.asarray(data['vector'], dtype=complex)/math.sqrt(data['vector_norm_squared'])
    target_vector = state.copy()
    for label, value in ordered:
        angle = -float(parameter*value/2)
        target_vector = (math.cos(angle)*target_vector
                         + 1j*math.sin(angle)*(algebra.tensor(label)@target_vector))
    indices = [sum(1 << (2*j+(0 if mask >> j & 1 else 1)) for j in range(8))
               for mask in range(256)]
    initial = np.zeros((1 << 16, 1), complex)
    target = np.zeros_like(initial)
    initial[0, 0] = 1
    target[indices, 0] = target_vector
    full_word, _, _, _, preparation = experiment.plan()
    probe = circuits.probe('X', 8, 0, 16)
    end, _, variables = cert.consume(full_word, len(full_word)-len(probe), probe)
    if end != len(full_word) or variables:
        raise ArithmeticError('source ground preparation suffix is not the declared readout')
    prefix = full_word[:-len(probe)]
    if any(max(event['targets']) >= 16 for event in prefix):
        raise ArithmeticError('ground preparation prefix uses an undeclared register')
    actual, report = native.complete(template, 16, initial, target,
                                     prefix_word=prefix,
                                     prefix_error=F(preparation['complete_ground_isometry_error']))
    protected = actual.success[indices, 0].reshape(64, 4).T
    offcode = np.ones(1 << 16, dtype=bool)
    offcode[indices] = False
    discarded_weight = float(np.linalg.norm(actual.success[offcode])**2)+actual.failure_probability
    decoded = protected@protected.conj().T+discarded_weight*np.eye(4)/4
    spectrum = np.linalg.eigvalsh(decoded)
    reference_spectrum = np.linalg.eigvalsh(np.asarray(data['rho'], float))
    if np.min(spectrum) <= 0:
        raise ArithmeticError('complete decoded native state is not faithful')
    decrease = float(np.sum(spectrum*np.log(spectrum))
                     - np.sum(reference_spectrum*np.log(reference_spectrum)))
    if decrease <= 4e-6:
        raise ArithmeticError('complete decoded native entropy diagnostic lost its sign')
    norm_sum = sum(abs(v) for v in terms.values())
    splitting = norm_sum**3*parameter**3/3
    error = (data['ground_vector_error']+splitting
             + F(report['complete_isometry_error_bound']))
    # Decode every failure and every off-code component into one fixed regional
    # state. Contractivity bounds the complete decoded output without conditioning.
    ideal_floor = data['actual_floor']-parameter*data['current_norm_upper']
    entropy_error = 2*error*F(5)  # -log eigenvalue norm <5 throughout this segment
    if error >= ideal_floor-F(1, 100) or entropy_error >= F(1, 10**7):
        raise ArithmeticError('complete source pulse leaves the entropy control domain')
    finite = regional.finite_pulse(parameter)
    final_entropy_upper = finite['actual_regional_entropy_change_upper']+entropy_error
    if final_entropy_upper >= F(-4, 10**6):
        raise ArithmeticError('native pulse lost the strict entropy/area witness')
    report.update({'parameter': str(parameter), 'physical_work_bits': 16,
                   'full_work_rows': 1 << 16, 'source_current_terms': {k: str(v) for k, v in terms.items()},
                   'exact_current_synthesis_checked': True,
                   'symmetric_product_error_upper': str(splitting),
                   'complete_actual_ground_output_halftrace_error_upper': str(error),
                   'decoded_entropy_error_upper': str(entropy_error),
                   'decoded_regional_entropy_change_upper': str(final_entropy_upper),
                   'all_failures_and_offcode_decoded_by_fixed_trace_preserving_map': True,
                   'failure_and_offcode_decoder_state': 'I_4/4',
                   'decoded_trace_residual': float(abs(np.trace(decoded)-1)),
                   'decoded_minimum_eigenvalue_diagnostic': float(np.min(spectrum)),
                   'decoded_entropy_decrease_diagnostic': decrease,
                   'selected_ground_preparation_reexecuted_here': True,
                   'ground_preparation_pulses': len(prefix),
                   'actual_input_is_blank': True,
                   'energy_reference_or_selected_ground_swapped': False,
                   'current_control_is_autonomous_H_B_evolution': False})
    return report


@lru_cache(maxsize=1)
def readout():
    """Two complete blank-to-pointer experiments distinguish the current pulse."""
    data, terms = regional.state(), current_terms()
    full_word, _, _, _, preparation = experiment.plan()
    probe = circuits.probe('X', 8, 0, 16)
    prefix = full_word[:-len(probe)]
    end, _, variables = cert.consume(full_word, len(prefix), probe)
    if end != len(full_word) or variables:
        raise ArithmeticError('ground prefix/readout interface changed')
    indices = [sum(1 << (2*j+(0 if mask >> j & 1 else 1)) for j in range(8))
               for mask in range(256)]
    initial = np.zeros((1 << 17, 1), complex)
    initial[0, 0] = 1
    state = np.asarray(data['vector'], complex)/math.sqrt(data['vector_norm_squared'])
    records = []
    for parameter in (F(0), F(1, 10000)):
        ordered = (list(terms.items())+list(reversed(list(terms.items())))) if parameter else []
        template = [event for label, value in ordered
                    for event in cert.pauli(''.join(axis+('X' if axis in 'XY' else 'I')
                                                    for axis in label), -parameter*value/2)]
        vector = state.copy()
        for label, value in ordered:
            angle = -float(parameter*value/2)
            vector = math.cos(angle)*vector+1j*math.sin(angle)*(algebra.tensor(label)@vector)
        target = np.zeros_like(initial)
        target[indices, 0] = (vector+vector[np.arange(256)^1])/2
        target[np.asarray(indices)+(1 << 16), 0] = (vector-vector[np.arange(256)^1])/2
        actual, report = native.complete(template+probe, 17, initial, target,
                                         prefix_word=prefix,
                                         prefix_error=F(preparation['complete_ground_isometry_error']))
        splitting = sum(abs(v) for v in terms.values())**3*parameter**3/3
        error = F(report['complete_isometry_error_bound'])+splitting+data['ground_vector_error']
        probabilities = [float(np.linalg.norm(actual.success[(1 << 16)*j:(1 << 16)*(j+1)])**2)
                         for j in range(2)]
        report.update({'parameter': str(parameter), 'full_work_rows': 1 << 17,
                       'pointer_success_probabilities': probabilities,
                       'reserved_six_level_pointer_probabilities': [0.0]*4,
                       'actual_ground_complete_output_TV_error_upper': str(error),
                       'input_is_blank': True, 'physical_reference_is_same_H_B_ground': True})
        records.append(report)
    parameter = F(1, 10000)
    derivative = regional.pairing(data['tangent'], regional.pauli('XI'))
    error = sum(F(r['actual_ground_complete_output_TV_error_upper']) for r in records)
    lower = (parameter*(derivative-4*data['current_norm_upper']*data['ground_vector_error'])/2
             -parameter**2*data['current_norm_upper']**2-error)
    contrast = records[1]['pointer_success_probabilities'][0]-records[0]['pointer_success_probabilities'][0]
    if not float(lower) > 17e-6 or contrast < float(lower)-2e-11:
        raise ArithmeticError('native complete pointer contrast lost its certified signal')
    return {'experiments': records, 'plus_probability_increase_lower': str(lower),
            'plus_probability_increase_diagnostic': contrast,
            'conditional_success_probabilities_used': False}
