"""Complete six-level preparation, bond pulse and local modular readout."""
from collections import Counter
from fractions import Fraction as F
import math

import numpy as np
from scipy.linalg import expm

from model.centered_native_preparation import tree
from model.certified_reference import certificate as cert, interval
from model.native_spin_observer import algebra, circuits, execute


def plan(parameter):
    if type(parameter) is not F or not -1 < parameter < 1:
        raise ValueError('exact source perturbation in (-1,1) required')
    # Exact selected ground in little-endian logical site order.
    logical = np.array([-1, -3, 3, 1], dtype=np.int64)
    inclusion = algebra.inclusion(2)
    values = tuple(int(x.real) for x in inclusion @ logical)
    design = tree.design(values, 4)
    preparation = []
    angular = F(0)
    for bit, raw, coefficients in design:
        a, b = tree.level_error(values, bit, raw, coefficients, 4)
        angular += a+b
        preparation.extend(tree.level_word(bit, coefficients, 4))
    cnots = (1 << 4)-2
    angular += cnots*(2*interval.pi_angle_error(F(math.pi/2), F(1, 2))
                       + interval.pi_angle_error(F(-math.pi), F(-1)))
    time = F(3, 4)
    template = circuits.exchange(time, 1+parameter)+circuits.probe('X', 2, 0, 4)
    measurement = list(circuits.emit(template))
    end, error, variables = cert.consume(measurement, 0, template)
    if end != len(measurement) or variables:
        raise ArithmeticError('pulse/probe angles have unbound input')
    word = preparation+measurement
    port = tree.port_bound(Counter((e['kind'], e['angle']) for e in word))
    bound = angular+error+port
    h = (sum(algebra.tensor(axis*2) for axis in 'XYZ')-np.eye(4))/4
    state = expm(-1j*float(time*(1+parameter))*h) @ (logical/math.sqrt(20))
    projectors = [(np.eye(4)+sign*algebra.tensor('XI'))/2 for sign in (1, -1)]
    target = np.vstack([inclusion @ p @ state[:, None] for p in projectors])
    initial = np.zeros((32, 1), complex)
    initial[0, 0] = 1
    return word, initial, target, {
        'parameter': str(parameter), 'time': str(time),
        'pulses': len(word), 'preparation_pulses': len(preparation),
        'exchange_and_readout_pulses': len(measurement),
        'complete_isometry_error_bound': str(bound),
        'preparation_angle_error': str(angular), 'pulse_readout_angle_error': str(error),
        'actual_port_generator_error': str(port),
        'word_sha256': cert.word_digest(word)}


def run(parameter):
    word, initial, target, result = plan(parameter)
    actual = execute.complete(word, 5, initial)
    distance = math.sqrt(float(np.linalg.norm(actual.success-target)**2)
                         + actual.failure_probability)
    if distance > float(F(result['complete_isometry_error_bound']))+2e-11:
        raise ArithmeticError('complete native readout differs from source target')
    probabilities = [float(np.linalg.norm(actual.success[16*j:16*(j+1)])**2)
                     for j in range(2)]
    ideal = [float(np.linalg.norm(target[16*j:16*(j+1)])**2) for j in range(2)]
    p = (1+0.6*math.cos(float(F(3, 4)*(1+parameter))))/2
    if max(abs(ideal[0]-p), abs(ideal[1]-(1-p))) > 2e-14:
        raise ArithmeticError('native encoding or complete outcome sign changed')
    if abs(sum(probabilities)+actual.failure_probability-1) > 2e-11:
        raise ArithmeticError('complete source output lost normalization')
    result.update({
        'source_pointer_probabilities': probabilities,
        'ideal_pointer_probabilities': ideal,
        'reserved_six_level_pointer_probabilities': [0.0]*4,
        'native_first_failure_weights': actual.failure_weights.tolist(),
        'native_first_failure_slots': actual.native_failure_slots,
        'complete_first_failure_probability': actual.failure_probability,
        'complete_output_distance_diagnostic': distance,
        'normalization_residual': actual.normalization_residual,
        'rounded_negative_count': actual.rounded_negative_count,
        'rounded_negative_mass': actual.rounded_negative_mass,
        'success_postselected_or_renormalized': False,
        'source_preparation_executed_from_blank': True})
    return result


def check():
    return [run(p) for p in (F(-1, 16), F(0), F(1, 16))]
