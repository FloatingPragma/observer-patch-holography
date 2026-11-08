"""Native charge creation followed by a complete randomized energy readout."""
from fractions import Fraction as F
from functools import lru_cache
import math

import numpy as np

from model.certified_reference import certificate as cert
from model.native_ordered_ground import experiment
from model.native_spin_observer import algebra, circuits
from model.spin_modular_geometry import native, regional
from . import energy


def energy_probe(label):
    if len(label) != 8 or any(axis not in 'IXYZ' for axis in label):
        raise ValueError('complete eight-spin Pauli required')
    basis = []
    for site, axis in enumerate(label):
        if axis == 'X': basis += circuits.spin_rotation('Y', cert._pi(-1, 4), site, 8)
        if axis == 'Y': basis += circuits.spin_rotation('X', cert._pi(1, 4), site, 8)
    parity = [event for site, axis in enumerate(label) if axis != 'I'
              for event in cert.cnot(2*site, 17)]
    return cert.inverse(basis)+parity+basis


def physical_pauli(label, state):
    physical = ''.join(axis+('X' if axis in 'XY' else 'I') for axis in label)
    flip = sum(1 << j for j, axis in enumerate(physical) if axis in 'XY')
    diagonal = sum(1 << j for j, axis in enumerate(physical) if axis in 'YZ')
    rows = np.arange(len(state), dtype=np.int64)
    signs = np.ones(len(state), dtype=np.int8)
    for j in range(16):
        if diagonal >> j & 1: signs *= 1-2*((rows >> j) & 1)
    phase = (1j)**physical.count('Y')
    result = np.empty_like(state)
    result[rows ^ flip] = phase*signs[:, None]*state
    return result


@lru_cache(maxsize=1)
def creation():
    full, _, _, _, preparation = experiment.plan()
    probe = circuits.probe('X', 8, 0, 16)
    prefix = full[:-len(probe)]
    end, _, variables = cert.consume(full, len(prefix), probe)
    if end != len(full) or variables:
        raise ArithmeticError('source preparation/readout boundary changed')
    z = [('I', (0,), cert._pi(-1, 2)), ('Z', (0,), cert._pi(1, 2))]
    data = regional.state()
    v = np.asarray(data['vector'], complex)/math.sqrt(data['vector_norm_squared'])
    v *= np.array([1-2*(mask & 1) for mask in range(256)])
    target = np.zeros((1 << 17, 1), complex)
    indices = [sum(1 << (2*j+(0 if mask >> j & 1 else 1)) for j in range(8))
               for mask in range(256)]
    target[indices, 0] = (v+v[np.arange(256)^1])/2
    target[np.asarray(indices)+(1 << 16), 0] = (v-v[np.arange(256)^1])/2
    blank = np.zeros_like(target)
    blank[0, 0] = 1
    actual, report = native.complete(z+probe, 17, blank, target,
                                    prefix_word=prefix,
                                    prefix_error=F(preparation['complete_ground_isometry_error']))
    report.update({'input_is_blank': True, 'complete_charge_record_labels': [1, -1],
                   'charge_probabilities': [float(np.linalg.norm(actual.success[(1 << 16)*j:(1 << 16)*(j+1)])**2)
                                            for j in range(2)],
                   'reserved_six_level_charge_pointer_probabilities': [0.0]*4,
                   'actual_ground_output_TV_error_upper': str(
                       F(report['complete_isometry_error_bound'])+data['ground_vector_error']),
                   'H_B_energy_reference_changed': False})
    return actual, report


@lru_cache(maxsize=1)
def protocol():
    created, prefix = creation()
    selector_template = [('Y', (bit,), cert._pi(-1, 4)) for bit in range(6)]
    blank = np.zeros((64, 1), complex)
    blank[0, 0] = 1
    selector, selector_report = native.complete(selector_template, 6, blank, np.ones((64, 1))/8)
    selector_probabilities = np.abs(selector.success[:, 0])**2
    offset, coefficients = energy.physical_terms()
    settings = list(coefficients.items())+[('I'*8, F(0))]*(64-len(coefficients))
    if len(settings) != 64 or len(coefficients) != 44:
        raise ArithmeticError('complete energy selector carrier changed')
    ground = energy.work()['complete_ground_certificate']['ground_energy_interval']
    reference = (F(ground[0])+F(ground[1]))/2
    reference_error = (F(ground[1])-F(ground[0]))/2
    initial = np.zeros((1 << 18, 1), complex)
    initial[:1 << 17] = created.success
    table, probabilities, moments = [], [0., 0.], [0., 0.]
    total = created.failure_probability+float(np.vdot(created.success, created.success).real)*selector.failure_probability
    for index, (label, coefficient) in enumerate(settings):
        transformed = physical_pauli(label, created.success)
        target = np.vstack([(created.success+transformed)/2, (created.success-transformed)/2])
        actual, row = native.complete(energy_probe(label), 18, initial, target)
        joint = [[float(np.linalg.norm(actual.success[(r*2+q)*(1 << 16):(r*2+q+1)*(1 << 16)])**2)
                  for r in range(2)] for q in range(2)]
        selector_weight = float(selector_probabilities[index])
        scores = [offset+64*coefficient*sign-reference for sign in (1, -1)]
        for q in range(2):
            probabilities[q] += selector_weight*sum(joint[q])
            moments[q] += selector_weight*sum(float(score)*weight for score, weight in zip(scores, joint[q], strict=True))
        total += selector_weight*(sum(sum(r) for r in joint)+actual.failure_probability)
        row.update({'setting': index, 'pauli': label, 'coefficient': str(coefficient),
                    'selector_success_probability': selector_weight,
                    'charge_energy_success_probabilities': joint,
                    'score_by_energy_pointer': [str(x) for x in scores],
                    'reserved_six_level_energy_pointer_probabilities': [0.0]*4,
                    'input_prefix_success_was_not_renormalized': True})
        table.append(row)
    if abs(total-1) > 2e-11:
        raise ArithmeticError('joint source charge-energy protocol lost total probability')
    error = (F(prefix['actual_ground_output_TV_error_upper'])
             + F(selector_report['complete_isometry_error_bound'])
             + max(F(row['complete_isometry_error_bound']) for row in table))
    # Each failure has score0; every success score has |score|<=33 in J units.
    score_bound = max(abs(F(score)) for row in table for score in row['score_by_energy_pointer'])
    if score_bound >= 33:
        raise ArithmeticError('declared common-clock score range changed')
    energy_error = 2*score_bound*error+reference_error
    by_charge = {row['charge']: row for row in energy.work()['channels']}
    summaries = []
    for index, q in enumerate((1, -1)):
        expected = by_charge[q]
        if not (float(expected['probability_lower']-error)-2e-11 <= probabilities[index]
                <= float(expected['probability_upper']+error)+2e-11):
            raise ArithmeticError('complete charge readout disagrees with its exact sector')
        if not (float(expected['unnormalized_excitation_lower']-energy_error)-2e-11 <= moments[index]
                <= float(expected['unnormalized_excitation_upper']+energy_error)+2e-11):
            raise ArithmeticError('source energy readout disagrees with its commutator work')
        summaries.append({'charge': q, 'unconditioned_success_probability': probabilities[index],
                          'unconditioned_excitation_moment': moments[index],
                          'conditional_mean_excitation_diagnostic': moments[index]/probabilities[index],
                          'conditional_ratio_is_not_a_particle_rest_mass': True})
    return {'creation': prefix, 'selector': selector_report, 'settings': table,
            'ground_energy_reference': str(reference), 'ground_energy_reference_error': str(reference_error),
            'complete_protocol_TV_error_upper': str(error),
            'unconditioned_energy_moment_error_upper': str(energy_error),
            'maximum_absolute_energy_score': str(score_bound),
            'complete_probability_residual': abs(total-1),
            'charge_energy_readouts': summaries,
            'all_failure_records_have_explicit_score': '0',
            'success_postselection_or_renormalization': False,
            'full_work_rows_per_energy_instrument': 1 << 18}


def control_cycle():
    z = [('I', (0,), cert._pi(-1, 2)), ('Z', (0,), cert._pi(1, 2))]
    inclusion = algebra.inclusion(1)/math.sqrt(2)
    _, report = native.complete(z+z, 2, inclusion, inclusion)
    report.update({'ideal_system_energy_change_for_every_input': '0',
                   'positive_primitive_pulse_count': 4,
                   'sum_absolute_ideal_angles': '2*pi',
                   'nonzero_control_effort_determines_nonzero_energy_work': False})
    return report
