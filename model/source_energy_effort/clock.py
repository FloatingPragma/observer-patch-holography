"""Source-clock action readout, including its physical relative phase."""
from fractions import Fraction as F
from functools import lru_cache
import math

import numpy as np

from model.native_joint_spin_spectroscopy.circuits import controlled_pauli
from model.native_ordered_ground import experiment
from model.native_spectral_detector.circuits import hadamard_template
from model.native_spin_observer import algebra, circuits
from model.certified_reference import certificate as cert
from model.spin_modular_geometry import native, regional
from . import energy


@lru_cache(maxsize=1)
def run():
    time = F(1, 4096)
    constant, coefficients = energy.physical_terms()
    bounds = energy.work()['complete_ground_certificate']['ground_energy_interval']
    reference, reference_error = (F(bounds[0])+F(bounds[1]))/2, (F(bounds[1])-F(bounds[0]))/2
    terms = [('I'*8, constant-reference)]+list(coefficients.items())
    ordered = terms+list(reversed(terms))
    full, _, _, _, preparation = experiment.plan()
    probe = circuits.probe('X', 8, 0, 16)
    prefix = full[:-len(probe)]
    end, _, variables = cert.consume(full, len(prefix), probe)
    if end != len(full) or variables:
        raise ArithmeticError('same-source preparation prefix changed')
    z = [('I', (0,), cert._pi(-1, 2)), ('Z', (0,), cert._pi(1, 2))]
    template = z+probe+hadamard_template(17)
    for label, coefficient in ordered:
        physical = ''.join(axis+('X' if axis in 'XY' else 'I') for axis in label)+'II'
        template += controlled_pauli(physical, -time*coefficient/2, 17)
    template += [('X', (17,), cert._pi(1, 4))]  # Z readout of -Y: positive sine response
    data = regional.state()
    v = np.asarray(data['vector'], complex)/math.sqrt(data['vector_norm_squared'])
    v *= np.array([1-2*(mask & 1) for mask in range(256)])
    channels = [(v+sign*v[np.arange(256)^1])/2 for sign in (1, -1)]
    evolved = np.column_stack(channels)
    for label, coefficient in ordered:
        angle = -float(time*coefficient/2)
        evolved = (math.cos(angle)*evolved
                   + 1j*math.sin(angle)*(algebra.tensor(label)@evolved))
    indices = np.array([sum(1 << (2*j+(0 if mask >> j & 1 else 1)) for j in range(8))
                        for mask in range(256)])
    target = np.zeros((1 << 18, 1), complex)
    for q in range(2):
        target[indices+(q << 16), 0] = (channels[q]+1j*evolved[:, q])/2
        target[indices+(q << 16)+(1 << 17), 0] = (1j*channels[q]+evolved[:, q])/2
    initial = np.zeros_like(target)
    initial[0, 0] = 1
    actual, report = native.complete(template, 18, initial, target,
                                     prefix_word=prefix,
                                     prefix_error=F(preparation['complete_ground_isometry_error']))
    norm_sum = sum(abs(coefficient) for _, coefficient in terms)
    splitting = norm_sum**3*time**3/3
    output_error = (F(report['complete_isometry_error_bound'])+data['ground_vector_error']
                    +splitting+time*reference_error)
    derivative_error = 2*output_error/time+time**2*F(28)**3/6
    joint = [[float(np.linalg.norm(actual.success[(r*2+q)*(1 << 16):(r*2+q+1)*(1 << 16)])**2)
              for r in range(2)] for q in range(2)]
    derivative = [(row[0]-row[1])/float(time) for row in joint]
    source = {row['charge']: row for row in energy.work()['channels']}
    for q, value in zip((1, -1), derivative, strict=True):
        row = source[q]
        if not (float(row['unnormalized_excitation_lower']-derivative_error)-2e-7 <= value
                <= float(row['unnormalized_excitation_upper']+derivative_error)+2e-7):
            raise ArithmeticError('controlled clock action rate disagrees with the work readout')
    report.update({'time': str(time), 'common_ground_energy_reference': str(reference),
                   'ground_reference_error': str(reference_error),
                   'controlled_scalar_phase_retained': True,
                   'symmetric_clock_product_error_upper': str(splitting),
                   'complete_true_clock_output_TV_error_upper': str(output_error),
                   'energy_derivative_error_upper': str(derivative_error),
                   'joint_charge_clock_success_probabilities': joint,
                   'unnormalized_mean_energy_from_sine_derivative': derivative,
                   'full_work_rows': 1 << 18, 'all_failure_clock_scores': '0',
                   'input_is_blank': True, 'native_physical_clock_reference_changed': False,
                   'native_controlled_clock_is_unread_C_E_dynamics': False})
    return report


def maintenance(time=F(7, 3)):
    if type(time) is not F or time < 0:
        raise ValueError('exact nonnegative common-clock duration required')
    rows = []
    for row in energy.work()['channels']:
        rows.append({'charge': row['charge'], 'duration': time,
                     'conditional_relative_action_lower': time*row['mean_excitation_lower'],
                     'conditional_relative_action_upper': time*row['mean_excitation_upper'],
                     'autonomous_system_power_exact': F(0),
                     'bare_protected_repair_energy_input_exact': F(0)})
    return {'source_energy_action_relation': 'relative_action(T)=T*mean_excitation_energy',
            'created_channel_action': rows,
            'arbitrary_C_E_repairs_between_exact_clock_intervals_change_the_state': False,
            'controller_reset_heat_or_apparatus_energy_derived': False,
            'nonzero_energy_requires_nonzero_continuing_maintenance_power': False}
