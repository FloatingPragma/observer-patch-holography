"""Actual source preparation, local parity creation and complete X readout."""
from collections import Counter
from fractions import Fraction as F

import numpy as np

from model.native_ordered_ground import experiment
from model.native_spin_observer import circuits, execute
from model.centered_native_preparation import tree
from model.certified_reference import certificate as cert


def run():
    word, initial, target, checked, preparation = experiment.plan()
    pointer = list(circuits.emit(circuits.probe('X', 8, 0, 16)))
    if word[-len(pointer):] != pointer:
        raise ArithmeticError('the actual prepared ground no longer ends in its X pointer')
    template = [('I', (0,), cert._pi(-1, 2)), ('Z', (0,), cert._pi(1, 2))]
    creation = list(circuits.emit(template))
    _, angle_error, variables = cert.consume(creation, 0, template)
    if variables:
        raise ArithmeticError('creation pulse contains an unbound angle')
    creation_error = angle_error+tree.port_bound(Counter((e['kind'], e['angle']) for e in creation))
    base_error = F(preparation['complete_ground_isometry_error'])
    total_error = base_error+creation_error
    if total_error >= F(1, 10**8):
        raise ArithmeticError('creation and detector exceed the complete error allowance')
    changed = word[:-len(pointer)]+creation+pointer
    base = execute.complete(word, 17, initial)
    created = execute.complete(changed, 17, initial)
    # Z anticommutes with X: the ideal two pointer branches exchange, while
    # the postmeasurement field vector receives Z. This controls back-action.
    expected = target[np.arange(1 << 17) ^ (1 << 16)].copy()
    expected *= np.array([1-2*(i & 1) for i in range(1 << 17)])[:, None]
    for output, wanted, error in ((base, target, base_error),
                                  (created, expected, total_error)):
        distance = np.sqrt(np.linalg.norm(output.success-wanted)**2+output.failure_probability)
        if distance > float(error)+execute.ADMISSION:
            raise ArithmeticError('actual complete source output misses its analytic comparison')
    p0 = float(np.linalg.norm(base.success[:1 << 16])**2)
    p1 = float(np.linalg.norm(created.success[:1 << 16])**2)
    if p0-p1 <= .67 or base.failure_probability > 1e-20 or created.failure_probability > 1e-20:
        raise ArithmeticError('complete executed readout lost its strict contrast or failure control')
    return {'complete_work_bits': 17, 'complete_work_rows': 1 << 17,
            'reference_program_pulses': len(word), 'created_program_pulses': len(changed),
            'creation_pulses': len(creation), 'X_pointer_pulses': len(pointer),
            'reference_word_sha256': base.word_sha256,
            'created_word_sha256': created.word_sha256,
            'reference_first_failure_slots': base.native_failure_slots,
            'created_first_failure_slots': created.native_failure_slots,
            'reference_complete_isometry_error': str(base_error),
            'created_complete_isometry_error': str(total_error),
            'all_failure_labels_and_zero_pointer_slots_retained': True,
            'failure_probability_diagnostic_below': '1/100000000000000000000',
            'both_actual_source_programs_executed': True,
            'success_renormalized': False,
            'analytic_complete_contrast_lower': str(F(67, 100)-base_error-total_error),
            'host_residual_admission': str(F(str(execute.ADMISSION))),
            'host_floating_execution_is_interval_proof': False,
            'ten_repair_lifetime_tape_executed': False}
