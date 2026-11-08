"""Complete decoded native core, with its exact encoded-source attachment."""
from fractions import Fraction as F
import numpy as np

from model.native_spin_observer import execute
from model.spin_covariant_localization import localization
from . import circuits, geometry


def input_and_target():
    data = localization.source_dictionary()
    modes = np.column_stack([np.asarray(v, float)/np.sqrt(float(n))
                             for v, n in zip(data['columns'], data['norms'], strict=True)])
    probe = modes@np.sqrt(np.asarray(data['weights'], float))
    dark = np.eye(256)[:, 0]-modes@(modes.T@np.eye(256)[:, 0])
    dark /= np.linalg.norm(dark)
    # A coherent source/reference input tests both the bright and dark blocks.
    field = np.column_stack((probe, dark))/np.sqrt(2)
    initial = np.zeros((4096, 2), complex)
    initial[:256] = field
    unitary = geometry.matrix()
    transformed = unitary@field
    target = np.zeros_like(initial)
    for label in (0, *range(8, 16)):
        selected = np.zeros_like(transformed)
        if label == 0:
            selected[8:] = transformed[8:]
        else:
            selected[label-8] = transformed[label-8]
        target[label*256:(label+1)*256] = unitary.T@selected
    return initial, target


def run():
    core, _, report = circuits.plan()
    initial, target = input_and_target()
    result = execute.complete(core, 12, initial)
    distance = float(np.sqrt(np.linalg.norm(result.success-target)**2+result.failure_probability))
    if distance > float(F(report['complete_source_halfdiamond_error']))+execute.ADMISSION:
        raise ArithmeticError('complete native coordinate core violates its certified source bound')
    return {'executed_work_bits': 12, 'all_quantum_output_rows': 4096,
            'retained_reference_columns': 2,
            'native_pulses': len(core), 'word_sha256': result.word_sha256,
            'complete_output_distance_diagnostic': distance,
            'all_first_failure_probability': result.failure_probability,
            'normalization_residual': result.normalization_residual,
            'first_failure_slots': result.native_failure_slots,
            'rounded_negative_count': result.rounded_negative_count,
            'rounded_negative_mass': result.rounded_negative_mass,
            'ideal_decoded_reference_input_bright_weight': '1/2',
            'ideal_decoded_reference_input_dark_weight': '1/2',
            'complete_record_probabilities': [float(np.linalg.norm(
                result.success[label*256:(label+1)*256])**2) for label in range(16)],
            'decoded_core_executed': True,
            'encoded_twenty_bit_wrapper_executed': False,
            'full_encoded_detector_has_analytic_all_input_bound': True}
