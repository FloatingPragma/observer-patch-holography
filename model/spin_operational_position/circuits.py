"""Literal admitted native words for the complete finite coordinate detector."""
from collections import Counter
from copy import deepcopy
from fractions import Fraction as F
from functools import lru_cache
import math

from model.centered_native_preparation import tree
from model.certified_reference import certificate as cert, interval as iv
from model.native_ordered_ground import source
from model.native_spin_observer import circuits as spin
from model.native_spectral_detector.circuits import hadamard_template
from model.observer.source_scalar import compiler
from . import geometry


def diagonal_phase_zero(bits):
    """Exactly I-2|0><0|, with every parity factor and final Gray CNOT."""
    if type(bits) is not int or bits < 1:
        raise ValueError('positive phase register size required')
    angle = cert._pi(1, 1 << bits)
    word = [('I', (0,), angle)]
    for pivot in range(bits):
        count = 1 << pivot
        for j in range(count):
            word.append(('Z', (pivot,), angle))
            if pivot:
                first = j ^ (j >> 1)
                following = ((j+1) % count) ^ (((j+1) % count) >> 1)
                toggle = first ^ following
                word += cert.cnot(toggle.bit_length()-1, pivot)
    return word


def multi_controlled_not(controls, target, pattern):
    """Flip the target exactly on the specified control word, using its phase."""
    if (type(controls) is not tuple or not controls or target in controls
            or len(set(controls)) != len(controls)
            or any(type(i) is not int or i < 0 for i in (*controls, target))
            or type(pattern) is not int or not 0 <= pattern < (1 << len(controls))):
        raise ValueError('distinct nonnegative controls/target and an exact control pattern required')
    # H (I-2|controls=pattern,target=1><...|) H is controlled X.
    before = sum((spin.xgate(bit) for j, bit in enumerate(controls) if pattern >> j & 1), [])
    before += spin.xgate(target)
    phase = spin.remap(diagonal_phase_zero(len(controls)+1), (*controls, target))
    basis = hadamard_template(target)
    return basis+before+phase+cert.inverse(before)+basis


def coordinate_label():
    """Twelve decoded bits: eight field, three coordinate, one validity flag."""
    # Valid coordinate anchors have high five field bits zero. All remaining
    # 248 field states get exactly the single dark label (flag=0, xyz=000).
    word = multi_controlled_not(tuple(range(3, 8)), 11, 0)
    for bit in range(3):
        word += multi_controlled_not((11, bit), 8+bit, 3)
    return word


def label_of(mask):
    if type(mask) is not int or not 0 <= mask < 256:
        raise ValueError('an eight-bit decoded field value required')
    return 8+mask if mask < 8 else 0


def inverse_word(word):
    return [compiler.pulse(event['kind'], -event['angle'], list(event['targets']))
            for event in reversed(word)]


def _template_word(template):
    word = list(spin.emit(template))
    end, error, variables = cert.consume(word, 0, template)
    if end != len(word) or variables:
        raise ArithmeticError('unaccounted fixed template angle')
    return word, error


def _tree_word(values):
    word = []
    error = F()
    for bit, raw, coefficients in tree.design(values, 8):
        a, b = tree.level_error(values, bit, raw, coefficients, 8)
        error += a+b
        word += list(tree.level_word(bit, coefficients, 8))
    fixed = 254*(2*iv.pi_angle_error(F(math.pi/2), F(1, 2))
                 +iv.pi_angle_error(F(-math.pi), F(-1)))
    return word, error+fixed


@lru_cache(maxsize=1)
def _plan():
    checked = geometry.certify()
    phase_word, phase_error = _template_word(diagonal_phase_zero(8))
    forward = []
    angle_error = F()
    tree_rows = []
    for normal in geometry.witness():
        preparation, error = _tree_word(normal)
        forward += inverse_word(preparation)+phase_word+preparation
        angle_error += 2*error+phase_error
        tree_rows.append({'pulses': len(preparation), 'tree_operator_error': str(error)})
    walsh, walsh_error = _template_word(sum((hadamard_template(i) for i in range(3)), []))
    forward += walsh
    angle_error += walsh_error
    label_word, label_error = _template_word(coordinate_label())
    core = forward+label_word+inverse_word(forward)
    angle_error = 2*angle_error+label_error
    # Actual spin encoding and Marshall frame, rather than a replacement
    # of the prepared ground by its positive-coordinate vector.
    marshall = sum(([('I', (2*x,), cert._pi(-1, 2)),
                      ('Z', (2*x,), cert._pi(1, 2))] for x in source.cube()[1]), [])
    before_template = cert.inverse(spin.encode(8))+marshall
    before, wrapper_error = _template_word(before_template)
    mapping = tuple(2*x for x in range(8))+tuple(range(16, 20))
    physical = before+[dict(event, targets=[mapping[j] for j in event['targets']])
                       for event in core]+inverse_word(before)
    angle_error += 2*wrapper_error
    port_error = tree.port_bound(Counter((event['kind'], event['angle']) for event in physical))
    source_error = angle_error+port_error
    detector_error = (source_error+F(checked['complete_detector_error_from_reflectors'])
                      +F(checked['complete_detector_error_from_true_ground']))
    if detector_error >= F(1, 10**6):
        raise ArithmeticError('complete source detector misses its declared tolerance')
    report = {'spin_sites': 8, 'full_spin_dimension': 256,
              'physical_work_bits': 20, 'decoded_core_bits': 12,
              'coordinate_labels': list(range(8, 16)), 'dark_label': 0,
              'unused_binary_labels_retained': list(range(1, 8)),
              'reserved_nonbinary_pointer_slots_retained': True,
              'physical_pulses': len(physical), 'decoded_core_pulses': len(core),
              'physical_word_sha256': cert.word_digest(physical),
              'decoded_core_word_sha256': cert.word_digest(core),
              'amplitude_trees': tree_rows, 'literal_angle_operator_error': str(angle_error),
              'literal_port_operator_error': str(port_error),
              'complete_source_halfdiamond_error': str(source_error),
              'complete_actual_coordinate_halfdiamond_error': str(detector_error),
              'full_twenty_bit_native_tape_executed': False,
              'declared_source_control_and_blank_pointer_services_reused': True,
              'probability_or_failure_renormalization': False,
              'dark_subspace_luders_backaction_retained': True}
    return tuple(core), tuple(physical), report


def plan():
    return deepcopy(_plan())
