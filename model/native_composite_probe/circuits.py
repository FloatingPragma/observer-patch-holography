"""Factorized coherent selection of all 24 signed cubic Pauli words."""
from fractions import Fraction as F
from functools import lru_cache
import math

from ..certified_reference import certificate as cert
from ..observer.source_scalar import compiler
from ..centered_native_preparation import tree

FIELD_MODES = (0, 1, 2, 3, 4, 5, 8, 9)
LOCAL_BITS = 8
BITS = 13
PAULI_SELECTORS = (8, 9, 10)
NEIGHBOR_SELECTORS = (11, 12)


def remap(template, targets):
    return [(kind, tuple(targets[j] for j in bits), angle) for kind, bits, angle in template]


def hadamard(bit):
    return [('H', (bit,), cert._pi(1, 2)), ('I', (bit,), cert._pi(-1, 2))]


def phase_projector(conditions, angle):
    """exp(i angle product_j |b_j><b_j|), including its scalar phase."""
    if (type(conditions) is not tuple or not conditions or len(conditions) > 4
            or len({bit for bit, _ in conditions}) != len(conditions)
            or any(type(bit) is not int or not 0 <= bit < BITS
                   or type(value) is not int or value not in (0, 1)
                   for bit, value in conditions)
            or type(angle) is not cert.PiMultiple):
        raise ValueError('one through four distinct literal projector bits required')
    result = []
    for mask in range(1 << len(conditions)):
        label = ''.join('Z' if mask >> j & 1 else 'I' for j in range(len(conditions)))
        sign = (-1)**sum(value for j, (_, value) in enumerate(conditions) if mask >> j & 1)
        result += remap(cert.pauli(label, cert.PiMultiple(angle.value*sign/(1 << len(conditions)))),
                        tuple(bit for bit, _ in conditions))
    return result


def flip(bit, controls=()):
    return hadamard(bit)+phase_projector(controls+((bit, 1),), cert._pi(1))+hadamard(bit)


def preparation_template():
    # The three Pauli bits are uniform. The neighbor bits have weights
    # (1,1,1,0)/3; the single variable angle has an independent certificate.
    result = sum((hadamard(bit) for bit in PAULI_SELECTORS), [])
    low, high = NEIGHBOR_SELECTORS
    result += [('Y', (high,), None), ('Y', (low,), cert._pi(-1, 8))]
    result += cert.cnot(high, low)
    result += [('Y', (low,), cert._pi(-1, 8))]+cert.cnot(high, low)
    return result


def selection_template():
    first, up_choice, down_choice = PAULI_SELECTORS
    low, high = NEIGHBOR_SELECTORS
    result = phase_projector(((0, 1),), cert._pi(1))
    # X followed by Z is iY. Reversing the order gives -iY.
    result += flip(1)+phase_projector(((first, 1), (1, 1)), cert._pi(1))
    for neighbor, up in enumerate((2, 4, 6)):
        controls = ((low, neighbor & 1), (high, (neighbor >> 1) & 1))
        result += phase_projector(controls+((up_choice, 1), (up, 1)), cert._pi(1))
        result += flip(up, controls)
        result += flip(up+1, controls)
        result += phase_projector(controls+((down_choice, 1), (up+1, 1)), cert._pi(1))
    return result


def literal(template, *, inverse=False):
    variable = float.fromhex('-0x1.3b2028082e8d4p-1')
    events = []
    for kind, targets, angle in template:
        number = (variable if angle is None else
                  float(angle.value)*math.pi if type(angle) is cert.PiMultiple else float(angle))
        events.append(compiler.pulse(kind, number, list(targets)))
    if inverse:
        events = [compiler.pulse(e['kind'], -e['angle'], e['targets']) for e in reversed(events)]
    return tuple(events)


@lru_cache(maxsize=1)
def word():
    prepare = literal(preparation_template())
    unprepare = tuple(compiler.pulse(e['kind'], -e['angle'], e['targets']) for e in reversed(prepare))
    return prepare+literal(selection_template())+unprepare


@lru_cache(maxsize=1)
def blocks():
    return tuple((instrument, targets) for _, instrument, targets in tree.blocks(word()))


@lru_cache(maxsize=1)
def certificate():
    prepare = literal(preparation_template())
    offset, fixed, variables = cert.consume(prepare, 0, preparation_template())
    if offset != len(prepare) or len(variables) != 1:
        raise ArithmeticError('unaccounted coherent selector preparation')
    variable_error = cert.two_vector_error(variables[0], 2, 1, False)
    selected = literal(selection_template())
    offset, select_error, unexpected = cert.consume(selected, 0, selection_template())
    if offset != len(selected) or unexpected:
        raise ArithmeticError('unaccounted composite SELECT gate')
    ports = cert.controls_bound(word(), {})
    error = 2*(fixed+variable_error)+select_error+ports
    if error >= F(1, 10**8):
        raise ArithmeticError('complete composite native control bound exceeded')
    return {'native_pulses': len(word()), 'native_blocks': len(blocks()),
            'native_word_sha256': cert.word_digest(word()),
            'complete_operator_and_untouched_reference_error': str(error),
            'literal_port_generator_error': str(ports),
            'fixed_angle_error': str(2*fixed+select_error),
            'coherent_neighbor_preparation_and_inverse_error': str(2*variable_error),
            'source_processor_dimension': 6, 'logical_pointer_labels': 32,
            'declared_native_pointer_outcomes': 7776,
            'all_first_export_failure_outcomes_retained': True,
            'ideal_success_Kraus_operator': 'C_0/3',
            'Pauli_terms': 24, 'sum_absolute_Pauli_coefficients': '3'}
