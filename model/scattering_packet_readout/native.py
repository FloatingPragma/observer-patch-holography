"""Complete small native control of the reflected collision echo and Y read."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import math

import numpy as np
from scipy.linalg import expm

from ..field_thermal.thermal import execute_with_exports, failure_probability
from ..full_fock_field.field import FieldGeometry, FullFock
from ..observer.source_scalar import compiler
from ..certified_reference import certificate as cert
from ..certified_reference import source as port_source


@lru_cache(maxsize=2)
def hamiltonian(interacting=True):
    geometry = FieldGeometry((1, 1, 2), speed=.125, interaction=1. if interacting else 0.)
    return FullFock(4).hamiltonian(geometry).toarray()


def pauli(label):
    matrices = {'I': np.eye(2), 'X': np.array([[0, 1], [1, 0]]),
                'Y': np.array([[0, -1j], [1j, 0]]), 'Z': np.diag([1, -1])}
    result = np.ones((1, 1), complex)
    for kind in reversed(label):
        result = np.kron(result, matrices[kind])
    return result


@lru_cache(maxsize=2)
def source_terms(interacting=True):
    matrix = hamiltonian(interacting)
    terms, reconstructed = [], np.zeros((16, 16), complex)
    for letters in product('IXYZ', repeat=4):
        label = ''.join(letters)
        basis = pauli(label)
        coefficient = np.trace(basis@matrix)/16
        if coefficient.imag != 0:
            raise ArithmeticError('the source Pauli coefficient is not real')
        if coefficient.real:
            exact = F(float(coefficient.real))
            terms.append((label, exact))
            reconstructed += float(exact)*basis
    if not np.array_equal(reconstructed, matrix):
        raise ArithmeticError('the dyadic full-Fock Pauli expansion changed the source')
    return tuple(terms)


def controlled_evolution(interacting, duration, steps):
    """Target exp(-i duration H), controlled by field-external bit four."""
    if type(duration) is not F or type(steps) is not int or not 1 <= steps <= 64:
        raise ValueError('exact rational duration and steps in 1..64 required')
    word = []
    for _ in range(steps):
        for label, coefficient in source_terms(interacting):
            angle = float(duration*coefficient/(2*steps))
            word.extend(compiler.pauli_rotation(label+'I', -angle))
            word.extend(compiler.pauli_rotation(label+'Z', angle))
    return word


def echo_word(interacting=True, duration=F(1, 4), steps=8):
    free = controlled_evolution(False, duration, steps)
    backward = compiler.inverse(controlled_evolution(interacting, duration, steps))
    # Chronological application gives e^(-i H0 T)e^(+2i H T)e^(-i H0 T).
    return [*free, *backward, *backward, *free]


def exact_echo(interacting=True, duration=F(1, 4)):
    t = float(duration)
    free = expm(-1j*t*hamiltonian(False))
    return free@expm(2j*t*hamiltonian(interacting))@free


def exact_control_template(interacting, duration, steps):
    result = []
    for _ in range(steps):
        for label, coefficient in source_terms(interacting):
            angle = duration*coefficient/(2*steps)
            result += cert.pauli(label+'I', -angle)
            result += cert.pauli(label+'Z', angle)
    return result


def literal_certificate(word, interacting, duration=F(1, 4), steps=8):
    """Independent exact target angles and twelve-port error enclosures."""
    free = exact_control_template(False, duration, steps)
    backward = cert.inverse(exact_control_template(interacting, duration, steps))
    expected = [('X', (0,), cert._pi(1, 2)), ('X', (1,), cert._pi(1, 2)),
                ('H', (4,), cert._pi(1, 2)), ('I', (4,), cert._pi(-1, 2))]
    expected += free+backward+backward+free+[('X', (4,), cert._pi(-1, 4))]
    offset, angle_error, variables = cert.consume(word, 0, expected)
    if offset != len(word) or variables:
        raise ArithmeticError('the native echo has an unaccounted source event')
    table = {}
    port_error = cert.controls_bound(word, table)
    port_source.complete_export_certificate()
    if angle_error+port_error >= F(1, 10**8):
        raise ArithmeticError('literal source complete-isometry error exceeded its bound')
    return {'native_word_sha256': cert.word_digest(word),
            'native_pulses': len(word), 'unique_literal_controls': len(table),
            'all_exact_target_angles_and_routing_checked': True,
            'literal_control_complete_isometry_error_upper': F(1, 10**8),
            'proper_code_exports_complete': True,
            'physical_hardware_accuracy_asserted': False}


def execute(interacting=True, duration=F(1, 4), steps=8, initial=None):
    """Every one of 32 field/pointer rows and every native export are kept."""
    if type(interacting) is not bool or type(duration) is not F or not 0 < duration <= 1:
        raise ValueError('boolean control and exact duration in (0,1] required')
    preparation = []
    if initial is None:
        field = np.zeros((16, 1), complex)
        field[0, 0] = 1
        preparation = [compiler.pulse('X', math.pi/2, [0]),
                       compiler.pulse('X', math.pi/2, [1])]
        prepared = np.zeros_like(field)
        prepared[3, 0] = -1
    else:
        field = np.asarray(initial, complex)
        if field.ndim != 2 or field.shape[0] != 16 or not np.isfinite(field).all():
            raise ValueError('finite full sixteen-row reference columns required')
        if abs(float(np.vdot(field, field).real)-1) > 1e-12:
            raise ValueError('unit total norm required')
        prepared = field
    pointer = [compiler.pulse('H', math.pi/2, [4]),
               compiler.pulse('I', -math.pi/2, [4])]
    readout = [compiler.pulse('X', -math.pi/4, [4])]
    word = [*preparation, *pointer, *echo_word(interacting, duration, steps), *readout]
    input_state = np.zeros((32, field.shape[1]), complex)
    input_state[:16] = field
    success, failures = execute_with_exports(word, 5, input_state)
    echo = exact_echo(interacting, duration)
    after = echo@prepared
    target = np.vstack(((prepared-1j*after)/2, (after-1j*prepared)/2))
    failure_mass = failure_probability(failures)
    error = math.sqrt(float(np.linalg.norm(success-target)**2)+failure_mass)
    gram = success.conj().T@success
    for events in failures:
        for event in events:
            gram += event.conj().T@event
    gram_error = float(np.linalg.norm(gram-field.conj().T@field))
    probability = [float(np.linalg.norm(success[:16])**2), float(np.linalg.norm(success[16:])**2), 0., 0., 0., 0.]
    reference_probability = (1+float(np.trace(prepared.conj().T@echo@prepared).imag))/2
    term_norm = sum(abs(x) for _, x in source_terms(interacting))
    product_error = 6*duration**2*term_norm**2/steps
    if (gram_error > 1e-8 or error > float(product_error)+1e-8 or
            abs(sum(probability)+failure_mass-1) > 1e-8):
        raise ArithmeticError('complete native Ramsey source control failed')
    return {'success': success, 'first_failure_events': failures,
            'target': target, 'word': word,
            'report': {'interacting': interacting, 'duration': str(duration), 'steps': steps,
                       'field_rows': 16, 'joint_quantum_rows': 32,
                       'reference_columns': field.shape[1], 'pointer_probabilities': probability,
                       'exact_Y_plus_probability': reference_probability,
                       'failure_probability': failure_mass,
                       'complete_vector_error': error, 'reference_gram_error': gram_error,
                       'product_operator_error_upper': str(product_error),
                       'source_resources': compiler.resources(word),
                       'complete_outputs_retained': True,
                       'analytic_scattering_phase_gate_used': False,
                       'fine_spacing_packet_experiment': False}}
