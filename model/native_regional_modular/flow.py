"""Compile and certify the full corrected regional modular generator."""
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import math

import numpy as np
from scipy.linalg import expm
import sympy as sp

from ..certified_reference import certificate as literal, source as native
from ..field_thermal.thermal import car_restriction, restrict
from ..observer.source_scalar import compiler
from ..plaquette_reference import preparation
from ..prepared_regional_entropy import entropy
from ..reference_recovery.recovery import execute_complete
from ..reference_observation.observation import frozen
from ..regional_modular_expansion import expansion


PARAMETER = F(1, 2)
STEPS = 32


def parameters(parameter, steps):
    if type(parameter) not in (int, F) or abs(parameter) > 1:
        raise ValueError('an exact dimensionless modular parameter in [-1,1] is required')
    if type(steps) is not int or not 1 <= steps <= 64:
        raise ValueError('one through 64 complete symmetric steps are required')
    return F(parameter), steps


@lru_cache(maxsize=1)
def generator():
    """Exact full P8, with a checked high-bit to little-endian convention."""
    polynomial = expansion.polynomial(expansion.plaquette())
    original = expansion.pauli_terms(polynomial)
    terms = tuple((label[::-1], F(value)) for label, value in sorted(original.items()))
    sigma = {'I': sp.eye(2), 'X': sp.Matrix([[0, 1], [1, 0]]),
             'Y': sp.Matrix([[0, -sp.I], [sp.I, 0]]), 'Z': sp.diag(1, -1)}
    reconstruction = sp.zeros(16)
    for label, coefficient in terms:
        reconstruction += sp.Rational(coefficient)*sp.kronecker_product(
            *(sigma[letter] for letter in reversed(label)))
    if reconstruction != polynomial or polynomial != polynomial.conjugate().T:
        raise ArithmeticError('complete Hermitian regional generator was not reconstructed')
    parity = sp.diag(*[(-1)**mask.bit_count() for mask in range(16)])
    number = sp.diag(*[mask.bit_count() for mask in range(16)])
    if polynomial*parity != parity*polynomial or polynomial*number != number*polynomial:
        raise ArithmeticError('the regional source generator lost its parity/number symmetry')
    if any(sum(letter in 'XY' for letter in label) % 2 for label, _ in terms):
        raise ArithmeticError('a split term is not an even regional operator')
    return polynomial, terms


def word(parameter=PARAMETER, steps=STEPS):
    parameter, steps = parameters(parameter, steps)
    _, terms = generator()
    half = parameter/(2*steps)
    sweep = list(terms)+list(reversed(terms))
    return [event for _ in range(steps) for label, coefficient in sweep
            for event in compiler.pauli_rotation(label, float(half*coefficient))]


def template(parameter=PARAMETER, steps=STEPS):
    parameter, steps = parameters(parameter, steps)
    _, terms = generator()
    half = parameter/(2*steps)
    return [event for _ in range(steps) for label, coefficient in (*terms, *reversed(terms))
            for event in literal.pauli(label, half*coefficient)]


def exact_bounds(prior, parameter=PARAMETER, steps=STEPS, emitted=None):
    """Certify the actual emitted source tape, including each abort branch."""
    parameter, steps = parameters(parameter, steps)
    emitted = word(parameter, steps) if emitted is None else emitted
    source = entropy.source_bounds(prior)
    native.specification()
    exports = native.complete_export_certificate()
    offset, angle, variables = literal.consume(emitted, 0, template(parameter, steps))
    if variables or offset != len(emitted):
        raise ArithmeticError('regional source word has an unaccounted pulse')
    table = {}
    controls = literal.controls_bound(emitted, table)
    _, terms = generator()
    norm_sum = sum((abs(value) for _, value in terms), F())
    # Equal derivatives through order two and a unitary third-derivative bound.
    splitting = norm_sum**3*abs(parameter)**3/(3*steps**2)
    modular = abs(parameter)*source['modular_polynomial_error']
    literal_error = angle+controls
    if literal_error >= F(1, 10**8):
        raise ArithmeticError('literal regional controls exceed the declared source-error budget')
    total = literal_error+splitting+modular
    counts = Counter(event['kind'] for event in emitted)
    paired = counts['J']
    failure_slots = 4*(len(emitted)-paired)+2*paired
    digest = hashlib.sha256()
    for kind, represented_angle in sorted(table):
        values = native.literal_controls(kind, represented_angle)
        digest.update((kind+'|'+represented_angle.hex()+'|'+','.join(map(str, values))+';').encode())
    return {'modular_parameter': parameter, 'symmetric_steps': steps,
            'exact_Pauli_terms': {label: value for label, value in terms},
            'Pauli_labels_are_little_endian': True,
            'generator_norm_sum_upper': norm_sum,
            'native_word_sha256': literal.word_digest(emitted),
            'integrated_controls_sha256': digest.hexdigest(),
            'fixed_and_variable_angle_error': angle,
            'literal_port_generator_error': controls,
            'literal_complete_output_error': literal_error,
            'ideal_symmetric_product_error': splitting,
            'actual_prepared_modular_generator_bridge': modular,
            'complete_half_diamond_error_upper': total,
            'actuator_abort_probability_upper': literal_error**2,
            'source_certificate': entropy.PRIOR,
            'inherited_source_hashes_checked': True,
            'inherited_native_control_proof_reexecuted': False,
            'resources': {'regional_work_bits': 4, 'exterior_field_bits_retained': 4,
                          'classical_preparation_record_bits_retained': 8,
                          'native_pulses': len(emitted), 'pair_code_pulses': paired,
                          'proper_code_load_events': len(emitted),
                          'proper_code_export_events': len(emitted),
                          'first_failure_Kraus_slots': failure_slots,
                          'all_zero_angle_requests_retained': True,
                          'extra_initialization_or_source_repreparation': False},
            'complete_export_certificate': exports}


def execute(parameter=PARAMETER, steps=STEPS):
    parameter, steps = parameters(parameter, steps)
    return _execute(parameter, steps)


@lru_cache(maxsize=4)
def _execute(parameter, steps):
    return execute_complete(word(parameter, steps), 4, np.eye(16, dtype=complex))


def diagnostic(run, parameter=PARAMETER, steps=STEPS):
    parameter, steps = parameters(parameter, steps)
    polynomial, _ = generator()
    target = expm(1j*float(parameter)*np.asarray(polynomial, complex))
    success_error = float(np.linalg.norm(run.success-target, 2))
    failure = max(0., float(np.linalg.norm(run.failure_gram, 2)))
    gram_error = float(np.linalg.norm(run.success.conj().T@run.success+
                                      run.failure_gram-np.eye(16), 2))
    return {'success_operator_error': success_error,
            'complete_operator_error_diagnostic': math.sqrt(success_error**2+failure),
            'complete_Gram_error': gram_error,
            'failure_effect_norm': failure,
            'all_input_columns_executed': run.success.shape[1],
            'raw_arrays_polarized_or_normalized': False}


@dataclass(frozen=True)
class Continued:
    """Continue every retained preparation column with the same local channel."""
    prepared: object
    regional_run: object
    columns: np.ndarray
    actuator_failure_probability: float

    def failure(self, event):
        """Every actuator complement map on all 256 records and exterior bits."""
        return tuple(np.kron(np.eye(16), branch)@self.prepared.columns
                     for branch in self.regional_run.failure(event))

    @property
    def failure_probability(self):
        return self.prepared.failure_probability+self.actuator_failure_probability

    @property
    def field_density(self):
        return self.columns@self.columns.conj().T


def continue_prepared(prepared, run=None):
    if not isinstance(prepared, preparation.Prepared):
        raise ValueError('the complete retained plaquette preparation is required')
    run = execute() if run is None else run
    if run.bits != 4 or run.success.shape != (16, 16):
        raise ValueError('the complete regional input action is required')
    columns = np.kron(np.eye(16), run.success)@prepared.columns
    regional = restrict(prepared.field_density, car_restriction(8, (0, 1, 2, 3)))
    failure = float(np.trace(regional@run.failure_gram).real)
    if failure < -1e-15:
        raise ArithmeticError('negative actuator failure probability')
    return Continued(prepared, run, frozen(columns), max(0., failure))


def prepared_diagnostic(prepared, continued):
    maps = car_restriction(8, (0, 1, 2, 3))
    before = restrict(prepared.field_density, maps)
    after = restrict(continued.field_density, maps)
    values, vectors = np.linalg.eigh(before)
    if min(values) <= 0:
        raise ArithmeticError('the source regional success reference is not faithful')
    actual = (vectors*np.exp(-1j*float(PARAMETER)*np.log(values)))@vectors.conj().T
    phase = np.vdot(actual, continued.regional_run.success)
    aligned = continued.regional_run.success*np.exp(-1j*np.angle(phase))
    gram = continued.columns.conj().T@continued.columns
    return {'actual_reference_modular_action_error_up_to_scalar': float(np.linalg.norm(aligned-actual, 2)),
            'regional_reference_stationarity_error': float(np.linalg.norm(after-before, ord='nuc')),
            'all256_classical_records_retained': continued.columns.shape[1] == 256,
            'complete_trace_error': abs(float(np.trace(gram).real)+continued.failure_probability-1),
            'preparation_abort_probability': prepared.failure_probability,
            'actuator_abort_probability': continued.actuator_failure_probability,
            'success_field_trace': float(np.trace(gram).real),
            'success_renormalized': False}
