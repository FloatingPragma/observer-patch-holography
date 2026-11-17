"""One local source word implements the existing native balanced mixer.

The supplied coherent proper-code services are unchanged. The balanced angle
is an operational design, not a consequence of A1--A3. Original binary arrays
are compared on all eight seam coordinates; their full spectator lift keeps
every one of the 512 native/padding addresses. No normalization or projection
is applied to the source execution.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
import hashlib
import math
from pathlib import Path

import numpy as np

from ..common_source import exact_output as exact, output_certificate as report
from ..native_interface import construction as interface
from ..native_quantum_history import construction as quantum
from ..observer.source_discrimination import circuit
from ..observer.source_scalar import compiler


ROOT = Path(__file__).resolve().parents[2]
KINDS = ('mixer', 'phase', 'y_analyzer')
SOURCE_PATHS = tuple(sorted(set((*report.SOURCE_PATHS, *interface.SOURCE_PATHS,
    'model/compiled_native_interface/__init__.py',
    'model/compiled_native_interface/__main__.py',
    'model/compiled_native_interface/construction.py',
    'docs/research/COHERENT_SOURCE_CLOCKS.md'))))
ERROR_GATE = F(1, 10**11)


def _kind(value):
    if type(value) is not str or value not in KINDS:
        raise ValueError('a declared mixer, phase, or y_analyzer source word is required')
    return value


def _reflection_word():
    """exp(-i*pi*|3><3|), retaining scalar phase and all Walsh factors.

    Each cyclic Gray walk accumulates Z parity on its greatest support bit
    and returns to its initial tensor chart. There are eight phase pulses
    and six phase-exact source CNOTs: 26 primitive pulses in all.
    """
    tape = [compiler.pulse('I', -math.pi/8, [0])]
    for pivot in range(3):
        count = 1 << pivot
        for index in range(count):
            gray = index ^ (index >> 1)
            mask = (1 << pivot) | gray
            sign = (-1)**((3 & mask).bit_count())
            tape.append(compiler.pulse('Z', -sign*math.pi/8, [pivot]))
            following = (index+1) % count
            difference = gray ^ (following ^ (following >> 1))
            if difference:
                tape.extend(compiler.cnot(difference.bit_length()-1, pivot))
    return tape


def word(kind='mixer'):
    """Fresh chronological source pulses; spectator qubits are untouched."""
    kind = _kind(kind)
    phase = [compiler.pulse('Z', math.pi/4, [0])]
    if kind == 'phase':
        return phase
    # CNOT maps seam addresses (2,1) to (2,3). Conjugating the reflection
    # at 3 by +/-pi/8 target-Y rotations gives H on that ordered pair.
    balanced = (compiler.cnot(0, 1)
                + [compiler.pulse('Y', math.pi/8, [0])]
                + _reflection_word()
                + [compiler.pulse('Y', -math.pi/8, [0])]
                + compiler.cnot(0, 1))
    return phase+balanced if kind == 'y_analyzer' else balanced


@lru_cache(maxsize=3)
def _executed(kind):
    matrix = compiler.execute(word(kind), 3)
    return matrix.dtype.str, matrix.tobytes()


def local_matrix(kind='mixer'):
    """Detach the raw source-executed 8x8 map, before any correction."""
    dtype, payload = _executed(_kind(kind))
    return np.frombuffer(payload, dtype=np.dtype(dtype)).reshape(8, 8)


def full_matrix(kind='mixer'):
    """Lift the executed map to all 512 coordinates, without compression.

    The six spectator bits are the high bits of circuit.encode. The tensor
    identity is fixed wire addressing, not an additional active operation.
    The target's declared extension also mixes padding pairs (8*s+2,8*s+1).
    It is not the identity on all padding.
    """
    matrix = np.kron(np.eye(64, dtype=int), local_matrix(kind))
    return np.frombuffer(matrix.tobytes(), dtype=matrix.dtype).reshape(512, 512)


def embedding():
    """The existing binary-shell embedding, with no new basis convention."""
    codes = tuple(circuit.encode(state) for state in quantum.STATES)
    if len(set(codes)) != 66:
        raise ArithmeticError('the binary source embedding lost injectivity')
    for left, right in interface.phase_pairs():
        a, b = codes[left], codes[right]
        if a % 8 != 2 or b % 8 != 1 or a//8 != b//8:
            raise ArithmeticError('native phase pairs no longer share one spectator address')
    result = np.zeros((512, 66), dtype=int)
    result[codes, np.arange(66)] = 1
    return result


def ideal_native_matrix():
    """Ideal 66D surrogate; it is not returned as source-executed data."""
    embedding()
    return interface.mixer()


# An exact component is a+b*sqrt(2); a matrix entry has real and imaginary
# pairs ((a,b),(c,d)). Only these small declared algebraic targets use it.
ZERO = (F(0), F(0))
ONE = (F(1), F(0))


def _add(a, b):
    return a[0]+b[0], a[1]+b[1]


def _multiply(a, b):
    return a[0]*b[0]+2*a[1]*b[1], a[0]*b[1]+a[1]*b[0]


def _complex_product(a, b):
    rr, ii = _multiply(a[0], b[0]), _multiply(a[1], b[1])
    return (rr[0]-ii[0], rr[1]-ii[1]), _add(_multiply(a[0], b[1]), _multiply(a[1], b[0]))


def _conjugate(a):
    return a[0], (-a[1][0], -a[1][1])


def target(kind='mixer'):
    """Independent complete exact 8D target in Q(sqrt(2))+i Q(sqrt(2))."""
    kind = _kind(kind)
    balanced = {(i, i): (ONE, ZERO) for i in range(8) if i not in (1, 2)}
    for row, column, sign in ((2, 2, 1), (2, 1, 1), (1, 2, 1), (1, 1, -1)):
        balanced[row, column] = ((F(0), F(sign, 2)), ZERO)
    phase = {i: ((F(0), F(1, 2)), (F(0), F((-1)**(i & 1), 2))) for i in range(8)}
    if kind == 'mixer':
        return balanced
    if kind == 'phase':
        return {(i, i): value for i, value in phase.items()}
    return {key: (_multiply(real, phase[key[1]][0]),
                  _multiply(real, phase[key[1]][1])) for key, (real, _) in balanced.items()}


def ideal_local_matrix(kind='mixer'):
    """Floating independent oracle only; never consumed by source execution."""
    result = np.zeros((8, 8), dtype=complex)
    for key, (real, imag) in target(kind).items():
        result[key] = complex(float(real[0])+math.sqrt(2)*float(real[1]),
                              float(imag[0])+math.sqrt(2)*float(imag[1]))
    return result


def exact_effect_audit():
    """Check complete target unitarity and X/Y pointer effects algebraically."""
    zero = (ZERO, ZERO)
    half = (F(1, 2), F(0))
    for kind in KINDS:
        matrix = target(kind)
        for row in range(8):
            for column in range(8):
                gram, effect = zero, zero
                for k in range(8):
                    term = _complex_product(_conjugate(matrix.get((k, row), zero)),
                                            matrix.get((k, column), zero))
                    gram = _add(gram[0], term[0]), _add(gram[1], term[1])
                    weight = F(0) if k == 1 else F(1) if k == 2 else F(1, 2)
                    effect = (_add(effect[0], tuple(weight*x for x in term[0])),
                              _add(effect[1], tuple(weight*x for x in term[1])))
                if gram != ((ONE, ZERO) if row == column else zero):
                    raise ArithmeticError('the complete exact source target is not unitary')
                expected = (half, ZERO) if row == column else zero
                if kind == 'phase' and row == column:
                    expected = ((F(0) if row == 1 else F(1) if row == 2 else F(1, 2), F(0)), ZERO)
                if kind == 'mixer' and (row, column) in ((1, 2), (2, 1)):
                    expected = half, ZERO
                if kind == 'y_analyzer' and (row, column) in ((1, 2), (2, 1)):
                    expected = ZERO, (F(1 if row == 1 else -1, 2), F(0))
                if effect != expected:
                    raise ArithmeticError('the exact complete native quadrature effect changed')
    return {'complete_unitarity_entries_checked': 3*64,
            'complete_pointer_effect_entries_checked': 3*64,
            'exact_field': 'Q(sqrt(2)) + i Q(sqrt(2))',
            'arbitrary_untouched_references': True}


def certify(matrix=None, *, kind='mixer'):
    """Original-binary, all-input/reference bounds against the exact target.

    eps >= ||A-U||F >= ||A-U||op implies half-diamond error
    eps*(2+eps)/2 and channel norm <=(1+eps)^2, since U is unitary.
    Two uses separated by exact CPTP/adaptive processing cost delta*(1+g).
    This retains numerical norm gain. No proof-object authentication is
    claimed: receipts freshly execute the source and compute this comparison.
    """
    kind = _kind(kind)
    matrix = local_matrix(kind) if matrix is None else matrix
    if type(matrix) is not np.ndarray:
        raise ValueError('an original base ndarray is required before array dispatch')
    matrix = exact._array(matrix, (8, 8))
    wanted = target(kind)
    norm = cross = ZERO
    reported = F(0)
    for row in range(8):
        for column in range(8):
            ar, ai = exact._component(matrix[row, column])
            real, imag = wanted.get((row, column), (ZERO, ZERO))
            reported += ar*ar+ai*ai
            norm = _add(norm, _add(_multiply(real, real), _multiply(imag, imag)))
            cross = _add(cross, (ar*real[0]+ai*imag[0], ar*real[1]+ai*imag[1]))
    if norm != (F(8), F(0)):
        raise ArithmeticError('the exact declared target norm changed')
    lower, upper = exact._radical_sum({F(1): cross[0], F(2): cross[1]})
    squared_lo, squared_hi = exact._nonnegative(reported+8-2*upper, reported+8-2*lower)
    epsilon = exact._sqrt_interval(squared_hi)[1]
    delta = epsilon*(2+epsilon)/2
    gain = (1+epsilon)**2
    return {'operator_error_upper': epsilon,
            'frobenius_error_lower': exact._sqrt_interval(squared_lo)[0],
            'frobenius_error_upper': epsilon,
            'half_diamond_upper': delta,
            'channel_norm_upper': gain,
            'two_call_half_diamond_upper': delta*(1+gain),
            'full_lift_dimension': 512,
            'untouched_reference_dimension_unrestricted': True}


def resources(kind='mixer'):
    """Existing compiler accounting; physical service error remains unbounded."""
    return compiler.resources(word(_kind(kind)))


def _outward(value, *, upper=True):
    scaled = value*10**12
    integer = -(-scaled.numerator//scaled.denominator) if upper else scaled.numerator//scaled.denominator
    return str(F(integer, 10**12))


def run():
    effect_audit = exact_effect_audit()
    transported = embedding()
    local = ideal_local_matrix()
    if not np.array_equal(np.kron(np.eye(64), local)@transported,
                          transported@ideal_native_matrix()):
        # The two independent float spellings of 1/sqrt(2) can differ by one
        # ulp; actual matching is exact from the checked ordered address pairs.
        if np.linalg.norm(np.kron(np.eye(64), local)@transported-
                          transported@ideal_native_matrix(), 'fro') > 1e-14:
            raise ArithmeticError('the whole binary encoding changed the mixer target')
    # Exact Walsh coefficients of the reflection sum to -pi at address3
    # and zero at every other address. This checks all spectator-independent
    # seam coordinates, not only the two native orientation columns.
    walsh = [sum((F((-1)**((3 & mask).bit_count()+(index & mask).bit_count()), 8)
                  for mask in range(8)), F(0)) for index in range(8)]
    if walsh != [F(index == 3) for index in range(8)]:
        raise ArithmeticError('the complete diagonal reflection lost its scalar or phase')
    rows = []
    for kind in KINDS:
        bound = certify(kind=kind)
        report._gate(bound['operator_error_upper'], ERROR_GATE, 'actual compiled native interface')
        counts = resources(kind)
        if counts['native_pulses'] != {'mixer': 34, 'phase': 1, 'y_analyzer': 35}[kind]:
            raise ArithmeticError('source schedule pulse count changed')
        rows.append({'operation': kind, 'primitive_pulses': counts['native_pulses'],
            'single_code_pulses': counts['single_code_pulses'],
            'pair_code_pulses': counts['pair_code_pulses'],
            'proper_code_service_events': counts['proper_code_service_events'],
            'identity_phase_pulses': counts['identity_phase_pulses'],
            'maximum_work_qubits': counts['maximum_work_qubits'],
            'full_spectator_lift_dimension': 512,
            'operator_error_upper': _outward(bound['operator_error_upper']),
            'half_diamond_upper': _outward(bound['half_diamond_upper']),
            'channel_norm_upper': _outward(bound['channel_norm_upper']),
            'two_call_half_diamond_upper': _outward(bound['two_call_half_diamond_upper']),
            'binary_array_bounds_include_all_eight_seam_coordinates': True,
            'physical_service_error_upper': None})
    return {'schema': 'oph.compiled_native_interface.v1', 'all_pass': True,
        'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS},
        'source_operations': rows,
        'exact_operator_checks': effect_audit,
        'encoding': {'binary_dimension': 66, 'full_carrier_dimension': 78,
            'full_code_dimension': 512, 'local_seam_dimension': 8,
            'spectator_labels': 64, 'native_orientation_pairs': 10,
            'ordered_seam_pair': [2, 1], 'six_high_spectator_bits_unchanged': True,
            'padding_target_is_identity': False},
        'imaginary_quadrature': {'additional_source_pulses_before_balanced_analyzer': 1,
            'native_effects': '(I +/- Y_orientation)/2',
            'conjugate_phase_mixture_plus_probabilities': ['1', '0'],
            'previous_real_effect_plus_probabilities': ['1/2', '1/2'],
            'reduced_orientation_qubit_tomography_on_known_active_support': True,
            'full_66_dimensional_tomography_claimed': False},
        'scope': {'one_existing_source_control_family': True,
            'independent_balanced_mixer_service_needed': False,
            'coherent_d2_d4_proper_code_services_supplied': True,
            'balanced_angle_source_selected': False,
            'native_committed_partial_injections_alone_sufficient': False,
            'full_source_map_and_padding_retained': True,
            'subsequent_floating_density_products_certified': False,
            'physical_service_error_upper': None,
            'new_empirical_observations': False,
            'new_Lean_theorem': False}}
