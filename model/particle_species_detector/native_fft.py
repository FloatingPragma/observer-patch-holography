"""Constructive native-source interface for the1024-mode ideal Fourier detector.

The executed array is the complete N1 Gaussian surrogate.  Literal native
words can leave that sector; their full-Fock, complete-export error is
bounded independently, never projected away and called exact execution.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import math

import numpy as np

from ..observer.source_scalar import compiler
from ..common_source.chiral_candidate import gaussian_source as gaussian
from ..certified_reference import certificate as cert, continuation as templates


SIDE = 8
MODES = 1024


@dataclass(frozen=True)
class Operation:
    kind: str
    first: int
    second: int
    phase_pi: F = F()


def address(coordinates, spin):
    if (len(coordinates) != 3 or any(type(x) is not int or not 0 <= x < SIDE for x in coordinates)
            or type(spin) is not int or spin not in (0, 1)):
        raise ValueError('cube coordinates0..7 and spin0 or1 required')
    x, y, z = coordinates
    return 2*((x*SIDE+y)*SIDE+z)+spin


@lru_cache(maxsize=1)
def plan():
    """Immutable three-axis DIT FFT, including every input bit reversal."""
    operations = []
    for axis in range(3):
        others = [j for j in range(3) if j != axis]
        for fixed in product(range(SIDE), repeat=2):
            for spin in (0, 1):
                line = []
                for coordinate in range(SIDE):
                    xyz = [0, 0, 0]; xyz[axis] = coordinate
                    for i, value in zip(others, fixed, strict=True): xyz[i] = value
                    line.append(address(xyz, spin))
                # Three-bit reversal exchanges1<->4 and3<->6.
                operations.extend((Operation('swap', line[1], line[4]), Operation('swap', line[3], line[6])))
                for width in (2, 4, 8):
                    half = width//2
                    for start in range(0, SIDE, width):
                        for j in range(half):
                            operations.append(Operation('butterfly', line[start+j], line[start+j+half], F(1)-F(2*j, width)))
    result = tuple(operations)
    if (sum(o.kind == 'butterfly' for o in result) != 4608 or
            sum(o.kind == 'swap' for o in result) != 768):
        raise ArithmeticError('complete FFT stage or bit-reversal count changed')
    return result


def adjacent_route(operation):
    """Exact CAR transposition/bring-together route as adjacent FSWAP indices."""
    p, q = operation.first, operation.second
    if not 0 <= p < q < MODES: raise ValueError('ordered valid Fourier mode addresses required')
    if operation.kind == 'swap':
        return tuple(range(q-1, p-1, -1))+tuple(range(p+1, q))
    if operation.kind == 'butterfly':
        return tuple(range(q-1, p, -1))
    raise ValueError('known Fourier operation required')


def apply(state, inverse=False):
    """Execute ideal Gaussian Fourier action on ALL1024 N1 rows/references."""
    if type(inverse) is not bool or np.ma.isMaskedArray(state): raise ValueError('unmasked input and Boolean inverse required')
    value = np.asarray(state)
    vector = value.ndim == 1
    if vector: value = value[:, None]
    if (value.ndim != 2 or value.shape[0] != MODES or not 1 <= value.shape[1] <= MODES or
            value.size > 2**20 or value.dtype.kind not in 'ficu' or not np.isfinite(value).all()):
        raise ValueError('all1024 finite N1 rows and bounded reference columns required')
    value = np.asarray(value, complex).copy()
    operations = tuple(reversed(plan())) if inverse else plan()
    for op in operations:
        p, q = op.first, op.second
        if op.kind == 'swap': value[[p, q]] = value[[q, p]]; continue
        a, b = value[p].copy(), value[q].copy()
        phase = np.exp(1j*math.pi*float(op.phase_pi))
        if inverse:
            value[p] = (a+b)/math.sqrt(2)
            value[q] = phase.conjugate()*(-a+b)/math.sqrt(2)
        else:
            value[p] = (a-phase*b)/math.sqrt(2)
            value[q] = (a+phase*b)/math.sqrt(2)
    return value[:, 0] if vector else value


def butterfly_word(phase_pi, inverse=False):
    if phase_pi not in (F(1), F(3, 4), F(1, 2), F(1, 4)):
        raise ValueError('an actual FFT8 butterfly phase required')
    phase = math.pi*float(phase_pi)
    word = [compiler.pulse('I', phase/2, [1]), compiler.pulse('Z', -phase/2, [1])]
    word += list(gaussian.local_word('hopping', -1j, math.pi/4))
    if inverse:
        word = [{**event, 'angle': -event['angle']} for event in reversed(word)]
    return tuple(word)


def butterfly_template(phase_pi, inverse=False):
    hop = (templates.phase(0, cert._pi(1, 2))+
           templates.real_hopping(cert._pi(1, 4))+
           templates.phase(0, cert._pi(-1, 2)))
    template = templates.phase(1, cert.PiMultiple(phase_pi))+hop
    return cert.inverse(template) if inverse else template


def word_bound(word, template):
    offset, error, variables = cert.consume(word, 0, template)
    if offset != len(word) or variables: raise ArithmeticError('unconsumed native Fourier source event')
    return error+cert.controls_bound(word, {})


def word_counts(word):
    return {'pulses': len(word),
            'first_export_abort_slots': sum(6-(1 << len(event['targets'])) for event in word)}


def _root(power):
    """zeta8^power in the exact integer basis1,zeta,zeta²,zeta³."""
    power %= 8
    return tuple((-1 if power >= 4 else 1) if j == power % 4 else 0 for j in range(4))


def _add(a, b): return tuple(x+y for x, y in zip(a, b, strict=True))
def _negative(a): return tuple(-x for x in a)


def _multiply(a, b):
    result = (0, 0, 0, 0)
    for j, x in enumerate(a):
        for k, y in enumerate(b):
            result = _add(result, tuple(x*y*z for z in _root(j+k)))
    return result


def _conjugate(a):
    result = (0, 0, 0, 0)
    for j, value in enumerate(a): result = _add(result, tuple(value*z for z in _root(-j)))
    return result


@lru_cache(maxsize=1)
def exact_fft8_identity():
    """All entries and Gram sums, without host trigonometry or residuals."""
    zero, one = (0, 0, 0, 0), _root(0)
    matrix = [[one if i == j else zero for j in range(8)] for i in range(8)]
    matrix[1], matrix[4] = matrix[4], matrix[1]
    matrix[3], matrix[6] = matrix[6], matrix[3]
    for width in (2, 4, 8):
        for start in range(0, 8, width):
            for j in range(width//2):
                p, q = start+j, start+j+width//2
                a, b = matrix[p], matrix[q]
                wb = [_multiply(_root(-8*j//width), entry) for entry in b]
                matrix[p] = [_add(x, y) for x, y in zip(a, wb, strict=True)]
                matrix[q] = [_add(x, _negative(y)) for x, y in zip(a, wb, strict=True)]
    if matrix != [[_root(-k*x) for x in range(8)] for k in range(8)]:
        raise ArithmeticError('exact normalized FFT8 entry identity failed')
    for k in range(8):
        for j in range(8):
            dot = zero
            for x in range(8): dot = _add(dot, _multiply(matrix[k][x], _conjugate(matrix[j][x])))
            if dot != (8*int(k == j), 0, 0, 0): raise ArithmeticError('exact Fourier orthogonality failed')
    return {'coefficient_ring': 'Z[zeta8]/(zeta8^4+1)',
            'exact_unnormalized_entries_checked': 64, 'exact_Gram_entries_checked': 64,
            'per_axis_normalization': '1/sqrt8',
            'three_axis_normalization': '1/sqrt512',
            'three_axis_identity': 'tensor three FFT8 transforms and identity on physical spin'}


@lru_cache(maxsize=1)
def certificate():
    swap = gaussian.local_word('fswap')
    swap_error = word_bound(swap, templates.swap())
    swap_counts = word_counts(swap)
    phases = sorted({op.phase_pi for op in plan() if op.kind == 'butterfly'})
    errors = {(phase, inverse): word_bound(butterfly_word(phase, inverse), butterfly_template(phase, inverse))
              for phase in phases for inverse in (False, True)}
    route_swaps = sum(len(adjacent_route(op))*(2 if op.kind == 'butterfly' else 1) for op in plan())
    counts = {phase: sum(op.kind == 'butterfly' and op.phase_pi == phase for op in plan()) for phase in phases}
    fft_errors = [route_swaps*swap_error+sum((counts[p]*errors[p, inverse] for p in phases), F())
                  for inverse in (False, True)]
    pulse_count = route_swaps*swap_counts['pulses']+sum(counts[p]*len(butterfly_word(p)) for p in phases)
    abort_slots = route_swaps*swap_counts['first_export_abort_slots']+sum(counts[p]*word_counts(butterfly_word(p))['first_export_abort_slots'] for p in phases)
    # Selected-mode excitation from a completely initialized vacuum. Y(-pi/2)
    # sends the blank bit to|1> with no unspecified relative scalar.
    excitation = [compiler.pulse('Y', -math.pi/2, [0])]
    excitation_error = word_bound(excitation, [('Y', (0,), cert._pi(-1, 2))])
    copy = compiler.cnot(0, 1)
    copy_error = word_bound(copy, cert.cnot(0, 1))
    # A fresh pointer is copied/exported for each mode. All binary strings
    # and abort histories remain; only ideal N1 inputs are one-hot.
    pointers_error = MODES*copy_error
    return {'geometry': [8, 8, 8], 'physical_spin_modes_per_site': 2,
            'exact_Fourier_factorization': exact_fft8_identity(),
            'N1_dimension': MODES, 'FFT_stages': 9, 'two_mode_butterflies': 4608,
            'explicit_bit_reversal_transpositions': 768,
            'adjacent_FSWAPs_per_Fourier_transform': route_swaps,
            'native_pulses_per_Fourier_transform': pulse_count,
            'first_export_abort_slots_per_Fourier_transform': abort_slots,
            'Fourier_sign_convention': 'U[k,x]=exp(-2pi i k dot x/8)/sqrt512, separately for each spin',
            'forward_complete_Fock_source_error': str(fft_errors[0]),
            'inverse_complete_Fock_source_error': str(fft_errors[1]),
            'selected_mode_excitation_source_error': str(excitation_error),
            'all_mode_occupation_pointer_source_error': str(pointers_error),
            'occupation_pointer_CNOTs': MODES, 'occupation_pointer_native_pulses': MODES*len(copy),
            'occupation_pointer_first_export_abort_slots': MODES*word_counts(copy)['first_export_abort_slots'],
            'occupation_pointer_complete_final_labels_each': 6,
            'preparation_plus_detection_complete_source_error': str(excitation_error+sum(fft_errors, F())+pointers_error),
            'initialization_premise': 'complete field vacuum and a fresh blank pointer for each occupation test',
            'ideal_occupation_POVM_on_N1': 'all1024 one-hot mode projectors; no output mode removed',
            'full_occupation_dictionary': 'all binary occupation strings plus every source/export abort history; no conditioning on one-hot survival',
            'original_H_propagation': 'separate exact model-law service, not synthesized by this Fourier interface',
            'full_Fock_native_state_array_executed': False,
            'native_map_assumed_exactly_Gaussian_or_number_preserving': False,
            'N1_surrogate_projection_of_native_residuals_used': False,
            'ideal_second_quantized_Fourier_identity': True,
            'source_error_scope': 'complete native outputs on arbitrary full-Fock inputs and untouched references, except selected excitation requires initialized vacuum'}
