"""Compile certified rank-one tests into the existing nine-bit source words."""
from __future__ import annotations
from fractions import Fraction as F
from functools import lru_cache
import math
import numpy as np

from . import intervals as iv
from model.native_effect_coverage import coverage, detector
from model.native_interface import tomography as tomo
from model.native_program import interpreter
from model.observer.source_scalar import compiler

RAW_GATE = F(1, 10**11)
PULSE_GATE = F(2, 10**11)
MAX_PULSES = 23579
MAX_RATIONAL_BITS = 64


def rational(value):
    if type(value) is not F or max(abs(value.numerator).bit_length(), value.denominator.bit_length()) > MAX_RATIONAL_BITS:
        raise ValueError('an exact Fraction with at most64 numerator/denominator bits required')
    return value


def target(vector, strength=F(1)):
    if type(vector) not in (list, tuple) or len(vector) != 66:
        raise ValueError('66 exact Gaussian-rational source coordinates required')
    rows = []
    for pair in vector:
        if type(pair) not in (list, tuple) or len(pair) != 2:
            raise ValueError('a real/imaginary rational pair required')
        rows.append(tuple(rational(x) for x in pair))
    strength = rational(strength)
    norm = sum((x*x+y*y for x, y in rows), F())
    if not norm or not 0 <= strength <= 1:
        raise ValueError('a nonzero vector and strength in [0,1] required')
    return tuple(rows), strength, norm


def _walsh(values):
    """The entire rational transform is exactly representable in binary64."""
    length = len(values)
    if length & (length-1) or any(type(v) is not int for v in values):
        raise ValueError('a power-of-two list of integer angle numerators required')
    return tuple(iv.binary_exact(sum(v*(-1)**((i & mask).bit_count())
                                     for i, v in enumerate(values)), iv.AQ*length)
                 for mask in range(length))


def _state_word(weights):
    if len(weights) != 512 or sum(weights, F()) != 1 or min(weights) < 0:
        raise ValueError('all512 exact nonnegative squared amplitudes required')
    word = []; parameters = []
    for bit in reversed(range(9)):
        controls = 8-bit; angles = []
        for prefix in range(1 << controls):
            block = weights[prefix*(1 << (bit+1)):(prefix+1)*(1 << (bit+1))]
            angles.append(-iv.split_angle(sum(block[:1 << bit], F()),
                                           sum(block[1 << bit:], F())))
        coefficients = _walsh(angles)
        for mask, coefficient in enumerate(coefficients):
            label = ['I']*9; label[bit] = 'Y'
            for control in range(controls):
                if mask >> control & 1:
                    label[bit+1+control] = 'Z'
            word.extend(compiler.pauli_rotation(label, coefficient))
        # The exact inverse Walsh identity attaches the emitted angles to
        # every conditional split, including all zero-weight prefixes.
        for prefix, angle in enumerate(angles):
            reconstructed = sum((F(c)*(-1)**((prefix & mask).bit_count())
                                 for mask, c in enumerate(coefficients)), F())
            if reconstructed != F(angle, iv.AQ):
                raise ArithmeticError('a state-tree angle changed during native compilation')
        parameters.extend(angles)
    if len(word) != 12291 or len(parameters) != 511:
        raise ArithmeticError('the complete nine-bit amplitude-tree schedule changed')
    return word


def _phase_word(phases):
    if len(phases) != 512:
        raise ValueError('the complete work-register phase list is required')
    coefficients = _walsh(phases); word = []
    for mask, coefficient in enumerate(coefficients):
        label = ['Z' if mask >> bit & 1 else 'I' for bit in range(9)]
        word.extend(compiler.pauli_rotation(label, coefficient))
    for index, phase in enumerate(phases):
        actual = sum((F(c)*(-1)**((index & mask).bit_count())
                      for mask, c in enumerate(coefficients)), F())
        if actual != F(phase, iv.AQ):
            raise ArithmeticError('the emitted diagonal word lost a phase')
    if len(word) != 11270:
        raise ArithmeticError('the complete nine-bit diagonal schedule changed')
    return word


def rank_one_word(vector, strength=F(1), output=0):
    """One admitted word; its selected native output effect targets lambda P_v.

    No source-chart projection occurs inside the word. The selected address
    is moved from literal work address0 to an actual native source code.
    """
    vector, strength, norm = target(vector, strength)
    if type(output) is not int or not 0 <= output < 66:
        raise ValueError('an original native output label required')
    codes = tomo.codes(); padding = next(i for i in range(512) if i not in codes)
    weights = [F()]*512; phases = [0]*512
    for code, (x, y) in zip(codes, vector, strict=True):
        weights[code] = strength*(x*x+y*y)/norm
        phases[code] = iv.phase_angle(x, y)
    weights[padding] = 1-strength
    preparation = _state_word(weights)+_phase_word(phases)
    # The inverse is executed as actual negative native requests. An adjoint
    # of a returned floating matrix is not substituted for this word.
    move_zero = tomo.route_word(tuple(('X', bit) for bit in range(9)
                                      if codes[output] >> bit & 1))
    word = compiler.inverse(preparation)+move_zero
    if len(word) > MAX_PULSES:
        raise ArithmeticError('the rank-one word exceeded its fixed source budget')
    interpreter.source_word(word, 9)
    return word


def pair_word(a=0, b=1):
    """A short native word for the unequal-weight (3|a>+4|b>)/5 projector."""
    tomo.pair(a, b)
    theta = iv.split_angle(F(9), F(16))
    local = compiler.cnot(0, 1)
    for mask in range(4):
        label = ['Y', 'Z' if mask & 1 else 'I', 'Z' if mask & 2 else 'I']
        coefficient = iv.binary_exact(theta*(-1 if mask & 1 else 1), 4*iv.AQ)
        local.extend(compiler.pauli_rotation(label, coefficient))
    local += compiler.cnot(0, 1)
    routing = tomo.route_word(tomo.route(a, b))
    word = routing+local+compiler.inverse(routing)
    interpreter.source_word(word, 9)
    return word


def _complex(real=(0, 0), imag=(0, 0)):
    return real, imag


def primitive_target(kind, angle):
    """Full-six-level ideal pulse intervals at its exact binary64 angle."""
    sine, cosine = iv.sincos(F(angle)); zero = iv.exact(0); one = iv.exact(1)
    out = [[_complex(one if i == j else zero, zero) for j in range(6)] for i in range(6)]
    if kind == 'J':
        for j in range(3): out[j][j] = _complex(cosine, sine)
    elif kind == 'I':
        for j in range(2): out[j][j] = _complex(cosine, sine)
    elif kind == 'Z':
        out[0][0] = _complex(cosine, sine); out[1][1] = _complex(cosine, iv.neg(sine))
    elif kind == 'X':
        out[0][0] = out[1][1] = _complex(cosine, zero)
        out[0][1] = out[1][0] = _complex(zero, sine)
    elif kind == 'Y':
        out[0][0] = out[1][1] = _complex(cosine, zero)
        out[0][1] = _complex(sine, zero); out[1][0] = _complex(iv.neg(sine), zero)
    elif kind == 'H':
        from math import isqrt
        root = isqrt(iv.Q*iv.Q//2); half_root = root, root+1
        part = iv.mul(sine, half_root)
        out[0][0] = _complex(cosine, part); out[1][1] = _complex(cosine, iv.neg(part))
        out[0][1] = out[1][0] = _complex(zero, part)
    else:
        raise ValueError('an existing native pulse kind required')
    return out


@lru_cache(maxsize=32768)
def certify_primitive(kind, angle):
    actual = compiler.native_pulse(kind, angle); wanted = primitive_target(kind, angle)
    squared = sum((iv.scalar_error(actual[i, j].real, wanted[i][j][0])**2+
                   iv.scalar_error(actual[i, j].imag, wanted[i][j][1])**2
                   for i in range(6) for j in range(6)), F())
    if squared > RAW_GATE**2:
        raise ArithmeticError('a stored full-six-level pulse exceeds its exact target gate')
    return True


def certify_word(word):
    interpreter.source_word(word, 9)
    pi_lower, pi_upper = iv.pi_interval()
    if max(abs(F(math.pi)-pi_lower), abs(F(math.pi)-pi_upper)) > F(1, 10**15):
        raise ArithmeticError('the fixed basis/parity pi payload exceeds its allowance')
    kinds = {(event['kind'], event['angle']) for event in word}
    for kind, angle in kinds:
        certify_primitive(kind, angle)
    return len(kinds)


def bounds(pulses, *, general):
    if type(pulses) is not int or not 1 <= pulses <= MAX_PULSES or type(general) is not bool:
        raise ValueError('a bounded nonempty native synthesis word required')
    x = pulses*PULSE_GATE
    target_error = F(11 if general else 1, iv.AQ)
    operator_error = x/(1-x)+target_error
    failure = pulses*PULSE_GATE**2/(1-x)**2
    complete = (operator_error*(2+operator_error)+failure)/2
    return {'pulses':pulses,'fixed_binary_angle_target_error_upper':str(target_error),
            'stored_primitive_operator_gate':str(PULSE_GATE),
            'full_word_operator_error_upper':str(operator_error),
            'all_first_processor_failure_mass_upper':str(failure),
            'complete_halfdiamond_upper':str(complete),
            'complete_map_gain_upper':str(1/(1-x)**2+failure)}


def work_row(word, address):
    """Diagnostic row of the actual emitted word, not a forward-word adjoint.

    Reversed adjoints of these very stored blocks evaluate the desired row.
    The uniform operator certificate covers exact products and all references;
    binary64 vector propagation is only an independent diagnostic.
    """
    interpreter.source_word(word, 9)
    if type(address) is not int or not 0 <= address < 512:
        raise ValueError('one actual work-register address required')
    vector = np.zeros(512, complex); vector[address] = 1
    for event in reversed(word):
        full = compiler.native_pulse(event['kind'], event['angle'])
        size = 1 << len(event['targets'])
        vector = compiler._apply(vector.reshape(512, 1), full[:size, :size].conj().T,
                                 event['targets'], 9).ravel()
    return vector.conj()


def selected_row(word, output=0):
    if type(output) is not int or not 0 <= output < 66:
        raise ValueError('one original native output label required')
    return work_row(word, tomo.codes()[output])


def coarse_output(word, source_reference):
    """Success-array diagnostic after tracing the tree's fine tape and source.

    Processor first failures have separate exact norm bounds and are not
    silently renormalized into these success blocks. The returned aggregate
    contains66 decoded initial labels and the terminal padding export. The
    fine instrument is exposed by fine_success_row and processor_failure.
    """
    if not isinstance(source_reference, np.ndarray) or source_reference.ndim != 2 or source_reference.shape[0] != 66:
        raise ValueError('source by untouched-reference amplitude columns required')
    work = np.zeros((512, source_reference.shape[1]), complex)
    work[list(tomo.codes())] = source_reference
    actual = compiler.execute(word, 9, work)
    blocks = [np.outer(actual[code], actual[code].conj()) for code in tomo.codes()]
    padding = sorted(set(range(512))-set(tomo.codes()))
    blocks.append(actual[padding].T@actual[padding].conj())
    return blocks


def fine_success_row(word, slots):
    """One fine native map: sqrt(weight)|terminal><actual_row|.

    The row is a numerical diagnostic of the unchanged stored-pulse product;
    the weight, complete native slot address and terminal label are exact.
    Zero native outcomes return an explicit zero row and remain addressable.
    """
    weight, mapping, records = detector.branch(tuple(slots))
    if not mapping:
        return weight, None, np.zeros(66, complex), records
    (initial, terminal), = mapping.items()
    return weight, terminal, selected_row(word, initial)[list(tomo.codes())], records


def padding_kraus_row(word, address):
    """Internal export-Kraus diagnostic; the observed abort flag sums all446."""
    if type(address) is not int or not 0 <= address < 512 or address in tomo.codes():
        raise ValueError('one of the446 internal padding Kraus rows required')
    return work_row(word, address)[list(tomo.codes())]


def processor_failure(word, pulse_index, complement, source_reference):
    """Query one retained first-failure branch, without conditioning or scaling.

    This is the same reset-to-zero export as interpreter.execute_word.
    Spectators and reference columns remain present at the terminal failure.
    """
    interpreter.source_word(word, 9)
    if type(pulse_index) is not int or not 0 <= pulse_index < len(word):
        raise ValueError('an actual emitted pulse address required')
    event = word[pulse_index]; size = 1 << len(event['targets'])
    if type(complement) is not int or not size <= complement < 6:
        raise ValueError('a retained processor-complement row required')
    source_reference = np.asarray(source_reference)
    if (source_reference.ndim != 2 or source_reference.shape[0] != 66
            or source_reference.dtype.kind not in 'iufc' or not np.isfinite(source_reference).all()):
        raise ValueError('source by untouched-reference amplitude columns required')
    work = np.zeros((512, source_reference.shape[1]), complex)
    work[list(tomo.codes())] = source_reference
    for index, current in enumerate(word[:pulse_index+1]):
        full = compiler.native_pulse(current['kind'], current['angle'])
        dimension = 1 << len(current['targets'])
        if index == pulse_index:
            block = np.zeros((dimension, dimension), complex)
            block[0] = full[complement, :dimension]
        else:
            block = full[:dimension, :dimension]
        work = compiler._apply(work, block, current['targets'], 9)
    return work


def outcome_catalogue(word):
    """Lazy complete alphabet; even zero processor and native outcomes exist."""
    return {'processor_first_failure_labels':sum(6-(1 << len(e['targets'])) for e in word),
            'terminal_padding_kraus_rows':446,'observed_terminal_padding_abort_flags':1,
            'native_tree_raw_words':25**11,
            'nonzero_native_tree_words':66*2**11,
            'native_raw_slots_per_attempt':25,'native_attempts':11,
            'early_failure_stops_word':True,'all_failure_outputs_retained':True,
            'native_detector_applied_only_after_complete_address_export':True,
            'padding_projection_inside_coherent_word':False,
            'fine_native_weight_and_partial_map_api':'fine_success_row / native_effect_coverage.detector.branch',
            'processor_failure_query_is_normalized_sampler':False}
