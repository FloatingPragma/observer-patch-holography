"""An independent certificate of the existing literal plaquette source words.

Tree coefficients are read from the tape, not accepted from atan2 or norm
rounding claims. Exact inverse Walsh sums and trig enclosures check every
conditional two-vector. Fixed basis/CNOT/blank/frame gates and twelve-port
coefficient rounding receive separate full-source operator error bounds.
"""
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
import hashlib

from ..plaquette_reference import preparation as p
from ..observer.source_scalar import compiler
from . import interval as iv
from . import source


@dataclass(frozen=True)
class PiMultiple:
    value: F

    def __neg__(self):
        return PiMultiple(-self.value)


def _pi(n, d=1):
    return PiMultiple(F(n, d))


def inverse(template):
    return [(kind, targets, None if angle is None else -angle)
            for kind, targets, angle in reversed(template)]


def cnot(control, target):
    return [('H', (target,), _pi(1, 2)), ('J', (control, target), _pi(-1)),
            ('H', (target,), _pi(1, 2))]


def pauli(label, angle=None):
    support = [i for i, v in enumerate(label) if v != 'I']
    if not support:
        return [('I', (0,), angle)]
    basis = []
    for i, letter in enumerate(label):
        if letter == 'X': basis.append(('Y', (i,), _pi(-1, 4)))
        elif letter == 'Y': basis.append(('X', (i,), _pi(1, 4)))
        elif letter not in 'IZ': raise ValueError('tensor-Pauli label required')
    pivot = support[-1]
    parity = [event for i in support[:-1] for event in cnot(i, pivot)]
    return inverse(basis)+parity+[('Z', (pivot,), angle)]+inverse(parity)+basis


def frame_template():
    blocks = ((('ZI', _pi(-1, 4)), ('IZ', _pi(1, 4))),
              (('XY', _pi(1, 4)), ('YX', _pi(-1, 4))),
              (('XX', _pi(1, 4)), ('YY', _pi(1, 4))))
    forward = []
    for site_index, site in enumerate(p.GEOMETRY.sites):
        for count, block in zip(reversed(site), blocks, strict=True):
            for _ in range(count % 4):
                for label, angle in block:
                    forward += [(k, tuple(2*site_index+j for j in t), a)
                                for k, t, a in pauli(label, angle)]
    return inverse(forward)


@lru_cache(maxsize=128)
def _angle_error(actual, desired):
    return (iv.pi_angle_error(actual, desired.value) if type(desired) is PiMultiple
            else abs(actual-desired))


def consume(word, offset, template):
    if offset+len(template) > len(word):
        raise ArithmeticError('literal source word ended inside a required gate')
    error, variables = F(), []
    for index, (kind, targets, desired) in enumerate(template, offset):
        event = word[index]
        if set(event) != {'kind', 'targets', 'angle'} or event['kind'] != kind or tuple(event['targets']) != targets:
            raise ArithmeticError('literal source word differs from independent gate routing')
        if type(event['angle']) is not float:
            raise ValueError('the retained word must contain literal binary64 angles')
        actual = F(event['angle'])
        if desired is None:
            variables.append(actual)
        else:
            error += _angle_error(actual, desired)
    return offset+len(template), error, variables


def inverse_walsh(values):
    result = list(values)
    step = 1
    while step < len(result):
        for start in range(0, len(result), 2*step):
            for j in range(start, start+step):
                a, b = result[j], result[j+step]
                result[j], result[j+step] = a+b, a-b
        step *= 2
    return result


def two_vector_error(angle, a, b, signed):
    """Bound the full two-vector error, including all signed final leaves."""
    if type(a) is not int or type(b) is not int or type(signed) is not bool:
        raise ValueError('exact integer tree data and signed-leaf flag required')
    if signed:
        total = a*a+b*b
        if total:
            x = iv.sqrt_fraction(F(a*a, total))*(1 if a >= 0 else -1)
            y = iv.sqrt_fraction(F(b*b, total))*(1 if b >= 0 else -1)
        else:
            x, y = iv.Interval.point(1), iv.Interval.point(0)
    else:
        if min(a, b) < 0:
            raise ValueError('nonnegative subtree weights required')
        total = a+b
        x, y = ((iv.sqrt_fraction(F(a, total)), iv.sqrt_fraction(F(b, total))) if total
                else (iv.Interval.point(1), iv.Interval.point(0)))
    sine, cosine = iv.sine_cosine(angle)
    # exp(i theta Y)|0>=(cos(theta),-sin(theta)); L1 bounds Euclidean norm.
    return (cosine-x).absolute_upper()+(-sine-y).absolute_upper()


def amplitude_certificate(word, values, *, probability=False, offset=0, bits=8):
    """Check every level against integer target data, without target angles."""
    if type(bits) is not int or not 7 <= bits <= 12:
        raise ValueError('the retained binary-tree compiler uses seven through twelve bits')
    if (type(probability) is not bool or len(values) != 1 << bits
            or any(type(x) is not int for x in values)
            or (probability and min(values) < 0) or not any(values)):
        raise ValueError('complete nonzero integer amplitude/weight data required')
    fixed_error, level_errors = F(), []
    weights = list(values) if probability else [v*v for v in values]
    prefix_sum = [0]
    for weight in weights: prefix_sum.append(prefix_sum[-1]+weight)
    nodes = 0
    for target in reversed(range(bits)):
        controls = bits-target-1
        coefficients = []
        for mask in range(1 << controls):
            label = ['I']*bits; label[target] = 'Y'
            for j in range(controls):
                if mask >> j & 1: label[target+j+1] = 'Z'
            offset, fixed, variable = consume(word, offset, pauli(label))
            if len(variable) != 1:
                raise ArithmeticError('each commuting Walsh factor must have one literal coefficient')
            coefficients.append(variable[0]); fixed_error += fixed
        angles = inverse_walsh(coefficients)
        error = F()
        for prefix, angle in enumerate(angles):
            start, half = prefix*(1 << (target+1)), 1 << target
            signed = target == 0 and not probability
            if signed:
                a, b = values[start], values[start+1]
            else:
                a = prefix_sum[start+half]-prefix_sum[start]
                b = prefix_sum[start+2*half]-prefix_sum[start+half]
            error = max(error, two_vector_error(angle, a, b, signed))
            nodes += 1
        level_errors.append(error)
    return {'offset': offset, 'fixed_angle_error': fixed_error,
            'tree_error': sum(level_errors, F()), 'level_errors': level_errors,
            'checked_conditional_nodes': nodes}


def blank_template(index):
    if type(index) is not int or not 0 <= index < 256:
        raise ValueError('all 256 retained records are the declared domain')
    return [('X', (j,), _pi(1, 2)) for j in range(8) if index >> j & 1]+[
        ('I', (0,), _pi(-index.bit_count(), 2))]


def controls_bound(word, table):
    counts = Counter((event['kind'], event['angle']) for event in word)
    result = F()
    for key, count in counts.items():
        if key not in table:
            controls = source.literal_controls(*key)
            table[key] = source.generator_error(key[0], F(key[1]), controls)
        result += count*table[key]
    return result


def word_digest(word):
    digest = hashlib.sha256()
    for event in word:
        digest.update((event['kind']+'|'+event['angle'].hex()+'|'+
                       ','.join(map(str, event['targets']))+';').encode())
    return digest.hexdigest()


def certify_preparation():
    data = p.design()
    source.specification()
    exports = source.complete_export_certificate()
    table = {}
    probability_word = p.probability_word(data)
    probability = amplitude_certificate(probability_word, data.column_norm_squared, probability=True)
    if probability['offset'] != len(probability_word):
        raise ArithmeticError('unaccounted probability preparation event')
    probability_source = controls_bound(probability_word, table)
    probability_error = probability_source+probability['fixed_angle_error']+probability['tree_error']
    rows, maximum = [], F()
    for index in range(256):
        word = p.conditional_word(index, data)
        offset, blank_error, variables = consume(word, 0, blank_template(index))
        if variables:
            raise ArithmeticError('blanking angle must be charged explicitly')
        target = [row[index] for row in data.numerator]
        amplitude = amplitude_certificate(word, target, offset=offset)
        offset, frame_error, variables = consume(word, amplitude['offset'], frame_template())
        if variables or offset != len(word):
            raise ArithmeticError('inverse frame contains an unaccounted pulse')
        port_error = controls_bound(word, table)
        total = blank_error+amplitude['tree_error']+amplitude['fixed_angle_error']+frame_error+port_error
        maximum = max(maximum, total)
        rows.append({'record': index, 'pulses': len(word), 'native_word_sha256': word_digest(word),
                     'tree_error': amplitude['tree_error'],
                     'fixed_basis_and_CNOT_error': amplitude['fixed_angle_error'],
                     'blanking_error': blank_error, 'inverse_frame_error': frame_error,
                     'rounded_port_generator_error': port_error,
                     'complete_conditional_output_error': total})
    # Orthogonal records turn conditional error into the maximum branch bound,
    # not the sum over 256 possible trajectories. All abort labels stay in the
    # common Stinespring output; telescoping is valid on arbitrary references.
    complete = probability_error+maximum
    if complete >= p.EXECUTION_ADMISSION:
        raise ArithmeticError('proved literal-source error exceeds the previous declared tolerance')
    controls_digest = hashlib.sha256()
    for kind, angle in sorted(table):
        values = source.literal_controls(kind, angle)
        controls_digest.update((kind+'|'+angle.hex()+'|'+','.join(map(str, values))+';').encode())
    return {'probability_pulses': len(probability_word),
            'probability_native_word_sha256': word_digest(probability_word),
            'probability_tree_error': probability['tree_error'],
            'probability_fixed_angle_error': probability['fixed_angle_error'],
            'probability_rounded_port_generator_error': probability_source,
            'probability_complete_output_error': probability_error,
            'conditional_records': rows,
            'max_conditional_complete_output_error': maximum,
            'complete_preparation_output_error': complete,
            'all_conditional_tree_nodes': 257*255,
            'unique_literal_pulse_controls': len(table),
            'literal_integrated_port_controls_sha256': controls_digest.hexdigest(),
            'maximum_pulse_generator_error': max(table.values()),
            'exports': exports}
