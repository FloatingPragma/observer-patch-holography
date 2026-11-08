"""Refinable complex Householder words through the existing source API."""
from fractions import Fraction as F

from model.centered_native_preparation.tree import walsh
from model.certified_reference import certificate as cert
from model.native_ordered_ground import angles
from model.native_spin_observer import precision


def elimination_denominator(real, imaginary, sign):
    """Exact d=1-s*x_j in H=I-(x-s e_j)(x-s e_j)*/d."""
    if (type(real) is not F or type(imaginary) is not F
            or type(sign) is not int or sign not in (-1, 1)
            or sign*real > F(1, 2)):
        raise ValueError('exact complex pivot and a stable signed target required')
    return 1-sign*real, -sign*imaginary


def _normal(value):
    if (type(value) is not tuple or len(value) < 2 or len(value) & (len(value)-1)
            or any(type(z) is not tuple or len(z) != 2
                   or any(type(x) is not int for x in z) for z in value)
            or not any(a or b for a, b in value)):
        raise ValueError('a nonzero complete complex integer vector is required')
    return value


def diagonal(phases):
    """Exactly synthesize a rational-phase diagonal through commuting parities."""
    if (type(phases) is not tuple or len(phases) < 2 or len(phases) & (len(phases)-1)
            or any(type(x) is not F for x in phases)):
        raise ValueError('a complete rational phase table is required')
    n = len(phases).bit_length()-1
    coefficients = tuple(x/len(phases) for x in walsh(phases))
    word = [('I', (0,), coefficients[0])]
    for pivot in range(n):
        count = 1 << pivot
        for j in range(count):
            gray = j ^ (j >> 1)
            word.append(('Z', (pivot,), coefficients[(1 << pivot)|gray]))
            if pivot:
                following = ((j+1) % count) ^ (((j+1) % count) >> 1)
                word += cert.cnot((gray ^ following).bit_length()-1, pivot)
    return word


def preparation(normal, angle_bits):
    """Prepare an arbitrary complex rational direction, with no zero oracle.

    Integer zero components are decidable. Every nonzero quadrant is handled
    by the existing terminating rational square-root/arctangent algorithm.
    No host float chooses an angle, precision or certificate.
    """
    normal = _normal(normal)
    if type(angle_bits) is not int or angle_bits < 1:
        raise ValueError('positive exact angle precision required')
    n = len(normal).bit_length()-1
    weights = tuple(a*a+b*b for a, b in normal)
    prefix = [0]
    for w in weights:
        prefix.append(prefix[-1]+w)
    word = []
    for bit in reversed(range(n)):
        raw = []
        for label in range(1 << (n-bit-1)):
            start, half = label*(1 << (bit+1)), 1 << bit
            a = prefix[start+half]-prefix[start]
            b = prefix[start+2*half]-prefix[start+half]
            raw.append(angles.node(a, b, angle_bits)['angle'] if a+b else F())
        coefficients = tuple(x/len(raw) for x in walsh(raw))
        count = len(raw)
        for j in range(count):
            gray = j ^ (j >> 1)
            word.append(('Y', (bit,), coefficients[gray]))
            if count > 1:
                following = ((j+1) % count) ^ (((j+1) % count) >> 1)
                control = bit+1+(gray ^ following).bit_length()-1
                word += cert.cnot(control, bit)
    lo, hi = precision.pi_interval(angle_bits+2)
    pi = (lo+hi)/2
    phases = []
    for a, b in normal:
        if not (a or b):
            phases.append(F())
            continue
        quadrant = -angles.node(a*a, b*b, angle_bits)['angle']
        phases.append((quadrant if b >= 0 else -quadrant) if a >= 0
                      else (pi-quadrant if b >= 0 else -pi+quadrant))
    word += diagonal(tuple(phases))
    bound = F(n+2, 1 << angle_bits)
    return word, {'bits': n, 'complex_normal_integer_components': len(normal),
                  'pulses': len(word), 'angle_bits': angle_bits,
                  'preparation_unitary_operator_error_upper': str(bound),
                  'host_float_used_in_angle_selection': False}


def rank_one_phase(normal, phase, angle_bits):
    """Approximate I+(exp(i phase)-1)|v><v| on every input column."""
    if type(phase) is not F:
        raise ValueError('a rational rank-one phase is required')
    preparation_word, record = preparation(normal, angle_bits)
    phases = (phase,)+(F(),)*(len(normal)-1)
    marked = diagonal(phases)
    # Words act in list order: inverse U, marked phase, U.
    word = cert.inverse(preparation_word)+marked+preparation_word
    error = 2*F(record['preparation_unitary_operator_error_upper'])
    return word, {'dimension': len(normal), 'pulses': len(word),
                  'rank_one_operator_error_upper': str(error),
                  'all_input_columns_and_untouched_reference': True}
