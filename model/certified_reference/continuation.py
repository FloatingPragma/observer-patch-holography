"""Independent native-template proof for the retained original field word."""
from fractions import Fraction as F

from ..full_fock_field.field import source_word
from ..full_fock_field.local_simulation import split_bound
from ..plaquette_reference import preparation as p
from .certificate import PiMultiple, _pi, pauli, consume, controls_bound, word_digest


def divide(angle, divisor):
    return PiMultiple(angle.value/divisor) if type(angle) is PiMultiple else angle/divisor


def phase(bit, angle):
    return [('I', (bit,), divide(angle, 2)), ('Z', (bit,), -divide(angle, 2))]


def real_hopping(angle):
    return pauli('XX', -divide(angle, 2))+pauli('YY', -divide(angle, 2))


def swap():
    return real_hopping(_pi(1, 2))+phase(0, _pi(1, 2))+phase(1, _pi(1, 2))


def hopping(z, duration):
    real, imag = F(z.real), F(z.imag)
    if real and imag:
        raise ValueError('the fixed Weyl/contact links lie on a real or imaginary axis')
    if real:
        angle, magnitude = (_pi(0) if real > 0 else _pi(1)), abs(real)
    elif imag:
        angle, magnitude = _pi(1 if imag > 0 else -1, 2), abs(imag)
    else:
        angle, magnitude = _pi(0), F()
    return phase(0, -angle)+real_hopping(duration*magnitude)+phase(0, angle)


def contact(duration):
    return (pauli('II', -duration/4)+pauli('ZI', duration/4)
            +pauli('IZ', duration/4)+pauli('ZZ', -duration/4))


def field_template(duration, steps):
    if type(duration) is not F or type(steps) is not int or steps < 1:
        raise ValueError('exact model duration and positive integral step count required')
    half = duration/(2*steps)
    tape = []
    terms = tuple(term for group in p.GEOMETRY.groups for term in group)
    def append(template, left):
        tape.extend((k, tuple(j+left for j in t), a) for k, t, a in template)
    def move(first, second, z):
        route = tuple(range(second-1, first, -1))
        for left in route: append(swap(), left)
        append(hopping(z, half), first)
        for left in reversed(route): append(swap(), left)
    for _ in range(steps):
        for first, second, z in terms: move(first, second, z)
        for _ in range(2):
            for site in range(p.GEOMETRY.volume): append(contact(half*F(p.GEOMETRY.interaction)), 2*site)
        for first, second, z in reversed(terms): move(first, second, z)
    return tape


def certify():
    parameter, duration, steps = F(1, 2), F(-1, 64), 2
    word = source_word(p.GEOMETRY, float(duration), steps)
    offset, angle_error, variables = consume(word, 0, field_template(duration, steps))
    if offset != len(word) or variables:
        raise ArithmeticError('continued original field word has an unaccounted event')
    table = {}
    generator_error = controls_bound(word, table)
    splitting = split_bound(p.GEOMETRY.sides, abs(duration), steps)
    return {'modular_parameter': parameter, 'model_duration': duration, 'steps': steps,
            'native_pulses': len(word), 'native_word_sha256': word_digest(word),
            'fixed_angle_rounding_error': angle_error,
            'rounded_port_generator_error': generator_error,
            'literal_source_complete_operator_error': angle_error+generator_error,
            'ideal_Strang_error': splitting,
            'full_action_complete_operator_error': splitting+angle_error+generator_error,
            'all_Fock_inputs_and_untouched_references': True,
            'all_export_abort_complements_retained': True}
