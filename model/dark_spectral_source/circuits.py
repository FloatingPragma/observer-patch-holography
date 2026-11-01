"""Sparse phase-exact dark-state circuit, compiled to the existing source.

No statevector or desired unitary is passed to the native processor.  The
small Pauli plan contains only rational multiples of pi and is compiled
through the existing one/pair-code source words, including scalar phases.
"""
from fractions import Fraction as F
from functools import lru_cache
import math

from ..observer.source_scalar import compiler
from ..certified_reference import certificate as cert
from ..certified_reference import continuation as templates
from ..full_fock_field.field import source_word
from ..common_spin_protection.spin import GEOMETRY


MODES = 12
POINTER = 12
BITS = 13


def label(*items, bits=MODES):
    value = ['I']*bits
    for bit, pauli in items:
        if type(bit) is not int or not 0 <= bit < bits or pauli not in 'XYZ' or value[bit] != 'I':
            raise ValueError('distinct valid Pauli support required')
        value[bit] = pauli
    return ''.join(value)


def exact_bit_flip(bit):
    return [(label(), F(-1, 2)), (label((bit, 'X')), F(1, 2))]


def exact_cnot(control, target):
    # exp(-i*pi/4*(I-Zc)) exp(i*pi/4*(Xt-Zc Xt)).
    return [(label(), F(-1, 4)), (label((control, 'Z')), F(1, 4)),
            (label((target, 'X')), F(1, 4)),
            (label((control, 'Z'), (target, 'X')), F(-1, 4))]


def split(first, second):
    # exp(i*pi/8*(Ya Xb-Xa Yb)) sends |10> to (|10>+|01>)/sqrt2.
    return [(label((first, 'Y'), (second, 'X')), F(1, 8)),
            (label((first, 'X'), (second, 'Y')), F(-1, 8))]


@lru_cache(maxsize=1)
def preparation_plan():
    plan = []
    for bit in (0, 2, 4):
        plan += exact_bit_flip(bit)
    for first, second in ((4, 8), (4, 6), (8, 10)):
        plan += split(first, second)
    for bit in (4, 6, 8, 10):
        plan += exact_cnot(bit, bit+1)
    for bit in (6, 10):
        plan += [(label(), F(-1, 2)), (label((bit, 'Z')), F(1, 2))]
    # Inverse of the actual second-quantized onsite frame R_x=(i sigma_x)^x.
    for site in reversed(range(6)):
        for _ in range(site % 4):
            plan += [(label((2*site, 'Y'), (2*site+1, 'Y')), F(-1, 4)),
                     (label((2*site, 'X'), (2*site+1, 'X')), F(-1, 4))]
    return tuple(plan)


def controlled(plan):
    out = []
    for pauli, angle in plan:
        if len(pauli) != MODES:
            raise ValueError('twelve field bits required before adding the pointer')
        out += [(pauli+'I', angle/2), (pauli+'Z', -angle/2)]
    return tuple(out)


def inverse(plan):
    return tuple((pauli, -angle) for pauli, angle in reversed(plan))


def emit(plan):
    return [event for pauli, multiple in plan
            for event in compiler.pauli_rotation(pauli, float(multiple)*math.pi)]


def plan_template(plan):
    return [event for pauli, multiple in plan
            for event in cert.pauli(pauli, cert.PiMultiple(multiple))]


def preparation_word(*, pointer=True):
    if type(pointer) is not bool:
        raise ValueError('explicit pointer choice required')
    if pointer:
        return [compiler.pulse('Y', -math.pi/4, [POINTER])]+emit(controlled(preparation_plan()))
    return emit(preparation_plan())


def preparation_template(*, pointer=True):
    if pointer:
        return [('Y', (POINTER,), cert._pi(-1, 4))]+plan_template(controlled(preparation_plan()))
    return plan_template(preparation_plan())


def analyzer_word(phase):
    if type(phase) is not F or phase not in (F(), F(1, 2)):
        raise ValueError('exact Ramsey phase0 or pi/2 required')
    return (emit(inverse(controlled(preparation_plan())))+
            [compiler.pulse('I', float(phase)*math.pi/2, [POINTER]),
             compiler.pulse('Z', -float(phase)*math.pi/2, [POINTER]),
             compiler.pulse('Y', math.pi/4, [POINTER])])


def analyzer_template(phase):
    return (plan_template(inverse(controlled(preparation_plan())))+
            [('I', (POINTER,), cert.PiMultiple(phase/2)),
             ('Z', (POINTER,), cert.PiMultiple(-phase/2)),
             ('Y', (POINTER,), cert._pi(1, 4))])


def field_template(duration, steps):
    """Independent literal routing for the unchanged six-site H, not K."""
    if type(duration) is not F or type(steps) is not int or steps < 1:
        raise ValueError('exact model duration and positive integer steps required')
    half = duration/(2*steps)
    tape = []
    terms = tuple(term for group in GEOMETRY.groups for term in group)
    def append(template, left):
        tape.extend((kind, tuple(bit+left for bit in targets), angle)
                    for kind, targets, angle in template)
    def move(first, second, z):
        route = tuple(range(second-1, first, -1))
        for left in route: append(templates.swap(), left)
        append(templates.hopping(z, half), first)
        for left in reversed(route): append(templates.swap(), left)
    for _ in range(steps):
        for first, second, z in terms: move(first, second, z)
        for _ in range(2):
            for site in range(GEOMETRY.volume): append(templates.contact(half), 2*site)
        for first, second, z in reversed(terms): move(first, second, z)
    return tape


def field_word(duration, steps):
    return source_word(GEOMETRY, float(duration), steps)


def certify_word(word, template):
    offset, angle, variables = cert.consume(word, 0, template)
    if offset != len(word) or variables:
        raise ArithmeticError('native word has an unaccounted event or angle')
    controls = cert.controls_bound(word, {})
    return {'pulses': len(word), 'word_sha256': cert.word_digest(word),
            'exact_angle_rounding_bound': angle,
            'literal_port_generator_bound': controls,
            'complete_source_operator_error': angle+controls}
