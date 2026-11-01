"""Unchanged cube source field words and a local fixed-number Ramsey pointer."""
from fractions import Fraction as F
import math

from ..cube_eta_protection.cube import GEOMETRY
from ..observer.source_scalar import compiler
from ..common_source.chiral_candidate import gaussian_source as gaussian
from ..interacting_field import source
from ..full_fock_field.field import source_word
from ..certified_reference import certificate as cert, continuation as templates
from . import execute


POINTER = 16
BITS = 17


def block(word, bits, targets): return execute.compile_block(execute.key(word), bits), targets


def controlled_phase_plan():
    return (('III', F(1, 4)), ('ZII', F(-1, 8)), ('IZI', F(-1, 8)),
            ('IIZ', F(-1, 4)), ('ZIZ', F(1, 8)), ('IZZ', F(1, 8)))


def emit(plan):
    return [event for label, angle in plan for event in compiler.pauli_rotation(label, float(angle)*math.pi)]


def control_word(inverse=False):
    plan = controlled_phase_plan()
    return emit(tuple((label, -angle) for label, angle in reversed(plan)) if inverse else plan)


def pointer_prepare():
    return [compiler.pulse('Y', -math.pi/4, [2])]+control_word()


def pointer_analyzer(phase):
    if type(phase) is not F or phase not in (F(), F(1, 2)): raise ValueError('phase0 or pi/2 required')
    return control_word(True)+[compiler.pulse('I', float(phase)*math.pi/2, [2]),
                              compiler.pulse('Z', -float(phase)*math.pi/2, [2]),
                              compiler.pulse('Y', math.pi/4, [2])]


def pointer_template(phase=None):
    plan = controlled_phase_plan()
    if phase is None:
        return [('Y', (2,), cert._pi(-1, 4))]+[event for label, angle in plan for event in cert.pauli(label, cert.PiMultiple(angle))]
    return ([event for label, angle in reversed(plan) for event in cert.pauli(label, cert.PiMultiple(-angle))]+
            [('I', (2,), cert.PiMultiple(phase/2)), ('Z', (2,), cert.PiMultiple(-phase/2)),
             ('Y', (2,), cert._pi(1, 4))])


def field_blocks(duration, steps):
    if type(duration) is not F or type(steps) is not int or steps < 1: raise ValueError('exact duration and positive steps required')
    half = float(duration)/(2*steps)
    swap = gaussian.local_gate('fswap').word
    terms = tuple(term for group in GEOMETRY.groups for term in group)
    hops = {z: gaussian.local_gate('hopping', z, half).word for _, _, z in terms}
    _, contact = source.contact(half, GEOMETRY.interaction)
    blocks = []
    def hop(p, q, z):
        route = tuple(range(q-1, p, -1))
        for left in route: blocks.append(block(swap, 2, (left, left+1)))
        blocks.append(block(hops[z], 2, (p, p+1)))
        for left in reversed(route): blocks.append(block(swap, 2, (left, left+1)))
    for _ in range(steps):
        for p, q, z in terms: hop(p, q, z)
        for _ in range(2):
            for site in range(GEOMETRY.volume): blocks.append(block(contact['word'], 2, (2*site, 2*site+1)))
        for p, q, z in reversed(terms): hop(p, q, z)
    return tuple(blocks)


def flatten(blocks):
    return tuple(dict(kind=kind, angle=angle, targets=[targets[j] for j in local_targets])
                 for instrument, targets in blocks for kind, angle, local_targets in instrument.tape)


def field_binding(duration, steps):
    blocks = field_blocks(duration, steps)
    actual = source_word(GEOMETRY, float(duration), steps)
    if flatten(blocks) != actual: raise ArithmeticError('contracted blocks differ from the unchanged full cube source word')
    return blocks
