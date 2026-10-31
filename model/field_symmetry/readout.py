"""Complete source preparation, evolution and readout of the eta tower.

The initialized coherent even-charge preparation is an explicit apparatus
input. No source state, failed export or terminal address is normalized away.
"""
from __future__ import annotations

from fractions import Fraction as F
import math
import numpy as np
from scipy.sparse.linalg import expm_multiply

from model.full_fock_field import field
from model.full_fock_field.local_simulation import split_bound
from model.field_thermal.thermal import execute_with_exports, failure_probability
from model.observer.source_scalar import compiler


PULSE_CEILING = F(1, 10**11)
OUTPUT_CEILING = 2e-9


def geometry(plane=0, interaction=1.):
    if type(plane) is not int or plane not in (0, 1, 2):
        raise ValueError('integer plane index zero, one or two required')
    return field.FieldGeometry(tuple(1 if axis == plane else 2 for axis in range(3)),
                               interaction=interaction)


def preparation(g):
    """Even onsite rotation, from eight initialized blank field bits."""
    if type(g) is not field.FieldGeometry or sorted(g.sides) != [1, 2, 2]:
        raise ValueError('one complete four-site open plaquette required')
    word = []
    for site, point in enumerate(g.sites):
        sign = (-1)**sum(point)
        # Conjugating the rotation by CNOT makes the complete four-state
        # macro even. A final CNOT alone prepares the right vacuum state
        # but changes parity on other inputs and is not this even service.
        word += compiler.cnot(2*site, 2*site+1)
        word += [compiler.pulse('Y', -sign*math.pi/4, [2*site])]
        word += compiler.cnot(2*site, 2*site+1)
    return word


def tower_state(g, phase=0.):
    """Independent product amplitudes in the original CAR chart."""
    if type(g) is not field.FieldGeometry or g.modes > 16:
        raise ValueError('bounded complete field geometry required')
    if type(phase) not in (int, float) or not math.isfinite(phase):
        raise ValueError('finite real phase required')
    result = np.zeros(1 << g.modes, complex)
    for subset in range(1 << g.volume):
        address, sign = 0, 1
        for site, point in enumerate(g.sites):
            if subset >> site & 1:
                address |= 3 << (2*site)
                sign *= (-1)**sum(point)
        result[address] = sign*np.exp(1j*phase*subset.bit_count())/math.sqrt(2**g.volume)
    return result


def _ideal_pulse(event):
    """Analytic primitive columns, independent of source coefficient solving."""
    kind, angle = event['kind'], event['angle']
    if kind == 'J':
        block = np.diag([np.exp(1j*angle)]*3+[1.])
    else:
        x = np.array([[0, 1], [1, 0]], complex)
        y = np.array([[0, -1j], [1j, 0]], complex)
        z = np.diag([1., -1.])
        generator = {'X': x, 'Y': y, 'Z': z, 'H': (x+z)/math.sqrt(2), 'I': np.eye(2)}[kind]
        block = math.cos(angle)*np.eye(2)+1j*math.sin(angle)*generator
    result = np.zeros((6, len(block)), complex)
    result[:len(block)] = block
    return result


def pulse_budget(word):
    """Check each distinct actual primitive, including its export complement.

    This is a fixed numerical admission ceiling, not an outward certificate
    for floating arithmetic or physical hardware. The exact telescoping
    inequality also covers possible numerical norm growth of accepted blocks.
    """
    maximum = 0.
    for kind, angle in set((e['kind'], e['angle']) for e in word):
        event = {'kind': kind, 'angle': angle}
        target = _ideal_pulse(event)
        actual = compiler.native_pulse(kind, angle)
        if (type(actual) is not np.ndarray or actual.shape != (6, 6)
                or not np.isfinite(actual).all()):
            raise ArithmeticError('nonfinite or malformed primitive response')
        error = float(np.linalg.norm(actual[:, :target.shape[1]]-target, 2))
        if not math.isfinite(error) or error > float(PULSE_CEILING):
            raise ArithmeticError('primitive complete-output admission failed')
        maximum = max(maximum, error)
    total = len(word)*PULSE_CEILING
    if total >= 1:
        raise ValueError('word exceeds numerical telescoping domain')
    return total/(1-total), maximum


def split_tower_contract(g):
    """Check the actual split schedule groups complete commuting spatial bonds."""
    bonds = 0
    for group in g.groups:
        if len(group) % 2:
            raise ArithmeticError('incomplete spatial bond in eta split schedule')
        for first, second in zip(group[::2], group[1::2], strict=True):
            p, q, _ = first
            r, s, _ = second
            if p//2 != r//2 or q//2 != s//2 or len({p, q, r, s}) != 4:
                raise ArithmeticError('split factors do not form a complete commuting bond')
            bonds += 1
    if g.speed and bonds != g.spatial_edges:
        raise ArithmeticError('eta schedule omitted a complete spatial bond')
    return bonds


def execute(plane=0, *, centered=False, interaction=1., duration=F(1), steps=1):
    """Literal source execution on all 256 states, retaining every failure."""
    if type(centered) is not bool:
        raise ValueError('explicit Boolean centered-generator choice required')
    if (type(duration) not in (int, F) or abs(duration) > 2
            or type(steps) is not int or not 1 <= steps <= 4):
        raise ValueError('bounded exact duration and integer split count required')
    g = geometry(plane, interaction)
    split_tower_contract(g)
    prep = preparation(g)
    evolution = list(field.source_word(g, float(duration), steps))
    if centered:
        for mode in range(g.modes):
            angle = g.interaction*float(duration)/4
            evolution += [compiler.pulse('I', angle, [mode]),
                          compiler.pulse('Z', -angle, [mode])]
    word = prep+evolution+compiler.inverse(prep)
    budget, primitive_error = pulse_budget(word)
    initial = np.zeros(1 << g.modes, complex); initial[0] = 1
    actual, failures = execute_with_exports(word, g.modes, initial)
    if (type(failures) is not tuple or len(failures) != len(word)
            or actual.shape != initial.shape or not np.isfinite(actual).all()):
        raise ArithmeticError('complete source output inventory changed')
    for event, outputs in zip(word, failures, strict=True):
        d = 1 << len(event['targets'])
        if type(outputs) is not tuple or len(outputs) != 6-d:
            raise ArithmeticError('export complement outcome omitted')
        for output in outputs:
            if (type(output) is not np.ndarray or output.shape != (len(initial)//d,)
                    or not np.isfinite(output).all()):
                raise ArithmeticError('export complement output malformed')
    phi = 0. if centered else -g.interaction*float(duration)
    # Decode the exact product directly, without the producer's compiler.
    target = np.zeros_like(actual)
    a, b = (1+np.exp(1j*phi))/2, (np.exp(1j*phi)-1)/2
    for subset in range(1 << g.volume):
        address, coefficient = 0, 1.+0j
        for site, point in enumerate(g.sites):
            if subset >> site & 1:
                address |= 3 << (2*site)
                coefficient *= (-1)**sum(point)*b
            else:
                coefficient *= a
        target[address] = coefficient
    lost = failure_probability(failures)
    error = math.sqrt(float(np.linalg.norm(actual-target))**2+lost)
    mass = float(np.vdot(actual, actual).real)+lost
    if (not np.isfinite([lost, error, mass]).all() or lost < 0
            or abs(mass-1) > OUTPUT_CEILING or error > float(budget)):
        raise ArithmeticError('complete eta source history failed admission')
    probabilities = np.abs(actual)**2
    read_probabilities = [float(probabilities[(np.arange(len(actual)) >> (2*x)) & 1 == 1].sum())
                          for x in range(g.volume)]
    general_split = split_bound(g.sides, duration, steps, F(g.speed), F(g.interaction), refined=True)
    summary = {
        'plane': plane, 'sides': list(g.sides), 'generator': 'K=H-U*N/2' if centered else 'H',
        'interaction': g.interaction, 'duration': str(duration), 'steps': steps,
        'terminal_probabilities': probabilities.tolist(), 'read_one_probabilities': read_probabilities,
        'ideal_read_one_probability': math.sin(phi/2)**2,
        'first_export_failure_probability': lost, 'complete_output_error': error,
        'complete_probability_residual': abs(mass-1), 'raw_primitive_residual': primitive_error,
        'numerical_primitive_ceiling': str(PULSE_CEILING),
        'numerical_dilation_budget': str(budget),
        'state_specific_ideal_split_error': '0',
        'arbitrary_input_ideal_split_bound': str(general_split),
        'resources': compiler.resources(word), 'initialized_field_bits': g.modes,
        'terminal_raw_outcomes': len(actual), 'first_failure_event_slots': len(failures),
        'first_failure_kraus_outputs': sum(len(event) for event in failures),
        'preparation_is_supplied_even_charge_coherence': True,
        'floating_arithmetic_or_hardware_certified': False,
    }
    return {'output': actual, 'failures': failures, 'target': target, 'report': summary}


def independent_evolution(plane=0, duration=F(1)):
    """Complete Hamiltonian oracle for the eta state; no split-word shortcut."""
    g = geometry(plane)
    space = field.FullFock(g.modes)
    actual = expm_multiply(-1j*float(duration)*space.hamiltonian(g), tower_state(g))
    return float(np.linalg.norm(actual-tower_state(g, -float(duration))))
