"""Explicit native word for the source-derived cube ground, with aborts kept."""
from fractions import Fraction as F
from functools import lru_cache
import math

import numpy as np

from model.interacting_vacuum import ground
from model.cube_eta_protection import cube
from model.observer.source_scalar import compiler
from model.reference_recovery.recovery import execute_complete
from model.reference_observation.observation import frozen
from model.protected_particle_readout import execute
from model.certified_reference import certificate as cert, interval as iv
from model.field_symmetry.symmetry import onsite_frame
from model.field_stress.energy import GaussianInteger as GI


BITS = 8


def signed_weights():
    """Exact signed square-root numerators of G; no eigensolver input."""
    weights, signs = [0]*256, [1]*256
    for mask, (a, b) in ground.reference().items():
        if a and b:
            raise ArithmeticError('ground determinant is not a pure rational/sqrt3 term')
        if mask & 0xAAAA:
            raise ArithmeticError('hidden down modes must start blank')
        compact = sum(((mask >> (2*j)) & 1) << j for j in range(8))
        weights[compact] = a*a+3*b*b
        signs[compact] = 1 if a+b >= 0 else -1
    if sum(weights) != 1296 or sum(bool(w) for w in weights) != 62:
        raise ArithmeticError('source-derived ground normalization or support changed')
    return tuple(weights), tuple(signs)


def tree_word():
    weights, signs = signed_weights()
    # These host square roots propose literal angles only. certificate() below
    # admits the tape using independent rational/transcendental enclosures.
    return compiler.compile_real_state(np.array([s*math.sqrt(w)/36
                                                 for w, s in zip(weights, signs, strict=True)]))


def signed_tree_certificate(word, weights, signs):
    """Check a complete8-bit tree against exact signed square-root weights."""
    if (type(weights) is not tuple or type(signs) is not tuple or
            len(weights) != 256 or len(signs) != 256 or
            any(type(w) is not int or w < 0 for w in weights) or
            any(type(s) is not int or s not in (-1, 1) for s in signs) or not sum(weights)):
        raise ValueError('256 exact nonnegative weights and +/-1 signs required')
    sums = [0]
    for weight in weights: sums.append(sums[-1]+weight)
    offset, fixed, level_errors, nodes = 0, F(), [], 0
    for target in reversed(range(BITS)):
        controls = BITS-target-1
        coefficients = []
        for mask in range(1 << controls):
            label = ['I']*BITS; label[target] = 'Y'
            for j in range(controls):
                if mask >> j & 1: label[target+j+1] = 'Z'
            offset, error, variables = cert.consume(word, offset, cert.pauli(label))
            if len(variables) != 1:
                raise ArithmeticError('one variable coefficient required per Walsh factor')
            fixed += error; coefficients.append(variables[0])
        maximum = F()
        for prefix, angle in enumerate(cert.inverse_walsh(coefficients)):
            start, half = prefix*(1 << (target+1)), 1 << target
            left = sums[start+half]-sums[start]
            right = sums[start+2*half]-sums[start+half]
            total = left+right
            if total:
                x, y = iv.sqrt_fraction(F(left, total)), iv.sqrt_fraction(F(right, total))
                if target == 0:
                    x, y = signs[start]*x, signs[start+1]*y
            else:
                x, y = iv.Interval.point(1), iv.Interval.point(0)
            sine, cosine = iv.sine_cosine(angle)
            error = (cosine-x).absolute_upper()+(-sine-y).absolute_upper()
            maximum = max(maximum, error); nodes += 1
        level_errors.append(maximum)
    if offset != len(word):
        raise ArithmeticError('unconsumed preparation event')
    return {'conditional_nodes': nodes, 'level_errors': level_errors,
            'signed_tree_error': sum(level_errors, F()), 'fixed_angle_error': fixed}


def inverse_frame_plan(site):
    """Inverse exterior-power frame R_x=(iX)^x(iY)^y(iZ)^z, including phases."""
    if type(site) is not tuple or len(site) != 3 or any(type(x) is not int or x not in (0, 1) for x in site):
        raise ValueError('a source cube coordinate required')
    xyz = ((('XX', F(1, 4)), ('YY', F(1, 4))),
           (('XY', F(1, 4)), ('YX', F(-1, 4))),
           (('ZI', F(-1, 4)), ('IZ', F(1, 4))))
    forward = [pair for count, axis in zip(reversed(site), reversed(xyz), strict=True)
               for _ in range(count) for pair in axis]
    return tuple((label, -angle) for label, angle in reversed(forward))


def emit(plan):
    return [event for label, angle in plan
            for event in compiler.pauli_rotation(label, float(angle)*math.pi)]


def template(plan):
    return [event for label, angle in plan
            for event in cert.pauli(label, cert.PiMultiple(angle))]


def frame_blocks():
    blocks = []
    for x, site in enumerate(cube.GEOMETRY.sites):
        word = emit(inverse_frame_plan(site))
        if word:
            blocks.append((execute.compile_block(execute.key(word), 2), (2*x, 2*x+1)))
    return tuple(blocks)


def exact_frame_certificate():
    """Exact full local Fock action of every inverse-frame Pauli word.

    For a quarter-turn, exp(i sign*pi P/4)=(I+i sign P)/sqrt(2).
    Every local word has an even number of factors, so multiplying its
    Gaussian-integer numerators leaves an integer power of two denominator.
    This verifies the vacuum and doublon phases as well as the spin block.
    """
    matrices={'I':np.eye(2,dtype=complex), 'X':np.array([[0,1],[1,0]],complex),
              'Y':np.array([[0,-1j],[1j,0]],complex),
              'Z':np.diag([1,-1]).astype(complex)}
    identity=GI.from_exact(np.eye(4,dtype=np.int64))
    for site in cube.GEOMETRY.sites:
        value=identity; plan=inverse_frame_plan(site)
        if len(plan)%2:raise ArithmeticError('inverse frame lost paired quarter-turns')
        for label, angle in plan:
            if abs(angle)!=F(1,4):raise ArithmeticError('exact frame angle changed')
            pauli=GI.from_exact(np.kron(matrices[label[1]],matrices[label[0]]))
            factor=identity+(pauli.times_i() if angle>0 else GI.zero(4)-pauli.times_i())
            value=factor@value
        expected=np.eye(4,dtype=complex)
        expected[1:3,1:3]=onsite_frame(site).conj().T
        if not (value-GI.from_exact((2**(len(plan)//2))*expected)).is_zero():
            raise ArithmeticError('inverse frame differs on the complete local Fock space')
    return {'exact_local_Fock_columns':32,
            'vacuum_single_and_doublon_phases_checked':True,
            'arithmetic':'Gaussian integers divided by powers of two'}


def embed_blank_modes(value):
    """Isometry on all256 compact rows; no number or success normalization."""
    a = np.asarray(value)
    if a.ndim == 1: a = a[:, None]
    if a.ndim != 2 or a.shape[0] != 256:
        raise ValueError('complete8-bit state or reference columns required')
    result = np.zeros((65536, a.shape[1]), complex)
    for mask in range(256):
        result[sum(((mask >> j) & 1) << (2*j) for j in range(8))] = a[mask]
    return result


def ideal_ground():
    vector = np.zeros(65536, complex)
    for mask, (a, b) in ground.reference().items():
        vector[mask] = (a+b*math.sqrt(3))/36
    return cube.operators().frame.matrix().conj().T@vector


@lru_cache(maxsize=1)
def runs():
    initial = np.zeros(256); initial[0] = 1
    tree = execute_complete(tree_word(), 8, initial)
    frame = execute.complete(frame_blocks(), 16, embed_blank_modes(tree.success))
    return tree, frame


def certificate():
    word = tree_word(); weights, signs = signed_weights()
    checked = signed_tree_certificate(word, weights, signs)
    tree_error = checked['signed_tree_error']+checked['fixed_angle_error']+cert.controls_bound(word, {})
    frame_error, pulses = F(), 0
    for site in cube.GEOMETRY.sites:
        plan = inverse_frame_plan(site); tape = emit(plan)
        offset, error, variable = cert.consume(tape, 0, template(plan))
        if offset != len(tape) or variable: raise ArithmeticError('unaccounted inverse-frame event')
        frame_error += error+cert.controls_bound(tape, {}); pulses += len(tape)
    bound = tree_error+frame_error
    if bound >= F(1, 10**10): raise ArithmeticError('ground preparation failed exact precision target')
    return {'amplitudes_from_exact_original_action_determinants': True,
            'inverse_frame':exact_frame_certificate(),
            'ground_energy': '-sqrt(3)/4', 'ground_number': 4,
            'nonzero_hidden_occupation_amplitudes': 62,
            'squared_amplitude_denominator': 1296,
            'signed_conditional_nodes_checked': checked['conditional_nodes'],
            'tree_pulses': len(word), 'frame_pulses': pulses,
            'tree_word_sha256': cert.word_digest(word),
            'complete_tree_error': str(tree_error), 'complete_frame_error': str(frame_error),
            'complete_preparation_state_error': str(bound),
            'all_abort_outputs_retained': True, 'initialized_blank_input_required': True,
            'physical_hardware_accuracy_established': False,
            'ground_selected_by_autonomous_repair_or_refinement': False}
