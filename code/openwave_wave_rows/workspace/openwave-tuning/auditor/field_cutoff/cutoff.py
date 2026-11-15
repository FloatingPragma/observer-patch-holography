"""State-dependent cyclic-register transfer for the actual plaquette source.

No continuum, vacuum or mass identification is assumed. Numerical controls
use the complete one-particle sector of the declared three-step member.
"""
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
import hashlib
import json
import math
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / 'openwave-tuning/common_dynamics/vector_field'))
import autonomous as source


def exponential_tail_upper(value, degree):
    """Exact rational upper bound on sum_{k>=degree} value**k/k!."""
    value = Fraction(value)
    if value < 0 or type(degree) is not int or degree < 0:
        raise ValueError('nonnegative magnitude and integer order required')
    k = max(degree, math.floor(value) + 1)
    prefix = sum((value**j / math.factorial(j) for j in range(degree, k)), Fraction())
    first = value**k / math.factorial(k)
    return prefix + first / (1 - value / (k + 1))


def transfer_bound(modulus, steps, beta_sum, initial_cap=0, matter_increment=1):
    """Norm bound, hence an upper bound on terminal complete-record TV.

    beta_sum is the sum of absolute face coefficients per macro. For an
    arbitrary schedule use its total coefficient sum / steps, or the proof.
    A positive step's support may move at most matter_increment per edge.
    """
    if type(modulus) is not int or modulus < 2 or modulus % 2:
        raise ValueError('even finite modulus required')
    if any(type(v) is not int or v < 0 for v in (steps, initial_cap, matter_increment)):
        raise ValueError('nonnegative integral support budgets required')
    beta_sum = Fraction(beta_sum)
    if beta_sum < 0:
        raise ValueError('sum of absolute coefficients required')
    distance = modulus // 2 - initial_cap - steps * matter_increment
    if distance <= 0:
        return Fraction(2)
    return min(Fraction(2), 2 * exponential_tail_upper(steps * beta_sum, distance))


def initial(modulus):
    return {(0, bytes(288)): np.ones(12, dtype=complex) / np.sqrt(12)}


def step(state, modulus, alpha, beta):
    """Independent complete sparse implementation; no amplitude pruning."""
    neighbours, edge_ids = source.mb.common.geometry_tables()
    graph, B, A, labels = source.data()
    matter = defaultdict(lambda: np.zeros(12, dtype=complex))
    for (x, field), amplitudes in state.items():
        reflected = amplitudes.sum() / 6 - amplitudes
        for event, amplitude in enumerate(reflected):
            direction, move = divmod(event, 2)
            y = int(neighbours[x, direction]) if move else x
            field_out = field
            if move:
                edge = int(edge_ids[x, direction])
                sign = 1 if x == graph.seams[edge][0] else -1
                f = bytearray(field)
                f[edge] = (f[edge] + sign) % modulus
                field_out = bytes(f)
            matter[y, field_out][event] += amplitude
    face = A[:, 0]
    pivot = int(np.flatnonzero(face)[0])
    sign = int(face[pivot])
    groups = defaultdict(lambda: np.zeros((modulus, 12), dtype=complex))
    for (x, field), amplitudes in matter.items():
        f = np.frombuffer(field, dtype=np.uint8).astype(np.int64)
        coordinate = sign * int(f[pivot]) % modulus
        canonical = bytes(((f - coordinate * face) % modulus).astype(np.uint8))
        groups[x, canonical][coordinate] += amplitudes
    t = np.arange(modulus)
    fourier = np.exp(2j * np.pi * np.outer(t, t) / modulus) / np.sqrt(modulus)
    kick = (fourier * np.exp(1j * beta * np.cos(2 * np.pi * t / modulus))) @ fourier.conj().T
    output = {}
    for (x, field), amplitudes in groups.items():
        out = kick @ amplitudes
        f = np.frombuffer(field, dtype=np.uint8).astype(np.int64)
        for coordinate in range(modulus):
            encoded = (f + coordinate * face) % modulus
            centered = (encoded + modulus // 2) % modulus - modulus // 2
            energy = int(centered @ centered) / 2
            output[x, bytes(encoded.astype(np.uint8))] = out[coordinate] * np.exp(-1j * alpha * energy)
    return output


def embedded(state, modulus):
    result = {}
    for (x, field), values in state.items():
        encoded = np.frombuffer(field, dtype=np.uint8).astype(np.int64)
        signed = (encoded + modulus // 2) % modulus - modulus // 2
        # Fixed-width signed encoding gives one common output alphabet.
        result[x, signed.astype('<i2').tobytes()] = values
    return result


def differences(first, second):
    norm2 = tv = 0.
    for key in first.keys() | second.keys():
        a = first.get(key, np.zeros(12, dtype=complex))
        b = second.get(key, np.zeros(12, dtype=complex))
        norm2 += float(np.vdot(a-b, a-b).real)
        tv += float(np.sum(abs(abs(a)**2 - abs(b)**2))) / 2
    return math.sqrt(norm2), tv


def run():
    paths = [Path(__file__), Path(source.__file__), Path(source.cell_complex.__file__),
             Path(source.mb.common.__file__), Path(source.mb.__file__)]
    pins = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    beta = Fraction.from_float(.41)
    states = {}
    norms = {}
    for modulus in (16, 32, 64):
        state = initial(modulus)
        for _ in range(3):
            state = step(state, modulus, .17, .41)
        norms[modulus] = sum(float(np.vdot(v,v).real) for v in state.values())
        assert abs(norms[modulus]-1) < 1e-11
        states[modulus] = embedded(state, modulus)
    original = source.initial()
    for _ in range(3):
        original = source.step(original, .17, .41)
    equality = differences(states[16], embedded(original, 16))
    assert max(equality) < 1e-11
    comparisons = []
    for modulus in (16, 32):
        error, tv = differences(states[modulus], states[64])
        bound = transfer_bound(modulus, 3, beta) + transfer_bound(64, 3, beta)
        assert error <= float(bound) + 1e-11
        comparisons.append({'modulus': modulus, 'reference_modulus': 64,
                            'state_vector_error': error, 'complete_position_field_event_TV': tv,
                            'rigorous_state_bound': float(bound), 'bound_rational': str(bound)})
    # Boundary vectors falsify a uniform operator-norm approximation.
    # A cyclic shift sends |m/2-1> to |-m/2>; the integer shift sends it to
    # |m/2>, orthogonal after the stated embedding.
    boundary_shift_norm_error = math.sqrt(2)
    cases = []
    for modulus in (16, 32, 64, 128):
        b = transfer_bound(modulus, 3, beta)
        cases.append({'modulus': modulus, 'distance': modulus//2-3,
                      'state_norm_bound_to_integer_rotor': float(b), 'rational': str(b),
                      'complete_terminal_record_TV_bound': min(1.,float(b))})
    # Complete 208-face coefficients require a substantially larger register.
    homogeneous = []
    for modulus in (16, 1024, 2048):
        b = transfer_bound(modulus, 3, 208 * beta)
        homogeneous.append({'modulus': modulus, 'bound': float(b)})
    assert pins == {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    result = {'schema': 'oph.local12.field_register_transfer.v1',
              'checks': {'actual_Z16_source_reproduced': True, 'complete_branch_norms': True,
                         'Z32_and_Z64_refinement_with_common_integer_output_labels': True,
                         'exact_rational_Dyson_tail': True,
                         'boundary_mutation_refutes_uniform_operator_transfer': True,
                         'homogeneous_face_count_retained': True, 'source_bytes_unchanged': True},
              'source_equality_norm_and_TV': equality, 'norms': norms,
              'comparisons': comparisons, 'one_face_bounds': cases,
              'homogeneous208_face_bounds': homogeneous,
              'boundary_shift_norm_error': boundary_shift_norm_error,
              'theorem': 'For initial support |F_e|<=E0, total matter support motion J per edge, and magnetic coefficient sum Lambda over the schedule, d=m/2-E0-J>0 gives norm error <=2 sum_{k>=d} Lambda^k/k! under the centered embedding. The same bound controls any common complete terminal instrument with arbitrary quantum references.',
              'limitations': ['This compares a cyclic register to an unbounded INTEGER compact-rotor source, not to a real-valued field, Maxwell theory or a vacuum.',
                              'The tested one-face coefficients are alpha=.17,beta=.41; all other207faces are zero. Homogeneous error uses all208coefficients.',
                              'Every field register and physical reset/readout implementation changes with modulus. No Z16 native-word cost or error is transferred to Z32/Z64.',
                              'Numerical residuals are floating controls; the stated analytic bound is evaluated with exact rational arithmetic.',
                              'Initial support/leakage, geometric refinement, clock, pulse precision and interacting continuum limits remain separate.'],
              'source_sha256': pins, 'official_rows_changed': []}
    (HERE/'FIELD_CUTOFF.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'checks':result['checks'],'comparisons':comparisons,
                      'one_face_bounds':[{k:v for k,v in r.items() if k!='rational'} for r in cases],
                      'homogeneous':homogeneous}))


if __name__ == '__main__':
    run()
