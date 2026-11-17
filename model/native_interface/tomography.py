"""Full native state identification from compiled quadratures and current tapes.

The architecture supplies coherent proper-code routing and the existing source
analyzers. The detector is the finite current-only configuration controller.
Every trial keeps native export abort and detector timeout; nothing is
conditioned on success. This is destructive repeated-preparation tomography,
not nondisturbing observation of an unknown state on one continuing carrier.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
from itertools import combinations
import json
import math
from pathlib import Path

import numpy as np

from tools import strict_json
from ..compiled_native_interface import construction as compiled
from ..common_source import exact_output as exact
from ..native_clock import configuration
from ..native_modular_reference import local_preparation
from ..observer.source_scalar import compiler

ROOT = Path(__file__).resolve().parents[2]
DIMENSION = 66
GATE = F(1, 10**11)
SOURCE_PATHS = tuple(sorted(set((*compiled.SOURCE_PATHS,
    *configuration.SOURCE_PATHS, *local_preparation.SOURCE_PATHS,
    'model/native_interface/tomography.py'))))


def codes():
    return tuple(compiled.circuit.encode(x) for x in compiled.quantum.STATES)


def pair(a, b):
    if any(type(v) is not int or not 0 <= v < DIMENSION for v in (a, b)) or a >= b:
        raise ValueError('an increasing pair of native basis indices required')
    return a, b


def route(a, b):
    """Exact affine bit operations taking codes[a],codes[b] to 2,1."""
    pair(a, b)
    x, y = codes()[a], codes()[b]
    difference = x ^ y
    pivot = (difference & -difference).bit_length()-1
    operations = [('X', bit) for bit in range(9) if x >> bit & 1]
    operations += [('C', pivot, bit) for bit in range(9)
                   if bit != pivot and difference >> bit & 1]
    if pivot:
        operations += [('C', pivot, 0), ('C', 0, pivot), ('C', pivot, 0)]
    operations += [('C', 0, 1), ('X', 1)]
    return tuple(operations)


def permutation(operations):
    result = np.arange(512)
    for op in operations:
        if op[0] == 'X':
            result ^= 1 << op[1]
        else:
            result ^= ((result >> op[1]) & 1) << op[2]
    return result


def route_word(operations):
    tape = []
    for op in operations:
        if op[0] == 'X':
            tape.extend((compiler.pulse('I', -math.pi/2, [op[1]]),
                         compiler.pulse('X', math.pi/2, [op[1]])))
        else:
            tape.extend(compiler.cnot(op[1], op[2]))
    return tape


def word(a, b, kind):
    pair(a, b)
    if type(kind) is not str or kind not in ('mixer', 'y_analyzer'):
        raise ValueError('a real or imaginary quadrature analyzer required')
    routing = route_word(route(a, b))
    return routing+compiled.word(kind)+compiler.inverse(routing)


def ideal(a, b, kind):
    """Independent target only; source execution always consumes word()."""
    word(a, b, kind)  # Validate the complete request.
    p = permutation(route(a, b))
    target = np.kron(np.eye(64), compiled.ideal_local_matrix(kind))
    return target[np.ix_(p, p)]


def _integer_error(actual, target):
    actual = exact._array(actual, target.shape)
    square = F()
    for actual_value, wanted in zip(actual.flat, target.flat, strict=True):
        real, imag = exact._component(actual_value)
        square += (real-int(wanted))**2+imag**2
    return exact._sqrt_interval(square)[1]


def design_certificate():
    """Exact effect differences generate every real Hermitian coordinate."""
    zero = (compiled.ZERO, compiled.ZERO)
    for kind in ('mixer', 'y_analyzer'):
        target = compiled.target(kind)
        for a in range(8):
            for b in range(8):
                plus = compiled._complex_product(compiled._conjugate(target.get((2, a), zero)),
                                                 target.get((2, b), zero))
                minus = compiled._complex_product(compiled._conjugate(target.get((1, a), zero)),
                                                  target.get((1, b), zero))
                observed = tuple(tuple(x-y for x, y in zip(p, m, strict=True))
                                 for p, m in zip(plus, minus, strict=True))
                expected = zero
                if (a, b) in ((2, 1), (1, 2)):
                    expected = ((compiled.ONE, compiled.ZERO) if kind == 'mixer' else
                                (compiled.ZERO, (F(-1 if a == 2 else 1), F())))
                if observed != expected:
                    raise ArithmeticError('the exact analyzer quadrature effect changed')
    return {'exact_effect_entries_checked': 128,
            'diagonal_coordinates': DIMENSION,
            'independent_real_and_imaginary_pairs': DIMENSION*(DIMENSION-1)//2,
            'real_hermitian_rank': DIMENSION**2}


def certificate():
    """One source certificate covers every route, including inverse pulses.

    Tensor placement preserves the operator norm. Each route
    consists solely of the freshly checked X and CNOT blocks. The analyzer
    already has an exact full-eight-address certificate. No numerical polar
    projection or normalization enters the source comparison. Inverse blocks
    are executed and certified separately, not assumed to be binary adjoints.
    The composition theorem concerns exact products of these certified block
    maps; subsequent floating state propagation remains diagnostic.
    """
    x_word = route_word((('X', 0),))
    cx_word = compiler.cnot(0, 1)
    x_target = np.array([[0, 1], [1, 0]], dtype=int)
    cx_target = np.eye(4, dtype=int)[:, [0, 3, 2, 1]]
    errors = [_integer_error(compiler.execute(x_word, 1), x_target),
              _integer_error(compiler.execute(cx_word, 2), cx_target),
              _integer_error(compiler.execute(compiler.inverse(x_word), 1), x_target),
              _integer_error(compiler.execute(compiler.inverse(cx_word), 2), cx_target)]
    errors += [compiled.certify(kind=k)['operator_error_upper']
               for k in ('mixer', 'y_analyzer')]
    if any(e > GATE for e in errors):
        raise ArithmeticError('source primitive exceeds its uniform envelope')
    compiled.exact_effect_audit()
    design = design_certificate()
    configuration.source_audit()
    counts, pulse_counts = [], []
    for a, b in combinations(range(DIMENSION), 2):
        routing = route(a, b)
        p = permutation(routing)
        if (p[codes()[a]], p[codes()[b]]) != (2, 1) or len(set(p)) != 512:
            raise ArithmeticError('the exact route does not isolate its two labels')
        counts.append(2*len(routing)+1)
        pulse_counts.append(2*len(route_word(routing))+len(compiled.word('y_analyzer')))
    # Telescoping operator error: product(1+eps_j)-1, with every raw gain.
    operator_error = (1+GATE)**max(counts)-1
    source_error = operator_error*(2+operator_error)/2
    detector = configuration.budget(mode='measure')
    return {'dimension': DIMENSION, 'settings': 1+DIMENSION*(DIMENSION-1),
            'exact_design': design,
            'exact_routes_checked': len(counts), 'all_route_addresses_checked': 512,
            'maximum_certified_blocks': max(counts),
            'maximum_source_pulses_per_trial': max(pulse_counts),
            'proper_code_service_events_per_trial': 2*max(pulse_counts),
            'native_attempts_per_trial': detector['native_attempts_max'],
            'source_operator_error_upper': str(operator_error),
            'source_half_diamond_upper': str(source_error),
            'detector_timeout_upper': str(detector['failure_upper']),
            'complete_record_half_l1_upper': str(source_error+detector['failure_upper']),
            'all_failures_retained': True,
            'same_readout_guarantee_for_C_and_E': True,
            'reference_safe_record_instrument': True,
            'composition_of_certified_block_maps': True,
            'floating_state_propagation_error_certified': False,
            'coherent_code_routing_and_repeatable_preparation_supplied': True,
            'physical_service_error_included': False}


def records(factor, a=None, b=None, kind=None, *, raw=True):
    """Complete label/timeout/export law, retaining untouched references.

    factor has axes (native basis, reference, purification). Output blocks
    act on the reference, with final work discarded as in configuration
    mode='measure'. Exact fine-tape identities justify the s_x probabilities.
    The returned arrays are numerical diagnostics, not interval enclosures.
    """
    if (type(factor) is not np.ndarray or factor.dtype.kind not in 'fc'
            or factor.ndim != 3 or factor.shape[0] != DIMENSION
            or min(factor.shape[1:]) < 1 or not np.isfinite(factor).all()
            or np.linalg.norm(factor)**2 > 1+1e-12 or type(raw) is not bool):
        raise ValueError('a finite subnormalized native/reference factor required')
    reference, rank = factor.shape[1:]
    vectors = np.zeros((512, reference*rank), complex)
    vectors[list(codes())] = factor.reshape(DIMENSION, -1)
    if kind is not None:
        tape = word(a, b, kind)
        vectors = compiler.execute(tape, 9, vectors) if raw else ideal(a, b, kind)@vectors
    elif a is not None or b is not None:
        raise ValueError('basis readout has no pair')
    tensor = vectors.reshape(512, reference, rank)
    blocks = np.einsum('irk,isk->irs', tensor, tensor.conj())
    success = configuration.success_probabilities(mode='measure') if raw else (F(1),)*DIMENSION
    result = {str(i): float(s)*blocks[c] for i, (c, s) in enumerate(zip(codes(), success, strict=True))}
    result['timeout'] = sum(((1-float(s))*blocks[c] for c, s in zip(codes(), success, strict=True)),
                            np.zeros((reference, reference), complex))
    padding = sorted(set(range(512))-set(codes()))
    result['export_abort'] = blocks[padding].sum(axis=0)
    return result


def reconstruct(diagonal, x, y):
    """Linear reconstruction; no PSD projection, truncation or renormalization."""
    keys = set(combinations(range(DIMENSION), 2))
    if (set(x) != keys or set(y) != keys or len(diagonal) != DIMENSION
            or not all(type(z) in (int, float) and math.isfinite(z)
                       for z in [*diagonal, *x.values(), *y.values()])):
        raise ValueError('all finite basis means and both quadratures required')
    result = np.diag(np.asarray(diagonal, dtype=complex))
    for a, b in keys:
        result[a, b] = (x[a, b]-1j*y[a, b])/2
        result[b, a] = result[a, b].conjugate()
    return result


def sampling_bound(samples, eta, service_error, *, dimension=DIMENSION):
    """Martingale/union Chebyshev bound for complete repeated trials.

    Each complete-tape score lies in [-1,1], with aborts scored zero.
    Given every previous record, its mean must differ by at most 2 delta
    from the same fixed target mean. Subtracting the conditional means
    gives martingale differences with conditional variance at most 1.
    Thus K/(N eta^2) bounds any mean's sampling failure. For the one basis
    setting the sum of all indicator variances is at most 1, so its joint
    contribution is 1/(N eta^2), rather than 66 separate unit-variance scores.
    A reference-uniform preparation/instrument bound supplies this premise
    even for a continuing carrier and correlated inputs. Independent copies
    are sufficient but unnecessary; arbitrary uncontrolled drift is excluded.
    """
    if (type(samples) is not int or samples <= 0 or type(dimension) is not int
            or dimension < 2 or any(type(v) is not F for v in (eta, service_error))
            or eta <= 0 or not 0 <= service_error <= 1):
        raise ValueError('positive samples/tolerance and a complete exact error required')
    count = 1+dimension*(dimension-1)
    error = eta+2*service_error
    frobenius = exact._sqrt_interval(F(dimension*(dimension+1), 2))[1]*error
    half_trace = exact._sqrt_interval(F(dimension))[1]*frobenius/2
    return {'simultaneous_failure_upper': min(F(1), F(count, samples)/eta**2),
            'operator_error_upper': frobenius, 'half_trace_error_upper': half_trace,
            'trials_used': samples*count}


def prepared_sampling_certificate(proof):
    """Use the actual bounded local preparation before every readout trial.

    Its complete ideal success block is q*sigma with q=5/33, independent of
    the incoming source and retaining its reference. All three preparation
    failures score zero. The target reconstructed matrix is therefore
    q*sigma, not a success-conditioned sample. Its modular flow is sigma's.
    Norm amplification by the raw readout blocks is charged before sampling.
    The preparation theorem is reference-uniform, hence supplies the same
    conditional-mean bound after arbitrary previously retained trial records.
    """
    preparation = local_preparation.budget()
    q = F(5, 33)
    gain = (1+F(proof['source_operator_error_upper']))**2
    error = (gain*preparation['complete_flagged_reference_halftrace_upper']
             +q*F(proof['complete_record_half_l1_upper']))
    sample = sampling_bound(10**14, F(1, 10000), error)
    return {'target': '(5/33)*local_preparation.reference()',
            'target_trace': str(q),
            'complete_trial_error_upper': str(error),
            'reference_uniform_conditioning_on_previous_records': True,
            'independent_initial_copies_required': False,
            'preparation_failures_discarded_from_sample_count': False,
            'modular_action_unchanged_by_target_scalar': True,
            'native_attempts_per_trial_upper': preparation['total_native_attempts_upper']
                +proof['native_attempts_per_trial'],
            'sampling': {k: str(v) if isinstance(v, F) else v for k,v in sample.items()},
            'physical_preparation_or_control_accuracy_certified': False}


def modular_bound(operator_error, estimate_floor, parameter):
    """Reference log stability, given a certified positive estimator floor.

    The floor is a premise to verify for the reconstructed Hermitian matrix,
    not a floating eigensolver diagnostic. Its normalization is immaterial to
    modular action. The bound follows from the resolvent identity and Duhamel.
    """
    if (any(type(v) is not F for v in (operator_error, estimate_floor, parameter))
            or not 0 <= operator_error < estimate_floor):
        raise ValueError('a certified eigenvalue floor exceeding the operator error required')
    return min(F(1), abs(parameter)*operator_error/(estimate_floor-operator_error))


def run():
    proof = certificate()
    # Deterministic complex state: synthetic diagnostic, not observed data.
    z = np.array([complex((i % 7)-3, (i % 5)-2) for i in range(DIMENSION)])
    z /= np.linalg.norm(z)
    rho = .7*np.outer(z, z.conj())+.3*np.eye(DIMENSION)/DIMENSION
    pairs = tuple(combinations(range(DIMENSION), 2))
    rebuilt = reconstruct(list(map(float, rho.diagonal().real)),
                          {p: float(2*rho[p].real) for p in pairs},
                          {p: float(-2*rho[p].imag) for p in pairs})
    residual = float(np.linalg.norm(rebuilt-rho, 'fro'))
    if residual > 1e-14:
        raise ArithmeticError('complex-state reconstruction failed')
    controls = []
    for a, b in ((0, 1), (7, 42), (64, 65)):
        factor = np.zeros((DIMENSION, 2, 2), complex)
        factor[a, 0, 0], factor[b, 1, 0] = .5, .5j
        factor[a, 1, 1], factor[b, 1, 1] = .5, .5j
        for kind in ('mixer', 'y_analyzer'):
            actual = records(factor, a, b, kind)
            target = records(factor, a, b, kind, raw=False)
            distance = sum(np.linalg.svd(actual[k]-target[k], compute_uv=False).sum()
                           for k in actual)/2
            if distance > float(F(proof['complete_record_half_l1_upper']))+1e-12:
                raise ArithmeticError('complete reference-output diagnostic exceeds bound')
            controls.append({'pair': [a, b], 'kind': kind,
                             'reference_dimension': 2, 'complete_error_diagnostic_pass': True})
    sample = sampling_bound(10**14, F(1, 10000), F(proof['complete_record_half_l1_upper']))
    return {'schema': 'oph.native_record_tomography.v1', 'all_pass': True,
            'certificate': proof, 'independent_reference_controls': controls,
            'reconstruction': {'exact_real_rank': DIMENSION**2,
                'without_phase_real_rank': DIMENSION*(DIMENSION+1)//2,
                'synthetic_complex_state_reconstructed': True,
                'native_only_diagonal_rank': DIMENSION,
                'reported_counts_are_laboratory_data': False},
            'conservative_sampling_example': {k: str(v) if isinstance(v, F) else v for k, v in sample.items()},
            'continuing_native_preparation': prepared_sampling_certificate(proof),
            'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', type=Path); group.add_argument('--check', type=Path)
    args = parser.parse_args(); result = run()
    if args.check:
        compiled.report.compare_receipt(strict_json.load(args.check), result)
        print('Complete native record tomography verified.')
    elif args.write:
        with args.write.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(result, indent=2, sort_keys=True)+'\n')
    else:
        print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
