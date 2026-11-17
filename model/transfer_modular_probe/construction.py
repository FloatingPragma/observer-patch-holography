"""Source-derived record-transfer spectrum, compiled on a declared coherent code.

The native Markov history and its reconstructed quantum representation are
different operational objects. This program executes the latter's modular
flow; it does not replace the former's dissipative repair by a unitary.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import permutations
import argparse
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from tools import strict_json
from ..common_source import exact_output as exact
from ..compiled_native_interface import construction as compiled
from ..native_clock import clock
from ..observer.scalar_response.experiment import canonical, compare_receipt
from ..observer.source_scalar import compiler
from ..primitive_observer import runtime, transfer


ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = tuple(sorted(set((*compiled.SOURCE_PATHS, *transfer.SOURCE_PATHS,
    'model/native_clock/clock.py', 'model/transfer_modular_probe/__init__.py',
    'model/transfer_modular_probe/__main__.py', 'model/transfer_modular_probe/construction.py',
    'Lean/Dynamics/PrivateInner.lean', 'Lean/SourceLaw/ReflectionPositivity.lean'))))
ERROR_GATE = F(1, 10**11)
SETTINGS = (0, 1, 4, 8)


def _setting(s):
    if type(s) is not int or s not in SETTINGS:
        raise ValueError('a declared integer modular parameter 0,1,4,8 is required')
    return s


def _mode(j):
    if type(j) is not int or j not in (1, 2):
        raise ValueError('one of the two nonconstant transfer modes is required')
    return j


@lru_cache(maxsize=1)
def reduced():
    return clock.lumping(transfer.native_counts())


def _det(matrix):
    return sum(((-1)**sum(p[i] > p[j] for i in range(3) for j in range(i+1, 3))
                *math.prod(matrix[i][p[i]] for i in range(3)) for p in permutations(range(3))), F())


def spectrum_audit():
    """Exact characteristic polynomial of the source's strongly lumped kernel."""
    p = reduced()['matrix']
    expected = ((F(14, 15), F(1, 15), F(0)), (F(1, 15), F(9, 10), F(1, 30)),
                (F(0), F(1, 6), F(5, 6)))
    if p != expected:
        raise ArithmeticError('the exact native record transfer changed')
    # Four evaluations identify the two monic cubics without floating roots.
    for x in map(F, (0, 1, 2, 3)):
        actual = _det(tuple(tuple((x if i == j else 0)-p[i][j] for j in range(3)) for i in range(3)))
        if actual != (x-1)*(x*x-F(5, 3)*x+F(311, 450)):
            raise ArithmeticError('native transfer characteristic polynomial mismatch')
    pi = reduced()['stationary']
    if any(pi[i]*p[i][j] != pi[j]*p[j][i] for i in range(3) for j in range(3)):
        raise ArithmeticError('the reduced transfer lost reversibility')
    return {'characteristic_polynomial': '(x-1)(x^2-5x/3+311/450)',
            'eigenvalues_in_code_order': ['1', '(25+sqrt(3))/30', '(25-sqrt(3))/30'],
            'reference_eigenvalues': ['3/8', '(25+sqrt(3))/80', '(25-sqrt(3))/80'],
            'stationary_class_weights': list(map(str, pi)), 'transfer_trace': '8/3',
            'basis': 'spectral basis of the symmetrized three-class transfer; fourth code address is padding'}


def _dyadic(x, bits=80):
    scale = 1 << bits
    return F((x*scale).numerator//(x*scale).denominator, scale)


def _log_interval(x):
    """Exact atanh-series enclosure for 0<x<=1; no binary log is a certificate."""
    if not F(1, 2) <= x <= 1:
        raise ValueError('log enclosure requires an exact argument in [1/2,1]')
    z = (1-x)/(1+x)
    value = -2*sum((z**(2*k+1)/F(2*k+1) for k in range(16)), F())
    tail = 2*z**33/(33*(1-z*z))
    return value-tail, value


@lru_cache(maxsize=1)
def energies():
    spectrum_audit()
    lo, hi = exact._sqrt_interval(F(3))
    lo, hi = _dyadic(lo), _dyadic(hi)+F(1, 1 << 80)
    result = [(F(0), F(0))]
    for a, b in (((25+lo)/30, (25+hi)/30), ((25-hi)/30, (25-lo)/30)):
        left, right = _log_interval(a)[0], _log_interval(b)[1]
        result.append((-right, -left))
    return tuple(result)


def _phase_box(interval):
    """Enclose exp(i*x) by a real-argument Taylor integral remainder."""
    lo, hi = interval
    x = _dyadic((lo+hi)/2)
    real = imag = F(0)
    for k in range(41):
        term = x**k/F(math.factorial(k))
        if k % 2:
            imag += (-1)**((k-1)//2)*term
        else:
            real += (-1)**(k//2)*term
    radius = abs(x)**41/F(math.factorial(41))+max(abs(x-lo), abs(x-hi))
    return (real-radius, real+radius), (imag-radius, imag+radius)


def phase_targets(s):
    s = _setting(s)
    return tuple(_phase_box((s*lo, s*hi)) for lo, hi in (*energies(), (F(0), F(0))))


def diagonal_word(s):
    s = _setting(s)
    phases = [s*float((lo+hi)/2) for lo, hi in energies()]+[0.]
    return compiler.compile_diagonal(phases)


def balanced_word(j):
    j = _mode(j)
    b = np.eye(4)
    b[0, 0] = b[0, j] = b[j, 0] = 1/math.sqrt(2)
    b[j, j] = -1/math.sqrt(2)
    return compiler.compile_real_orthogonal(b)


def _target_boxes(kind, *, s=0, j=1):
    zero = ((F(0), F(0)), (F(0), F(0)))
    result = [[zero for _ in range(4)] for _ in range(4)]
    for i in range(4):
        result[i][i] = ((F(1), F(1)), (F(0), F(0)))
    if kind == 'phase':
        for i, box in enumerate(phase_targets(s)):
            result[i][i] = box
    elif kind == 'balanced':
        j = _mode(j); lo, hi = exact._sqrt_interval(F(1, 2))
        for a, b in ((0, 0), (0, j), (j, 0)):
            result[a][b] = ((lo, hi), (F(0), F(0)))
        result[j][j] = ((-hi, -lo), (F(0), F(0)))
    elif kind == 'quadrature':
        j = _mode(j)
        result[j][j] = ((F(0), F(0)), (F(-1), F(-1)))
    else:
        raise ValueError('unknown complete target')
    return result


def certify(matrix, kind, *, s=0, j=1):
    if type(matrix) is not np.ndarray:
        raise ValueError('an original base ndarray is required')
    matrix = exact._array(matrix, (4, 4))
    boxes = _target_boxes(kind, s=s, j=j)
    squared = F(0)
    for a in range(4):
        for b in range(4):
            for value, (lo, hi) in zip(exact._component(matrix[a, b]), boxes[a][b], strict=True):
                squared += max(abs(value-lo), abs(value-hi))**2
    eps = exact._sqrt_interval(squared)[1]
    return eps


def _source(word, kind, **target):
    raw = compiler.execute(word, 2)
    eps = certify(raw, kind, **target)
    if eps > ERROR_GATE:
        raise ArithmeticError('raw compiled map exceeds its complete-target error envelope')
    resources = compiler.resources(word)
    additive = ('native_pulses', 'proper_code_service_events', 'proper_code_load_events',
                'proper_code_export_events', 'single_code_pulses', 'pair_code_pulses',
                'identity_phase_pulses', 'zero_angle_requests')
    return raw, tuple(sorted((key, resources[key]) for key in additive))


def _unitary(raw, resources, name):
    return runtime.Instrument(kraus=((raw,),), resources=resources, binding={'operation': name})


def program(j, intervals=(4,), *, quadrature='x'):
    """Continue one state through all reads; never reset after measurement."""
    j = _mode(j)
    if (type(intervals) is not tuple or not 1 <= len(intervals) <= 2
            or quadrature not in ('x', 'y') or (len(intervals) == 2 and quadrature != 'x')):
        raise ValueError('one x/y read or two continuing x reads required')
    for s in intervals:
        _setting(s)
    bword = balanced_word(j)
    b, bc = _source(bword, 'balanced', j=j)
    inverse, ic = _source(compiler.inverse(bword), 'balanced', j=j)
    events = [runtime.Event('coherent_mode_preparation', instruments=(_unitary(b, bc, 'balanced'),))]
    for i, s in enumerate(intervals):
        u, uc = _source(diagonal_word(s), 'phase', s=s)
        events.append(runtime.Event(f'modular_segment_{i}', instruments=(_unitary(u, uc, 'transfer_modular'),)))
        if quadrature == 'y':
            phases = [0.]*4; phases[j] = -math.pi/2
            qword = compiler.compile_diagonal(phases)
            q, qc = _source(qword, 'quadrature', j=j)
            # Certify the inverse by conjugating its raw map back to the same target.
            qi = compiler.execute(compiler.inverse(qword), 2)
            if certify(qi.conj().T, 'quadrature', j=j) > ERROR_GATE:
                raise ArithmeticError('inverse quadrature source error')
            events.append(runtime.Event(f'quadrature_phase_{i}', instruments=(_unitary(q, qc, 'quadrature'),)))
        events.append(runtime.Event(f'analysis_basis_{i}', instruments=(_unitary(b, bc, 'balanced'),)))
        projectors = tuple(np.diag([int(k == a) for k in range(4)]) for a in range(4))
        read = runtime.Instrument(kraus=tuple((p,) for p in projectors), visible=True,
            resources=(('four_outcome_reads', 1),),
            binding={'pointer': 'complete computational read; basis and undo are separate events',
                     'quadrature': quadrature, 'all_padding_outcomes_retained': True})
        events.append(runtime.Event(f'retained_{quadrature}_read_{i}', instruments=(read,)))
        events.append(runtime.Event(f'undo_analysis_basis_{i}', instruments=(_unitary(inverse, ic, 'inverse_balanced'),)))
        if quadrature == 'y':
            events.append(runtime.Event(f'undo_quadrature_phase_{i}', instruments=(_unitary(qi, qc, 'inverse_quadrature'),)))
    initial = np.zeros((4, 4), complex); initial[0, 0] = 1
    return runtime.Program(f'transfer_modular_mode_{j}', 4, initial, tuple(events),
        binding={'modular_intervals': intervals, 'mode': j,
                 'algebra': 'reconstructed three-mode code plus explicit padding',
                 'native_dissipative_history_is_not_replaced': True})


def _sequence_allowance(gates):
    delta = ERROR_GATE*(2+ERROR_GATE)/2
    gain = (1+ERROR_GATE)**2
    return delta*sum((gain**k for k in range(gates)), F())


def factor_execute(p, *, stop_after=None, checkpoint=None):
    """Evaluate this rank-one Kraus protocol without cancellation in tiny outputs.

    Uses the validated shared Program, original raw arrays and original |0>
    input. Every Kraus group has exactly one member here. No state is projected,
    normalized or discarded. Checkpoints are scoped to this pure-state protocol.
    """
    p = runtime.validate_program(p)
    initial = np.zeros((4, 4), complex); initial[0, 0] = 1
    if p.work_dimension != 4 or not np.array_equal(p.initial_density, initial):
        raise ValueError('the declared four-code blank preparation is required')
    if any(e.control is not None or e.sequence or e.repeat != 1 or len(e.instruments) != 1
           or e.instruments[0].schur is not None
           or any(len(g) != 1 for g in e.instruments[0].kraus) for e in p.events):
        raise ValueError('only this uncontrolled single-Kraus-per-outcome protocol is supported')
    end = len(p.events) if stop_after is None else stop_after
    if type(end) is not int or not 0 <= end <= len(p.events):
        raise ValueError('valid event cursor required')
    cursor = 0
    branches = [((), np.array([1, 0, 0, 0], complex), p.preparation_resources)]
    if checkpoint is not None:
        cursor = checkpoint['cursor']
        if (checkpoint['program'] != p.digest or type(cursor) is not int or not 0 <= cursor <= end):
            raise ValueError('checkpoint must belong to this program and cursor range')
        # These small finite checkpoints are checked against a fresh prefix,
        # so an altered record, amplitude or resource cannot become a new input.
        prefix = factor_execute(p, stop_after=cursor)
        if checkpoint != prefix.checkpoint:
            raise ValueError('checkpoint differs from its complete source prefix')
        branches = [(tuple(r), np.array([complex(a, b) for a, b in v]), tuple(map(tuple, cost)))
                    for r, v, cost in checkpoint['branches']]
    for event in p.events[cursor:end]:
        instrument = event.instruments[0]
        branches = [(record+((k,) if instrument.visible else ()), group[0]@vector,
                     runtime.add_resources(cost, instrument.resources))
                    for record, vector, cost in branches for k, group in enumerate(instrument.kraus)]
    cp = {'program': p.digest, 'cursor': end,
          'branches': [[list(r), [[float(z.real), float(z.imag)] for z in v], [list(x) for x in cost]]
                       for r, v, cost in branches]}
    output = tuple(runtime.Branch(r, np.outer(v, v.conj()), cost) for r, v, cost in branches)
    return SimpleNamespace(branches=output, checkpoint=cp,
                           resources=dict(output[0].resources))


def correlation_join(rows):
    """Attach the same two frequencies to actual native current-bit moments.

    The marked transfer insertions are not treated as detector effects. The
    new X/Y probes measure spectral phases, whose weighted sum is the analytic
    continuation of this particular connected Euclidean correlation.
    """
    kernel = transfer.native_counts()
    weights = (1/7260+math.sqrt(3)/16335, 1/7260-math.sqrt(3)/16335)
    e = tuple(float(sum(x)/2) for x in energies()[1:])
    euclidean = []
    for n in (0, 1, 4, 8):
        law = transfer.exact_record_law(kernel, n)
        covariance = law[1, 1][0]-F(25, 1089)
        predicted = sum(w*math.exp(-n*v) for w, v in zip(weights, e, strict=True))
        if abs(predicted-float(covariance)) > 1e-15:
            raise ArithmeticError('spectral transfer modes fail the exact native record correlation')
        euclidean.append({'unread_native_attempts': n, 'exact_current_bit_covariance': str(covariance),
                          'spectral_value_diagnostic': predicted})
    modular = []
    for s in SETTINGS:
        read = {(r['mode'], r['quadrature']): 2*r['record_weights'][0]-1
                for r in rows if r['modular_parameter'] == s}
        observed = sum(w*complex(read[j, 'x'], read[j, 'y']) for j, w in enumerate(weights, 1))
        predicted = sum(w*complex(math.cos(s*v), math.sin(s*v)) for w, v in zip(weights, e, strict=True))
        if abs(observed-predicted) > 1e-14:
            raise ArithmeticError('source spectral probe fails the continued correlation')
        modular.append({'modular_parameter': s, 'source_phase_weighted_sum': [observed.real, observed.imag],
                        'analytic_value_diagnostic': [predicted.real, predicted.imag]})
    return {'weights_in_mode_order': ['1/7260+sqrt(3)/16335', '1/7260-sqrt(3)/16335'],
            'euclidean_native_record_checks': euclidean, 'modular_source_probe_checks': modular,
            'native_record_insertions_are_not_POVM_effects': True,
            'analytic_continuation_does_not_identify_native_waiting_with_modular_time': True}


def run():
    audit = spectrum_audit()
    rows = []
    for j in (1, 2):
        energy = float(sum(energies()[j])/2)
        for s in SETTINGS:
            for quadrature in ('x', 'y'):
                p = program(j, (s,), quadrature=quadrature)
                h = factor_execute(p)
                target = (1+(math.cos(s*energy) if quadrature == 'x' else math.sin(s*energy)))/2
                expected = [target if k == 0 else 1-target if k == j else 0. for k in range(4)]
                residual = max(abs(b.probability-v) for b, v in zip(h.branches, expected, strict=True))
                if residual > 1e-11:
                    raise ArithmeticError('complete modular source readout disagrees with its spectrum')
                rows.append({'mode': j, 'modular_parameter': s, 'quadrature': quadrature,
                    'complete_records': [list(b.record) for b in h.branches],
                    'record_weights': [b.probability for b in h.branches],
                    'ideal_probability_diagnostic': expected,
                    'probability_residual_diagnostic': residual,
                    'gate_comparison_half_diamond_upper': '1/10000000000',
                    'resources_per_branch': h.resources})
    continuing = []
    for j in (1, 2):
        p = program(j, (4, 1)); h = factor_execute(p)
        prefix = factor_execute(p, stop_after=5)
        resumed = factor_execute(p, checkpoint=strict_json.loads(json.dumps(prefix.checkpoint)))
        q = [(1+math.cos(s*float(sum(energies()[j])/2)))/2 for s in (4, 1)]
        expected = {(0, 0): q[0]*q[1], (0, j): q[0]*(1-q[1]),
                    (j, 0): (1-q[0])*(1-q[1]), (j, j): (1-q[0])*q[1]}
        residual = max(abs(b.probability-expected.get(b.record, 0.)) for b in h.branches)
        checkpoint = sum(float(np.linalg.norm(a.density-b.density, 'fro'))
                         for a, b in zip(h.branches, resumed.branches, strict=True))
        if max(residual, checkpoint) > 1e-11:
            raise ArithmeticError('continuing modular record or checkpoint mismatch')
        continuing.append({'mode': j, 'modular_intervals': [4, 1],
            'complete_records': [list(b.record) for b in h.branches],
            'record_weights': [b.probability for b in h.branches],
            'probability_residual_diagnostic': residual, 'checkpoint_distance_diagnostic': checkpoint,
            'gate_comparison_half_diamond_upper': '1/10000000000',
            'resources_per_branch': h.resources})
    if _sequence_allowance(7) > F(1, 10**10):
        raise ArithmeticError('portable complete sequence allowance is insufficient')
    return {'schema': 'oph.transfer_modular_probe.v1', 'spectrum': audit,
        'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS},
        'energy_intervals': [[str(x) for x in interval] for interval in energies()],
        'single_readouts': rows, 'continuing_histories': continuing,
        'record_correlation_join': correlation_join(rows),
        'scope': {'coherent_source_services_supplied': True, 'new_empirical_data': False,
                  'reference_is_reconstructed_transfer_not_native_Palm': True,
                  'complete_effective_code_modular_flow': True,
                  'physical_service_error_bounded': False,
                  'compiled_gate_certificates_include_all_four_columns_and_arbitrary_references': True,
                  'runtime_state_arithmetic_is_diagnostic': True,
                  'native_waiting_equals_modular_evolution': False}}


def main():
    parser = argparse.ArgumentParser(__doc__); parser.add_argument('--check', type=Path)
    args = parser.parse_args(); result = run()
    if args.check:
        compare_receipt(result, strict_json.loads(args.check.read_text()))
        print('Transfer-derived modular probes replayed with all outcomes and continuing states.')
    else:
        print(canonical(result), end='')
