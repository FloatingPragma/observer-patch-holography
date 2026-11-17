"""Exact quantum-memory limits for persistently exploring native observers.

Before each fresh fair coin, a classical controller chooses the law and seam
from previous committed records. Every seam has probability at least eta/30.
No quantum operation intervenes between native events. All histories, including
zero raw outcomes, are retained; t means t actual events, not padded halts.

Equal records translate both configurations by the same current, preserving
z=x-y. Their one-edge agreement probabilities are
  C: 1-min(abs(z_a-z_b),2)/2; E: 1[z_a=z_b].
The existing all-pair census gives uniform-seam maxima 5/6 and 7/10. Splitting
any admitted selector into eta times the uniform law and its remainder gives
r_C=1-eta/6, r_E=1-3*eta/10. The C bound also covers history/seam-dependent
C/E switching. Summing surviving common-history weights inductively gives
H_xy(t)<=r**t, without an iid assumption about the controller.

Each complete history is still sqrt(w_h) U_h P_h. Consequently, for every
input and untouched reference, the trace distance between Gamma_t(rho) and
Gamma_t(Delta(rho)) is at most (D-1)*r**t/2. Positivity gives
||rho_xy||_1<=sqrt(p_x*p_y), and sum_{x!=y}sqrt(p_x*p_y)<=D-1.
The comparison channel is entanglement breaking. Any final decoder of a
d-dimensional encoded Bell input therefore has error at least
1-1/d-(D-1)*r**t/2. This does not assert finite-time entanglement breaking,
loss of classical history memory, a physical clock, or a defect energy.

The executable finite-memory DP below contracts exact history-pair weights,
not amplitudes across distinct committed records. It does not execute new
source pulses. Existing source instruments supply the complete raw grammar.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import math
from pathlib import Path

import numpy as np

from . import record_recovery as recovery
from .source_discrimination import circuit


ROOT = Path(__file__).resolve().parents[2]


def _fraction(value, *, positive=False):
    if (type(value) is not F or type(value.numerator) is not int
            or type(value.denominator) is not int or value.denominator <= 0
            or math.gcd(value.numerator, value.denominator) != 1
            or not 0 <= value <= 1 or positive and value == 0):
        raise ValueError('an original canonical Fraction in the declared probability range required')
    # Fraction is mathematically immutable but its private slots remain
    # writable. Never expose a default or retain caller-owned rational slots.
    return F(value.numerator, value.denominator)


def _binary(value):
    if type(value) is not bool:
        raise ValueError('an original Boolean binary-domain switch required')
    return value


def actions():
    """Fixed complete action menu; the law choice itself is archived."""
    return tuple((law, edge) for law in ('C', 'E') for edge in recovery.seams())


@dataclass(frozen=True)
class RawBranch:
    action: int
    coin: int
    raw: int
    inputs: tuple[int, ...]
    outputs: tuple[int, ...]


@lru_cache(maxsize=2)
def _branches(binary):
    states = recovery.carrier(binary)
    index = {state: i for i, state in enumerate(states)}
    rows = []
    for action, (law, edge) in enumerate(actions()):
        for coin in (0, 1):
            groups = [[] for _ in range(1 << circuit.RECORD_BITS[law])]
            decoded = circuit.record_values(law, coin)
            inverse = {value: raw for raw, value in decoded.items()}
            for i, state in enumerate(states):
                out, record = recovery.recorded_step(law, state, edge, coin)
                groups[inverse[record[-1]]].append((i, index[out]))
            for raw, pairs in enumerate(groups):
                inputs, outputs = tuple(zip(*pairs)) if pairs else ((), ())
                if len(set(outputs)) != len(outputs):
                    raise ArithmeticError('native recorded branch ceased to be a partial injection')
                rows.append(RawBranch(action, coin, raw, inputs, outputs))
    return tuple(rows)


def branch_maps(*, binary=False):
    """Every raw branch, including aborts, with amplitude 1/sqrt(2)."""
    # Even frozen dataclasses can be changed through __dict__/object.__setattr__.
    # Public rows must not alias the owned cached source maps used next time.
    return tuple(RawBranch(row.action, row.coin, row.raw, row.inputs, row.outputs)
                 for row in _branches(_binary(binary)))


@lru_cache(maxsize=1)
def _certificate():
    states, edges = recovery.carrier(), recovery.seams()
    counts = np.zeros((60, 78, 78), dtype=np.int64)
    rows = branch_maps()
    for row in rows:
        counts[row.action][np.ix_(row.inputs, row.inputs)] += 1
        for i, out_i in zip(row.inputs, row.outputs, strict=True):
            for j, out_j in zip(row.inputs, row.outputs, strict=True):
                if tuple(a-b for a, b in zip(states[i], states[j], strict=True)) != tuple(
                        a-b for a, b in zip(states[out_i], states[out_j], strict=True)):
                    raise ArithmeticError('equal committed records changed signed difference')
    maxima = {}
    for law_index, law in enumerate(('C', 'E')):
        total = np.zeros((78, 78), dtype=np.int64)
        for e, (a, b) in enumerate(edges):
            actual = counts[30*law_index+e]
            for i, x in enumerate(states):
                for j, y in enumerate(states):
                    delta = (x[a]-y[a])-(x[b]-y[b])
                    expected = 2-min(abs(delta), 2) if law == 'C' else 2*int(delta == 0)
                    if actual[i, j] != expected:
                        raise ArithmeticError('complete native branches violate the agreement formula')
            total += actual
        maximum = max(int(total[i, j]) for i in range(78) for j in range(78) if i != j)
        if maximum != (50 if law == 'C' else 42):
            raise ArithmeticError('native all-pair agreement maximum changed')
        if not np.array_equal(np.diag(total), np.full(78, 60)):
            raise ArithmeticError('complete raw branch probabilities fail normalization')
        maxima[law] = str(F(maximum, 60))
    return {'states': 78, 'ordered_pairs': 78**2, 'seams': 30,
            'law_pair_seam_checks': 2*30*78**2,
            'raw_branches': len(rows), 'zero_raw_branches': sum(not row.inputs for row in rows),
            'primitive_input_addresses': sum(len(row.inputs) for row in rows),
            'uniform_offdiagonal_agreement_maxima': maxima,
            'complete_branch_normalization': 'sum_coin,raw K* K=I for every law and seam',
            'equal_record_difference': 'output_x-output_y=input_x-input_y'}


def exact_certificate():
    result = _certificate()
    return {**result, 'uniform_offdiagonal_agreement_maxima':
            dict(result['uniform_offdiagonal_agreement_maxima'])}


def bounds(eta, attempts, *, law_scope='switching', binary=False, logical_dimension=2):
    """Exact bounds for continuing histories with at least `attempts` events.

    The theorem concerns all history-dependent policies satisfying the stated
    lower bound, not only the finite-memory examples. Scalar rates are sharp:
    put the selector's remainder on an edge invisible to a maximizing pair.
    """
    eta = _fraction(eta, positive=True)
    attempts = recovery.attempts(attempts)
    size = len(recovery.carrier(_binary(binary)))
    if type(law_scope) is not str or law_scope not in ('C', 'E', 'switching'):
        raise ValueError('law scope C, E or switching required')
    if type(logical_dimension) is not int or not 2 <= logical_dimension <= size:
        raise ValueError('an original encoded dimension between 2 and the carrier size required')
    rate = 1-eta*(F(3, 10) if law_scope == 'E' else F(1, 6))
    survival = rate**attempts
    distance = min(F(1), F(size-1, 2)*survival)
    return {'eta': eta, 'attempts': attempts, 'dimension': size, 'law_scope': law_scope,
            'logical_dimension': logical_dimension, 'one_step_agreement_upper': rate,
            'pair_agreement_upper': survival,
            'complete_reference_trace_distance_upper': distance,
            'full_input_optimal_entanglement_fidelity_upper': F(1, size)+F(size-1, size)*survival,
            'encoded_final_recovery_trace_error_lower': max(F(0), 1-F(1, logical_dimension)-distance)}


@dataclass(frozen=True)
class Controller:
    """An explicit finite classical controller, independent of the input.

    probabilities[m][action] chooses a law/seam BEFORE the fresh fair coin.
    updates[m][action][coin][raw] includes every unused raw address too.
    Initial memory is fixed; the complete record archive is never overwritten.
    """
    eta: F
    probabilities: tuple
    updates: tuple
    initial_memory: int = 0


def validate_controller(controller):
    if type(controller) is not Controller:
        raise ValueError('an explicit finite Controller required')
    eta = _fraction(controller.eta, positive=True)
    probabilities, updates = controller.probabilities, controller.updates
    if type(probabilities) is not tuple or not probabilities or type(updates) is not tuple:
        raise ValueError('immutable complete controller tables required')
    size = len(probabilities)
    if (len(updates) != size or type(controller.initial_memory) is not int
            or not 0 <= controller.initial_memory < size):
        raise ValueError('controller memory address outside its declared table')
    detached_probabilities = []
    for m, row in enumerate(probabilities):
        if type(row) is not tuple or len(row) != 60:
            raise ValueError('the complete law/seam probability menu is required')
        detached_probabilities.append(tuple(_fraction(probability) for probability in row))
        if sum(row, F(0)) != 1 or any(row[e]+row[e+30] < eta/30 for e in range(30)):
            raise ValueError('selector violates normalization or persistent exploration')
        if type(updates[m]) is not tuple or len(updates[m]) != 60:
            raise ValueError('complete record-memory updates required')
        for a, (law, _) in enumerate(actions()):
            coins = updates[m][a]
            if type(coins) is not tuple or len(coins) != 2:
                raise ValueError('both fresh-coin update menus required')
            for records in coins:
                if type(records) is not tuple or len(records) != 1 << circuit.RECORD_BITS[law]:
                    raise ValueError('every raw record including aborts is required')
                if any(type(value) is not int or not 0 <= value < size for value in records):
                    raise ValueError('record-memory update leaves the declared memory')
    return Controller(eta, tuple(detached_probabilities), updates, controller.initial_memory)


def example_controller(eta=F(1, 2)):
    """Two-state record-dependent law switching with different preferred seams."""
    eta = _fraction(eta, positive=True)
    probabilities, updates = [], []
    for m in (0, 1):
        row = [F(0)]*60
        for e in range(30):
            row[30*m+e] = eta/30 + (1-eta if e == m else 0)
        probabilities.append(tuple(row))
        update = []
        for law, _ in actions():
            update.append(tuple(tuple(
                (m+coin+int(circuit.record_values(law, coin).get(raw, 0) != 0)) % 2
                for raw in range(1 << circuit.RECORD_BITS[law])) for coin in (0, 1)))
        updates.append(tuple(update))
    return validate_controller(Controller(eta, tuple(probabilities), tuple(updates)))


def finite_controller_kernel(controller, attempts, *, binary=False):
    """Exact H_xy=sum_{histories compatible with x,y} probability(history).

    Backward DP merges only sums with identical memory and current input pair.
    No classical histories are coherently added or normalized conditionally.
    Integer numerators over one common denominator avoid floating underflow.
    """
    controller = validate_controller(controller)
    attempts = recovery.attempts(attempts)
    rows = branch_maps(binary=_binary(binary))
    dimension = len(recovery.carrier(binary))
    memories = len(controller.probabilities)
    denominator = 2*math.lcm(*(p.denominator for row in controller.probabilities for p in row))
    weights = tuple(tuple(int(p*denominator/2) for p in row) for row in controller.probabilities)
    previous = np.ones((memories, dimension, dimension), dtype=object)
    indexed = tuple((row, np.ix_(row.inputs, row.inputs), np.ix_(row.outputs, row.outputs))
                    for row in rows if row.inputs)
    for _ in range(attempts):
        following = np.zeros_like(previous)
        for m in range(memories):
            for row, source, target in indexed:
                weight = weights[m][row.action]
                if weight:
                    next_memory = controller.updates[m][row.action][row.coin][row.raw]
                    following[m][source] += weight*previous[next_memory][target]
        previous = following
    divisor = denominator**attempts
    result = tuple(tuple(F(int(value), divisor) for value in row)
                   for row in previous[controller.initial_memory])
    if any(result[i][i] != 1 for i in range(dimension)):
        raise ArithmeticError('complete controller output probabilities changed')
    e_only = all(not any(row[:30]) for row in controller.probabilities)
    upper = bounds(controller.eta, attempts, law_scope='E' if e_only else 'switching',
                   binary=binary)['pair_agreement_upper']
    if any(not 0 <= result[i][j] <= (1 if i == j else upper)
           for i in range(dimension) for j in range(dimension)):
        raise ArithmeticError('adaptive common-history kernel violates the exact theorem')
    return result


def source_bindings():
    from .record_recovery_receipt import SOURCE_PATHS
    paths = set(SOURCE_PATHS) | {'model/observer/adaptive_record_recovery.py',
                               'model/observer/test_adaptive_record_recovery.py'}
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}


def _json_safe(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {key: _json_safe(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_safe(item) for item in value]
    return value


def report():
    controller = example_controller()
    kernel = finite_controller_kernel(controller, 3)
    counts = Counter(value for i, row in enumerate(kernel) for j, value in enumerate(row) if i != j)
    return _json_safe({'schema': 'oph.observer.adaptive_record_recovery.v1',
        'native_certificate': exact_certificate(),
        'grammar': {'selector': 'p(seam|complete past)>=eta/30; law may depend on past and selected seam',
                    'coin': 'fresh fair independent coin after law/seam selection',
                    'intermediate_work_operations': 'native C/E only; classical control and retention',
                    'clock': 'actual native event count; no halted-branch padding',
                    'reference': 'arbitrary untouched reference; no input commutation assumption'},
        'bounds_at_eta_half': {law: bounds(F(1, 2), t, law_scope=law)
                              for law, t in (('C', 58), ('E', 31))},
        'executed_exact_controller': {'attempts': 3, 'memory_states': 2,
            'probabilities': controller.probabilities, 'updates': controller.updates,
            'initial_memory': controller.initial_memory, 'eta': controller.eta,
            'offdiagonal_kernel_multiplicities': {str(k): v for k, v in sorted(counts.items())},
            'optimal_full_input_entanglement_fidelity': sum(map(sum, kernel), F(0))/78**2,
            'all_input_probability_sums': '1 on every diagonal',
            'archive_bit_upper_per_event': 10,
            'archive_scope': 'law1+seam5+coin1+raw3; all zero/abort addresses retained',
            'source_pulse_execution': False},
        'scope': {'finite_horizon_entanglement_breaking_claimed': False,
                  'classical_history_memory_excluded': False,
                  'active_quantum_recovery_excluded_by_theorem': False,
                  'native_energy_or_physical_duration_derived': False,
                  'whole_program_or_refinement_no_go': False}})
