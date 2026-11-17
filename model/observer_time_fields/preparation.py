"""Native labelled tracial preparation and complete stopped field histories.

The target is selected by coarse native auxiliary coins with bounded rejection.
A full78 current-adaptive reset prepares that target; only its two local bits
are disclosed. The same auxiliary supplies every later seam selector. Its
coarse effect I/2 permits arbitrary auxiliary/work/reference correlations.

For arbitrary reference R, successful target t has the exact block
  (1-tau)/66 |t><t| tensor sum_x s[t,x] rho_R[x,x].
Every ideal-minus-actual success block is positive. The complete half-diamond
distance from labelled tracial replacement is exactly
  max_x (1-(1-tau)/66 sum_t s[t,x]).
Preparation failures retain their flag and R, explicitly discarding work at
that comparison boundary. Private replay retains actual terminal carriers.
Neither successful targets nor their reference are renormalized. The public
local bits are intentionally correlated with the prepared work.

The continuation implements the existing incident-attempt field experiment,
including first-bit-controlled native feedback. Remote repairs continue;
zero-current incident attempts count; selector/wait timeouts are outputs.
Only the existing local-field cylinder dictionary joins the effective model:
no hidden-terminal-state equality, federation, or physical-time claim.

Every finite successful continuation is a CP subinstrument of its completed
native continuation. Each scheduled tick completes with state-independent
probability z=p(1-tau_s)(1-r**H)/(1-r), r=(1-p)(1-tau_s), p=3/10.
Thus the exact complete native-continuation half-diamond error is
max_x[1-(1-f_x)z**reads], also with the declared feedback. Its field-record TV
bound transfers through the existing effective dictionary; no conditioning
or per-read multiplication of the preparation error is used.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
import hashlib
from itertools import product
import json
import math
from pathlib import Path
import random
from typing import NamedTuple

import numpy as np
from tools import strict_json

from ..native_clock import configuration as reset
from ..observer import native
from ..repair_scheduler import coins
from . import core
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('preparation_results.json')
SOURCE_PATHS = tuple(sorted(set((*reset.SOURCE_PATHS, *core.SOURCE_PATHS,
    'model/repair_scheduler/coins.py', 'model/observer_time_fields/preparation.py'))))
SEED_SCANS, RETRIES, TARGET_ROUNDS, SELECTOR_ROUNDS, CAP = 4, 32, 32, 8, 64
SCHEMA = 'oph.observer_time_fields.native_labelled_preparation.v1'
SCOPE = {
    'carrier': 'elementary twelve-port total-two native carrier; full78 input, binary66 success',
    'quantum_domain': 'arbitrary joint auxiliary/work/reference input; auxiliary discarded',
    'public_preparation_output': 'two actual initial local bits or an explicit failure flag',
    'private_boundary': 'discard target, resolved original label, all auxiliary outputs and fine tapes',
    'preparation_failure_work': 'discarded at comparison interface; retained only by private replay',
    'continuation_output': 'local field word, incident read clocks and all stopping flags',
    'comparison': 'complete unconditioned half-diamond preparation; complete field-record TV afterwards',
    'remaining_services': 'two native carriers, joint seam-current access, adaptive seam control and classical memory',
    'initial_bit_or_stationary_preparation_service_required': False,
    'auxiliary_independence_or_reset_required': False,
    'physical_time_or_full_federation_claimed': False,
    'thermodynamic_erasure_of_private_environment_claimed': False,
}


def _parameters(seed_scans, retries, target_rounds, selector_rounds, cap, read_ports, feedback):
    if any(type(v) is not int or not 1 <= v <= maximum for v, maximum in
           ((seed_scans, 128), (retries, 128), (target_rounds, 32),
            (selector_rounds, 8), (cap, 64))):
        raise ValueError('original bounded integer controller parameters required')
    if (type(read_ports) is not tuple or not 1 <= len(read_ports) <= 4
            or any(type(v) is not int or v not in (0, 1) for v in read_ports)
            or type(feedback) is not bool):
        raise ValueError('one through four local reads and an original Boolean feedback flag required')


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                     allow_nan=False).encode()).hexdigest()


def _rational_hash(values):
    digest = hashlib.sha256()
    for row in values:
        digest.update((' '.join(map(str, row))+'\n').encode())
    return digest.hexdigest()


def _up(value, bits=80):
    value, scale = F(value), 1 << bits
    return F(-((-value.numerator*scale)//value.denominator), scale)


def _down(value, bits=80):
    value, scale = F(value), 1 << bits
    return F(value.numerator*scale//value.denominator, scale)


class PreparationTable(NamedTuple):
    targets: tuple
    inputs: tuple
    probabilities: tuple  # rows x, columns t, including target-selector success
    selector_timeout: F
    source_digest: str

    @property
    def failure_by_input(self):
        return tuple(1-sum(row, F()) for row in self.probabilities)


def preparation_table(seed_scans=SEED_SCANS, retries=RETRIES, target_rounds=TARGET_ROUNDS):
    _parameters(seed_scans, retries, target_rounds, 1, 1, (0,), False)
    table = _preparation_table(seed_scans, retries, target_rounds)
    # Fraction has mutable protected slots. A read-only outer tuple alone
    # would still expose cache-owned numbers to callers.
    return PreparationTable(table.targets, table.inputs,
        tuple(tuple(F(p.numerator, p.denominator) for p in row) for row in table.probabilities),
        F(table.selector_timeout.numerator, table.selector_timeout.denominator), table.source_digest)


@lru_cache(maxsize=12)
def _preparation_table(seed_scans, retries, target_rounds):
    targets, inputs = tuple(native.states(binary=True)), tuple(native.states())
    selector = coins.selector_certificate(66, target_rounds)
    alpha_per_target = F(selector['accepted_probability_per_address'])
    columns = tuple(reset.success_probabilities(seed_scans, retries, full=True,
                     mode='reset', target=target) for target in targets)
    rows = tuple(tuple(alpha_per_target*column[x] for column in columns)
                 for x in range(len(inputs)))
    tau = F(selector['timeout_probability'])
    bound = reset.budget(seed_scans, retries, full=True, mode='reset', target=targets[0])['failure_upper']
    if (len(rows) != 78 or len(columns) != 66
            or any(not 0 <= p <= alpha_per_target for row in rows for p in row)
            or any(not tau <= 1-sum(row, F()) <= tau+(1-tau)*bound for row in rows)):
        raise ArithmeticError('labelled native reset lost its complete all-input domination bound')
    return PreparationTable(targets, inputs, rows, tau, _rational_hash(rows))


def preparation_certificate(seed_scans=SEED_SCANS, retries=RETRIES, target_rounds=TARGET_ROUNDS):
    table = preparation_table(seed_scans, retries, target_rounds)
    worst = max(table.failure_by_input)
    bound = reset.budget(seed_scans, retries, full=True, mode='reset', target=table.targets[0])
    public_marginal = {''.join(map(str, label)): str(F(sum(tuple(t[p] for p in (0, 1)) == label
                                    for t in table.targets), 66))
                       for label in product((0, 1), repeat=2)}
    return {
        'input_columns': 78, 'native_targets': 66,
        'source_success_table_sha256': table.source_digest,
        'target_selector': coins.selector_certificate(66, target_rounds),
        'reset_native_attempts_max': bound['native_attempts_max'],
        'reset_failure_upper': str(bound['failure_upper']),
        'complete_halfdiamond_upper': str(table.selector_timeout+(1-table.selector_timeout)*bound['failure_upper']),
        'exact_complete_halfdiamond_upper_80bit_enclosure': str(_up(worst)),
        'exact_complete_halfdiamond_lower_80bit_enclosure': str(_down(worst)),
        'exact_complete_halfdiamond_sha256': _rational_hash(((worst,),)),
        'ideal_initial_label_probabilities': public_marginal,
        'successful_targets_exactly_uniform_for_every_input': False,
        'successful_reference_claimed_unchanged': False,
        'retained_local_label_independent_of_work': False,
        'success_fine_word': 'scalar |target><resolved original configuration|',
        'comparison': 'positive success deficits plus disjoint failure flags; no success renormalization',
    }


def apply_preparation(density, seed_scans=SEED_SCANS, retries=RETRIES, target_rounds=TARGET_ROUNDS):
    """Numerical execution of the proved complete reduced quantum instrument.

    Success has one work-basis block per target; failures retain R and discard
    work. Floating arrays are diagnostics. Exact probabilities above supply
    the reference-uniform proof, without a floating positivity premise.
    """
    table = preparation_table(seed_scans, retries, target_rounds)
    if (type(density) is not np.ndarray or density.dtype.kind not in 'fc'
            or density.ndim != 2 or density.shape[0] != density.shape[1]
            or density.shape[0] not in (78, 156, 234) or not np.isfinite(density).all()):
        raise ValueError('a finite full78 operator with one through three reference levels required')
    r = density.shape[0]//78
    blocks = density.reshape(78, r, 78, r)
    diagonal = np.stack([blocks[x, :, x, :] for x in range(78)])
    marginal = diagonal.sum(axis=0)
    weights = np.array([[float(p) for p in row] for row in table.probabilities])
    success = np.einsum('xt,xij->tij', weights, diagonal)
    selector_failure = float(table.selector_timeout)*marginal
    reset_failure = sum((float(f-table.selector_timeout)*diagonal[x]
                         for x, f in enumerate(table.failure_by_input)), np.zeros((r, r), complex))
    return {'success': success, 'target_timeout': selector_failure,
            'reset_timeout': reset_failure, 'reference': marginal}


class StopKernel(NamedTuple):
    success: tuple
    selector_timeout: tuple
    incident_timeout: tuple
    denominator: int
    source_digest: str


def _left(matrix, values):
    """Integer sparse-left product; no floating probability or fitted kernel."""
    width = len(values[0])
    sparse = tuple(tuple((column, coefficient) for column, coefficient in enumerate(row)
                         if coefficient) for row in matrix)
    return tuple(tuple(sum(coefficient*values[column][j]
                           for column, coefficient in row)
                       for j in range(width)) for row in sparse)


def stopped_kernel(cap=CAP, selector_rounds=SELECTOR_ROUNDS):
    _parameters(1, 1, 1, selector_rounds, cap, (0,), False)
    return _stopped_kernel(cap, selector_rounds)


@lru_cache(maxsize=8)
def _stopped_kernel(cap, selector_rounds):
    data = core.source()
    coin_den = 16**selector_rounds
    attempt_den = data.denominator*coin_den
    eye = tuple(tuple(int(i == j) for j in range(66)) for i in range(66))
    power, succeeded, rejected = eye, tuple((0,)*66 for _ in range(66)), tuple((0,)*66 for _ in range(66))
    denominator = 1
    for _ in range(cap):
        incident = _left(data.incident_counts, power)
        remote = _left(data.remote_counts, power)
        succeeded = tuple(tuple(attempt_den*old+(coin_den-1)*new for old, new in zip(a, b, strict=True))
                          for a, b in zip(succeeded, incident, strict=True))
        rejected = tuple(tuple(attempt_den*old+data.denominator*new for old, new in zip(a, b, strict=True))
                         for a, b in zip(rejected, power, strict=True))
        power = tuple(tuple((coin_den-1)*v for v in row) for row in remote)
        denominator *= attempt_den
    if any(sum(succeeded[i][j]+rejected[i][j]+power[i][j] for i in range(66)) != denominator
           for j in range(66)):
        raise ArithmeticError('finite selector/incident/remote stopping lost a complete output')
    return StopKernel(succeeded, rejected, power, denominator,
                      _rational_hash((*succeeded, *rejected, *power)))


class FiniteFieldLaw(NamedTuple):
    records: tuple
    counts: tuple
    denominator: int
    read_ports: tuple
    feedback: bool
    initial_index: int
    source_digest: str
    controller_parameters: tuple  # seed scans, retries, target/selector rounds, cap

    def probabilities(self):
        _validate_field_law(self)
        return {key: F(sum(values), self.denominator)
                for key, values in zip(self.records, self.counts, strict=True)}


def record_alphabet(read_ports=core.READ_PORTS):
    _parameters(1, 1, 1, 1, 1, read_ports, False)
    records = [('target_timeout',), ('reset_timeout',)]
    for length in range(len(read_ports)):
        for initial, word, flag in product(product((0, 1), repeat=2),
                                          product((0, 1), repeat=length),
                                          ('selector_timeout', 'incident_timeout')):
            records.append((flag, initial, word))
    records.extend(('complete', initial, word) for initial, word in
                   product(product((0, 1), repeat=2), product((0, 1), repeat=len(read_ports))))
    return tuple(records)


def _vector(matrix, vector):
    return tuple(sum(coefficient*value for coefficient, value in zip(row, vector, strict=True)
                     if coefficient and value) for row in matrix)


def finite_field_law(initial_index=0, *, feedback=False, seed_scans=SEED_SCANS,
                     retries=RETRIES, target_rounds=TARGET_ROUNDS,
                     selector_rounds=SELECTOR_ROUNDS, cap=CAP, read_ports=core.READ_PORTS):
    _parameters(seed_scans, retries, target_rounds, selector_rounds, cap, read_ports, feedback)
    if type(initial_index) is not int or not 0 <= initial_index < 78:
        raise ValueError('an original full78 input index required')
    return _finite_field_law(initial_index, feedback, seed_scans, retries, target_rounds,
                             selector_rounds, cap, read_ports)


@lru_cache(maxsize=12)
def _finite_field_law(initial_index, feedback, seed_scans, retries, target_rounds,
                      selector_rounds, cap, read_ports):
    table, kernel, data = preparation_table(seed_scans, retries, target_rounds), stopped_kernel(cap, selector_rounds), core.source()
    if table.targets != data.states:
        raise ArithmeticError('native preparation and field basis order disagree')
    weights = table.probabilities[initial_index]
    prep_failure = 1-sum(weights, F())
    denominator = math.lcm(*(p.denominator for p in (*weights, table.selector_timeout, prep_failure)))
    starting = tuple(int(p*denominator) for p in weights)
    terminals = {('target_timeout',): ((int(table.selector_timeout*denominator),), denominator),
                 ('reset_timeout',): ((int((prep_failure-table.selector_timeout)*denominator),), denominator)}
    branches = {(label, ()): tuple(v if tuple(state[p] for p in (0, 1)) == label else 0
                                   for v, state in zip(starting, table.targets, strict=True))
                for label in product((0, 1), repeat=2)}
    for cursor, port in enumerate(read_ports):
        following = {}
        denominator *= kernel.denominator
        feedback_den = 2 if feedback and cursor == 0 else 1
        for (initial, word), state in branches.items():
            for flag, operator in (('selector_timeout', kernel.selector_timeout),
                                   ('incident_timeout', kernel.incident_timeout)):
                terminals[flag, initial, word] = (_vector(operator, state), denominator)
            evolved = _vector(kernel.success, state)
            for bit in (0, 1):
                selected = tuple(v if target[port] == bit else 0
                                 for v, target in zip(evolved, table.targets, strict=True))
                if feedback_den == 2:
                    selected = (_vector(data.internal_counts, selected) if bit else
                                tuple(2*v for v in selected))
                following[initial, (*word, bit)] = selected
        denominator *= feedback_den
        branches = following
    terminals.update({('complete', initial, word): (state, denominator)
                      for (initial, word), state in branches.items()})
    records = record_alphabet(read_ports)
    counts = tuple(tuple(v*(denominator//d) for v in values) for values, d in
                   (terminals[key] for key in records))
    if (set(terminals) != set(records) or sum(map(sum, counts)) != denominator
            or any(v < 0 for row in counts for v in row)):
        raise ArithmeticError('complete prepared field histories lost probability mass')
    return FiniteFieldLaw(records, counts, denominator, read_ports, feedback, initial_index,
                          _rational_hash(counts),
                          (seed_scans, retries, target_rounds, selector_rounds, cap))


def _same_exact(left, right):
    if type(left) is not type(right):
        return False
    if isinstance(right, tuple):
        return len(left) == len(right) and all(_same_exact(a, b) for a, b in zip(left, right, strict=True))
    return left == right


def _validate_field_law(law):
    """Bind every terminal population to its source, not a supplied digest."""
    if (type(law) is not FiniteFieldLaw or type(law.controller_parameters) is not tuple
            or len(law.controller_parameters) != 5):
        raise ValueError('an exact source-derived finite prepared field law required')
    scans, retries, targets, selectors, cap = law.controller_parameters
    canonical = finite_field_law(law.initial_index, feedback=law.feedback, seed_scans=scans,
        retries=retries, target_rounds=targets, selector_rounds=selectors,
        cap=cap, read_ports=law.read_ports)
    if not _same_exact(law, canonical):
        raise ValueError('finite field law differs from complete source-derived terminal populations')


def resources(*, feedback=False, seed_scans=SEED_SCANS, retries=RETRIES,
              target_rounds=TARGET_ROUNDS, selector_rounds=SELECTOR_ROUNDS, cap=CAP,
              read_ports=core.READ_PORTS):
    _parameters(seed_scans, retries, target_rounds, selector_rounds, cap, read_ports, feedback)
    table = preparation_table(seed_scans, retries, target_rounds)
    rb = reset.budget(seed_scans, retries, full=True, mode='reset', target=table.targets[0])
    prep = table.selector_timeout+(1-table.selector_timeout)*rb['failure_upper']
    scheduler, waits = F(1, 16**selector_rounds), F(7, 10)**cap
    continuation = len(read_ports)*(cap*scheduler+waits)
    exact_error = _complete_error(seed_scans, retries, target_rounds, selector_rounds, cap, len(read_ports))
    return {
        'preparation_native_attempts_max': rb['native_attempts_max']+7*target_rounds,
        'field_native_attempts_max': len(read_ports)*cap*(1+5*selector_rounds)+int(feedback),
        'all_native_attempts_max': rb['native_attempts_max']+7*target_rounds+len(read_ports)*cap*(1+5*selector_rounds)+int(feedback),
        'maximum_scheduled_incident_ticks': len(read_ports),
        'maximum_feedback_incident_ticks': int(feedback),
        'maximum_auxiliary_native_attempts': 7*target_rounds+len(read_ports)*cap*5*selector_rounds,
        'preparation_complete_halfdiamond_upper': str(prep),
        'field_selector_timeout_union_upper': str(len(read_ports)*cap*scheduler),
        'field_incident_timeout_union_upper': str(len(read_ports)*waits),
        'continuation_complete_record_TV_upper': str(continuation),
        'complete_field_record_TV_union_upper': str(prep+continuation),
        'complete_field_record_TV_upper': str(_up(exact_error)),
        'exact_native_continuation_halfdiamond_enclosure': [str(_down(exact_error)), str(_up(exact_error))],
        'exact_native_continuation_halfdiamond_sha256': _rational_hash(((exact_error,),)),
        'finite_tick_completion_probability': str(_completion_probability(cap, selector_rounds)),
        'finite_tick_probability_state_independent': True,
        'reference_bound_multiplied_by_history_length': False,
        'failed_trial_costs_included': True, 'coherent_services': 0,
        'record_storage': 'bounded native fine replay retained privately; initial bits and field word public',
        'excluded_costs': 'classical arithmetic, carrier manufacture, seam/current access and physical wall-clock duration',
    }


def _completion_probability(cap, selector_rounds):
    alpha = 1-F(1, 16**selector_rounds)
    remote = F(7, 10)*alpha
    return F(3, 10)*alpha*(1-remote**cap)/(1-remote)


@lru_cache(maxsize=16)
def _complete_error(seed_scans, retries, target_rounds, selector_rounds, cap, reads):
    failure = max(preparation_table(seed_scans, retries, target_rounds).failure_by_input)
    return 1-(1-failure)*_completion_probability(cap, selector_rounds)**reads


@dataclass(frozen=True)
class History:
    law: str
    flag: str
    initial_label: tuple | None
    word: tuple
    read_clock: tuple
    terminal_work: tuple
    terminal_auxiliary: tuple
    native_records: tuple  # private chronological (carrier, phase, seam, coin, fine label)
    reset_reason: str | None

    @property
    def native_attempts(self):
        return len(self.native_records)

    @property
    def public_record(self):
        if self.flag in ('target_timeout', 'reset_timeout'):
            return (self.flag,)
        return self.flag, self.initial_label, self.word


def execute_history(initial, auxiliary_initial, bits, *, law='C', feedback=False,
                    seed_scans=SEED_SCANS, retries=RETRIES, target_rounds=TARGET_ROUNDS,
                    selector_rounds=SELECTOR_ROUNDS, cap=CAP, read_ports=core.READ_PORTS):
    """Microscopic native C/E trajectory, with one continuing auxiliary.

    The iterator samples native placement outcomes in this basis-state
    simulator. Controllers read actual source records, never hidden states.
    All target and seam decisions read auxiliary coins.
    """
    _parameters(seed_scans, retries, target_rounds, selector_rounds, cap, read_ports, feedback)
    if type(law) is not str or law not in ('C', 'E'):
        raise ValueError('original native C/E law required')
    state = reset._configuration(initial, binary=False)
    auxiliary = reset._configuration(auxiliary_initial, binary=False)
    stream, records = iter(bits), []
    label, word, clocks, ticks, reset_reason = None, (), (), 0, None

    def bit():
        try:
            value = next(stream)
        except StopIteration as exc:
            raise ValueError('native placement tape ended before a complete terminal outcome') from exc
        if type(value) is not int or value not in (0, 1):
            raise ValueError('original integer native coin outcomes required')
        return value

    def work_step(seam, phase):
        nonlocal state
        c = bit()
        final, (observed_coin, current) = native.native(state, seam, c)
        if observed_coin != c:
            raise ArithmeticError('native placement record changed')
        mark = current if law == 'C' else state[seam[0]]-state[seam[1]]
        records.append(('work', phase, seam, observed_coin, mark))
        state = final
        return current

    def select(size, rounds, phase):
        nonlocal auxiliary
        width = (size-1).bit_length()
        seam = core.source().seams[0]
        for _ in range(rounds):
            word_bits = []
            for _ in range(width):
                c = bit()
                final, (observed_coin, current) = native.native(auxiliary, seam, c)
                if observed_coin != c:
                    raise ArithmeticError('native auxiliary coin record changed')
                mark = current if law == 'C' else auxiliary[seam[0]]-auxiliary[seam[1]]
                records.append(('auxiliary', phase, seam, observed_coin, mark))
                auxiliary = final
                word_bits.append(observed_coin)
            address = coins.decode_coin_word(tuple(word_bits), size)
            if address is not None:
                return address
        return None

    def finish(flag):
        history = History(law, flag, label, word, clocks, state, auxiliary, tuple(records), reset_reason)
        if history.native_attempts > resources(feedback=feedback, seed_scans=seed_scans, retries=retries,
                target_rounds=target_rounds, selector_rounds=selector_rounds, cap=cap,
                read_ports=read_ports)['all_native_attempts_max']:
            raise ArithmeticError('microscopic execution exceeded its declared native budget')
        return history

    target = select(66, target_rounds, 'target')
    if target is None:
        return finish('target_timeout')
    target_state = tuple(native.states(binary=True))[target]
    result = reset.execute(state, stream, seed_scans, retries, full=True,
                           mode='reset', target=target_state, law=law)
    records.extend(('work', 'preparation', seam, c, mark) for seam, c, mark in result.records)
    state, reset_reason = result.output, result.outcome.reason
    if not result.outcome.success:
        return finish('reset_timeout')
    if result.outcome.final != target_state or state != target_state:
        raise ArithmeticError('the current-only reset did not produce its selected target')
    # Read the controller's conclusive final label, not the hidden state.
    label = tuple(result.outcome.final[p] for p in (0, 1))
    local = label
    for cursor, port in enumerate(read_ports):
        for _ in range(cap):
            address = select(30, selector_rounds, 'schedule')
            if address is None:
                return finish('selector_timeout')
            seam = core.source().seams[address]
            current = work_step(seam, 'scheduled')
            local = core.update_field(local, seam, current)
            if set(seam).intersection((0, 1)):
                ticks += 1
                word, clocks = (*word, local[port]), (*clocks, ticks)
                break
        else:
            return finish('incident_timeout')
        if feedback and cursor == 0 and word[0]:
            current = work_step((0, 1), 'feedback')
            local = core.update_field(local, (0, 1), current)
            ticks += 1
    if clocks != core.read_clock(word, feedback=feedback):
        raise ArithmeticError('record-controlled feedback was not charged as an incident tick')
    return finish('complete')


def replay_history(history, initial, auxiliary_initial, **kwargs):
    """Replay every source record and controller choice; reject altered fine labels."""
    if type(history) is not History:
        raise ValueError('an immutable microscopic History required')
    actual = execute_history(initial, auxiliary_initial, (row[3] for row in history.native_records),
                             law=history.law, **kwargs)
    if actual != history:
        raise ValueError('native trajectory or its controller/public output was altered')
    return True


def correlated_coin_audit(law):
    """Actual source arrays on one entangled auxiliary, with adaptive work action."""
    maps = coins.coin_maps(law, core.source().seams[0], binary=False)
    d, minor = 78, 4  # minor is work-qubit tensor reference-qubit
    vector = np.zeros(d*minor, complex)
    vector[0*minor], vector[1*minor+3] = .6, .8j
    branches = {(): np.outer(vector, vector.conj())}
    maximum = 0.
    for _ in range(3):
        following = {}
        for word, density in branches.items():
            blocks = density.reshape(d, minor, d, minor)
            marginal = np.trace(blocks, axis1=0, axis2=2)
            for coin in (0, 1):
                result = np.zeros_like(blocks)
                for c, _, mapping in maps:
                    if c == coin:
                        source, target = np.array(mapping, dtype=int).T
                        result[target[:, None], :, target[None, :], :] += blocks[source[:, None], :, source[None, :], :]/2
                reduced = np.trace(result, axis1=0, axis2=2)
                maximum = max(maximum, float(np.linalg.norm(reduced-marginal/2, 'nuc')/2))
                u = np.array([[1, 1], [1, -1]], complex)/np.sqrt(2) if coin else np.array([[0, 1], [1, 0]], complex)
                operation = np.kron(u, np.eye(2))
                result = np.einsum('ab,ibjc,dc->iajd', operation, result, operation.conj())
                following[(*word, coin)] = result.reshape(d*minor, d*minor)
        branches = following
    if maximum > 64*np.finfo(float).eps:
        raise ArithmeticError('correlated reused native auxiliary failed its reduced-effect diagnostic')
    return {'law': law, 'sequential_uses': 3, 'complete_coin_words': len(branches),
            'initial_auxiliary_work_reference_entangled': True,
            'adaptive_work_operations_between_uses': True,
            'array_reduced_effect_audit_pass': True,
            'array_roundoff_allowance': str(F(64, 2**52)),
            'floating_diagnostics_not_accuracy_certificate': True,
            'exact_certificate': 'partial-injection coarse effects I/2 on every full78 and binary66 source column'}


def quantum_preparation_audit():
    table = preparation_table()
    failures = table.failure_by_input
    worst = max(range(78), key=failures.__getitem__)
    other = (worst+17) % 78
    psi = np.zeros(156, complex)
    psi[2*worst], psi[2*other+1] = .6, .8j
    output = apply_preparation(np.outer(psi, psi.conj()))
    ideal = output['reference']/66
    deficits = np.array([ideal-block for block in output['success']])
    minimum = min(float(np.linalg.eigvalsh(block).min()) for block in deficits)
    failure = float(np.trace(output['target_timeout']+output['reset_timeout']).real)
    distance = (sum(np.linalg.norm(block, 'nuc') for block in deficits)
                +np.linalg.norm(output['target_timeout'], 'nuc')
                +np.linalg.norm(output['reset_timeout'], 'nuc'))/2
    mass = output['success'].sum(axis=0)+output['target_timeout']+output['reset_timeout']
    allowance = 512*np.finfo(float).eps
    if minimum < -allowance or abs(distance-failure) > allowance or np.linalg.norm(mass-output['reference'], 'nuc') > allowance:
        raise ArithmeticError('complete reference preparation failed its positive-deficit diagnostic')
    exact_failure = F(9, 25)*failures[worst]+F(16, 25)*failures[other]
    return {'reference_dimension': 2, 'entangled_input_columns': [worst, other],
            'exact_complete_halftrace_distance_enclosure': [str(_down(exact_failure)), str(_up(exact_failure))],
            'array_positive_deficit_and_mass_audit_pass': True,
            'array_roundoff_allowance': str(F(512, 2**52)),
            'floating_diagnostics_not_accuracy_certificate': True,
            'success_not_renormalized': True}


def field_certificate(law):
    if type(law) is not FiniteFieldLaw:
        raise ValueError('an exact finite prepared field law required')
    probabilities = law.probabilities()
    ideal = core.field_law(law.read_ports, feedback=law.feedback).probabilities()
    variation = sum((abs(p-(ideal.get((key[1], key[2]), F()) if key[0] == 'complete' else 0))
                     for key, p in probabilities.items()), F())/2
    complete_mass = sum((p for key, p in probabilities.items() if key[0] == 'complete'), F())
    if (any(p > ideal[key[1], key[2]] for key, p in probabilities.items() if key[0] == 'complete')
            or variation != 1-complete_mass):
        raise ArithmeticError('successful prepared field records lost ideal domination')
    flags = {flag: sum((p for key, p in probabilities.items() if key[0] == flag), F())
             for flag in ('target_timeout', 'reset_timeout', 'selector_timeout', 'incident_timeout', 'complete')}
    return {'initial_full78_column': law.initial_index, 'feedback': law.feedback,
            'records_including_impossible_zero_branches': len(law.records),
            'complete_terminal_population_sha256': law.source_digest,
            'mass_exactly_one': sum(probabilities.values(), F()) == 1,
            'complete_success_probability_lower': str(_down(complete_mass)),
            'complete_record_TV_to_existing_ideal_upper': str(_up(variation)),
            'complete_record_TV_equals_total_failure_probability': True,
            'terminal_flag_probability_enclosures': {flag: [str(_down(p)), str(_up(p))] for flag, p in flags.items()},
            'successful_record_probability_enclosures': [
                {'initial_bits': list(key[1]), 'word': list(key[2]),
                 'read_clock': list(core.read_clock(key[2], feedback=law.feedback)),
                 'probability': [str(_down(p)), str(_up(p))]}
                for key, p in probabilities.items() if key[0] == 'complete'],
            'preparation_failed_work_retained_in_comparison': False,
            'field_timeout_terminal_work_populations_retained': True}


def _bits(seed):
    generator = random.Random(seed)
    while True:
        yield generator.randrange(2)


def microscopic_certificate():
    rows = []
    inputs = tuple(native.states())
    for law, feedback in product(('C', 'E'), (False, True)):
        histories = []
        for index, state in enumerate(inputs):
            auxiliary = inputs[(index*17+3) % 78]
            history = execute_history(state, auxiliary, _bits(20261025+index), law=law, feedback=feedback)
            replay_history(history, state, auxiliary, feedback=feedback)
            if history.flag == 'complete' and history.terminal_work not in core.source().states:
                raise ArithmeticError('continuing field execution left the binary native sector')
            histories.append(history)
        rows.append({'law': law, 'feedback': feedback, 'full78_inputs_executed': len(histories),
                     'terminal_flags': {flag: sum(h.flag == flag for h in histories) for flag in
                         ('target_timeout', 'reset_timeout', 'selector_timeout', 'incident_timeout', 'complete')},
                     'native_attempts_total': sum(h.native_attempts for h in histories),
                     'native_attempts_max_observed': max(h.native_attempts for h in histories),
                     'private_replay_sha256': _digest([h.__dict__ for h in histories]),
                     'same_auxiliary_reused_without_reset': True,
                     'all_public_fields_computed_from_currents_and_conclusive_initial_bits': True})
    return rows


def run():
    source_audit, coin_audit = reset.source_audit(), coins.effect_certificate()
    prepared = preparation_certificate()
    exact = [finite_field_law(feedback=flag) for flag in (False, True)]
    fields = [field_certificate(law) for law in exact]
    bound = F(resources(feedback=True)['complete_field_record_TV_upper'])
    ideal_separation = core.total_variation(core.field_law(), core.field_law(feedback=True))
    observed_separation = sum((abs(p-q) for p, q in zip(exact[0].probabilities().values(),
                                      exact[1].probabilities().values(), strict=True)), F())/2
    if (any(F(row['complete_record_TV_to_existing_ideal_upper']) > bound for row in fields)
            or observed_separation < ideal_separation-2*bound or ideal_separation != F(5, 99)):
        raise ArithmeticError('native preparation/continuation did not preserve the existing feedback prediction')
    return {'schema': SCHEMA, 'all_pass': True,
            'parameters': {'seed_scans': SEED_SCANS, 'retries': RETRIES, 'target_rounds': TARGET_ROUNDS,
                           'selector_rounds': SELECTOR_ROUNDS, 'incident_attempt_cap': CAP,
                           'read_ports': list(core.READ_PORTS)},
            'source_sha256': bound_pins({p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS}),
            'scope': dict(SCOPE), 'source_kraus_audit': source_audit, 'native_coin_effect_audit': coin_audit,
            'preparation': prepared, 'resources': resources(feedback=True),
            'reference_execution': quantum_preparation_audit(),
            'correlated_auxiliary_execution': [correlated_coin_audit(law) for law in ('C', 'E')],
            'finite_field_laws': fields,
            'feedback_separation': {'existing_ideal_TV': str(ideal_separation),
                'prepared_finite_TV_lower': str(_down(observed_separation)),
                'uniform_all_input_TV_separation_lower': str(ideal_separation-2*bound),
                'same_coarse_field_record_alphabet_used': True},
            'microscopic_native_executions': microscopic_certificate(),
            'effective_dictionary': {'existing_local_field_and_feedback_dictionary_reused': True,
                'no_stationary_initial_state_or_initial_bit_service_assumed': True,
                'only_public_field_cylinder_laws_compared': True,
                'hidden_native_terminal_state_equals_effective_terminal_state_claimed': False}}


def compare_receipt(supplied, actual=None):
    """Fresh replay is mandatory; an alleged actual receipt is not an authority."""
    canonical = run()
    expected = json.dumps(canonical, sort_keys=True, separators=(',', ':'), allow_nan=False)
    for value in (supplied,) if actual is None else (supplied, actual):
        try:
            text = json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)
        except (ValueError, TypeError) as exc:
            raise ValueError('noncanonical native preparation receipt') from exc
        if text != expected:
            raise ValueError('native preparation receipt differs from complete source replay')
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group()
    action.add_argument('--write', action='store_true')
    action.add_argument('--check', nargs='?', const=str(OUTPUT))
    parser.add_argument('--output', type=Path, default=OUTPUT)
    args = parser.parse_args()
    actual = run()
    if args.check:
        supplied = strict_json.load(Path(args.check))
        if json.dumps(supplied, sort_keys=True, allow_nan=False) != json.dumps(actual, sort_keys=True, allow_nan=False):
            raise SystemExit('native labelled preparation receipt differs')
        print('native labelled preparation receipt verified')
    elif args.write:
        with args.output.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(actual, indent=2, sort_keys=True, allow_nan=False)+'\n')
        print(f'wrote {args.output}')
    else:
        print(json.dumps(actual, indent=2, sort_keys=True, allow_nan=False))


if __name__ == '__main__':
    main()
