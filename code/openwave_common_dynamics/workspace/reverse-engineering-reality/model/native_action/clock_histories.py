"""Seeded complete source histories stopped at incident repair opportunities.

The coupling samples a uniform union seam. Erasing slots absent from one
screen leaves that screen's uniform native scheduler. These reproducible
examples are not frequency estimates or a complete-law comparison certificate.
"""
from __future__ import annotations

import base64
from functools import lru_cache
import hashlib
import json
import random
import struct
from typing import NamedTuple

from tools import strict_json
from ..observer import native
from . import screen_compatibility as screen


LEVELS = (3, 5)
MAX_CAP = 500000
SOURCE_PATHS = tuple(sorted(set((*screen.SOURCE_PATHS,
    'model/native_action/clock_histories.py', 'tools/strict_json.py'))))


class UnionSource(NamedTuple):
    ports: tuple
    left_edges: tuple
    right_only: tuple
    # left edge index, right edge, canonical edge, left/right reversal flags.
    shared: tuple
    matched_addresses: tuple
    observed_labels: tuple
    observed_addresses: tuple
    incident_edge_indices: tuple
    source_digest: str

    @property
    def slots(self):
        return len(self.left_edges)+len(self.right_only)


@lru_cache(maxsize=1)
def source():
    """All seams of both actual source screens; no prescribed tape replaces them."""
    balls = tuple(screen.root_ball(level, 9) for level in LEVELS)
    screen.identify(*balls)
    addresses = tuple(dict(ball.addresses) for ball in balls)
    inverse = tuple({p: label for label, p in table.items()} for table in addresses)
    graphs = tuple(screen.screen(level) for level in LEVELS)
    edges = tuple(tuple(sorted(tuple(sorted((int(a), int(b))))
        for a, b in zip(graph.seam_a, graph.seam_b, strict=True))) for graph in graphs)
    if any(len(set(menu)) != graph.seams for menu, graph in zip(edges, graphs, strict=True)):
        raise ArithmeticError('the actual source seam menu is incomplete or duplicated')
    right_set, shared, paired = set(edges[1]), [], set()
    incident = []
    for index, (a, b) in enumerate(edges[0]):
        if a not in inverse[0] or b not in inverse[0]:
            continue
        ca, cb = inverse[0][a], inverse[0][b]
        canonical = tuple(sorted((ca, cb)))
        right_oriented = tuple(addresses[1][p] for p in canonical)
        other = tuple(sorted(right_oriented))
        if other not in right_set or other in paired:
            raise ArithmeticError('the source port map lost or duplicated a shared seam')
        paired.add(other)
        shared.append((index, other, canonical, ca > cb, right_oriented[0] > right_oriented[1]))
        if canonical[0] < 12 or canonical[1] < 12:
            incident.append(index)
    right_only = tuple(edge for edge in edges[1] if edge not in paired)
    root_edges = tuple(row[2] for row in shared if row[0] in incident)
    observed = tuple(sorted(set(range(12)) | {p for edge in root_edges for p in edge}))
    for side, menu in enumerate(edges):
        root = {addresses[side][p] for p in range(12)}
        actual = {edge for edge in menu if root.intersection(edge)}
        mapped = {tuple(sorted(addresses[side][p] for p in edge)) for edge in root_edges}
        if actual != mapped:
            raise ArithmeticError('a root-incident source seam is absent from the shared clock')
    if len(root_edges) != 33 or len(observed) != 15 or (0, 1) not in root_edges:
        raise ArithmeticError('the declared root clock and feedback census changed')
    digest = hashlib.sha256(b'oph.union-screen-clock.v1')
    for menu in edges:
        for edge in menu:
            digest.update(struct.pack('>II', *edge))
    digest.update(json.dumps(shared, separators=(',', ':')).encode())
    for path in SOURCE_PATHS:
        digest.update(path.encode())
        digest.update(hashlib.sha256((screen.ROOT/path).read_bytes()).digest())
    return UnionSource(tuple(int(graph.ports) for graph in graphs), edges[0], right_only,
        tuple(shared), tuple(tuple(sorted(table.items())) for table in addresses), observed,
        tuple(tuple(table[p] for p in observed) for table in addresses), tuple(incident), digest.hexdigest())


class Checkpoint(NamedTuple):
    source_digest: str
    seed: int
    ticks: int
    cap: int
    feedback: bool
    slots: int
    scheduled_hits: int
    uniform_attempts: tuple
    feedback_attempts: tuple
    field_values: tuple
    initial_field_labels: tuple
    initial_totals: tuple
    initial_state_sha256: tuple
    states: tuple
    rng_state: tuple
    records: tuple


def _parameters(seed, ticks, cap, feedback, stop):
    if (type(seed) is not int or not 0 <= seed < 2**32
            or type(ticks) is not int or not 1 <= ticks <= 4
            or type(cap) is not int or not 1 <= cap <= MAX_CAP
            or type(feedback) is not bool
            or type(stop) is not int or not 0 <= stop <= cap):
        raise ValueError('exact bounded seed, one to four hits, cap, Boolean feedback and slot limit required')


def _initial(seed, ticks, cap, feedback):
    data = source()
    rng = random.Random(seed)
    states = [bytearray(int(rng.randrange(3) == 0) for _ in range(ports)) for ports in data.ports]
    # Replace both matched restrictions by one common Bernoulli product draw.
    # The independently drawn outside ports remain untouched.
    left, right = map(dict, data.matched_addresses)
    for label in sorted(left):
        value = int(rng.randrange(3) == 0)
        states[0][left[label]] = states[1][right[label]] = value
    values = tuple(state[table[0]] for state, table in zip(states, (left, right), strict=True))
    return Checkpoint(data.source_digest, seed, ticks, cap, feedback, 0, 0, (0, 0), (0, 0),
        values, values, tuple(map(sum, states)), tuple(hashlib.sha256(state).hexdigest() for state in states),
        tuple(map(bytes, states)), rng.getstate(), ((), ()))


def _native_step(state, seam, coin):
    a, b = seam
    before = state[a], state[b]
    after, (_, current) = native.native(before, (0, 1), coin)
    if sum(after) != sum(before):
        raise ArithmeticError('a native source branch changed total occupation')
    state[a], state[b] = after
    return current


def _continue(checkpoint, stop):
    data = source()
    states = [bytearray(state) for state in checkpoint.states]
    rng = random.Random()
    rng.setstate(checkpoint.rng_state)
    shared = {row[0]: row[1:] for row in data.shared}
    incident = set(data.incident_edge_indices)
    addresses = tuple(dict(row) for row in data.matched_addresses)
    attempts, repairs = list(checkpoint.uniform_attempts), list(checkpoint.feedback_attempts)
    fields = list(checkpoint.field_values)
    records = [list(row) for row in checkpoint.records]
    slots, hits = checkpoint.slots, checkpoint.scheduled_hits

    def record(side, kind, canonical, coin, current):
        fields[side] += current*(int(canonical[1] == 0)-int(canonical[0] == 0))
        if fields[side] != states[side][addresses[side][0]]:
            raise ArithmeticError('incident-current field counter disagrees with the actual source')
        # The controller uses only this counter, never the hidden source state.
        records[side].append((kind, slots, hits, len(records[side])+1,
                              *canonical, coin, current, fields[side]))

    while slots < stop and hits < checkpoint.ticks:
        index, coin = rng.randrange(data.slots), rng.randrange(2)
        slots += 1
        if index >= len(data.left_edges):
            _native_step(states[1], data.right_only[index-len(data.left_edges)], coin)
            attempts[1] += 1
            continue
        edge = data.left_edges[index]
        attempts[0] += 1
        if index not in shared:
            _native_step(states[0], edge, coin)
            continue
        other, canonical, reverse_left, reverse_right = shared[index]
        currents = (
            _native_step(states[0], edge, coin ^ reverse_left)*(1-2*reverse_left),
            _native_step(states[1], other, coin ^ reverse_right)*(1-2*reverse_right))
        attempts[1] += 1
        if index not in incident:
            continue
        hits += 1
        for side, current in enumerate(currents):
            record(side, 'scheduled', canonical, coin, current)
        if hits == 1 and checkpoint.feedback:
            # Draw even when neither side requests feedback; later scheduling
            # must not depend on an observer's field value or hidden state.
            feedback_coin = rng.randrange(2)
            for side, table in enumerate(addresses):
                if fields[side] == 1:
                    oriented = table[0], table[1]
                    reverse = oriented[0] > oriented[1]
                    current = _native_step(states[side], tuple(sorted(oriented)),
                                           feedback_coin ^ reverse)*(1-2*reverse)
                    repairs[side] += 1
                    record(side, 'feedback', (0, 1), feedback_coin, current)
    if tuple(map(sum, states)) != checkpoint.initial_totals:
        raise ArithmeticError('continuing full-screen execution lost particle number')
    return checkpoint._replace(slots=slots, scheduled_hits=hits, uniform_attempts=tuple(attempts),
        feedback_attempts=tuple(repairs), field_values=tuple(fields), states=tuple(map(bytes, states)),
        rng_state=rng.getstate(), records=tuple(tuple(row) for row in records))


@lru_cache(maxsize=12)
def _compute(seed, ticks, cap, feedback, stop):
    return _continue(_initial(seed, ticks, cap, feedback), stop)


def _builtins(value):
    return type(value) in (str, int, bool, bytes, type(None)) or (
        type(value) is tuple and all(_builtins(item) for item in value))


def _exact_equal(first, second):
    if type(first) is not type(second):
        return False
    if type(first) is tuple:
        return len(first) == len(second) and all(_exact_equal(a, b) for a, b in zip(first, second))
    return first == second


def validate_checkpoint(checkpoint):
    """A digest is not ancestry: compare with the actual seeded source prefix."""
    if type(checkpoint) is not Checkpoint or not all(_builtins(value) for value in checkpoint):
        raise ValueError('an immutable typed source checkpoint is required')
    _parameters(checkpoint.seed, checkpoint.ticks, checkpoint.cap, checkpoint.feedback, checkpoint.slots)
    expected = _compute(checkpoint.seed, checkpoint.ticks, checkpoint.cap,
                        checkpoint.feedback, checkpoint.slots)
    if not all(_exact_equal(a, b) for a, b in zip(checkpoint, expected, strict=True)):
        raise ValueError('checkpoint differs from its seeded source replay')
    return expected


def execute(seed, *, ticks=4, cap=MAX_CAP, feedback=True, stop_after_slots=None, checkpoint=None):
    stop = cap if stop_after_slots is None else stop_after_slots
    _parameters(seed, ticks, cap, feedback, stop)
    if checkpoint is None:
        return _compute(seed, ticks, cap, feedback, stop)
    checked = validate_checkpoint(checkpoint)
    if (checked.seed, checked.ticks, checked.cap, checked.feedback) != (seed, ticks, cap, feedback):
        raise ValueError('continuation cannot replace its source preparation or controller')
    if stop < checked.slots:
        raise ValueError('continuation cannot rewind its native history')
    return _continue(checked, stop)


def _payload(checkpoint):
    result = checkpoint._asdict()
    result['states'] = [base64.b64encode(state).decode('ascii') for state in checkpoint.states]
    return result


def checkpoint_to_json(checkpoint):
    return json.dumps(_payload(validate_checkpoint(checkpoint)), separators=(',', ':'),
                      sort_keys=True, allow_nan=False)


def checkpoint_from_json(text):
    if type(text) is not str or len(text) > 1000000:
        raise ValueError('a bounded serialized full-source checkpoint is required')
    data = strict_json.loads(text)
    if type(data) is not dict or set(data) != set(Checkpoint._fields):
        raise ValueError('the complete checkpoint field menu is required')
    _parameters(data['seed'], data['ticks'], data['cap'], data['feedback'], data['slots'])
    expected = _compute(data['seed'], data['ticks'], data['cap'], data['feedback'], data['slots'])
    # Comparing canonical JSON preserves int/Boolean distinctions and rejects
    # all changed states, PRNG positions, records, counters and source bindings.
    if json.dumps(data, separators=(',', ':'), sort_keys=True, allow_nan=False) != json.dumps(
            _payload(expected), separators=(',', ':'), sort_keys=True, allow_nan=False):
        raise ValueError('serialized checkpoint differs from seeded source replay')
    return expected


def summarize(checkpoint):
    checkpoint = validate_checkpoint(checkpoint)
    data = source()
    status = 'complete' if checkpoint.scheduled_hits == checkpoint.ticks else (
        'timeout' if checkpoint.slots == checkpoint.cap else 'prefix')
    rows = []
    for side, level in enumerate(LEVELS):
        rows.append({'level': level, 'ports': data.ports[side],
            'supplied_bernoulli_preparations': 1,
            'initial_field_label': checkpoint.initial_field_labels[side],
            'initial_total_particles': checkpoint.initial_totals[side],
            'initial_state_sha256': checkpoint.initial_state_sha256[side],
            'terminal_state_sha256': hashlib.sha256(checkpoint.states[side]).hexdigest(),
            'incident_records': [[row[0], *row[2:]] for row in checkpoint.records[side]],
            'record_columns': ['kind', 'scheduled_hit', 'incident_tick',
                               'canonical_a', 'canonical_b', 'canonical_coin', 'signed_current', 'field_0'],
            'scheduled_field_word': [row[-1] for row in checkpoint.records[side] if row[0] == 'scheduled'],
            'terminal_observed_bits': [checkpoint.states[side][p] for p in data.observed_addresses[side]],
            'uniform_native_attempts': checkpoint.uniform_attempts[side],
            'feedback_native_attempts': checkpoint.feedback_attempts[side],
            'total_native_attempts': checkpoint.uniform_attempts[side]+checkpoint.feedback_attempts[side],
            'identity_comparison_slots': checkpoint.slots-checkpoint.uniform_attempts[side],
            'total_incident_clock_ticks': len(checkpoint.records[side]),
            'complete_source_state_bytes': len(checkpoint.states[side])})
    return {'seed': checkpoint.seed, 'requested_scheduled_hits': checkpoint.ticks,
        'external_comparison_slot_cap': checkpoint.cap, 'comparison_slots_executed': checkpoint.slots,
        'scheduled_hits_completed': checkpoint.scheduled_hits, 'status': status,
        'feedback_enabled': checkpoint.feedback, 'source_digest': data.source_digest,
        'canonical_observed_ports': list(data.observed_labels), 'screens': rows,
        'coupled_incident_records_equal': rows[0]['incident_records'] == rows[1]['incident_records'],
        'coupled_terminal_observed_bits_equal': rows[0]['terminal_observed_bits'] == rows[1]['terminal_observed_bits'],
        'deterministic_seeded_example_not_a_sampled_confidence_bound': True,
        'full_quantum_reference_state_executed': False}


def run():
    data = source()
    rows = []
    for ticks in (2, 4):
        for seed in (0, 1, 2, 3):
            full = execute(seed, ticks=ticks)
            prefix = execute(seed, ticks=ticks, stop_after_slots=full.slots//2)
            encoded = checkpoint_to_json(prefix)
            resumed = execute(seed, ticks=ticks, checkpoint=checkpoint_from_json(encoded))
            if resumed != full:
                raise ArithmeticError('serialized continuation differs from the full source execution')
            row = summarize(full)
            row['checkpoint_comparison_slot'] = prefix.slots
            row['serialized_full_state_checkpoint_bytes'] = len(encoded.encode())
            row['checkpoint_sha256'] = hashlib.sha256(encoded.encode()).hexdigest()
            row['serialized_continuation_exactly_reproduces_states_records_and_rng'] = True
            rows.append(row)
    return {'source_levels': list(LEVELS), 'global_ports': list(data.ports),
        'global_seams': [len(data.left_edges), len(data.right_only)+len(data.shared)],
        'matched_ball_ports': len(data.matched_addresses[0]), 'shared_seams': len(data.shared),
        'comparison_union_seams': data.slots, 'scheduled_root_incident_seams': len(data.incident_edge_indices),
        'clock_record_support_ports': len(data.observed_labels),
        'full_source_scheduler_not_a_selected_seam_tape': True,
        'marginal_scheduler': 'erase identity comparison slots: each native seam has equal weight',
        'preparation': 'Bernoulli density 1/3, common matched ball, independent outside',
        'sampling': 'Python Random exact-integer choices; fixed seeds are reproducible examples only',
        'comparison_cap_is_an_external_resource_cutoff_not_an_incident_clock_reading': True,
        'classical_diagonal_C_and_E_source_laws_agree': True,
        'forced_timeout_control': summarize(execute(1, ticks=2, cap=1)),
        'examples': rows}
