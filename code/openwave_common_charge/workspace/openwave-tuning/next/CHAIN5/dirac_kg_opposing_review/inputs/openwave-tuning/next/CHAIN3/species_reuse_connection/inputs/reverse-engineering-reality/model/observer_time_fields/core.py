"""Exact native fields on an observer's incident-attempt clock.

This is the diagonal, binary total-two sector. The original native rule
supplies every transition. Remote repairs continue while the observer waits.
The completed column kernel is K=B(I-A)^-1. Its symmetric positive transfer
T=(I-A)^-1/2 B (I-A)^-1/2 has exactly the same stationary local-field cylinder
laws, including the declared record-controlled internal-seam intervention.
The dictionary requires each insertion to commute with I-A. It does not
identify native event evolution with the unitary analytic continuation.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
import hashlib
from itertools import combinations, product
import json
import math
from typing import NamedTuple

import numpy as np
import sympy as sp

from ..native_quantum_history import construction as quantum
from ..observer import native
from ..relational_event_clock import kernel as movement


SOURCE_PATHS = tuple(sorted(set((*quantum.SOURCE_PATHS, *movement.SOURCE_PATHS,
    'model/observer/native.py', 'model/observer_time_fields/__init__.py',
    'model/observer_time_fields/core.py'))))
READ_PORTS = (0, 1, 0, 1)


def _digest(value):
    return hashlib.sha256(json.dumps(value, separators=(',', ':'), sort_keys=True).encode()).hexdigest()


def _tuples(matrix):
    return tuple(tuple(int(value) for value in row) for row in matrix.tolist())


class SourceSplit(NamedTuple):
    states: tuple
    seams: tuple
    region: tuple
    denominator: int
    remote_counts: tuple
    incident_counts: tuple
    internal_counts: tuple
    triangle: tuple
    triangle_counts: tuple
    incident_edges: tuple
    remote_edges: tuple
    tick_probability: F
    source_digest: str


@lru_cache(maxsize=1)
def source():
    """Original native enumeration, with immutable integer column matrices."""
    states, seams = tuple(quantum.STATES), tuple(tuple(e) for e in quantum.mixing.seams())
    region = seams[0]
    if region != (0, 1) or len(states) != 66 or len(seams) != 30:
        raise ArithmeticError('the retained native two-port instance changed')
    incident = tuple(i for i, edge in enumerate(seams) if set(edge).intersection(region))
    remote = tuple(i for i, edge in enumerate(seams) if not set(edge).intersection(region))
    triangle = next(((*region, v) for v in range(12) if v not in region
                     and all(tuple(sorted(edge)) in seams for edge in combinations((*region, v), 2))), None)
    if triangle is None:
        raise ArithmeticError('an incident source triangle is required for the strict spectral certificate')
    triangle_edges = {tuple(sorted(edge)) for edge in combinations(triangle, 2)}
    a, b, internal, tri = (np.zeros((66, 66), dtype=np.int64) for _ in range(4))
    for column, state in enumerate(states):
        for edge_index, edge in enumerate(seams):
            for coin in (0, 1):
                final, record = native.native(state, edge, coin)
                if record != (coin, state[edge[0]]-final[edge[0]]):
                    raise ArithmeticError('native signed-current grammar changed')
                row = quantum.INDEX[final]
                (b if edge_index in incident else a)[row, column] += 1
                if edge == region:
                    internal[row, column] += 1
                if edge in triangle_edges:
                    tri[row, column] += 1
                before = tuple(state[p] for p in region)
                after = update_field(before, edge, record[1], region=region)
                if after != tuple(final[p] for p in region):
                    raise ArithmeticError('record-only local field readback failed on a native branch')
    denominator = 2*len(seams)
    for matrix, mass in ((a, 2*len(remote)), (b, 2*len(incident)), (internal, 2), (tri, 6)):
        if (not np.array_equal(matrix, matrix.T) or np.any(matrix < 0)
                or np.any(matrix.sum(axis=0) != mass)):
            raise ArithmeticError('the exact native subtransfer lost symmetry, positivity of entries or mass')
    if not np.array_equal((tri-3*np.eye(66, dtype=np.int64))@(tri-6*np.eye(66, dtype=np.int64)),
                          np.zeros((66, 66), dtype=np.int64)):
        raise ArithmeticError('binary triangle spectrum is not contained in {3,6}')
    # Each unnormalized seam count is I plus its endpoint-swap involution.
    # This independent operator check supplies PSD for the remaining seams.
    for edge_index, edge in enumerate(seams):
        swap = np.zeros((66, 66), dtype=np.int64)
        for column, state in enumerate(states):
            target = list(state)
            target[edge[0]], target[edge[1]] = target[edge[1]], target[edge[0]]
            swap[quantum.INDEX[tuple(target)], column] = 1
        counts = np.zeros((66, 66), dtype=np.int64)
        for column, state in enumerate(states):
            for coin in (0, 1):
                final, _ = native.native(state, edge, coin)
                counts[quantum.INDEX[final], column] += 1
        if (not np.array_equal(counts, np.eye(66, dtype=np.int64)+swap)
                or not np.array_equal(swap@swap, np.eye(66, dtype=np.int64))):
            raise ArithmeticError('a native seam is not its certified positive exchange average')
    arrays = tuple(_tuples(matrix) for matrix in (a, b, internal, tri))
    digest = _digest([states, seams, region, denominator, *arrays])
    return SourceSplit(states=states, seams=seams, region=region, denominator=denominator,
        remote_counts=arrays[0], incident_counts=arrays[1], internal_counts=arrays[2],
        triangle=triangle, triangle_counts=arrays[3], incident_edges=incident,
        remote_edges=remote, tick_probability=F(len(incident), len(seams)), source_digest=digest)


def update_field(field, seam, current, *, region=(0, 1)):
    """Only disclosed local loads and a signed incident record enter here."""
    if (type(region) is not tuple or len(region) != 2
            or any(type(v) is not int or not 0 <= v < 12 for v in region)
            or len(set(region)) != len(region)
            or type(field) is not tuple or len(field) != len(region)
            or any(type(v) is not int or v not in (0, 1) for v in field)
            or type(seam) is not tuple or len(seam) != 2
            or any(type(v) is not int or not 0 <= v < 12 for v in seam)
            or seam not in tuple(quantum.mixing.seams())
            or type(current) is not int or current not in (-1, 0, 1)):
        raise ValueError('binary local fields and an exact signed native seam record required')
    result = tuple(value+current*(int(port == seam[1])-int(port == seam[0]))
                   for port, value in zip(region, field, strict=True))
    if any(v not in (0, 1) for v in result):
        raise ValueError('the current record is inconsistent with the disclosed binary field')
    return result


class IntegerKernel(NamedTuple):
    counts: tuple
    denominator: int
    source_digest: str


@lru_cache(maxsize=1)
def completed():
    """Exact B# (60 I-A#)^-1, stored with one integer denominator."""
    data = source()
    d = data.denominator*sp.eye(66)-sp.Matrix(data.remote_counts)
    inverse = d.inv(method='DM')
    if d*inverse != sp.eye(66):
        raise ArithmeticError('native no-tick elimination failed exact inversion')
    matrix = sp.Matrix(data.incident_counts)*inverse
    denominator = int(sp.ilcm(*(value.q for value in matrix)))
    counts = tuple(tuple(int(value*denominator) for value in row) for row in matrix.tolist())
    if (any(value < 0 for row in counts for value in row)
            or any(sum(row) != denominator for row in counts)
            or any(sum(row[i] for row in counts) != denominator for i in range(66))):
        raise ArithmeticError('the completed native event kernel is not bistochastic')
    return IntegerKernel(counts, denominator, data.source_digest)


def spectral_certificate():
    """Exact inequalities for the source-derived transfer and its logarithm.

    An incident triangle contributes at least 3 I to the integer B count.
    All other seam counts are PSD. A Gershgorin bound on the full symmetric
    native kernel improves D<=I. No floating eigensystem supplies a premise.
    """
    data = source()
    total = tuple(tuple(a+b for a, b in zip(ar, br, strict=True))
                  for ar, br in zip(data.remote_counts, data.incident_counts, strict=True))
    native_lower = min(F(2*total[i][i]-data.denominator, data.denominator) for i in range(66))
    d_upper = min(F(1), 1-native_lower+data.tick_probability)
    b_lower = F(3, data.denominator)
    lower = b_lower/d_upper
    if not 0 < lower < 1 or data.tick_probability != F(3, 10):
        raise ArithmeticError('the exact local-clock spectral domain changed')
    return {'triangle': list(data.triangle), 'triangle_count_annihilator': '(C-3I)(C-6I)=0',
            'incident_transfer_lower': b_lower, 'native_attempt_transfer_lower': native_lower,
            'D_lower': data.tick_probability, 'D_upper': d_upper,
            'effective_transfer_lower': lower, 'effective_transfer_upper': F(1),
            'full_effective_hamiltonian_defined': True,
            'kernel': 'K=B(I-A)^-1', 'effective_transfer': 'T=(I-A)^-1/2 B (I-A)^-1/2',
            'Hamiltonian': 'H=-log(T)', 'field_dictionary_requires_commutation_with_D': True}


def logarithm_remainder(degree):
    """Operator tail for sum_(j=1..m)(I-T)^j/j, without evaluating T."""
    if type(degree) is not int or not 1 <= degree <= 4096:
        raise ValueError('an integer logarithm degree from one through 4096 required')
    lower = spectral_certificate()['effective_transfer_lower']
    return (1-lower)**(degree+1)/((degree+1)*lower)


def _counts_apply(matrix, state):
    return tuple(sum(value*coefficient for value, coefficient in zip(row, state, strict=True)
                     if value and coefficient) for row in matrix)


class FieldLaw(NamedTuple):
    read_ports: tuple
    feedback: bool
    records: tuple
    counts: tuple
    denominator: int
    source_digest: str

    def probabilities(self):
        return {record: F(sum(state), self.denominator)
                for record, state in zip(self.records, self.counts, strict=True)}


def _read_ports(read_ports):
    if (type(read_ports) is not tuple or not 1 <= len(read_ports) <= 6
            or any(type(p) is not int or p not in source().region for p in read_ports)):
        raise ValueError('one through six reads at the declared two observer ports required')
    return read_ports


def _law_from_kernel(read_ports, feedback, kernel):
    data = source()
    # The actual uniform preparation is disclosed only on the two local
    # ports. Each branch is unnormalized, with common initial denominator66.
    branches = {(label, ()): tuple(int(tuple(state[p] for p in data.region) == label)
                                  for state in data.states)
                for label in product((0, 1), repeat=2)}
    denominator = 66
    for cursor, port in enumerate(read_ports):
        future = {}
        for (initial, word), state in branches.items():
            evolved = _counts_apply(kernel.counts, state)
            for bit in (0, 1):
                selected = tuple(value if native_state[port] == bit else 0
                                 for value, native_state in zip(evolved, data.states, strict=True))
                if feedback and cursor == 0:
                    selected = (_counts_apply(data.internal_counts, selected) if bit else
                                tuple(2*value for value in selected))
                future[initial, (*word, bit)] = selected
        denominator *= kernel.denominator*(2 if feedback and cursor == 0 else 1)
        branches = future
        if sum(sum(state) for state in branches.values()) != denominator:
            raise ArithmeticError('the complete continued local-field law lost probability mass')
    return FieldLaw(read_ports, feedback, tuple(branches), tuple(branches.values()),
                    denominator, kernel.source_digest)


@lru_cache(maxsize=12, typed=True)
def _field_law(read_ports, feedback):
    return _law_from_kernel(read_ports, feedback, completed())


def field_law(read_ports=READ_PORTS, *, feedback=False):
    """All initial local labels, read words and full terminal populations.

    Every scheduled interval ends at one actual incident attempt. If the
    first field bit is one, feedback applies one additional native repair
    on the internal source seam before the next scheduled interval. That
    intervention itself is another incident tick, not a free clock pause.
    """
    read_ports = _read_ports(read_ports)
    if type(feedback) is not bool:
        raise ValueError('an explicit Boolean native-feedback policy required')
    return _field_law(read_ports, feedback)


def read_clock(word, *, feedback=False):
    if (type(word) is not tuple or not 1 <= len(word) <= 6
            or any(type(v) is not int or v not in (0, 1) for v in word)
            or type(feedback) is not bool):
        raise ValueError('a nonempty binary read word and Boolean feedback policy required')
    return tuple(index+1+int(feedback and index > 0 and word[0] == 1) for index in range(len(word)))


def total_variation(left, right):
    if type(left) is not FieldLaw or type(right) is not FieldLaw or left.records != right.records:
        raise ValueError('two complete field laws on the same record menu required')
    a, b = left.probabilities(), right.probabilities()
    return sum((abs(a[key]-b[key]) for key in a), F())/2


def _partitions(items):
    if not items:
        yield ()
        return
    first, *rest = items
    for partition in _partitions(rest):
        yield ((first,), *partition)
        for index in range(len(partition)):
            yield (*partition[:index], (first, *partition[index]), *partition[index+1:])


def field_statistics(law):
    if type(law) is not FieldLaw:
        raise ValueError('an exact continued field law required')
    probabilities = law.probabilities()
    count = len(law.read_ports)
    moments = {subset: sum((p for (_, word), p in probabilities.items()
                           if all(word[index] for index in subset)), F())
               for length in range(count+1) for subset in combinations(range(count), length)}
    cumulant = sum(((-1)**(len(partition)-1)*math.factorial(len(partition)-1)
                    *math.prod(moments[tuple(sorted(block))] for block in partition)
                    for partition in _partitions(list(range(count)))), F())
    feedback_probability = moments[0,] if law.feedback else F(0)
    return {'mass': moments[()], 'means': tuple(moments[i,] for i in range(count)),
            'all_one_probability': moments[tuple(range(count))],
            'joint_cumulant': cumulant, 'feedback_probability': feedback_probability,
            'expected_native_attempts': F(count)/source().tick_probability+feedback_probability,
            'expected_incident_ticks': F(count)+feedback_probability}


def dictionary_certificate():
    data = source()
    d = data.denominator*np.eye(66, dtype=np.int64)-np.array(data.remote_counts, dtype=np.int64)
    for port in data.region:
        field = np.diag([state[port] for state in data.states])
        if np.any(d@field-field@d):
            raise ArithmeticError('an observer-local field no longer commutes with remote evolution')
    intervention = np.array(data.internal_counts, dtype=np.int64)
    if np.any(d@intervention-intervention@d):
        raise ArithmeticError('the internal native repair no longer commutes with remote evolution')
    outside = next(port for port in range(12) if port not in data.region)
    field = np.diag([state[outside] for state in data.states])
    commutator = d@field-field@d
    witnesses = np.argwhere(commutator != 0)
    if not len(witnesses):
        raise ArithmeticError('the nonlocal dictionary negative control disappeared')
    row, column = map(int, witnesses[0])
    return {'local_field_ports': list(data.region), 'all_local_field_projectors_commute_with_D': True,
            'internal_native_intervention_commutes_with_D': True,
            'all_record_adaptive_cylinders_in_this_dictionary_preserved': True,
            'stationary_uniform_boundary_required': True,
            'outside_port': outside, 'outside_commutator_entry': [row, column],
            'outside_commutator_exact': F(int(commutator[row, column]), data.denominator)}


def movement_clock_control():
    """The existing actual-change clock has an exact minus-one eigenmode."""
    port = movement.CLOCK_PORT
    source_clock = movement.source('C')
    for mark, _, _, _, mapping in source_clock.groups:
        for initial, final in mapping:
            before, after = quantum.STATES[initial][port], quantum.STATES[final][port]
            if (after == before) != (mark == 0):
                raise ArithmeticError('actual-change tick no longer flips the binary clock field')
    # Reuse the source's exact finite all-input survival proof to justify
    # completion; the sign identity then passes through the entire CP series.
    tail = movement.tail_certificate((2,))[0]['maximum_timeout_probability']
    if F(tail) >= 1:
        raise ArithmeticError('the movement-clock completed map lacks a uniform escape bound')
    return {'clock_port': port, 'eigenobservable': '1-2*n_clock',
            'completed_adjoint_eigenvalue': '-1', 'two_attempt_timeout_upper': tail,
            'same_reading_does_not_define_a_positive_Euclidean_transfer': True}


def holding_control():
    """Independent identity holds cancel from the complete local clock."""
    data = source()
    a = sp.Matrix(data.remote_counts)/data.denominator
    b = sp.Matrix(data.incident_counts)/data.denominator
    activity = sp.Rational(1, 2)
    held_a = (1-activity)*sp.eye(66)+activity*a
    held_b = activity*b
    d = sp.eye(66)-a
    if sp.eye(66)-held_a != activity*d or held_b != activity*b:
        raise ArithmeticError('independent holding no longer rescales both stopped factors equally')
    # D_a=aD cancels exactly in both B_a D_a^-1 and
    # D_a^-1/2 B_a D_a^-1/2. No refitted Hamiltonian is needed.
    return {'tested_activity': F(activity), 'D_activity_equals_activity_times_D': True,
            'complete_kernel_unchanged': True, 'symmetric_transfer_and_Hamiltonian_unchanged': True,
            'all_compatible_field_and_feedback_laws_unchanged': True,
            'finite_attempt_caps_and_attempt_costs_unchanged': False}


def scalar_remote_control():
    """Deleting real remote motion changes continuing field predictions."""
    data = source()
    wrong = IntegerKernel(data.incident_counts, 2*len(data.incident_edges), 'scalar-remote-control')
    law = _law_from_kernel(READ_PORTS, False, wrong)
    difference = total_variation(field_law(), law)
    if difference <= 0:
        raise ArithmeticError('deleting remote field evolution stopped changing the continued record law')
    return {'incorrect_remote_map': '(1-p)I', 'incorrect_completed_map': 'B/p',
            'four_read_record_total_variation': difference,
            'actual_remote_evolution_cannot_be_dropped': True}


def summary():
    """Compact exact evidence; terminal populations stay available via API."""
    data = source()
    baseline, feedback = field_law(), field_law(feedback=True)
    difference = total_variation(baseline, feedback)
    if difference <= 0:
        raise ArithmeticError('record-controlled native feedback no longer changes the continuing field law')
    return {'source_digest': data.source_digest, 'states': len(data.states),
            'region': list(data.region), 'incident_seams': len(data.incident_edges),
            'remote_seams': len(data.remote_edges), 'tick_probability': data.tick_probability,
            'spectral_certificate': spectral_certificate(), 'dictionary': dictionary_certificate(),
            'read_ports': list(READ_PORTS), 'complete_record_count': len(baseline.records),
            'initial_local_label_retained': True,
            'baseline': field_statistics(baseline), 'feedback': field_statistics(feedback),
            'baseline_feedback_total_variation': difference,
            'complete_baseline_terminal_counts_sha256': _digest([baseline.records, baseline.counts, baseline.denominator]),
            'complete_feedback_terminal_counts_sha256': _digest([feedback.records, feedback.counts, feedback.denominator]),
            'feedback_is_one_extra_native_incident_attempt': True,
            'feedback_read_clock_if_first_bit_one': list(read_clock((1, 0, 0, 0), feedback=True)),
            'movement_clock_control': movement_clock_control(),
            'holding_control': holding_control(), 'scalar_remote_control': scalar_remote_control(),
            'diagonal_source_sector_only': True,
            'native_event_map_is_unitary_continuation': False,
            'physical_clock_duration_or_relativistic_field_identified': False}
