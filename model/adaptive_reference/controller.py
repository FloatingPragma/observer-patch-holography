"""A reversible, record-driven classical controller with no record erasure.

The two bits are working memory. Seed records and the subsequently appended
native records are different registers and are never recycled or discarded
by the execution. Stationarity refers to the reduced field/working-memory
state, after tracing that growing archive.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
from itertools import product

from ..observer import native


LABELS = tuple(('native', coin, raw) for coin in (0, 1) for raw in range(8)) + (('outside',),)
PORTS = (0, 1, 5)


def bits(value):
    if (type(value) is not tuple or len(value) != 2
            or any(type(x) is not int or x not in (0, 1) for x in value)):
        raise ValueError('a pair of classical bits is required')
    return value


def label(value):
    if type(value) is not tuple or value not in LABELS:
        raise ValueError('one complete native/outside record is required')
    if any(type(x) is not int for x in value[1:]):
        raise ValueError('record indices must be integers')
    return value


def selector(memory):
    c, _ = bits(memory)
    return ('C', (0, 1)) if c == 0 else ('E', (0, 5))


def feedback(record):
    record = label(record)
    return 0 if record == ('outside',) else record[1] ^ (record[2] & 1)


def advance(memory, record):
    c, d = bits(memory)
    return d ^ feedback(record), c


def retreat(memory, record):
    u, v = bits(memory)
    return v, u ^ feedback(record)


def permutation_certificate():
    memories = tuple(product(range(2), repeat=2))
    for r in LABELS:
        outputs = [advance(m, r) for m in memories]
        if len(set(outputs)) != 4 or any(retreat(advance(m, r), r) != m for m in memories):
            raise ArithmeticError('working-memory update is not reversible given the record')
        # All incoming terms at (u,v) use exactly the selector v, independent
        # of r. That stronger property, not mere reversibility, proves balance.
        for m in memories:
            if retreat(m, r)[0] != m[1]:
                raise ArithmeticError('incoming row mixes different instruments')
    return {'working_bits': 2, 'working_states': 4, 'complete_record_slots': 17,
            'forward': '(c,d),r -> (d XOR (coin XOR (raw mod2)),c); outside uses f=0',
            'inverse_given_record': '(u,v),r -> (v,u XOR f(r))',
            'all_68_memory_record_pairs_checked': True,
            'incoming_selector_at_output_u_v': 'v for every retained record',
            'controller_reads_field_loads_or_hidden_state': False,
            'working_memory_erasure_or_reset_per_tick': False,
            'supplied_classical_controller_law_not_quantum_source_compilation': True}


def seed_certificate():
    """Two exact native coins, with a stated unchanged spectator basis input.

    This is a simulator/native-law attachment. It does not claim preparation
    of the spectator as an interacting field state by an omitted pulse word.
    """
    loads = tuple(int(j in (2, 3)) for j in range(12))
    if loads not in native.states(binary=True):
        raise ArithmeticError('seed input is not an admitted native binary basis state')
    rows = []
    for c, d in product(range(2), repeat=2):
        first, r1 = native.native(loads, (0, 1), c)
        second, r2 = native.native(first, (0, 1), d)
        if second != loads or r1 != (c, 0) or r2 != (d, 0):
            raise ArithmeticError('native seed disturbed its spectator or biased its records')
        rows.append({'memory': [c, d], 'seed_records': [[c, 0], [d, 0]],
                     'probability': F(1, 4)})
    return {'input': list(loads), 'input_status': 'declared deterministic native simulator basis initialization',
            'seam': [0, 1], 'native_events': 2, 'uniform_memory': rows,
            'exact_native_coin_probability': F(1, 2),
            'full_C_raw_slots_per_event_including_zero': 8,
            'full_two_event_raw_archive_slots': 64,
            'nonzero_two_event_archive_slots': 4,
            'spectator_unchanged_on_every_branch': True,
            'seed_records_retained_and_correlated_with_working_memory': True,
            'seed_angle_or_floating_tolerance_assumed': False,
            'interacting_field_preparation_of_seed_spectator_claimed': False}


@dataclass(frozen=True)
class Checkpoint:
    memory: tuple[int, int]
    history: tuple[tuple, ...] = ()

    def __post_init__(self):
        bits(self.memory)
        if type(self.history) is not tuple or any(label(r) != r for r in self.history):
            raise ValueError('a complete immutable native history is required')

    def commit(self, record):
        record = label(record)
        return Checkpoint(advance(self.memory, record), self.history + (record,))

    def recover_initial_memory(self):
        result = self.memory
        for record in reversed(self.history):
            result = retreat(result, record)
        return result

    def export(self):
        return {'memory': list(self.memory), 'history': [list(r) for r in self.history]}

    @staticmethod
    def restore(value):
        if type(value) is not dict or set(value) != {'memory', 'history'}:
            raise ValueError('an exact controller checkpoint is required')
        return Checkpoint(tuple(value['memory']), tuple(tuple(r) for r in value['history']))
