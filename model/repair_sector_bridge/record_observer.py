"""Bounded distance records reconstructed from addressed native currents.

A known binary preparation and the uniform-seam scheduler are declared.
The controller receives no working loads. Its tracked state is reconstructed
from the preparation label and the complete ordered classical record stream.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path

from .bridge import certificate, lump_kernel, normalized, stochastic


@lru_cache(maxsize=1)
def tables():
    from model.observer import native, quantum_mixing
    states = native.states(binary=True)
    seams = tuple(quantum_mixing.seams())
    adjacency = {i: set() for i in range(12)}
    for a, b in seams:
        adjacency[a].add(b)
        adjacency[b].add(a)
    distances = []
    for state in states:
        a, b = (i for i, value in enumerate(state) if value)
        frontier, seen = {a}, {a}
        for distance in range(12):
            if b in frontier:
                distances.append(distance)
                break
            frontier = {v for u in frontier for v in adjacency[u]} - seen
            seen |= frontier
        else:
            raise ArithmeticError('disconnected carrier')
    return states, seams, tuple(distances)


LETTERS = ((0, 0), (0, 1), (1, 0), (1, -1))


def encode(address, coin, current):
    if (type(address) is not int or not 0 <= address < 30
            or type(coin) is not int or type(current) is not int
            or (coin, current) not in LETTERS):
        raise ValueError('native binary seam/current letter required')
    return 4 * address + LETTERS.index((coin, current))


@dataclass(frozen=True)
class Observer:
    initial_index: int
    state_index: int
    capacity: int
    tape: tuple[int, ...] = ()
    status: str = 'active'

    @classmethod
    def prepare(cls, initial_index, capacity):
        if type(initial_index) is not int or not 0 <= initial_index < 66:
            raise ValueError('known binary preparation label in 0..65 required')
        if type(capacity) is not int or capacity < 0:
            raise ValueError('nonnegative finite record capacity required')
        return cls(initial_index, initial_index, capacity,
                   status='full' if capacity == 0 else 'active')

    def consume(self, sequence, letter):
        """Read only an addressed coin/current letter and its sequence number.

An error halts tracking. A full controller requires the source protocol to
halt as well; its stored distance refers to the last admitted event.
"""
        if self.status != 'active':
            return self
        if type(sequence) is not int or sequence != len(self.tape):
            return replace(self, status='sequence_error')
        if type(letter) is not int or not 0 <= letter < 120:
            return replace(self, status='invalid_record')
        states, seams, _ = tables()
        address, kind = divmod(letter, 4)
        coin, current = LETTERS[kind]
        a, b = seams[address]
        memory = states[self.state_index]
        total = memory[a] + memory[b]
        low, high = total // 2, (total + 1) // 2
        expected_a = low if coin == 0 else high
        if current != memory[a] - expected_a:
            return replace(self, status='invalid_record')
        updated = list(memory)
        updated[a] -= current
        updated[b] += current
        index = states.index(tuple(updated))
        tape = self.tape + (letter,)
        return replace(self, state_index=index, tape=tape,
                       status='full' if len(tape) == self.capacity else 'active')

    def distance(self):
        if self.status not in ('active', 'full'):
            return None
        return tables()[2][self.state_index]


def memory_bound(capacity):
    if type(capacity) is not int or capacity < 0:
        raise ValueError('nonnegative finite record capacity required')
    counter = capacity.bit_length()  # ceil(log2(capacity+1))
    retained = 7 * capacity + 7 + 7 + counter + 3 + 2
    return {
        'record_capacity': capacity,
        'addressed_record_alphabet': 120,
        'bits_per_addressed_record': 7,
        'retained_native_tape_bits': 7 * capacity,
        'preparation_label_bits': 7,
        'tracked_state_bits': 7,
        'counter_bits': counter,
        'status_register_bits': 3,
        'current_distance_bits': 2,
        'retained_register_bits_upper': retained,
        'incoming_packet_workspace_bits_upper': 7 + counter,
        'total_including_incoming_packet_bits_upper': retained + 7 + counter,
        'optional_full_distance_tape_bits': 2 * capacity,
        'scope': 'logical classical registers; immutable program/geometry and physical implementation costs are separate',
    }


def dobrushin(kernel):
    kernel = stochastic(kernel)
    return max(sum(abs(a - b) for a, b in zip(left, right)) / 2
               for left in kernel for right in kernel)


def mixing_horizon(kernel, weights, tolerance, max_steps=1000):
    """Exact worst-initial-state TV, without sampling or a spectral fit."""
    kernel = stochastic(kernel)
    pi = normalized(weights, len(kernel))
    if not certificate(kernel, weights)['stationary']:
        raise ValueError('stationary target required')
    if type(tolerance) not in (int, F) or not 0 < tolerance < 1:
        raise ValueError('exact tolerance in (0,1) required')
    if type(max_steps) is not int or max_steps < 0:
        raise ValueError('nonnegative finite step limit required')
    power = tuple(tuple(F(i == j) for j in range(len(kernel))) for i in range(len(kernel)))
    previous = None
    for step in range(max_steps + 1):
        errors = tuple(sum(abs(row[j] - pi[j]) for j in range(len(kernel))) / 2 for row in power)
        worst = max(errors)
        if worst <= tolerance:
            return {'attempts': step, 'worst_TV_exact': str(worst),
                    'previous_worst_TV_exact': str(previous),
                    'per_initial_class_TV_exact': [str(x) for x in errors],
                    'tolerance': str(tolerance)}
        previous = worst
        power = tuple(tuple(sum(power[i][k] * kernel[k][j] for k in range(len(kernel)))
                            for j in range(len(kernel))) for i in range(len(kernel)))
    raise ValueError('requested tolerance not certified within the declared step limit')


def report():
    from model.observer import native
    states, seams, distances = tables()
    index = {state: i for i, state in enumerate(states)}
    kernel = [[F(0) for _ in states] for _ in states]
    checked = 0
    for initial, state in enumerate(states):
        for address, seam in enumerate(seams):
            for coin in (0, 1):
                output, (_, current) = native.native(state, seam, coin)
                letter = encode(address, coin, current)
                observer = Observer.prepare(initial, 1).consume(0, letter)
                if observer.status != 'full' or states[observer.state_index] != output:
                    raise ArithmeticError('native current reconstruction differs')
                if observer.distance() != distances[index[output]]:
                    raise ArithmeticError('record-only distance differs')
                kernel[initial][index[output]] += F(1, 60)
                checked += 1
    labels, quotient = lump_kernel(kernel, distances)
    weights = [distances.count(label) for label in labels]
    horizon = mixing_horizon(quotient, weights, F(1, 100))
    return {
        'schema': 'oph.native-current-distance-observer.v1',
        'exhaustive_native_branch_checks': checked,
        'known_preparation_required': True,
        'working_source_loads_exposed_to_controller': False,
        'scheduler': 'independent uniform choice of thirty addressed seams and fair native coin',
        'distance_labels': list(labels),
        'quotient_kernel': [[str(x) for x in row] for row in quotient],
        'stationary_distance_law': [str(x) for x in normalized(weights, len(weights))],
        'dobrushin_one_attempt_contraction': str(dobrushin(quotient)),
        'one_percent_relaxation': horizon,
        'memory_at_one_percent_horizon': memory_bound(horizon['attempts']),
        'source_clock': 'native attempt count, with source halted at declared tape capacity',
        'physical_temperature_or_particle_identity_selected': False,
    }


def receipt():
    import hashlib
    from model.particle_common_contract.provenance import import_closure
    root = Path(__file__).resolve().parents[2]
    seeds = {
        'model/repair_sector_bridge/record_observer.py',
        'model/repair_sector_bridge/test_record_observer.py',
        'model/repair_sector_bridge/bridge.py',
        'model/repair_sector_bridge/test_bridge.py',
        'model/repair_sector_bridge/PROOF.md',
        'model/repair_sector_bridge/RECORD_OBSERVER.md',
        'model/repair_sector_bridge/README.md',
    }
    value = report()
    value['input_sha256'] = {
        name: hashlib.sha256((root / name).read_bytes()).hexdigest()
        for name in sorted(import_closure(seeds))
    }
    return value


def main():
    import argparse
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    value = receipt()
    path = Path(__file__).with_name('record_observer_results.json')
    encoded = json.dumps(value, indent=2, sort_keys=True) + '\n'
    if args.write:
        path.write_text(encoded)
    elif args.check:
        if not path.is_file() or path.read_text() != encoded:
            raise SystemExit('Current-record observer receipt differs from source replay')
        print('Current-record observer: 3960 branches, 39 attempts, source hashes verified')
    else:
        print(encoded, end='')


if __name__ == '__main__':
    main()
