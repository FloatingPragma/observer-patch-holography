"""A literal-native current-tape winding meter, with finite reset and timeout.

No coherent pulse, target matrix, analog phase or ideal pointer is used. The
meter's one-load carrier starts at an arbitrary unknown port. Its controller
uses only actual native coin/current records and the fixed carrier seams.
"""
from collections import deque
from fractions import Fraction as F
from itertools import cycle
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
from native_kernel import native
CARRIER_SEAMS = tuple(tuple(e) for e in json.loads((HERE/'source_packet.json').read_text())['seams'])

EDGES = tuple(CARRIER_SEAMS)
CYCLE = (0, 1, 8, 6, 2, 10)


def route(start, finish):
    previous = {start: None}
    todo = deque([start])
    while todo:
        v = todo.popleft()
        if v == finish:
            break
        for edge in EDGES:
            if v not in edge:
                continue
            w = edge[0] if edge[1] == v else edge[1]
            if w not in previous:
                previous[w] = v
                todo.append(w)
    path = []
    current = finish
    while current != start:
        path.append(current)
        current = previous[current]
    return tuple(reversed(path))


def execute(initial_port, signed_turns, coins, retries=16):
    """Literal complete-path executor; final hidden loads are not policy inputs."""
    if type(initial_port) is not int or not 0 <= initial_port < 12:
        raise ValueError('a valid one-load source column is required')
    if type(retries) is not int or not 1 <= retries <= 128:
        raise ValueError('bounded positive retry count required')
    signed_turns = tuple(signed_turns)
    if any(type(x) is not int for x in signed_turns):
        raise ValueError('integer native tape entries required')
    maximum_steps = 6*sum(abs(x) for x in signed_turns)
    if maximum_steps > 4096:
        return {'success': False, 'reason': 'meter_capacity', 'records': [],
                'final_loads': [int(v == initial_port) for v in range(12)]}
    loads = tuple(int(v == initial_port) for v in range(12))
    records = []
    iterator = iter(coins)
    known_port = None

    def attempt(edge, stage, requested_direction=None):
        nonlocal loads
        coin = next(iterator)
        if type(coin) is not int or coin not in (0, 1):
            raise ValueError('literal native fair-coin alphabet required')
        loads, (_, current) = native(loads, edge, coin)
        records.append({'stage': stage, 'seam': list(edge), 'coin': coin,
                        'current': current,
                        'requested_cycle_direction': requested_direction})
        return current

    def finish(success, reason, step_counter=0):
        return {'success': success, 'reason': reason, 'records': records,
                'final_loads': list(loads), 'retained_signed_cycle_steps': step_counter,
                'winding': step_counter//6 if success else None,
                'fine_history_probability': str(F(1, 2**len(records)))}

    # A nonzero current identifies the endpoint reached, with no hidden read.
    for edge in EDGES:
        for _ in range(retries):
            current = attempt(edge, 'discover')
            if current:
                known_port = edge[1] if current > 0 else edge[0]
                break
        if known_port is not None:
            break
    if known_port is None:
        return finish(False, 'discovery_timeout')

    def move_to(destination, stage, direction=None):
        nonlocal known_port
        edge = tuple(sorted((known_port, destination)))
        assert edge in EDGES
        expected = 1 if known_port == edge[0] else -1
        for _ in range(retries):
            current = attempt(edge, stage, direction)
            if current:
                assert current == expected
                known_port = destination
                return True
        return False

    for destination in route(known_port, CYCLE[0]):
        if not move_to(destination, 'reset_route'):
            return finish(False, 'routing_timeout')
    cursor = 0
    signed_steps = 0
    for turns in signed_turns:
        direction = 1 if turns > 0 else -1
        for _ in range(6*abs(turns)):
            next_cursor = (cursor+direction) % 6
            if not move_to(CYCLE[next_cursor], 'measure', direction):
                return finish(False, 'measurement_timeout', signed_steps)
            cursor = next_cursor
            signed_steps += direction
    assert cursor == 0 and known_port == CYCLE[0] and sum(loads) == 1
    assert signed_steps == 6*sum(signed_turns)
    return finish(True, 'complete_closed_pointer_word', signed_steps)


def reset_success(initial_port, retries):
    """Exact law from first nonzero-current block and subsequent routing."""
    hold = F(1, 2**retries)
    success = F()
    earlier_incident_blocks = 0
    for a, b in EDGES:
        if initial_port not in (a, b):
            continue
        destination = b if initial_port == a else a
        path_length = len(route(destination, CYCLE[0]))
        success += hold**earlier_incident_blocks*(1-hold)**(1+path_length)
        earlier_incident_blocks += 1
    assert earlier_incident_blocks == 5
    return success


def execute_class(initial_port, regional_class, coins, retries=16):
    """Finite controller input; unused class labels produce explicit refusal."""
    if type(regional_class) is not int or regional_class not in range(6):
        raise ValueError('six-state retained regional counter required')
    lift = {0: 0, 1: 1, 5: -1}
    if regional_class not in lift:
        return {'success': False, 'reason': 'counter_outside_primary_source',
                'counter_class': regional_class, 'records': [],
                'final_loads': [int(v == initial_port) for v in range(12)],
                'winding': None}
    return execute(initial_port, (lift[regional_class],), coins, retries)


def independent_reset_distribution(initial_port, retries):
    """Second computation from the actual two-outcome native transition."""
    loads = tuple(int(v == initial_port) for v in range(12))
    unresolved = F(1)
    success = F()
    hold = F(1, 2**retries)
    for edge in EDGES:
        for _ in range(retries):
            next_unresolved = F()
            for coin in (0, 1):
                updated, (_, current) = native(loads, edge, coin)
                weight = unresolved/2
                if current:
                    destination = updated.index(1)
                    success += weight*(1-hold)**len(route(destination, CYCLE[0]))
                else:
                    assert updated == loads
                    next_unresolved += weight
            unresolved = next_unresolved
    assert unresolved == hold**5
    return success, unresolved


def favorable_coins(initial_port, signed_turns, retries):
    """A deterministic complete replay tape, not a probability sample."""
    loads = tuple(int(v == initial_port) for v in range(12))
    known = None
    for edge in EDGES:
        for _ in range(retries):
            possible = [(c, *native(loads, edge, c)) for c in (0, 1)]
            c, updated, (_, current) = next((x for x in possible if x[2][1]), possible[0])
            yield c
            loads = updated
            if current:
                known = loads.index(1)
                break
        if known is not None:
            break
    destinations = list(route(known, CYCLE[0]))
    cursor = 0
    for turns in signed_turns:
        direction = 1 if turns > 0 else -1
        for _ in range(6*abs(turns)):
            cursor = (cursor+direction) % 6
            destinations.append(CYCLE[cursor])
    for destination in destinations:
        edge = tuple(sorted((loads.index(1), destination)))
        coin = 0 if loads[edge[0]] else 1
        yield coin
        loads, (_, current) = native(loads, edge, coin)
        assert current and loads[destination] == 1


def bounded_source_counter():
    """Primary source: one actual load, known G, zero initial current tape.

    The pair (G,current load) has 144 states across the twelve declared
    preparations. Its charge image has 133 classes; the zero class does not
    replace the retained background label G in the source process.
    """
    lift = {0: 0, 1: 1, 5: -1}
    events = 0
    cuts = 0
    classes = set()
    for origin in range(12):
        background = tuple(int(v == origin) for v in range(12))
        for occupied in range(12):
            loads = tuple(int(v == occupied) for v in range(12))
            degree = tuple(g-n for g, n in zip(background, loads))
            counter = tuple(q % 6 for q in degree)
            classes.add(degree)
            assert tuple(lift[x] for x in counter) == degree
            for region in range(1 << 12):
                actual = sum(degree[v] for v in range(12) if region >> v & 1)
                measured_class = sum(counter[v] for v in range(12) if region >> v & 1) % 6
                assert lift[measured_class] == actual
                cuts += 1
            for seam in EDGES:
                for coin in (0, 1):
                    updated, (_, current) = native(loads, seam, coin)
                    assert sum(updated) == 1 and set(updated) <= {0, 1}
                    updated_counter = list(counter)
                    updated_counter[seam[0]] = (updated_counter[seam[0]]+current) % 6
                    updated_counter[seam[1]] = (updated_counter[seam[1]]-current) % 6
                    assert tuple(lift[x] for x in updated_counter) == tuple(g-n for g, n in zip(background, updated))
                    events += 1
    assert len(classes) == 133 and events == 8640 and cuts == 589824
    # The finite decoder must not be silently extended to larger charges.
    assert 7 % 6 == 1 and lift[7 % 6] != 7
    return {'source': 'unchanged one-load native carrier with retained known G=n0 and zero initial tape',
            'source_states_including_background': 144, 'distinct_charge_classes': 133,
            'native_events_checked': events, 'regional_cuts_checked': cuts,
            'counter_initialization': 'all twelve counters zero at actual known one-load preparation',
            'counter_update': 'c_a=(c_a+j) mod6; c_b=(c_b-j) mod6 using actual native current',
            'regional_decoder': 'sum_{v in A} c_v mod6 lifted by 0->0,1->1,5->-1',
            'source_scope_exact': 'single load remains single under every native coin branch; G is retained and fixed',
            'larger_charge_alias_control': {'native_integer': 7, 'class': 1, 'bounded_decoder': 1,
                                           'outside_primary_source_domain': True},
            'bounded_meter_command': 'one signed full six-edge traversal for decoded +/-1; zero traversals for0',
            'maximum_measurement_hops': 6}


if __name__ == '__main__':
    raise SystemExit('Run run.py for the portable primary experiment.')
