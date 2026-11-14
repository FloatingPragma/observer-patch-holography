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
ROOT = HERE.parents[4]
RER = ROOT/'reverse-engineering-reality'
sys.path[:0] = [str(RER), str(ROOT/'openwave-tuning/candidate')]
from model.observer.native import native
from oph_openwave.geometry import CARRIER_SEAMS

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


def run(retries=16):
    assert all(tuple(sorted((a, b))) in EDGES for a, b in zip(CYCLE, CYCLE[1:]+CYCLE[:1]))
    assert max(len(route(v, 0)) for v in range(12)) == 3
    reset = [reset_success(v, retries) for v in range(12)]
    assert all(independent_reset_distribution(v, retries)[0] == reset[v] for v in range(12))
    hold = F(1, 2**retries)
    assert all(1-p <= hold**5+3*hold for p in reset)
    bounded = bounded_source_counter()
    bounded_timeout = 1-min(reset)*(1-hold)**6
    bounded_examples = []
    for regional_class, winding in ((0, 0), (1, 1), (5, -1)):
        for initial in (0, 5, 11):
            replay = execute_class(initial, regional_class,
                                   favorable_coins(initial, (winding,), retries), retries)
            assert replay['success'] and replay['winding'] == winding
            bounded_examples.append({'regional_counter_class': regional_class,
                                     'unknown_initial_pointer_port': initial,
                                     'complete_execution': replay})
    rejected_counter_labels = [execute_class(0, label, cycle((0, 1)), retries) for label in (2, 3, 4)]
    assert all(not row['success'] and row['winding'] is None for row in rejected_counter_labels)
    source_path = HERE.parent/'spatial/NATIVE_DUAL_WINDING.json'
    source = json.loads(source_path.read_text())
    source_hash = hashlib.sha256(source_path.read_bytes()).hexdigest()
    examples = []
    worst_timeout = F()
    cases = 0
    for seed in source['seed_records']:
        for vertex in range(12):
            turns = [int(tape)*(1 if vertex == a else -1)
                     for (a, b), tape in zip(EDGES, seed['tape']) if vertex in (a, b)]
            commanded_steps = 6*sum(abs(x) for x in turns)
            expected = seed['degree'][vertex]
            assert sum(turns) == expected
            for initial in range(12):
                complete_success = reset[initial]*(1-hold)**commanded_steps
                worst_timeout = max(worst_timeout, 1-complete_success)
                cases += 1
            if len(examples) < 6 and (expected or sum(abs(x) for x in turns)):
                initial = (vertex+7) % 12
                replay = execute(initial, turns, favorable_coins(initial, turns, retries), retries)
                assert replay['success'] and replay['winding'] == expected
                examples.append({'source_deficit': seed['deficit_port'],
                                 'source_excess': seed['excess_port'], 'contour_vertex': vertex,
                                 'signed_turns': turns, 'unknown_initial_pointer_port': initial,
                                 'complete_execution': replay})
    # A timeout is a retained failure, never a guessed zero winding.
    fixed_initial = 0
    def all_holds():
        for edge in EDGES:
            for _ in range(retries):
                yield 1 if fixed_initial == edge[0] else 0
    timeout = execute(fixed_initial, (1,), all_holds(), retries)
    assert not timeout['success'] and timeout['winding'] is None
    assert timeout['reason'] == 'discovery_timeout'
    # The meter cannot turn a signed half-label into a native source word.
    try:
        execute(0, (F(1, 2),), cycle((0, 1)), retries)
    except ValueError:
        half_rejected = True
    else:
        half_rejected = False
    assert half_rejected
    assert source_hash == hashlib.sha256(source_path.read_bytes()).hexdigest()
    files = [Path(__file__), source_path, RER/'model/observer/native.py',
             ROOT/'openwave-tuning/candidate/oph_openwave/geometry.py']
    return {'schema': 'oph.row1.literal_native_cycle_meter.v1',
            'checks': {'actual_six_edge_carrier_cycle': True,
                       'unknown_one_load_reset_has_exact_complete_law': True,
                       'independent_actual_branch_reset_reconstruction': True,
                       'all_spatial_seed_contours_and_all_pointer_inputs_covered': cases == 132*12*12,
                       'successful_native_pointer_winding_equals_spatial_degree': True,
                       'timeout_retained_without_decoded_charge': True,
                       'fractional_command_rejected': half_rejected,
                       'bounded_mod6_class_meter_exact_on_primary_source': True,
                       'source_receipt_unchanged': True},
            'carrier_cycle': list(CYCLE)+[CYCLE[0]], 'retries_per_commanded_edge': retries,
            'exact_reset_success_by_unknown_initial_port': list(map(str, reset)),
            'discovery_timeout': str(hold**5),
            'per_commanded_step_success': str(1-hold),
            'complete_success_formula': 'p_reset(initial_port)*(1-2^-R)^(6 sum_e |B[v,e] F_e|)',
            'complete_timeout_formula': '1 minus complete_success; all early stops retain their tape and current pointer state',
            'worst_complete_timeout_for132seeds': str(worst_timeout),
            'worst_complete_timeout_float': float(worst_timeout),
            'bounded_primary_source_counter': bounded,
            'compressed_primary_meter_worst_timeout': str(bounded_timeout),
            'compressed_primary_meter_worst_timeout_float': float(bounded_timeout),
            'compressed_primary_meter_attempt_cap': 39*retries,
            'compressed_primary_meter_examples': bounded_examples,
            'retained_counter_complement': rejected_counter_labels,
            'cases': cases, 'reset_attempt_cap': 33*retries,
            'measurement_attempt_cap': '6R sum_e |B[v,e] F_e|',
            'all_path_semantics': {'executor': 'execute', 'coin_probability': '1/2 independently per executed attempt',
                'fine_history_probability': '2^-number_of_executed_attempts',
                'records': ['stage', 'seam', 'coin', 'current', 'requested cycle direction'],
                'observed_controller_inputs': 'committed source tape F, known carrier geometry, own prior native current records; no hidden pointer load inspection',
                'complete_alphabet': 'success with integer winding, discovery/routing/measurement/capacity timeout with full current state and no decoded winding'},
            'replay_examples': examples, 'timeout_example': timeout,
            'scope': 'A classical native observer measures the canonical spatial class represented by its actual integer current tape. The pointer is apparatus; its cycle is not a new independently evolved field. The phase-field/cohomology interpretation is the separately checked spatial adapter.',
            'declared_resources': ['one additional native carrier in the one-load sector, at arbitrary unknown initial port',
                'finite observer tape and current-adaptive seam scheduler',
                'access to the committed source current records and its fixed dual-cell dictionary; bounded primary meter uses only twelve six-class counters plus retained background'],
            'not_used': ['coherent proper-code transfer', 'source M6 analog pulses', 'ideal projective pointer', 'known pointer-position preparation', 'discarded no-move outcomes', 'postselection without timeout record'],
            'source_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
            'official_rows_changed': []}


if __name__ == '__main__':
    result = run()
    (HERE/'NATIVE_CYCLE_METER.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({key: result[key] for key in ['checks', 'cases', 'worst_complete_timeout_float', 'reset_attempt_cap']}))
