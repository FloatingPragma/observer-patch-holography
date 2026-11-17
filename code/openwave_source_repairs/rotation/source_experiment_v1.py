#!/usr/bin/env python3
"""Complete finite classical source audit. Python standard library only."""
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import argparse
import gzip
import hashlib
import itertools
import json
import math
import random
import time

from native_endpoint import native

HERE = Path(__file__).resolve().parent
SEAMS = tuple(map(tuple, json.loads((HERE / 'seams.json').read_text())))
REF = (12, 13)
SEARCH_CAP, SELECTOR_CAP, HISTORY_LENGTH = 256, 16, 8
SIGNS = (-1, 1)
CALLS = 0


def canonical(obj):
    return json.dumps(obj, sort_keys=True, separators=(',', ':')) + '\n'


def frac(value):
    return str(F(value))


def sha(data):
    return hashlib.sha256(data).hexdigest()


def preparation(sign, port):
    return tuple(sign * int(k == port) for k in range(14))


def read(record, seam, current):
    result = list(record)
    a, b = seam
    result[a] -= current
    result[b] += current
    return tuple(result)


def body(q, edge, coin, law):
    global CALLS
    if law == 'native':
        actual, (_, j) = native(q, edge, coin)
        CALLS += 1
    else:
        a, b = edge
        j = q[a] - q[b] if law == 'swap' else 0
        actual = read(q, edge, j)
    observed = read(q, edge, j)
    assert actual == observed and sum(actual) == sum(q)
    return observed, j


def position(q):
    occupied = [k for k, value in enumerate(q) if value]
    assert len(occupied) == 1 and occupied[0] < 12
    return occupied[0]


def tv(a, b):
    return sum(abs(a.get(k, F(0)) - b.get(k, F(0))) for k in set(a) | set(b)) / 2


def native_columns():
    rows = []
    discovery = {}
    for sign, port in itertools.product(SIGNS, range(12)):
        q = preparation(sign, port)
        moves = []
        zeros = []
        for edge_id, edge in enumerate(SEAMS):
            for coin in (0, 1):
                out, current = body(q, edge, coin, 'native')
                row = [sign, port, edge_id, coin, current, position(out)]
                rows.append(row)
                (moves if current else zeros).append(row)
        assert len(moves) == 5 and len(zeros) == 55
        assert len({r[-1] for r in moves}) == 5
        discovery[(sign, port)] = moves
    # Vacuum and blank reference columns are retained separately, not excluded.
    vacuum = (0,) * 14
    for edge_id, edge in enumerate(SEAMS + (REF,)):
        for coin in (0, 1):
            out, j = body(vacuum, edge, coin, 'native')
            assert out == vacuum and j == 0
            rows.append([0, -1, edge_id, coin, j, -1])
    return rows, discovery


def discovery_law():
    # Five fair reference coins produce 32 labels. Two are retained retries.
    selector_fail = F(1, 16) ** SELECTOR_CAP
    accepted = 1 - selector_fail
    unresolved = accepted * F(11, 12)
    move_rows, fail_rows = [], []
    for t in range(1, SEARCH_CAP + 1):
        reach = unresolved ** (t - 1)
        move_rows.append([t, frac(reach * accepted / 60)])
        fail_rows.append([t, frac(reach * selector_fail)])
    one_neighbour = sum(F(p) for _, p in move_rows)
    selector_total = sum(F(p) for _, p in fail_rows)
    timeout = unresolved ** SEARCH_CAP
    assert 5 * one_neighbour + selector_total + timeout == 1
    return {'selector_cap': SELECTOR_CAP, 'search_cap': SEARCH_CAP,
            'per_directed_first_move_by_time': move_rows,
            'selector_failure_by_time': fail_rows,
            'per_neighbour_mass': frac(one_neighbour),
            'discovery_mass': frac(5 * one_neighbour),
            'selector_failure_mass': frac(selector_total),
            'search_timeout_mass': frac(timeout),
            'discovery_mass_decimal': float(5 * one_neighbour),
            'selector_failure_mass_decimal': float(selector_total),
            'search_timeout_mass_decimal': float(timeout),
            'complete_probability': '1',
            'history_encoding': 'SOURCE_OUTCOMES.md: exact finite prefix DAG'}


def suffixes(discovery):
    transition_rows = []
    tests = []
    history_rows = []
    for (sign, origin), moves in discovery.items():
        for first in moves:
            edge = SEAMS[first[2]]
            arrival = first[-1]
            laws = {}
            for law in ('native', 'hold', 'swap'):
                table = {}
                for start in edge:
                    for coin in (0, 1):
                        out, current = body(preparation(sign, start), edge, coin, law)
                        table[(start, coin)] = (position(out), current)
                        transition_rows.append([sign, origin, arrival, law, start,
                                                coin, current, position(out)])
                distributions = {}
                for start in edge:
                    by_time = [Counter() for _ in range(HISTORY_LENGTH)]
                    for bits in itertools.product((0, 1), repeat=HISTORY_LENGTH):
                        current_port = start
                        currents = []
                        ports = []
                        for time_index, coin in enumerate(bits):
                            current_port, j = table[(current_port, coin)]
                            currents.append(j)
                            ports.append(current_port)
                            by_time[time_index][current_port] += 1
                        history_rows.append([sign, origin, arrival, law, start,
                                             ''.join(map(str, bits)), currents, ports])
                    distributions[start] = [{k: F(v, 2 ** HISTORY_LENGTH)
                                              for k, v in counts.items()}
                                             for counts in by_time]
                distances = [tv(distributions[edge[0]][t], distributions[edge[1]][t])
                             for t in range(HISTORY_LENGTH)]
                expected = F(0) if law == 'native' else F(1)
                assert distances == [expected] * HISTORY_LENGTH
                if law == 'native':
                    stationary = {edge[0]: F(1, 2), edge[1]: F(1, 2)}
                    assert all(p == stationary for d in distributions.values() for p in d)
                laws[law] = list(map(frac, distances))
            tests.append({'sign': sign, 'departure': origin, 'arrival': arrival,
                          'distance_between_opposite_starts_by_time': laws,
                          'native_vacuum_probability': '0'})
    return transition_rows, tests, history_rows


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def norm(a):
    return math.sqrt(dot(a, a))


def rotate(x, axis, theta):
    c, s = math.cos(theta), math.sin(theta)
    product = cross(axis, x)
    return tuple(c*x[k] + s*product[k] + (1-c)*dot(axis, x)*axis[k] for k in range(3))


def geometry():
    phi = (1 + math.sqrt(5)) / 2
    points = [(-1,phi,0),(1,phi,0),(-1,-phi,0),(1,-phi,0),
              (0,-1,phi),(0,1,phi),(0,-1,-phi),(0,1,-phi),
              (phi,0,-1),(phi,0,1),(-phi,0,-1),(-phi,0,1)]
    assert tuple((a,b) for a in range(12) for b in range(a+1,12)
                 if abs(norm(sub(points[a], points[b])) - 2) < 1e-12) == SEAMS
    records = []
    for edge in SEAMS:
        a, b = edge
        midpoint = tuple((x+y)/2 for x,y in zip(points[a], points[b]))
        axis = tuple(x/norm(midpoint) for x in midpoint)
        arrow = sub(points[b], points[a])
        rotated = [rotate(x, axis, math.pi) for x in points]
        permutation = []
        for x in rotated:
            distances = [norm(sub(x, p)) for p in points]
            assert min(distances) < 2e-14
            permutation.append(distances.index(min(distances)))
        assert permutation[a] == b and permutation[b] == a
        assert sorted(permutation) == list(range(12))
        assert all(permutation[permutation[i]] == i for i in range(12))
        assert {tuple(sorted((permutation[u], permutation[v]))) for u,v in SEAMS} == set(SEAMS)
        # These are actual 3D proper rotations, not a reader relabelling rule.
        matrix_columns = [rotate(tuple(int(i == k) for i in range(3)), axis, math.pi)
                          for k in range(3)]
        determinant = dot(matrix_columns[0], cross(matrix_columns[1], matrix_columns[2]))
        assert abs(determinant - 1) < 2e-14
        uniform = {a: F(1,2), b: F(1,2)}
        rotated_uniform = {permutation[k]: p for k,p in uniform.items()}
        native_half_distance = tv(uniform, rotated_uniform)
        polar = {a: F(1)}
        control_distance = tv(polar, {permutation[a]: F(1)})
        plant = {a: F(3,4), b: F(1,4)}
        plant_distance = tv(plant, {permutation[k]: p for k,p in plant.items()})
        assert native_half_distance == 0 and control_distance == 1 and plant_distance == F(1,2)
        half_arrow = rotate(arrow, axis, math.pi)
        full_arrow = rotate(arrow, axis, 2*math.pi)
        assert norm(sub(half_arrow, tuple(-x for x in arrow))) < 2e-14
        assert norm(sub(full_arrow, arrow)) < 2e-14
        finer = []
        for quarter in (1, 2, 3):
            angle = quarter * math.pi / 4
            support = [rotate(points[k], axis, angle) for k in edge]
            distances = [norm(sub(x, points[k])) for x in support for k in edge]
            assert min(distances) > .7
            finer.append({'angle_in_pi': frac(F(quarter,4)),
                          'field_TV': '1', 'minimum_support_separation': min(distances),
                          'frame_displacement': norm(sub(rotate(arrow,axis,angle),arrow))})
        records.append({'edge': list(edge), 'midpoint': midpoint, 'axis': axis,
                        'half_turn_port_permutation': permutation,
                        'proper_determinant': determinant,
                        'native_field_half_turn_TV': frac(native_half_distance),
                        'polar_control_half_turn_TV': frac(control_distance),
                        'asymmetric_plant_half_turn_TV': frac(plant_distance),
                        'frame_half_turn_displacement': norm(sub(half_arrow,arrow)),
                        'frame_full_turn_displacement': norm(sub(full_arrow,arrow)),
                        'finer_turns': finer})
    different = [tv({a:F(1,2),b:F(1,2)}, {c:F(1,2),d:F(1,2)})
                 for (a,b),(c,d) in itertools.combinations(SEAMS,2)]
    assert min(different) == F(1,2)
    return {'points': points, 'rotations': records,
            'distinct_bond_minimum_full_reader_TV': frac(min(different)),
            'pairs_of_distinct_bonds': len(different)}


class Source:
    """Opaque to the controller: only source.event returns native records."""
    def __init__(self, q, rng, law):
        self._q = q
        self.rng = rng
        self.law = law

    def event(self, edge):
        coin = self.rng.randrange(2)
        law = 'native' if edge == REF else self.law
        self._q, current = body(self._q, edge, coin, law)
        return coin, current

    def audit(self, estimate):
        assert estimate == self._q


def trajectory(sign, start, seed, law):
    initial = preparation(sign, start)
    source = Source(initial, random.Random(seed), law)
    estimated = initial
    records = []
    frame = None
    status = 'search_timeout'
    for t in range(1, SEARCH_CAP + 1):
        ticket = None
        for retry in range(SELECTOR_CAP):
            bits = []
            for _ in range(5):
                coin, j = source.event(REF)
                assert j == 0
                bits.append(coin)
                records.append(['selector', t, retry, coin, j])
            value = sum(bit << k for k,bit in enumerate(bits))
            if value < 30:
                ticket = value
                break
        if ticket is None:
            status = 'selector_timeout'
            break
        edge = SEAMS[ticket]
        departure = position(estimated)
        coin, j = source.event(edge)
        estimated = read(estimated, edge, j)
        source.audit(estimated)
        records.append(['search', t, ticket, coin, j, position(estimated)])
        if j:
            frame = (departure, position(estimated))
            status = 'discovered'
            break
    if frame:
        edge = tuple(sorted(frame))
        for t in range(64):
            coin, j = source.event(edge)
            estimated = read(estimated, edge, j)
            source.audit(estimated)
            assert position(estimated) in edge and sum(estimated) == sign
            records.append(['holding', t, SEAMS.index(edge), coin, j, position(estimated)])
    return {'sign': sign, 'start': start, 'seed': seed, 'law': law,
            'status': status, 'frame': frame, 'records': records}


def execute():
    global CALLS
    CALLS = 0
    native_rows, discovery = native_columns()
    law = discovery_law()
    transitions, tests, histories = suffixes(discovery)
    spatial = geometry()
    trajectories = [trajectory(sign, port, seed + 101*port + sign, candidate_law)
                    for seed in (730701, 730709) for sign in SIGNS for port in range(12)
                    for candidate_law in ('native', 'hold')]
    assert all(x['status'] == ('discovered' if x['law'] == 'native' else 'search_timeout')
               for x in trajectories)
    result = {'candidate': 'A-native-discovered-apolar-bond-v1',
              'signed_starting_positions': 24, 'vacuum_included': True,
              'signed_ordered_discovery_contexts': len(tests),
              'body_native_column_rows': len(native_rows),
              'suffix_transition_rows': len(transitions),
              'exact_suffix_histories': len(histories),
              'history_length': HISTORY_LENGTH, 'literal_native_calls': CALLS,
              'enumerated_history_transitions_are_not_extra_native_calls': True,
              'discovery': law, 'source_tests': tests, 'spatial': spatial,
              'controller_runs': [{'sign': x['sign'], 'start': x['start'],
                                   'seed': x['seed'], 'law': x['law'],
                                   'status': x['status'], 'frame': x['frame'],
                                   'event_count': len(x['records'])} for x in trajectories],
              'source_antipodal_identification': True,
              'coherent_memory_used': False,
              'source_axis_preserved': True,
              'spinor_phase_or_exchange_statistic_claimed': False,
              'native_rigid_rotation_apparatus_implemented': False,
              'spatial_coordinates_declared': True,
              'official_row_award': False,
              'g6_other_group_confirmed': False,
              'full_physical_spin_row_solved': False}
    records = {'native_columns': native_rows, 'suffix_transition_table': transitions,
               'suffix_histories': histories, 'controller_trajectories': trajectories}
    return result, records


def freeze_check():
    frozen = json.loads((HERE / 'PREREG.json').read_text())
    for name, digest in frozen['files'].items():
        assert sha((HERE / name).read_bytes()) == digest, name


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--record', action='store_true')
    args = parser.parse_args()
    freeze_check()
    started = time.perf_counter()
    result, records = execute()
    data = canonical(records).encode()
    result['records_sha256'] = sha(data)
    elapsed = time.perf_counter() - started
    if args.record:
        with (HERE / 'RECEIPT.json').open('x') as handle:
            handle.write(canonical(result))
        with (HERE / 'RECORDS.json.gz').open('xb') as handle:
            handle.write(gzip.compress(data, mtime=0))
        with (HERE / 'RUN.json').open('x') as handle:
            handle.write(canonical({'elapsed_seconds': elapsed, 'exit_code': 0}))
    else:
        assert json.loads((HERE / 'RECEIPT.json').read_text()) == result
        assert gzip.decompress((HERE / 'RECORDS.json.gz').read_bytes()) == data
    print(canonical({'candidate': result['candidate'], 'literal_native_calls': CALLS,
                     'exact_suffix_histories': result['exact_suffix_histories'],
                     'contexts': len(result['source_tests']),
                     'native_opposite_start_TV': '0', 'hold_and_swap_TV': '1',
                     'discovery_probability': result['discovery']['discovery_mass_decimal'],
                     'native_half_turn_field_TV': '0', 'asymmetric_plant_TV': '1/2',
                     'full_physical_spin_row_solved': False, 'seconds': elapsed}), end='')


if __name__ == '__main__':
    main()
