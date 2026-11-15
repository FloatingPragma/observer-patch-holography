"""Standalone native-history and spatial-winding adversarial reproducer.

Only the Python standard library is used. The frozen input packet supplies
geometry and complete already-executed development histories, not fitted
physical targets. Running this file does not award an OpenWave cell.
"""
import argparse
import cmath
from collections import deque
from fractions import Fraction
import hashlib
import itertools
import json
import math
from pathlib import Path

try:
    from .native_kernel import native
except ImportError:
    from native_kernel import native

HERE = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dual_boundaries(edges):
    """Reconstruct oriented pentagon cycles by vertex links, without a producer import."""
    edge_lookup = {tuple(edge): i for i, edge in enumerate(edges)}
    triangles = [f for f in itertools.combinations(range(12), 3)
                 if all(tuple(e) in edge_lookup for e in itertools.combinations(f, 2))]
    require(len(triangles) == 20, 'icosahedral triangle count')
    dual_edges = []
    for a, b in edges:
        adjacent = [i for i, f in enumerate(triangles) if a in f and b in f]
        require(len(adjacent) == 2, 'nonmanifold primal seam')
        dual_edges.append(tuple(adjacent))
    raw = []
    for port in range(12):
        incident = [i for i, face in enumerate(triangles) if port in face]
        require(len(incident) == 5, 'port link is not pentagonal')
        start = min(incident)
        walk = [start]
        while len(walk) < 5:
            candidates = [f for f in incident if f not in walk
                          and len(set(triangles[f]) & set(triangles[walk[-1]])) == 2]
            require(candidates, 'broken port link')
            walk.append(min(candidates))
        require(len(set(triangles[walk[-1]]) & set(triangles[start])) == 2,
                'port link does not close')
        raw.append(walk)

    def signature(walk):
        result = {}
        for first, second in zip(walk, walk[1:] + walk[:1]):
            pair = tuple(sorted((first, second)))
            edge = dual_edges.index(pair)
            result[edge] = 1 if first < second else -1
        return result

    oriented = {0: raw[0]}
    queue = deque([0])
    while queue:
        port = queue.popleft()
        sig = signature(oriented[port])
        for edge, (a, b) in enumerate(edges):
            if port not in (a, b):
                continue
            other = b if port == a else a
            if other in oriented:
                require(signature(oriented[other])[edge] == -sig[edge],
                        'inconsistent surface orientation')
            else:
                candidate = raw[other]
                if signature(candidate)[edge] == sig[edge]:
                    candidate = list(reversed(candidate))
                oriented[other] = candidate
                queue.append(other)
    require(len(oriented) == 12, 'disconnected surface')
    # Positive primal current is translated to a signed dual edge phase.
    dual_sign = [signature(oriented[a])[edge] for edge, (a, _) in enumerate(edges)]
    for port, walk in oriented.items():
        sig = signature(walk)
        for edge in sig:
            a, b = edges[edge]
            require(sig[edge] * dual_sign[edge] == int(port == a) - int(port == b),
                    'topology/current orientation map failed')
    require(12 - 30 + 20 == 2, 'Euler characteristic')
    return dual_edges, oriented, dual_sign


def curve(first, second, turns, deform):
    """A continuous representative with shared vertex phases and interior deformation."""
    gamma_a = .31 * math.sin(first + 1) if deform else 0.
    gamma_b = .31 * math.sin(second + 1) if deform else 0.
    amplitude = .23 * math.cos(first + second + 1) if deform else 0.
    derivative_bound = abs(gamma_b - gamma_a) + 2 * math.pi * (abs(turns) + abs(amplitude))
    subdivisions = max(16, int(math.ceil(8 * derivative_bound / math.pi)))
    values = [cmath.exp(1j * (gamma_a + (gamma_b - gamma_a) * k / subdivisions
                              + 2 * math.pi * turns * k / subdivisions
                              + amplitude * math.sin(2 * math.pi * k / subdivisions)))
              for k in range(subdivisions + 1)]
    require(abs(values[0] - cmath.exp(1j * gamma_a)) < 1e-12, 'initial gluing')
    require(abs(values[-1] - cmath.exp(1j * gamma_b)) < 1e-12, 'terminal gluing')
    return values


def winding(values):
    require(abs(values[0] - values[-1]) < 1e-10, 'open phase contour')
    require(min(map(abs, values)) > 1e-10, 'zero phase contour')
    angles = [cmath.phase(b * a.conjugate()) for a, b in zip(values, values[1:])]
    require(max(map(abs, angles), default=0.) < math.pi / 2, 'unresolved phase increments')
    return math.fsum(angles) / (2 * math.pi)


def evaluate(packet):
    edges = [tuple(e) for e in packet['seams']]
    dual_edges, contours, signs = dual_boundaries(edges)
    native_checks = 0
    for x, y, coin in itertools.product(range(5), range(5), (0, 1)):
        initial = (x, y) + (0,) * 10
        actual, record = native(initial, (0, 1), coin)
        # Independent balanced-endpoint form of the source law.
        total = x + y
        left = total // 2 + coin * (total % 2)
        require(actual == (left, total - left) + (0,) * 10 and record == (coin, x - left),
                'native kernel differs from balanced endpoint law')
        native_checks += 1
    cached = {}

    def degrees(tape, deform):
        answer = []
        for port in range(12):
            walk = contours[port]
            samples = []
            for first, second in zip(walk, walk[1:] + walk[:1]):
                pair = tuple(sorted((first, second)))
                edge = dual_edges.index(pair)
                key = (edge, tape[edge], deform)
                if key not in cached:
                    cached[key] = curve(pair[0], pair[1], signs[edge] * tape[edge], deform)
                values = cached[key] if first < second else list(reversed(cached[key]))
                require(not samples or abs(samples[-1] - values[0]) < 1e-10,
                        'shared vertex phase mismatch')
                samples += values if not samples else values[1:]
            answer.append(winding(samples))
        return answer

    histories = events = noops = 0
    max_error = 0.
    complete_outputs = []
    for setting, seed in enumerate(packet['spatial_result']['seed_records']):
        expected_length = len(seed['actual_history'])
        actual_coins = set()
        total_probability = Fraction(0)
        for history in seed['complete_prescribed_schedule_histories']:
            state = tuple(seed['initial_G'])
            tape = [0] * 30
            coins = []
            require(len(history['events']) == expected_length, 'history length')
            for position,event in enumerate(history['events']):
                require(event['seam']==seed['actual_history'][position]['seam'], 'history changed declared seam schedule')
                index, coin = event['seam_index'], event['coin']
                require(tuple(event['seam']) == edges[index], 'event seam')
                state, (_, current) = native(state, edges[index], coin)
                require(current == event['current'] and list(state) == event['output_loads'],
                        'incompatible native event')
                tape[index] += current
                events += 1
                noops += int(current == 0)
                coins.append(coin)
            require(tape == history['final_tape'] and list(state) == history['final_loads'],
                    'final tape or load mismatch')
            probability = Fraction(history['probability'])
            require(probability == Fraction(1, 2 ** expected_length), 'history probability')
            total_probability += probability
            require(tuple(coins) not in actual_coins, 'duplicate coin history')
            actual_coins.add(tuple(coins))
            target = [g - n for g, n in zip(seed['initial_G'], state)]
            ordinary = degrees(tuple(tape), False)
            deformed = degrees(tuple(tape), True)
            error = max(abs(x - y) for values in (ordinary, deformed) for x, y in zip(values, target))
            max_error = max(max_error, error)
            require(error < 1e-10, 'spatial winding differs from actual relative charge')
            complete_outputs.append({'setting': setting, 'coins': coins,
                                     'probability': str(probability),
                                     'relative_charge': [-x for x in target],
                                     'spatial_winding': ordinary,
                                     'deformed_spatial_winding': deformed})
            histories += 1
        require(total_probability == 1 and len(actual_coins) == 2 ** expected_length,
                'incomplete native history distribution')
    require(histories == packet['spatial_result']['complete_coin_histories'], 'history count')
    reject = {}
    for name, action in [('half_turn_gluing', lambda: curve(0, 1, .5, False)),
                         ('zero_crossing', lambda: winding([1., 0., 1.])),
                         ('nonclosed_loop', lambda: winding([1., 1j, -1.]))]:
        try:
            action()
        except ValueError as error:
            reject[name] = str(error)
        else:
            raise ValueError('mutation survived: ' + name)
    # All-zero corrupted memory would erase a genuine move; the readout must
    # then disagree with its source charge instead of silently repairing it.
    one = [0] * 30
    one[0] = 1
    correct, corrupted = degrees(tuple(one), True), degrees((0,) * 30, True)
    require(max(abs(a - b) for a, b in zip(correct, corrupted)) > .9,
            'forgotten current mutation survived')
    return {'schema': 'oph.row1.portable_spatial_review.v1',
            'checks': {'native_kernel_balanced_form': native_checks == 50,
                       'complete_native_histories': histories == 456,
                       'native_record_spatial_degree_identity': max_error < 1e-10,
                       'nonlinear_homotopy_and_shared_vertex_gauge_invariance': max_error < 1e-10,
                       'no_operation_events_retained': noops > 0,
                       'half_zero_open_and_forgotten_current_mutations_rejected': len(reject) == 3},
            'method': 'Independent dual pentagon rotation system; full copied native kernel; all raw histories; phase-only contour integration with nonconstant vertex phases and nonlinear edge deformations.',
            'geometry': {'dual_vertices': 20, 'dual_edges': 30,
                         'fixed_punctures': 12, 'first_Betti_number': 11},
            'histories': histories, 'events': events, 'noops': noops,
            'maximum_phase_degree_error': max_error, 'rejecting_controls': reject,
            'complete_outputs': complete_outputs,
            'scope': 'Relative charge q=n-G equals minus spatial phase winding on actual native histories. G is the retained preparation. Phase interpolation is a declared observer map, no external texture dynamics. Nonlinear homotopies retain its topological class. Electric units and official admission are not established by this replay.',
            'full_criterion_met': False, 'official_rows_changed': []}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads((HERE / 'manifest.json').read_text())
    for name, pin in manifest['files_sha256'].items():
        require(sha(HERE / name) == pin, 'bundle source mismatch: ' + name)
    packet = json.loads((HERE / 'source_packet.json').read_text())
    result = evaluate(packet)
    require(all(result['checks'].values()), 'portable spatial review failed')
    result['input_packet_sha256'] = sha(HERE / 'source_packet.json')
    result['executed_script_sha256'] = sha(Path(__file__))
    try:
        import openwave
        result['OpenWave_version'] = openwave.__version__
    except ImportError:
        result['OpenWave_version'] = None
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: result[key] for key in ('checks', 'histories', 'events', 'noops', 'maximum_phase_degree_error')}, indent=2))


if __name__ == '__main__':
    main()
