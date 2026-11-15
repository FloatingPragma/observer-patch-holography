"""Independent stdlib oracle and adversarial audit of the row-1 law controls.

Only this directory receives outputs. The source was read before this oracle
was coded; independence here means a separate implementation, not a blind test.
No reserved acceptance configuration is read or executed.
"""
import ast
import cmath
import copy
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
import importlib.util
import itertools
import json
import math
from pathlib import Path
import random
import sys

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
RER = ROOT / 'reverse-engineering-reality'
NATIVE = RER / 'code/openwave_row1_native'
TARGET = RER / 'code/openwave_row1_protocol/law_controls.py'
RECEIPT = TARGET.with_name('LAW_CONTROLS.json')
SOURCE_COMMIT = '19288e2b53b680e244e8ec5f4c5d45e000b2ef51'
MANIFEST_PIN = 'd5e692b0597dbbeafabb88e6f7d417da8f7d7d1fd39a334c08fa73ffe40e8d90'
EPS = 1e-10
K = 32
TURN_BOUND = 3


def demand(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def all_file_pins(directory):
    return {str(p.relative_to(ROOT)): sha(p) for p in sorted(directory.rglob('*'))
            if p.is_file() and '__pycache__' not in p.parts}


def shortest_words(edges):
    """Exhaust simple edge words up to length three, then minimize (length, word)."""
    answer = []
    multiplicities = Counter()
    for origin in range(12):
        candidates = {v: [] for v in range(12) if v != origin}

        def extend(vertices, word):
            if len(word) == 3:
                return
            for e, (a, b) in enumerate(edges):
                if vertices[-1] not in (a, b):
                    continue
                nxt = b if vertices[-1] == a else a
                if nxt in vertices:
                    continue
                extended = word + (e,)
                candidates[nxt].append(extended)
                extend(vertices + (nxt,), extended)

        extend((origin,), ())
        for destination, words in candidates.items():
            best = min(words, key=lambda w: (len(w), w))
            multiplicities[(len(best), sum(len(w) == len(best) for w in words))] += 1
            answer.append((origin, destination, best))
    demand(multiplicities == {(1, 1): 60, (2, 2): 60, (3, 10): 12},
           'independent shortest-path multiplicities')
    return answer, {f'distance_{d}_shortest_paths_{m}': n
                    for (d, m), n in sorted(multiplicities.items())}


def geometry_from_oriented_triangles(edges):
    """Orient primal triangles, then cross each outgoing edge to turn around a port."""
    lookup = {e: i for i, e in enumerate(edges)}
    faces = [t for t in itertools.combinations(range(12), 3)
             if all(e in lookup for e in itertools.combinations(t, 2))]
    adjacent = [[i for i, f in enumerate(faces) if a in f and b in f] for a, b in edges]
    demand(len(faces) == 20 and all(len(pair) == 2 for pair in adjacent), 'closed carrier')

    def coefficients(face):
        return {lookup[tuple(sorted((a, b)))]: (1 if a < b else -1)
                for a, b in zip(face, face[1:] + face[:1])}

    orientation = {0: faces[0]}
    pending = [0]
    while pending:
        f = pending.pop()
        for e, sign in coefficients(orientation[f]).items():
            g = next(i for i in adjacent[e] if i != f)
            if g not in orientation:
                raw = faces[g]
                orientation[g] = raw if coefficients(raw)[e] == -sign else raw[::-1]
                pending.append(g)
            demand(coefficients(orientation[g])[e] == -sign, 'triangle orientation conflict')
    contours = {}
    for v in range(12):
        start = min(i for i, f in enumerate(faces) if v in f)
        walk = [start]
        while True:
            f = orientation[walk[-1]]
            successor = f[(f.index(v) + 1) % 3]
            e = lookup[tuple(sorted((v, successor)))]
            nxt = next(i for i in adjacent[e] if i != walk[-1])
            if nxt == start:
                break
            demand(nxt not in walk, 'non-simple pentagon')
            walk.append(nxt)
        demand(len(walk) == 5, 'pentagon length')
        contours[v] = walk
    dual = [tuple(pair) for pair in adjacent]
    signs = [coefficients(orientation[pair[0]])[e] for e, pair in enumerate(adjacent)]
    for v, walk in contours.items():
        actual = [0] * len(edges)
        for a, b in zip(walk, walk[1:] + walk[:1]):
            e = dual.index(tuple(sorted((a, b))))
            actual[e] = signs[e] * (1 if a < b else -1)
        demand(actual == [int(v == a) - int(v == b) for a, b in edges], 'oriented contour incidence')
    return dual, contours, signs


def boundary(tape, edges):
    return tuple(sum((tape[e] for e, (a, _) in enumerate(edges) if v == a), Q())
                 - sum((tape[e] for e, (_, b) in enumerate(edges) if v == b), Q())
                 for v in range(12))


def singleton_load(state):
    # Different predicate from the producer's sum/integrality predicate.
    return tuple(sorted(state)) == (Q(0),) * 11 + (Q(1),)


@lru_cache(maxsize=None)
def region_properties(q):
    # Every subset, including the empty and whole carrier. No extremal shortcut.
    denominator = math.lcm(*(x.denominator for x in q))
    values = [0]
    for rational in q:
        charge = int(rational * denominator)
        values += [x + charge for x in values]
    demand(len(values) == 4096, 'regional enumeration')
    return (all(x in (-denominator, 0, denominator) for x in values),
            all(-denominator <= x <= denominator for x in values),
            Q(min(values), denominator), Q(max(values), denominator))


@lru_cache(maxsize=None)
def encoded_edge(a, b, turns, deformation):
    """Independent analytic evaluation on the predeclared uniform mesh."""
    alpha = 0.31 * math.sin(a + 1) * deformation
    beta = 0.31 * math.sin(b + 1) * deformation
    ripple = 0.23 * math.cos(a + b + 1) * deformation
    values = []
    for i in range(K + 1):
        t = i / K
        phase = (1 - t) * alpha + t * beta + math.tau * float(turns) * t + ripple * math.sin(math.tau * t)
        values.append(complex(math.cos(phase), math.sin(phase)))
    return tuple(values), abs(values[-1] - complex(math.cos(beta), math.sin(beta)))


def phase_degree(samples):
    demand(len(samples) > 1 and abs(samples[-1] - samples[0]) < EPS, 'open curve')
    demand(all(abs(z) > EPS for z in samples), 'zero crossing')
    increments = [math.atan2((b / a).imag, (b / a).real) for a, b in zip(samples, samples[1:])]
    demand(all(abs(a) < math.pi / 2 for a in increments), 'unresolved sample increments')
    return math.fsum(increments) / math.tau


@lru_cache(maxsize=None)
def own_phase(tape, geometry_key, deformation):
    dual, walks, signs = geometry_key
    curves = [encoded_edge(a, b, signs[e] * tape[e], deformation) for e, (a, b) in enumerate(dual)]
    defect = max(gap for _, gap in curves)
    if defect > EPS:
        return None, defect
    degrees = []
    for walk in walks:
        samples = []
        for a, b in zip(walk, walk[1:] + walk[:1]):
            values = curves[dual.index(tuple(sorted((a, b))))][0]
            if a > b:
                values = values[::-1]
            demand(not samples or abs(samples[-1] - values[0]) < EPS, 'vertex discontinuity')
            samples.extend(values if not samples else values[1:])
        degrees.append(phase_degree(samples))
    return tuple(degrees), defect


def independent_oracle(edges, schedules, geometry):
    geometry_key = (tuple(geometry[0]), tuple(tuple(geometry[1][v]) for v in range(12)), tuple(geometry[2]))
    rng = random.Random(104729)
    histories = []
    arms = {}
    native_counts = Counter()
    prefix_counts = Counter()
    for law, repeats in (('native', 1), ('half_current', 1), ('random_integer', 3), ('omitted_record', 1)):
        counts = dict.fromkeys(('histories', 'continuity', 'closed_curves', 'degree_equals_current',
                               'one_unit_every_step', 'regional_bound_every_step',
                               'missing_record_against_loads', 'histories_with_omission'), 0)
        max_defect = max_error = 0.
        for repeat in range(repeats):
            for origin, destination, path in schedules:
                probability = Q(0)
                for coins in itertools.product((0, 1), repeat=len(path)):
                    probability += Q(1, 2 ** len(path))
                    start = tuple(Q(v == origin) for v in range(12))
                    state, tape, held = list(start), [Q(0)] * 30, [Q(0)] * 30
                    omitted = False
                    steps = []
                    for e, c in zip(path, coins):
                        a, b = edges[e]
                        before = tuple(state)
                        if law in ('native', 'omitted_record'):
                            lower, odd = divmod(state[a] + state[b], 2)
                            next_a = lower + (odd if c else 0)
                            current = state[a] - next_a
                        elif law == 'half_current':
                            current = (state[a] - state[b]) * Q(1, 2)
                        else:
                            current = Q(rng.randrange(3) - 1)
                        state[a] -= current
                        state[b] += current
                        tape[e] += current
                        if law == 'omitted_record' and current != 0 and not omitted:
                            omitted = True
                        else:
                            held[e] += current
                        q = tuple(g - n for g, n in zip(start, state))
                        regions, real_bound, _, _ = region_properties(q)
                        gates = {'continuity': boundary(tape, edges) == q,
                                 'one_unit': singleton_load(state), 'regional_bound': regions}
                        prefix_counts[law + '_steps'] += 1
                        prefix_counts[law + '_continuity'] += gates['continuity']
                        prefix_counts[law + '_integral_regional_domain'] += regions
                        prefix_counts[law + '_real_regional_bound'] += real_bound
                        prefix_counts[law + '_source_record_agreement'] += boundary(held, edges) == q
                        steps.append(dict(seam=e, coin=c, current=current, before=before,
                                          after=tuple(state), tape=tuple(tape), held=tuple(held), gates=gates))
                        if law == 'native':
                            native_counts['events'] += 1
                            native_counts['holds'] += current == 0
                            native_counts['moves'] += current != 0
                    if law == 'native':
                        native_counts['arrivals'] += state[destination] == 1
                    period = boundary(held, edges)
                    q = tuple(g - n for g, n in zip(start, state))
                    normal, defect = own_phase(tuple(held), geometry_key, 0.)
                    max_defect = max(max_defect, defect)
                    counts['histories'] += 1
                    counts['continuity'] += all(s['gates']['continuity'] for s in steps)
                    counts['one_unit_every_step'] += all(s['gates']['one_unit'] for s in steps)
                    counts['regional_bound_every_step'] += all(s['gates']['regional_bound'] for s in steps)
                    counts['missing_record_against_loads'] += period != q
                    counts['histories_with_omission'] += omitted
                    if normal is not None:
                        deformed, _ = own_phase(tuple(held), geometry_key, 1.)
                        error = max(abs(float(p) - w) for values in (normal, deformed) for p, w in zip(period, values))
                        max_error = max(max_error, error)
                        counts['closed_curves'] += 1
                        counts['degree_equals_current'] += error < EPS
                    histories.append(dict(law=law, repeat=repeat, origin=origin, destination=destination,
                                          path=path, coins=coins, start=start, state=tuple(state),
                                          tape=tuple(tape), held=tuple(held), omitted=omitted,
                                          steps=steps, period=period, q=q, defect=defect))
                demand(probability == 1, 'complete fair-coin law per schedule')
        arms[law] = dict(counts, maximum_phase_error=max_error, maximum_gluing_defect=max_defect)
    demand(native_counts == {'events': 888, 'holds': 564, 'moves': 324, 'arrivals': 132}, 'native census')
    return arms, histories, dict(native_counts), dict(prefix_counts)


def load_target():
    spec = importlib.util.spec_from_file_location('law_controls_audited', TARGET)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def rejected(action):
    try:
        action()
    except (ValueError, KeyError, TypeError, AssertionError) as error:
        return str(error)
    raise ValueError('adversarial control survived')


CHECKS = {
    'native_histories_enumerated_from_law', 'native_whole_unit_at_every_step',
    'native_all_regional_bounds_at_every_step', 'native_fixed_sample_contours',
    'half_law_preserves_continuity_but_fails_closure',
    'integer_fragmenting_law_passes_topology_but_fails_unit_gate',
    'integer_fragmenting_law_can_exceed_regional_domain',
    'omission_detected_against_loads_for_every_affected_history',
    'damaged_integer_tape_still_has_valid_spatial_degree',
}


def validate_receipt(receipt, expected_arms, pins, histories):
    demand(receipt['schema'] == 'oph.row1.executed_law_controls.v1', 'receipt schema')
    demand(receipt['development_only'] is True and receipt['official_rows_changed'] == [], 'scope escalation')
    demand(receipt['samples_per_edge'] == K and receipt['declared_absolute_turn_bound'] == TURN_BOUND, 'sampling declaration')
    demand(set(receipt['checks']) == CHECKS, 'missing or additional check name')
    demand(all(v is True for v in receipt['checks'].values()), 'failed or nonboolean check')
    demand(receipt['files_sha256'] == pins, 'source pin mismatch or omission')
    demand(set(receipt['arms']) == set(expected_arms), 'missing or additional law arm')
    for arm, expected in expected_arms.items():
        actual = receipt['arms'][arm]
        demand(set(actual) == set(expected), 'missing arm statistic')
        for key, value in expected.items():
            if key.startswith('maximum_'):
                demand(isinstance(actual[key], (int, float)) and math.isfinite(actual[key])
                       and abs(actual[key] - value) < EPS, f'wrong numeric residual {arm}/{key}')
            else:
                demand(type(actual[key]) is int and actual[key] == value, f'wrong count {arm}/{key}')
    demand(set(receipt['witnesses']) == set(expected_arms), 'missing law witness')
    for arm, witness in receipt['witnesses'].items():
        selected = next(h for h in histories if h['law'] == arm and (arm == 'native' or h['defect'] > EPS
                        or h['period'] != h['q'] or not all(s['gates']['one_unit'] for s in h['steps'])))
        demand(witness['origin'] == selected['origin'] and witness['destination'] == selected['destination']
               and witness['schedule'] == list(selected['path']), 'wrong witness schedule')
        events = [{'seam_index': s['seam'], 'coin': s['coin'], 'current': str(s['current']),
                   'before': list(map(str, s['before'])), 'after': list(map(str, s['after']))}
                  for s in selected['steps']]
        demand(witness['events'] == events, 'witness current/state/event omission')
        demand(witness['prefix_gates'] == [s['gates'] for s in selected['steps']], 'witness prefix gate')
        demand(witness['held_period'] == list(map(str, selected['period']))
               and witness['load_difference'] == list(map(str, selected['q'])), 'witness continuity')
        demand(abs(witness['terminal_gluing_defect'] - selected['defect']) < EPS, 'witness gluing defect')


def compare_implementation(module, edges, schedules, geometry, arms, histories):
    demand(list(module.schedules(edges)) == schedules, 'actual scheduler differs from independent edge-word search')
    actual_geometry = module.spatial.dual_boundaries(edges)
    demand(actual_geometry[0] == geometry[0] and actual_geometry[2] == geometry[2], 'dual edges/signs differ')
    # Walk starting points are immaterial; both use the same orientation.
    for v, walk in actual_geometry[1].items():
        own = geometry[1][v]
        offset = own.index(walk[0])
        demand(walk == own[offset:] + own[:offset], 'contour rotation differs')
    rng = random.Random(104729)
    tested_steps = 0
    phase_calls = 0
    for h in histories:
        actual = module.run_history(h['origin'], h['path'], h['coins'], edges, h['law'], rng)
        start, state, tape, held, gates, events, omitted = actual
        demand(tuple(start) == h['start'] and tuple(state) == h['state'] and tuple(tape) == h['tape']
               and tuple(held) == h['held'] and omitted == h['omitted'], 'actual history mismatch')
        demand(gates == [s['gates'] for s in h['steps']], 'actual per-step gate mismatch')
        for event, step in zip(events, h['steps']):
            demand(event == {'seam_index': step['seam'], 'coin': step['coin'], 'current': str(step['current']),
                             'before': list(map(str, step['before'])), 'after': list(map(str, step['after']))},
                   'actual event mismatch')
            demand(module.one_unit(step['after']) == singleton_load(step['after']), 'unit gate disagreement')
            q = tuple(g - n for g, n in zip(h['start'], step['after']))
            demand(module.regional_bound(q) == region_properties(q)[0], 'regional shortcut disagreement')
            tested_steps += 1
        for deformed in (False, True):
            values, defect = module.phase_readout(h['held'], actual_geometry, deformed)
            phase_calls += 1
            if h['defect'] > EPS:
                demand(values is None and defect > EPS, 'fractional history accepted')
            else:
                demand(values is not None and max(abs(float(p) - x) for p, x in zip(h['period'], values)) < EPS,
                       'actual phase curves disagree with independently reconstructed period')
    demand(tested_steps == 5328 and phase_calls == 5472, 'comparison coverage')
    return {'histories': len(histories), 'prefixes': tested_steps, 'phase_readouts': phase_calls,
            'terminal_pentagon_degrees_on_closed_arms': 2 * 12 * (456 + 1368 + 456)}


def curve_controls(module, geometry, edges):
    analytic_samples = 0
    max_formula_error = 0.
    for a, b in geometry[0]:
        for turns in (Q(-3), Q(-2), Q(-1), Q(-1, 2), Q(0), Q(1, 2), Q(1), Q(2), Q(3)):
            for deformation in (False, True):
                actual, defect = module.curve_fixed(a, b, turns, deformation)
                expected, expected_defect = encoded_edge(a, b, turns, float(deformation))
                demand(len(actual) == K + 1, 'adaptive or incomplete fixed sampling')
                demand(max(abs(abs(z) - 1) for z in actual) < EPS, 'off-circle curve')
                error = max(abs(x - y) for x, y in zip(actual, expected))
                max_formula_error = max(max_formula_error, error)
                demand(error < EPS and abs(defect - expected_defect) < EPS, 'curve formula mismatch')
                analytic_samples += len(actual)
    controls = {}
    controls['fixed_turn_32_alias_rejected'] = rejected(lambda: module.curve_fixed(0, 1, 32, False))
    controls['fixed_turn_4_outside_domain_rejected'] = rejected(lambda: module.curve_fixed(0, 1, 4, False))
    for name, values in (('zero', [1, 0, 1]), ('open', [1, 1j, -1]), ('large_increment', [1, -1, 1])):
        controls[name] = rejected(lambda values=values: module.spatial.winding(values))
    # An aliased smooth 32-turn loop passes the sampled-increment gate without a bandwidth contract.
    alias = [cmath.exp(1j * math.tau * 32 * k / K) for k in range(K + 1)]
    alias_reading = module.spatial.winding(alias)
    demand(abs(alias_reading) < EPS, 'alias demonstration')
    ordinary_counts = [len(module.curve_fixed(0, 1, j, False)[0]) - 1 for j in (0, 1, 2, 3)]
    adaptive_counts = [len(module.spatial.curve(0, 1, j, False)) - 1 for j in (0, 1, 2, 3)]
    demand(ordinary_counts == [32] * 4 and len(set(adaptive_counts)) > 1, 'fixed versus adaptive density')
    module.curve_fixed.cache_clear()
    original_density = module.SAMPLES_PER_EDGE
    module.SAMPLES_PER_EDGE = 4
    try:
        controls['insufficient_density'] = rejected(lambda: module.curve_fixed(0, 1, 1, False))
    finally:
        module.SAMPLES_PER_EDGE = original_density
        module.curve_fixed.cache_clear()

    # Free homotopy between F and F plus a dual coboundary.
    dual, contours, signs = geometry
    potential = [int(v == 0) for v in range(20)]
    change = tuple(Q(signs[e] * (potential[b] - potential[a])) for e, (a, b) in enumerate(dual))
    demand(boundary(change, edges) == (Q(0),) * 12 and any(change), 'nontrivial zero-charge coboundary')
    max_homotopy_degree = max_gluing = 0.
    moving_vertex = False
    for s in (0., .25, .5, .75, 1.):
        curves = []
        for a, b in dual:
            values = [cmath.exp(1j * math.tau * s * ((1 - k / K) * potential[a] + k / K * potential[b]))
                      for k in range(K + 1)]
            curves.append(values)
        moving_vertex |= abs(curves[0][0] - 1) > EPS
        for walk in contours.values():
            samples = []
            for a, b in zip(walk, walk[1:] + walk[:1]):
                values = curves[dual.index(tuple(sorted((a, b))))]
                values = values if a < b else values[::-1]
                if samples:
                    max_gluing = max(max_gluing, abs(samples[-1] - values[0]))
                samples += values if not samples else values[1:]
            max_homotopy_degree = max(max_homotopy_degree, abs(phase_degree(samples)))
    actual_end, _ = module.phase_readout(change, geometry)
    demand(moving_vertex and max_homotopy_degree < EPS and max_gluing < EPS
           and max(map(abs, actual_end)) < EPS, 'free homotopy')
    # The native rule halves an even stack; the protected unit is exactly one register load.
    splits = []
    for coin in (0, 1):
        output, record = module.native((2, 0) + (0,) * 10, (0, 1), coin)
        demand(output == (1, 1) + (0,) * 10 and record == (coin, 1), 'native even-stack fission')
        splits.append({'coin': coin, 'after': list(output), 'current': record[1]})
    return {'sample_points_compared_to_analytic_curve': analytic_samples,
            'maximum_formula_error': max_formula_error, 'rejecting_controls': controls,
            'fixed_subdivisions_turns_0_1_2_3': ordinary_counts,
            'frozen_adaptive_subdivisions_turns_0_1_2_3': adaptive_counts,
            'unchecked_32_turn_sampled_readout': alias_reading,
            'global_derivative_bound_per_subdivision': (.62 + math.tau * (3 + .23)) / K,
            'free_homotopy': {'stages': 5, 'vertex_phases_move': moving_vertex,
                              'maximum_degree': max_homotopy_degree, 'maximum_junction_gap': max_gluing},
            'even_integer_stack_split': splits}


def main():
    before = all_file_pins(NATIVE)
    demand(sha(NATIVE / 'manifest.json') == MANIFEST_PIN, 'frozen manifest pin')
    manifest = json.loads((NATIVE / 'manifest.json').read_text())
    for name, expected in manifest['files_sha256'].items():
        demand(sha(NATIVE / name) == expected, 'frozen member pin: ' + name)
    source_pin = sha(TARGET)
    packet = json.loads((NATIVE / 'source_packet.json').read_text())
    edges = [tuple(e) for e in packet['seams']]
    schedules, multiplicities = shortest_words(edges)
    geometry = geometry_from_oriented_triangles(edges)
    arms, histories, census, prefixes = independent_oracle(edges, schedules, geometry)
    # The entire independent result exists before importing any producer function.
    oracle = {'schema': 'oph.row1.independent_law_oracle.v1', 'arms': arms, 'native_events': census,
              'prefix_counts': prefixes, 'shortest_schedule_multiplicities': multiplicities}
    (HERE / 'ORACLE.json').write_text(json.dumps(oracle, indent=2) + '\n')
    print('Independent oracle complete:', json.dumps({'native_events': census, 'arms': arms}), flush=True)
    module = load_target()
    coverage = compare_implementation(module, edges, schedules, geometry, arms, histories)
    curve = curve_controls(module, geometry, edges)
    pins = {str(p.relative_to(RER / 'code')): sha(p)
            for p in (TARGET, NATIVE / 'native_kernel.py', NATIVE / 'native_charge_winding.py',
                      NATIVE / 'source_packet.json', NATIVE / 'manifest.json')}
    demand(RECEIPT.exists(), 'producer receipt absent; rerun after the author execution completes')
    receipt = json.loads(RECEIPT.read_text())
    validate_receipt(receipt, arms, pins, histories)
    replay = module.execute()
    validate_receipt(replay, arms, pins, histories)
    demand(replay == receipt, 'actual execute() receipt differs from saved receipt')
    mutations = {}

    def mutate(name, operation):
        candidate = copy.deepcopy(receipt)
        operation(candidate)
        mutations[name] = rejected(lambda: validate_receipt(candidate, arms, pins, histories))

    mutate('empty_checks', lambda d: d.update(checks={}))
    mutate('missing_check', lambda d: d['checks'].pop('native_whole_unit_at_every_step'))
    mutate('missing_source_pin', lambda d: d['files_sha256'].pop('openwave_row1_native/source_packet.json'))
    mutate('all_true_forged_counts', lambda d: d['arms']['random_integer'].update(one_unit_every_step=1368))
    mutate('missing_mutation_arm', lambda d: d['arms'].pop('half_current'))
    mutate('nan_residual', lambda d: d['arms']['native'].update(maximum_phase_error=float('nan')))
    mutate('missing_witness_event', lambda d: d['witnesses']['half_current']['events'].pop())
    mutate('scope_escalation', lambda d: d.update(official_rows_changed=[1]))
    mutate('adaptive_density_claim', lambda d: d.update(samples_per_edge=8))
    mutate('missing_prefix_check', lambda d: d['witnesses']['native']['prefix_gates'][0].pop('continuity'))
    # Mutation of executable behavior, not a forged receipt alone.
    authentic = module.phase_readout
    module.phase_readout = lambda tape, geometry, deformed=False: ([0.] * 12, 0.)
    try:
        mutations['constant_zero_phase_implementation'] = rejected(module.execute)
    finally:
        module.phase_readout = authentic
    authentic_regional = module.regional_bound
    module.regional_bound = lambda q: True
    try:
        mutations['literal_true_regional_implementation'] = rejected(module.execute)
    finally:
        module.regional_bound = authentic_regional
    tree = ast.parse(TARGET.read_text())
    checks_nodes = [n.value for n in ast.walk(tree) if isinstance(n, ast.Assign)
                    and any(isinstance(t, ast.Name) and t.id == 'checks' for t in n.targets)]
    demand(len(checks_nodes) == 1 and isinstance(checks_nodes[0], ast.Dict)
           and not any(isinstance(v, ast.Constant) for v in checks_nodes[0].values), 'literal receipt checks')
    demand(sha(TARGET) == source_pin, 'source changed during audit')
    demand(all_file_pins(NATIVE) == before, 'frozen package changed during audit')
    result = {
        'schema': 'oph.row1.independent_law_gate_audit.v1',
        'verdict': 'PASS_WITH_DECLARED_SCOPE', 'development_only': True, 'official_rows_changed': [],
        'source_commit': SOURCE_COMMIT, 'source_manifest_sha256': MANIFEST_PIN,
        'source_sha256': pins, 'producer_receipt_sha256': sha(RECEIPT), 'auditor_sha256': sha(Path(__file__)),
        'independence': 'Separate stdlib implementation; source and BR6/BR9 findings read before coding; oracle completed before producer import.',
        'native_census': census, 'schedule_multiplicities': multiplicities,
        'arms': arms, 'prefix_counts': prefixes, 'actual_implementation_comparison': coverage,
        'exhaustive_regions_per_distinct_prefix_charge': 4096,
        'distinct_prefix_charges': region_properties.cache_info().currsize,
        'curve_controls': curve, 'rejected_receipt_and_implementation_mutations': mutations,
        'checks': {'all_132_shortest_schedules_match': len(schedules) == 132,
                   'all_456_native_coin_histories_match': arms['native']['histories'] == 456,
                   'all_native_prefixes_have_exact_continuity': prefixes['native_continuity'] == 888,
                   'all_native_prefixes_have_integral_regional_domain': prefixes['native_integral_regional_domain'] == 888,
                   'actual_curves_tested': coverage['phase_readouts'] == 5472,
                   'producer_receipt_matches_independent_oracle_and_replay': replay == receipt,
                   'malicious_receipts_and_implementations_rejected': len(mutations) == 12,
                   'fixed_package_unchanged': all_file_pins(NATIVE) == before},
        'interpretation': {
            'carrier': 'Exactly one declared twelve-port icosahedral carrier; no glued-carrier transfer.',
            'reference': 'BF=G-n requires retained preparation G and zero initial tape.',
            'topology': 'Degree identifies the integer tape class; random integer fragmented laws also pass topology.',
            'unit': 'One unit is one integer register load; even stacks split under the same native law.',
            'homotopy': 'Coboundary-related tapes admit free homotopy with moving shared vertex phases; fixed-vertex homotopy is not asserted.',
            'observation': 'The curve is a declared encoding of current records, not an independent dynamical field sensor; the pointer reads counters.',
            'omission': 'A missing moving record is detected against source loads; a valid damaged integer tape has no internal omission detector.',
            'sampler': '32 subdivisions and absolute turn bound 3 are predeclared development-domain assumptions; arbitrary long histories require another contract.',
            'status': 'No official row validation, prospective acceptance execution, physical electric-unit identification, mass, spin or force result.'},
        'prior_findings_sha256': {str(p.relative_to(ROOT)): sha(p) for p in (
            ROOT / 'openwave-push/bridgework/BR6_dual_winding/AUDIT.md',
            ROOT / 'openwave-push/bridgework/BR6_dual_winding/REVIEW_FOR_CODEX.md',
            ROOT / 'openwave-push/bridgework/BR9_precedent_replication/AUDIT.md',
            ROOT / 'openwave-push/bridgework/BR9_precedent_replication/REVIEW_FOR_CODEX.md')},
    }
    demand(all(result['checks'].values()), 'audit check failed')
    (HERE / 'AUDIT_RESULT.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'verdict': result['verdict'], 'checks': result['checks'], 'mutations': mutations}, indent=2))


if __name__ == '__main__':
    main()
