#!/usr/bin/env python3
"""Independent exact audit. Writes only this directory; standard library only."""
from __future__ import annotations
import ast
from collections import Counter, defaultdict, deque
from copy import deepcopy
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import itertools
import json
from pathlib import Path
import random
import sys

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
PRODUCER = ROOT / 'reverse-engineering-reality/code/openwave_row1_integrity'
RECEIPTS = HERE.parent / 'row1_integrity'
EXPECTED_MANIFEST = '1f8c91d4820fa2fe5f439f0a8e344645b4c29b11e3541ac32b7c52a71f1d73c4'
sys.path.insert(0, str(PRODUCER))
import integrity as p
from observer_only import admit_record_only, admit_record_only_naive

# Declared geometry, read from its pinned packet. No producer transition,
# routing, sweep distribution, reset law, or validator computes expected laws.
EDGES = tuple(tuple(e) for e in json.loads((PRODUCER/'source_packet.json').read_text())['seams'])
CYCLE = (0, 1, 8, 6, 2, 10)
SUCCESS = 'complete_closed_pointer_word'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def save(name, value):
    (HERE/name).write_text(json.dumps(value, indent=2, sort_keys=True)+'\n')


def pins():
    assert sha(PRODUCER/'MANIFEST.json') == EXPECTED_MANIFEST
    manifest = json.loads((PRODUCER/'MANIFEST.json').read_text())
    inputs = json.loads((PRODUCER/'INPUTS.json').read_text())
    files = {str(PRODUCER.relative_to(ROOT)/'MANIFEST.json'): EXPECTED_MANIFEST}
    for rel, wanted in manifest['files_sha256'].items():
        file = PRODUCER/rel
        assert sha(file) == wanted, rel
        files[str(file.relative_to(ROOT))] = wanted
    for group, root in [('runtime_files_sha256', PRODUCER), ('external_files_sha256', ROOT)]:
        for rel, wanted in inputs[group].items():
            file = root/rel
            assert sha(file) == wanted, rel
            files[str(file.relative_to(ROOT))] = wanted
    result = json.loads((RECEIPTS/'FROZEN_RESULT.json').read_text())
    assert sha(PRODUCER/'RESULTS.json') == manifest['result_file_sha256']
    assert sha(RECEIPTS/'FROZEN_RESULT.json') == manifest['result_file_sha256']
    assert sha(RECEIPTS/'PORTABLE_MANIFEST.json') == EXPECTED_MANIFEST
    scientific = deepcopy(result)
    del scientific['elapsed_seconds']
    assert digest(scientific) == manifest['result_scientific_sha256']
    for file in [PRODUCER/'RESULTS.json', *(RECEIPTS/f for f in
                 ['FROZEN_RESULT.json', 'PORTABLE_MANIFEST.json', 'SUMMARY.json'])]:
        files[str(file.relative_to(ROOT))] = sha(file)
    assert result['frozen_inputs_before'] == result['frozen_inputs_after']
    assert len(result['frozen_inputs_before']['files']) == 54
    return {'producer_manifest_sha256': EXPECTED_MANIFEST, 'files_sha256': files,
            'frozen_runtime_and_external_entries': sum(len(inputs[g]) for g in ('runtime_files_sha256', 'external_files_sha256')),
            'scientific_result_sha256': digest(scientific)}, result


def native(loads, edge, coin):
    """Independent balancing form: split the pair total, not producer j(d,c)."""
    a, b = edge
    pair_total = loads[a] + loads[b]
    after = list(loads)
    after[a] = pair_total // 2 + coin * (pair_total % 2)
    after[b] = pair_total - after[a]
    return tuple(after), loads[a] - after[a]


@lru_cache(None)
def route(port):
    """Queue complete paths, independently of predecessor-map implementation."""
    todo = deque([(port,)])
    seen = {port}
    while todo:
        path = todo.popleft()
        if path[-1] == 0:
            return path[1:]
        for edge in EDGES:
            if path[-1] not in edge:
                continue
            nxt = edge[0] + edge[1] - path[-1]
            if nxt not in seen:
                seen.add(nxt)
                todo.append(path + (nxt,))
    raise AssertionError('Disconnected geometry')


def preparations():
    return [tuple(int(i == a) for i in range(12)) for a in range(12)] + [
        tuple(int(i == a)+int(i == b) for i in range(12))
        for a in range(12) for b in range(a, 12)]


def geometry_and_kernel():
    assert EDGES == p.EDGES == p.strict.SEAMS == p.sweep.SEAMS
    assert len(EDGES) == 30 and len(set(EDGES)) == 30
    assert Counter(v for edge in EDGES for v in edge) == Counter({v: 5 for v in range(12)})
    assert max(len(route(v)) for v in range(12)) == 3
    assert all(route(v) == p.meter.route(v, 0) == p.strict.route_to_zero(v) for v in range(12))
    branches = 0
    for G in preparations():
        for edge in EDGES:
            for coin in (0, 1):
                after, j = native(G, edge, coin)
                actual, (_, actual_j) = p.native(G, edge, coin)
                strict_after, strict_j = p.strict.native(G, edge, coin)
                sweep_after, d, sweep_j = p.sweep.native_event(G, edge, coin)
                assert after == actual == strict_after == sweep_after
                assert j == actual_j == strict_j == sweep_j
                assert d == G[edge[0]]-G[edge[1]] and sum(after) == sum(G) and min(after) >= 0
                branches += 1
    def defs(file):
        return {n.name: ast.dump(n, include_attributes=False)
                for n in ast.parse(file.read_text()).body if isinstance(n, ast.FunctionDef)}
    a = defs(PRODUCER/'native_cycle_meter.py')
    b = defs(PRODUCER/'provenance/PRIMARY_native_cycle_meter.original.py')
    names = ['execute', 'route', 'reset_success', 'favorable_coins', 'independent_reset_distribution']
    assert all(a[k] == b[k] for k in names)
    return {'independent_native_branches': branches, 'preparations': 90,
            'connected': True, 'maximum_route_length': 3, 'original_pointer_AST_functions': names}


def measure_state(w):
    sign = 1 if w > 0 else -1
    targets = tuple(CYCLE[(sign*k) % 6] for k in range(1, 6*abs(w)+1))
    if not targets:
        return (SUCCESS, w, 0)
    return (0, 'measure', targets, 0, 0)


def start_state(initial):
    return (initial, 'discover', tuple(range(30)), 0, 0)


def step(state, w, retries, coin):
    """One controller event; terminal tuples have length 3, active tuples 5."""
    port, stage, targets, spent, signed = state
    edge = EDGES[targets[0]] if stage == 'discover' else tuple(sorted((port, targets[0])))
    before = tuple(int(v == port) for v in range(12))
    after, j = native(before, edge, coin)
    record = {'stage': stage, 'seam': list(edge), 'coin': coin, 'current': j,
              'requested_cycle_direction': (1 if w > 0 else -1) if stage == 'measure' else None}
    if j:
        port = after.index(1)
        if stage == 'discover':
            path = route(port)
            nxt = (port, 'reset_route', path, 0, 0) if path else measure_state(w)
        elif stage == 'reset_route':
            nxt = (port, stage, targets[1:], 0, 0) if len(targets) > 1 else measure_state(w)
        else:
            signed += 1 if w > 0 else -1
            nxt = (port, stage, targets[1:], 0, signed) if len(targets) > 1 else (SUCCESS, w, signed)
    elif spent + 1 < retries:
        nxt = (port, stage, targets, spent + 1, signed)
    elif stage == 'discover' and len(targets) > 1:
        nxt = (port, stage, targets[1:], 0, signed)
    else:
        reason = {'discover': 'discovery_timeout', 'reset_route': 'routing_timeout',
                  'measure': 'measurement_timeout'}[stage]
        nxt = (reason, None, signed)
    return nxt, record, after


@lru_cache(None)
def pointer_law(initial, w, retries):
    """Event-by-event absorbing Markov recurrence, using no block formulas."""
    live = {start_state(initial): F(1)}
    law = defaultdict(F)
    attempts = 0
    while live:
        attempts += 1
        nxt = defaultdict(F)
        for state, mass in live.items():
            for coin in (0, 1):
                successor, _, _ = step(state, w, retries, coin)
                if len(successor) == 3:
                    law[(*successor, attempts)] += mass / 2
                else:
                    nxt[successor] += mass / 2
        live = dict(nxt)
    assert sum(law.values()) == 1
    return dict(law)


def law_serial(law):
    return [list(k)+[str(v)] for k, v in sorted(law.items(), key=lambda item: repr(item[0]))]


def pointer_laws():
    rows = []
    for retries in (1, 2, 3):
        for initial in range(12):
            for w in range(-2, 3):
                law = pointer_law(initial, w, retries)
                assert law == p.pointer_distribution(initial, w, retries)
                h = F(1, 2**retries)
                incident = [e[0]+e[1]-initial for e in EDGES if initial in e]
                reset = sum(h**k * (1-h)**(1+len(route(y))) for k, y in enumerate(incident))
                expected_success = reset * (1-h)**(6*abs(w))
                by_reason = defaultdict(F)
                for key, mass in law.items():
                    by_reason[key[0]] += mass
                    assert key[-1] <= (33+6*abs(w))*retries
                    if key[0] != SUCCESS:
                        assert key[1] is None
                assert by_reason[SUCCESS] == expected_success
                assert by_reason['discovery_timeout'] == h**5
                rows.append({'initial': initial, 'command': w, 'retries': retries,
                    'mass': '1', 'outcomes': {k: str(v) for k, v in by_reason.items()},
                    'maximum_attempts': max(k[-1] for k in law),
                    'expected_attempts': str(sum(k[-1]*v for k, v in law.items())),
                    'law_sha256': digest(law_serial(law))})
    # Arbitrary unknown classical priors are mixtures, without conditioning
    # on success. This deliberately nonuniform rational prior exercises it.
    prior = [F(x+1, 78) for x in range(12)]
    mixed = defaultdict(F)
    for x, weight in enumerate(prior):
        for key, mass in pointer_law(x, -2, 2).items():
            mixed[key] += weight*mass
    assert sum(mixed.values()) == 1
    effects = []
    for retries in (1, 2, 3):
        for w in range(-2, 3):
            instrument = {}
            for x in range(12):
                for outcome, likelihood in pointer_law(x, w, retries).items():
                    instrument.setdefault(outcome, [F()] * 12)[x] += likelihood
            diagonal_sum = [sum(effect[x] for effect in instrument.values()) for x in range(12)]
            assert diagonal_sum == [F(1)] * 12
            effects.append({'command': w, 'retries': retries, 'outcome_effect_count': len(instrument),
                'effect_sum_diagonal': [str(v) for v in diagonal_sum],
                'effects_sha256': digest([list(k)+[[str(v) for v in effect]] for k, effect in
                                          sorted(instrument.items(), key=lambda item: repr(item[0]))])})
    return {'complete_exact_laws': len(rows), 'recurrence': 'per native event; all coin mass retained',
            'rows': rows, 'nonuniform_unknown_initial_mixture_mass': '1',
            'mixture_prior': [str(x) for x in prior], 'mixture_sha256': digest(law_serial(mixed)),
            'classical_effect_instruments': effects, 'every_effect_sum_equals_identity': True}


def reference_execution(initial, w, retries, failure=None, partial=0):
    """Generate actual-input tapes from the independent controller, no frames."""
    state = start_state(initial)
    records = []
    after = tuple(int(v == initial) for v in range(12))
    while len(state) != 3:
        choices = [step(state, w, retries, c) for c in (0, 1)]
        port, stage, targets, spent, signed = state
        hold = next((c for c, item in enumerate(choices) if not item[1]['current']), 0)
        move = next((c for c, item in enumerate(choices) if item[1]['current']), None)
        fail_here = stage == failure and (stage != 'measure' or abs(signed) >= partial)
        # Discovery directly to zero has no routing stage. Hold such a
        # discovery edge to reach a later incident seam when testing routing.
        avoid_zero = stage == 'discover' and failure == 'reset_route' and move is not None and choices[move][2][0] == 1
        coin = hold if fail_here or avoid_zero or move is None else move
        state, record, after = choices[coin]
        records.append(record)
    reason, winding, signed = state
    return {'success': reason == SUCCESS, 'reason': reason, 'records': records,
            'final_loads': list(after), 'retained_signed_cycle_steps': signed,
            'winding': winding, 'fine_history_probability': str(F(1, 2**len(records)))}


def tape_only_columns(meter, w, retries):
    """Independent controller replay, excluding all terminal diagnostic loads."""
    compatible = []
    for initial in range(12):
        state = start_state(initial)
        valid = True
        for record in meter['records']:
            if len(state) == 3:
                valid = False
                break
            if type(record.get('coin')) is not int or record['coin'] not in (0, 1):
                valid = False
                break
            state, expected, _ = step(state, w, retries, record['coin'])
            if record != expected:
                valid = False
                break
        if valid and len(state) == 3 and state == (meter['reason'], meter['winding'], meter['retained_signed_cycle_steps']):
            compatible.append(initial)
    return compatible


def actual_executions():
    rows = []
    examples = []
    outcomes = Counter()
    for initial in range(12):
        for w in range(-2, 3):
            for failure, partial in [(None, 0), ('discover', 0), ('reset_route', 0)] + (
                    [('measure', 0), ('measure', 3)] if w else []):
                expected = reference_execution(initial, w, 2, failure, partial)
                actual = p.meter.execute(initial, (w,), (r['coin'] for r in expected['records']), 2)
                assert actual == expected
                errors, _ = p.strict.validate_pointer_tape(actual, initial, w, 2)
                assert not errors, errors
                columns = tape_only_columns(actual, w, 2)
                assert initial in columns
                outcomes[actual['reason']] += 1
                rows.append({'initial': initial, 'command': w, 'failure_target': failure,
                    'partial_target': partial, 'outcome': actual['reason'], 'attempts': len(actual['records']),
                    'tape_only_columns': columns, 'meter_sha256': digest(actual)})
                examples.append((initial, w, actual, columns))
    return {'executions': len(rows), 'outcomes': dict(outcomes), 'rows': rows}, examples


def source_joint_counts(G, order):
    """Integer word counts, preserving actual L0 state through all 30 steps."""
    prefix = Counter({G: 1})
    for e in (0, 5):
        nxt = Counter()
        for n, count in prefix.items():
            for coin in (0, 1):
                nxt[native(n, EDGES[e], coin)[0]] += count
        prefix = nxt
    counts = {(n, n): count for n, count in prefix.items()}
    for e in order:
        nxt = Counter()
        for (n0, n), count in counts.items():
            for coin in (0, 1):
                nxt[(n0, native(n, EDGES[e], coin)[0])] += count
        counts = nxt
    assert sum(counts.values()) == 2**32
    return counts


def joint_laws(frozen):
    starts = preparations()[:12] + [tuple(int(v == 0)+int(v == b) for v in range(12)) for b in range(4)]
    rows = []
    expected = frozen['complete_bounded_laws']['rows']
    for G in starts:
        for name, order in [('identity', range(30)), ('reverse', range(29, -1, -1))]:
            counts = source_joint_counts(G, order)
            for region in [(0, 1, 2), tuple(range(3, 12))]:
                source = defaultdict(F)
                for (n0, n1), count in counts.items():
                    w0 = sum(G[v]-n0[v] for v in region)
                    w1 = sum(G[v]-n1[v] for v in region)
                    source[(w0, w1)] += F(count, 2**32)
                for initial in range(12):
                    joint = defaultdict(F)
                    for (w0, w1), mass in source.items():
                        for (reason, winding, signed, attempts), pmass in pointer_law(initial, w1, 1).items():
                            joint[(w0, w1, reason, -winding if winding is not None else None,
                                   signed, 32+attempts)] += mass*pmass
                    assert sum(joint.values()) == 1
                    reasons = defaultdict(F)
                    for key, mass in joint.items():
                        reasons[key[2]] += mass
                        assert key[-1] <= 32+33+6*sum(G)
                    row = {'G': list(G), 'order': name, 'region': list(region), 'pointer_initial_column': initial,
                        'mass': '1', 'support': len(joint), 'outcome_mass': {k: str(v) for k, v in reasons.items()},
                        'expected_total_attempts': str(sum(k[-1]*mass for k, mass in joint.items())),
                        'maximum_attempts': max(k[-1] for k in joint),
                        'complete_joint_law_sha256': digest(law_serial(joint))}
                    assert row == expected[len(rows)], (len(rows), row, expected[len(rows)])
                    rows.append(row)
    assert len(rows) == 768
    return {'independently_regenerated_joint_maps': len(rows), 'every_frozen_row_exactly_equal': True,
            'mass_error': '0', 'source_coin_words_per_setting': 2**32,
            'independent_source_recurrences': 32, 'row_list_sha256': digest(rows),
            'all_failures_retained': True}


def independent_source(G, prefix, order, coins):
    n = G
    tapes = [[], []]
    latches = []
    for phase, events in enumerate((prefix, [(e, c, 'E') for e, c in zip(order, coins)])):
        for e, c, kind in events:
            a, b = EDGES[e]
            d = n[a]-n[b]
            n, j = native(n, (a, b), c)
            rec = {'kind': kind, 'seam': e, 'coin': c, 'current': j}
            if kind == 'E':
                rec['contrast'] = d
            tapes[phase].append(rec)
        latches.append(n)
    return {'source_records': tapes[0], 'sweep_records': tapes[1], 'sweep_order': list(order)}, latches


def full_sources():
    rng = random.Random(601114)
    stale = None
    rows = []
    for k, G in enumerate(preparations()):
        prefix = [(rng.randrange(30), rng.randrange(2), rng.choice(('C', 'E')))
                  for _ in range((k*17) % 129)]
        order = list(range(30)); rng.shuffle(order)
        coins = [rng.randrange(2) for _ in range(30)]
        expected, states = independent_source(G, prefix, order, coins)
        received, audit = p.source_execute(G, prefix, coins, order)
        assert received == expected
        assert audit['source_L0'] == list(states[0]) and audit['source_L1'] == list(states[1])
        latches = p.observe_source(G, received, sum(G))
        for i, key in enumerate(('L0', 'L1')):
            assert latches[key]['reconstruction'] == list(states[i])
            assert latches[key]['counters'] == [(g-n) % 6 for g, n in zip(G, states[i])]
        region = tuple(v for v in range(12) if (v+k) % 3 == 0)
        w = sum(G[v]-states[1][v] for v in region)
        expected_meter = reference_execution(k % 12, w, 2)
        m = p.meter.execute(k % 12, (w,), (r['coin'] for r in expected_meter['records']), 2)
        assert m == expected_meter
        claim = p.assemble_claim(G, received, region, m, 2)
        verdict = admit_record_only(claim, G, region, 2, producer=p)
        assert verdict['admitted'] and verdict['readout'] == -w
        assert verdict['compatible_pointer_initials'] == tape_only_columns(m, w, 2)
        assert claim['combined_attempts'] <= len(prefix)+30+(33+6*sum(G))*2
        polluted = deepcopy(claim)
        polluted.update({'source_L0': [999]*12, 'source_L1': None, 'source_after': 'poison',
                         'audit_only': {'final_n': [999]*12, 'pointer_initial': -999}})
        assert admit_record_only(polluted, G, region, 2, producer=p) == verdict
        w0 = sum(G[v]-states[0][v] for v in region)
        if stale is None and w0 != w:
            bad = deepcopy(claim); bad['charge_read_at'] = 'L0'
            assert not admit_record_only(bad, G, region, 2, producer=p)['admitted']
            stale = {'G': list(G), 'region': list(region), 'L0': list(states[0]), 'L1': list(states[1]),
                     'w0': w0, 'w1': w, 'readout': verdict['readout'], 'stale_claim_rejected': True}
        rows.append({'G': list(G), 'prefix_length': len(prefix), 'source_tape_sha256': digest(received),
                     'attempts': claim['combined_attempts'], 'charge': verdict['readout']})
    assert stale is not None
    return {'native_source_to_L1_pointer_examples': len(rows), 'all_90_preparations': True,
            'poisoned_source_metadata_invariance_cases': len(rows), 'stale_L0_witness': stale, 'rows': rows}


def make_claim(G, received, region, initial=0, retries=2):
    latches = p.observe_source(G, received, sum(G))
    w = sum(G[v]-latches['L1']['reconstruction'][v] for v in region)
    ptr = reference_execution(initial, w, retries)
    return p.assemble_claim(G, received, region, ptr, retries)


def facade_controls(examples):
    # Fixed closing sweep gives one charge for each region; use independent
    # source words until every bounded command needed by the pointer fixtures
    # has a legitimate full-instrument source envelope.
    sources = {}
    rng = random.Random(114601)
    for _ in range(500):
        G = preparations()[rng.randrange(90)]
        order = list(range(30)); rng.shuffle(order)
        received, states = independent_source(G, [], order, [rng.randrange(2) for _ in range(30)])
        for bits in range(4096):
            region = tuple(v for v in range(12) if bits >> v & 1)
            w = sum(G[v]-states[1][v] for v in region)
            sources.setdefault(w, (G, received, region))
        if set(sources) == set(range(-2, 3)):
            break
    assert set(sources) == set(range(-2, 3))
    naive_crosschecks = []
    for total in (1, 2):
        G = (total,)+(0,)*11
        received, _ = independent_source(G, [(0, 0, 'C')], list(range(30)), [0]*30)
        for initial in (0, 9):
            for failure, partial in [(None, 0), ('discover', 0), ('reset_route', 0), ('measure', 0), ('measure', 3)]:
                meter = reference_execution(initial, total, 2, failure, partial)
                claim = p.assemble_claim(G, received, (0,), meter, 2)
                optimized = admit_record_only(claim, G, (0,), 2, producer=p)
                assert optimized['admitted']
                assert admit_record_only_naive(claim, G, (0,), 2, producer=p) == optimized
                naive_crosschecks.append({'total': total, 'initial': initial, 'outcome': meter['reason'],
                                          'signed_steps': meter['retained_signed_cycle_steps']})
    invariant_cases = 0
    outcomes = Counter()
    for initial, w, meter, columns in examples:
        G, received, region = sources[w]
        claim = p.assemble_claim(G, received, region, meter, 2)
        want = admit_record_only(claim, G, region, 2, producer=p)
        assert want['admitted'] and want['compatible_pointer_initials'] == columns
        weight = F(1, 2**claim['combined_attempts'])
        assert want['compatible_initial_indicator'] == [int(x in columns) for x in range(12)]
        assert want['pointer_initial_classical_effect'] == [str(weight if x in columns else 0) for x in range(12)]
        for value in (None, 'malformed private snapshot', [int(v == (initial+1) % 12) for v in range(12)]):
            changed = deepcopy(claim)
            if value is None:
                del changed['meter']['final_loads']
            else:
                changed['meter']['final_loads'] = value
            changed['pointer_initial'] = 999
            assert admit_record_only(changed, G, region, 2, producer=p) == want
            invariant_cases += 1
        outcomes[meter['reason']] += 1
    # Concrete R=1 ambiguous discovery timeout: no final state oracle.
    G = (1,)+(0,)*11
    received, _ = independent_source(G, [], list(range(30)), [0]*30)
    w = sum(G[v]-p.observe_source(G, received, 1)['L1']['reconstruction'][v] for v in (0,))
    meter = reference_execution(0, w, 1, 'discover')
    claim = p.assemble_claim(G, received, (0,), meter, 1)
    assert tape_only_columns(meter, w, 1) == [0, 9]
    raw0 = p.admit(claim, G, (0,), 1)
    changed = deepcopy(claim); changed['meter']['final_loads'] = [int(v == 9) for v in range(12)]
    raw9 = p.admit(changed, G, (0,), 1)
    assert raw0['compatible_pointer_initials'] == [0] and raw9['compatible_pointer_initials'] == [9]
    fixed = admit_record_only(claim, G, (0,), 1, producer=p)
    assert fixed['compatible_pointer_initials'] == [0, 9] and fixed['readout'] is None
    assert admit_record_only_naive(claim, G, (0,), 1, producer=p) == fixed
    naive_crosschecks.append({'total': 1, 'initial': 'unknown_0_or_9', 'outcome': 'discovery_timeout', 'retries': 1})
    effect = [F(v) for v in fixed['pointer_initial_classical_effect']]
    weight = F(fixed['fine_record_weight_conditional_on_compatible_initial'])
    prior = [F(x+1, 78) for x in range(12)]
    assert sum(prior[x]*effect[x] for x in range(12)) == weight*(prior[0]+prior[9])
    assert admit_record_only(changed, G, (0,), 1, producer=p) == fixed
    del changed['meter']['final_loads']
    raw_missing = p.admit(changed, G, (0,), 1)
    assert not raw_missing['admitted']
    assert admit_record_only(changed, G, (0,), 1, producer=p) == fixed
    reject = {}
    G, received, region = sources[1]
    good = make_claim(G, received, region)
    for name in ['echo', 'bad_E', 'bad_current', 'wrong_total', 'stale_latch', 'bad_weight', 'bad_law', 'missing_sweep']:
        bad = deepcopy(good)
        if name == 'echo':
            bad['meter']['records'] = []; bad['meter']['fine_history_probability'] = '1'
        if name == 'bad_E': bad['sweep_records'][0]['contrast'] += 2
        if name == 'bad_current': bad['sweep_records'][0]['current'] += 1
        if name == 'wrong_total': bad['declared_total'] = 3-sum(G)
        if name == 'stale_latch': bad['charge_read_at'] = 'L0'
        if name == 'bad_weight': bad['fine_joint_probability'] = '1'
        if name == 'bad_law': bad['conditional_success_by_unknown_initial_pointer']['0'] = '1'
        if name == 'missing_sweep': bad['sweep_records'].pop()
        result = admit_record_only(bad, G, region, 2, producer=p)
        assert not result['admitted'], name
        reject[name] = result['reasons']
    timeout = reference_execution(0, 1, 2, 'measure', 3)
    bad = p.assemble_claim(G, received, region, timeout, 2)
    bad['decoded_relative_charge'] = 0
    assert not admit_record_only(bad, G, region, 2, producer=p)['admitted']
    reject['guessed_zero_after_partial_timeout'] = True
    return {'facade_verified_complete_instruments': len(examples), 'outcomes': dict(outcomes),
            'optimized_vs_naive_union_crosschecks': naive_crosschecks,
            'final_snapshot_removal_or_edit_invariance_checks': invariant_cases,
            'rejected_controls': reject, 'pointer_snapshot_defect': {
                'status': 'CONFIRMED_AND_REPAIRED_BY_ADDITIVE_FACADE',
                'identical_observed_pointer_tape_sha256': digest(meter['records']),
                'actual_possible_initials': [0, 9], 'original_admit_final_at_0': raw0,
                'original_admit_final_at_9': raw9, 'original_admit_without_final_snapshot': raw_missing,
                'facade_verdict': fixed, 'explicit_nonuniform_prior_example': [str(x) for x in prior],
                'observed_tape_marginal_under_example_prior': str(sum(prior[x]*effect[x] for x in range(12))),
                'claim': claim}}


def integrity_controls():
    rng = random.Random(61114)
    lost = 0
    for a, b in EDGES:
        for origin in (a, b):
            G = tuple(int(v == origin) for v in range(12))
            coin = 0 if origin == a else 1
            order = list(range(30)); rng.shuffle(order)
            records, audit = p.source_execute(G, [(EDGES.index((a, b)), coin, 'C')],
                [rng.randrange(2) for _ in range(30)], order, prefix_channel=lambda records: [])
            try:
                p.observe_source(G, records, 1)
            except ValueError:
                lost += 1
            else:
                raise AssertionError('Lost moving current escaped faithful closing records')
            delta = [audit['source_L0'][v]-G[v] for v in range(12)]
            assert sum(delta[x] != delta[y] for x, y in EDGES) == 9
    # Every possible equal-total source/reconstruction pair: incidence
    # contrasts distinguish equality. This verifies the finite graph lemma,
    # independently of particular prefixes or sweep words.
    pairs = 0
    states = preparations()
    for n in states:
        for estimate in states:
            if sum(n) != sum(estimate):
                continue
            delta = [a-b for a, b in zip(n, estimate)]
            assert all(delta[a] == delta[b] for a, b in EDGES) == (n == estimate)
            pairs += 1
    cancellation = []
    for name, vertices in [('opposite_currents', (0, 1, 0)), ('closed_six_cycle', CYCLE+(0,))]:
        G = (1,)+(0,)*11
        prefix = []
        for a, b in zip(vertices, vertices[1:]):
            e = tuple(sorted((a, b)))
            prefix.append((EDGES.index(e), 0 if a == e[0] else 1, 'C'))
        received, audit = p.source_execute(G, prefix, [0]*30, prefix_channel=lambda records: [])
        claim = make_claim(G, received, (0,))
        verdict = admit_record_only(claim, G, (0,), 2, producer=p)
        assert verdict['admitted']
        assert claim['latches']['L0']['reconstruction'] == audit['source_L0']
        assert claim['latches']['L1']['reconstruction'] == audit['source_L1']
        actual_attempts = audit['actual_source_attempts']+len(claim['meter']['records'])
        assert actual_attempts-claim['combined_attempts'] == len(prefix)
        cancellation.append({'name': name, 'accepted': True, 'lost_prefix_events': len(prefix),
            'load_reconstruction_correct': True, 'received_attempts': claim['combined_attempts'],
            'actual_emitted_attempts': actual_attempts,
            'received_fine_weight': claim['fine_joint_probability'],
            'actual_emitted_fine_weight': str(F(1, 2**actual_attempts)),
            'historic_authenticity_certified': False})
    malicious = []
    for total in (1, 2):
        for name, order in [('identity', list(range(30))), ('reverse', list(range(29, -1, -1))),
                            ('independent_shuffle', rng.sample(list(range(30)), 30))]:
            G = (total,)+(0,)*11
            def forge(record, reconstruction):
                a, b = record['seam']
                after, j = native(reconstruction, (a, b), record['coin'])
                record['contrast'] = reconstruction[a]-reconstruction[b]
                record['current'] = j
                return record
            received, audit = p.source_execute(G, [(0, 0, 'C')], [0]*30, order,
                prefix_channel=lambda records: [], sweep_channel=forge)
            reconstructed = p.observe_source(G, received, total)['L1']['reconstruction']
            differing = [v for v in range(12) if reconstructed[v] != audit['source_L1'][v]]
            region = (differing[0],) if differing else (0,)
            claim = make_claim(G, received, region)
            result = admit_record_only(claim, G, region, 2, producer=p)
            assert result['admitted'] and claim['latches']['L0']['reconstruction'] != audit['source_L0']
            malicious.append({'total': total, 'order': name, 'admitted': True,
                'reconstruction_L0': claim['latches']['L0']['reconstruction'], 'actual_L0': audit['source_L0'],
                'reconstruction_L1': reconstructed, 'actual_L1': audit['source_L1'],
                'region': list(region), 'admitted_charge': result['readout'],
                'actual_charge': sum(audit['source_L1'][v]-G[v] for v in region),
                'wrong_L1_charge': bool(differing),
                'faithful_closing_E_premise_violated': True})
    assert any(row['wrong_L1_charge'] for row in malicious)
    # AU5 adaptation is explicit: a synthetic final_n authenticates nothing
    # beyond received replay. Full closing-E admission is a separate check.
    G = (2,)+(0,)*11
    pointer = reference_execution(0, 0, 2)
    adapted = p._au5_adapter(G, [], (0,), pointer, 0, 2, [0]*12)
    assert p.strict.validate(adapted)[0] and adapted['audit_only']['final_n'] == list(G)
    missing = deepcopy(adapted); del missing['audit_only']['final_n']
    assert not p.strict.validate(missing)[0]
    return {'lost_moving_current_cases': lost, 'all_equal_total_state_pairs': pairs,
        'cancelling_and_closed_current_controls': cancellation, 'forged_E_counterexamples': malicious,
        'AU5_synthetic_final_n_is_replay_only': True, 'AU5_requires_final_n': True,
        'constant_discrepancy': 'delta=(1,...,1) is contrast-invisible but changes total by 12; algebra only, no out-of-domain execution',
        'required_premises': ['true known G', 'conserved source total one or two', 'faithful closing E channel',
                              'connected sweep', 'declared fair-coin native policy']}


def main():
    before, frozen = pins()
    save('PINNED_INPUTS.json', before)
    print('Pins checked; independent native and geometry checks.', flush=True)
    geometry = geometry_and_kernel()
    print('Independent event-level pointer probability laws.', flush=True)
    pointers = pointer_laws()
    actual, examples = actual_executions()
    print('Independent full source-to-pointer examples and facade controls.', flush=True)
    sources = full_sources()
    facade = facade_controls(examples)
    controls = integrity_controls()
    print('Regenerating all 768 exact joint-map digests independently.', flush=True)
    joint = joint_laws(frozen)
    after, _ = pins()
    assert before == after
    result = {'schema': 'oph.row1.integrity.independent.audit.v1',
        'verdict': 'VALIDATED_WITH_RECORD_ONLY_FACADE',
        'producer_result': 'Exact per-column laws and conditional source integrity validated; raw admission narrows timeout candidates using diagnostic pointer final_loads.',
        'public_admission': 'observer_only.admit_record_only; frozen producer and AU5 files unchanged',
        'geometry_and_kernel': geometry, 'independent_pointer_laws': pointers,
        'actual_native_pointer_executions': actual, 'full_sources': sources,
        'observer_only_facade': facade, 'integrity_controls': controls,
        'complete_joint_laws': joint, 'all_pins_unchanged': True,
        'producer_manifest_sha256': EXPECTED_MANIFEST,
        'reserved_or_load_four_executed': False, 'official_rows_changed': [],
        'boundary': 'Native finite source charge reconstruction at L1 under explicit premises. No historical authentication, autonomous counter sensor, observed particle identification, or official row award.'}
    save('RESULT.json', result)
    save('SUMMARY.json', {'verdict': result['verdict'], 'native_branches': geometry['independent_native_branches'],
        'independent_pointer_laws': pointers['complete_exact_laws'], 'actual_pointer_executions': actual['executions'],
        'full_source_examples': sources['native_source_to_L1_pointer_examples'],
        'facade_full_instruments': facade['facade_verified_complete_instruments'],
        'final_snapshot_invariance_checks': facade['final_snapshot_removal_or_edit_invariance_checks'],
        'exact_joint_maps': joint['independently_regenerated_joint_maps'],
        'source_discrepancy_pairs': controls['all_equal_total_state_pairs'],
        'lost_moving_current_cases': controls['lost_moving_current_cases'],
        'forged_E_counterexamples': len(controls['forged_E_counterexamples']),
        'probability_mass_error': '0', 'all_pins_unchanged': True,
        'unresolved_producer_admission_caveat': 'Use observer_only.py facade for record-only candidate sets.'})
    print((HERE/'SUMMARY.json').read_text(), flush=True)


if __name__ == '__main__':
    main()
