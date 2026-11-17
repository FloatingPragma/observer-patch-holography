"""Native source records, closing E sweep, finite counters and actual pointer.

Source state exists only in source_execute(). Observer admission receives the
known preparation and records; it has no source-terminal-load input.
"""
from __future__ import annotations
from collections import defaultdict
from fractions import Fraction
import closing_sweep as sweep
import strict_tape_validator as strict
import native_cycle_meter as meter
from native_kernel import native

EDGES = tuple(meter.EDGES)
INDEX = {edge: i for i, edge in enumerate(EDGES)}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def preparation(G):
    require(isinstance(G, (list, tuple)) and len(G) == 12 and
            all(type(x) is int and x >= 0 for x in G) and sum(G) in (1, 2),
            'I0: twelve nonnegative integers with total one or two required')
    return tuple(G)


def region_set(region):
    require(isinstance(region, (list, tuple)) and
            all(type(v) is int and 0 <= v < 12 for v in region) and
            len(set(region)) == len(region), 'I0: region must be a set of ports')
    return tuple(region)


def lift(residue, total):
    values = {q % 6: q for q in range(-total, total + 1)}
    require(type(residue) is int and residue in values, 'I2: counter outside admitted load interval')
    return values[residue]


def source_execute(G, prefix, sweep_coins, order=None, *, prefix_channel=None, sweep_channel=None):
    """Execute the source. Channels are adversarial diagnostic controls only.

    The returned audit dictionary is separate from the observer's records.
    Frozen closing_sweep executes the thirty actual E events unchanged.
    """
    G = preparation(G)
    state = G
    emitted = []
    for position, (seam_index, coin, kind) in enumerate(prefix):
        require(type(seam_index) is int and 0 <= seam_index < 30 and
                type(coin) is int and coin in (0, 1) and kind in ('C', 'E'),
                'I0: source event alphabet')
        a, b = EDGES[seam_index]
        contrast = state[a] - state[b]
        state, (_, current) = native(state, (a, b), coin)
        record = {'kind': kind, 'seam': seam_index, 'coin': coin, 'current': current}
        if kind == 'E':
            record['contrast'] = contrast
        emitted.append(record)
    retained = [dict(r) for r in emitted]
    if prefix_channel is not None:
        retained = prefix_channel(retained)
    order = list(range(30)) if order is None else list(order)
    result = sweep.closing_sweep(G, retained, state, order, sweep_coins, channel=sweep_channel)
    records = []
    for rec in result['records']:
        records.append({'kind': 'E', 'seam': INDEX[tuple(rec['seam'])], 'coin': rec['coin'],
                        'current': rec['current'], 'contrast': rec['contrast']})
    received = {'source_records': retained, 'sweep_records': records, 'sweep_order': order}
    audit = {'initial_G': list(G), 'source_L0': list(state), 'source_L1': result['source_after'],
             'emitted_source_records': emitted, 'frozen_sweep_report': result,
             'actual_source_attempts': len(prefix) + 30}
    return received, audit


def observe_source(G, received, total):
    """Only G, declared total and received records enter this verification."""
    G = preparation(G)
    require(type(total) is int and total == sum(G), 'I0: declared total differs from known preparation')
    require(isinstance(received, dict), 'I0: record envelope required')
    prefix = received.get('source_records')
    closing = received.get('sweep_records')
    order = received.get('sweep_order')
    require(isinstance(prefix, list) and isinstance(closing, list), 'I0: both record tapes required')
    require(isinstance(order, list) and all(type(e) is int for e in order) and
            sorted(order) == list(range(30)), 'I1: closing sweep must address every seam exactly once')
    require(len(closing) == 30, 'I1: thirty closing E records required')
    require(sweep.connected_and_spanning([EDGES[e] for e in order]), 'I1: disconnected closing sweep')
    reconstruction = list(G)
    counters = [0] * 12
    latches = {}
    for phase, records in (('L0', prefix), ('L1', closing)):
        for k, rec in enumerate(records):
            require(isinstance(rec, dict), 'I0: event record required')
            e, c, j = rec.get('seam'), rec.get('coin'), rec.get('current')
            require(type(e) is int and 0 <= e < 30 and type(c) is int and c in (0, 1) and
                    type(j) is int and rec.get('kind') in ('C', 'E'), 'I0: event alphabet')
            a, b = EDGES[e]
            predicted = reconstruction[a] - reconstruction[b]
            if phase == 'L1':
                require(rec['kind'] == 'E' and e == order[k], 'I1: sweep event is absent or unconnected')
            if rec['kind'] == 'E':
                d = rec.get('contrast')
                require(type(d) is int, 'I0: E contrast must be an integer')
                require(j == sweep.current_of(d, c), 'I3: non-native E current')
                require(d == predicted, 'I3: closing or prefix E contrast disagrees with reconstruction')
            else:
                require(j == sweep.current_of(predicted, c), 'I3: non-native C current')
            # Apply the RECEIVED current. No source load is accessed here.
            reconstruction[a] -= j
            reconstruction[b] += j
            counters[a] = (counters[a] + j) % 6
            counters[b] = (counters[b] - j) % 6
            require(min(reconstruction) >= 0 and sum(reconstruction) == total,
                    'I3: reconstructed history outside source sector')
        latches[phase] = {'counters': list(counters), 'reconstruction': list(reconstruction),
                          'record_count': len(prefix) + (30 if phase == 'L1' else 0)}
    return latches


def laws(total, commanded, retries):
    require(type(retries) is int and 1 <= retries <= 128, 'I0: bounded retry count required')
    return {str(x): str(strict.reset_success(x, retries) *
                       (1 - Fraction(1, 2 ** retries)) ** (6 * abs(commanded))) for x in range(12)}


def assemble_claim(G, received, region, pointer_result, retries):
    """Construct the observer envelope exclusively from received records."""
    G = preparation(G)
    region = region_set(region)
    latches = observe_source(G, received, sum(G))
    classes = {l: sum(latches[l]['counters'][v] for v in region) % 6 for l in ('L0', 'L1')}
    command = lift(classes['L1'], sum(G))
    attempts = len(received['source_records']) + len(received['sweep_records']) + len(pointer_result['records'])
    success = pointer_result['success'] is True
    return {**received, 'schema': 'oph.row1.integrity.instrument.v1', 'declared_total': sum(G),
            'latches': latches, 'latched_counter_classes': classes, 'charge_read_at': 'L1',
            'meter': pointer_result, 'counter_class': classes['L1'],
            'decoded_spatial_degree': pointer_result['winding'] if success else None,
            'decoded_relative_charge': -pointer_result['winding'] if success else None,
            'conditional_success_by_unknown_initial_pointer': laws(sum(G), command, retries),
            'combined_attempts': attempts, 'fine_joint_probability': str(Fraction(1, 2 ** attempts)),
            'source_attempts': len(received['source_records']) + 30,
            'sweep_attempts': 30, 'pointer_attempt_cap': (30 + 3 + 6 * sum(G)) * retries}


def execute(G, prefix, region, pointer_initial, pointer_coins, sweep_coins, retries=2, order=None):
    """A full realization. pointer_coins can be an iterator or a callable(w)."""
    received, audit = source_execute(G, prefix, sweep_coins, order)
    latches = observe_source(G, received, sum(G))
    command = lift(sum(latches['L1']['counters'][v] for v in region) % 6, sum(G))
    coins = pointer_coins(command) if callable(pointer_coins) else pointer_coins
    pointer_result = meter.execute(pointer_initial, (command,), coins, retries)
    claim = assemble_claim(G, received, region, pointer_result, retries)
    audit['pointer_initial'] = pointer_initial
    return claim, audit


def _au5_adapter(G, records, region, pointer_result, initial, retries, counters):
    """Frozen two-load schema, with synthetic final_n from G and records only.

    No executor audit field is copied. This is a schema adapter, not a claim
    that AU5's terminal-load comparison independently authenticates a source.
    """
    final = sweep.reconstruct(G, records)
    residue = sum(counters[v] for v in region) % 6
    command = lift(residue, 2)
    success = pointer_result['success'] is True
    attempts = len(records) + len(pointer_result['records'])
    p = strict.reset_success(initial, retries) * (1 - Fraction(1, 2 ** retries)) ** (6 * abs(command))
    result = {'source_records': records, 'counter_class': residue, 'meter': pointer_result,
              'decoded_spatial_degree': pointer_result['winding'] if success else None,
              'decoded_relative_charge': -pointer_result['winding'] if success else None,
              'conditional_meter_success': str(p), 'conditional_meter_timeout': str(1-p),
              'combined_attempts': attempts, 'fine_joint_probability': str(Fraction(1, 2 ** attempts)),
              'audit_only': {'initial_G': list(G), 'final_n': list(final),
                             'expected_w': sum(G[v] - final[v] for v in region)}}
    return strict.with_inputs(result, G, [r['seam'] for r in records], [r['coin'] for r in records],
                              region, initial, retries)


def admit(claim, G, region, retries):
    """Observer admission. No source state or executor audit argument exists.

    Pointer initial position is unknown: replay every one-load column and
    retain exactly the compatible columns. Native tape validation is mandatory.
    """
    reasons = []
    try:
        G = preparation(G)
        region = region_set(region)
        require(isinstance(claim, dict), 'I0: claim envelope required')
        latches = observe_source(G, claim, claim.get('declared_total'))
        classes = {l: sum(latches[l]['counters'][v] for v in region) % 6 for l in ('L0', 'L1')}
        command = lift(classes['L1'], sum(G))
        require(claim.get('latches') == latches and claim.get('latched_counter_classes') == classes and
                claim.get('counter_class') == classes['L1'], 'I2: counters or latches disagree with received tape')
        require(claim.get('charge_read_at') == 'L1', 'I2: pointer readout must refer to L1')
        expected_laws = laws(sum(G), command, retries)
        pointer = claim.get('meter')
        require(isinstance(pointer, dict) and isinstance(pointer.get('records'), list), 'I0: pointer tape required')
        candidates = []
        pointer_reasons = []
        for initial in range(12):
            errors, replay = strict.validate_pointer_tape(pointer, initial, command, retries)
            if not errors:
                candidates.append(initial)
            else:
                pointer_reasons.extend(errors)
        if not candidates:
            return {'admitted': False, 'reasons': ['I4: no unknown initial pointer column validates'] +
                    sorted(set(pointer_reasons)), 'compatible_pointer_initials': []}
        if sum(G) == 2:
            records = claim['source_records'] + claim['sweep_records']
            for initial in candidates:
                ok, errors = strict.validate(_au5_adapter(G, records, region, pointer, initial,
                                                          retries, latches['L1']['counters']))
                require(ok, 'I4: frozen AU5 adapter: ' + '; '.join(errors))
        success = pointer['success'] is True
        want = pointer['winding'] if success else None
        require(claim.get('decoded_spatial_degree') == want and
                claim.get('decoded_relative_charge') == (-want if success else None),
                'I5: pointer outcome or timeout readout mismatch')
        attempts = len(claim['source_records']) + 30 + len(pointer['records'])
        require(type(claim.get('combined_attempts')) is int and claim['combined_attempts'] == attempts and
                claim.get('fine_joint_probability') == str(Fraction(1, 2 ** attempts)) and
                claim.get('source_attempts') == len(claim['source_records']) + 30 and
                claim.get('sweep_attempts') == 30 and
                claim.get('pointer_attempt_cap') == (30 + 3 + 6 * sum(G)) * retries,
                'I6: attempts, sweep clock or fine weights disagree')
        require(claim.get('conditional_success_by_unknown_initial_pointer') == expected_laws,
                'I7: conditional pointer laws disagree')
        return {'admitted': True, 'reasons': [], 'compatible_pointer_initials': candidates,
                'charge_read_at': 'L1', 'readout': -want if success else None,
                'outcome': pointer['reason'], 'record_attempts': attempts}
    except (ValueError, TypeError, KeyError, IndexError, OverflowError) as error:
        reasons.append(str(error) if str(error).startswith(('I', 'R')) else 'I0: malformed record: ' + str(error))
        return {'admitted': False, 'reasons': reasons, 'compatible_pointer_initials': []}


def sweep_distribution(initial, order=None):
    """Exact weighted finite-state recursion; no literal 2^30-word enumeration."""
    initial = preparation(initial)
    order = list(range(30)) if order is None else order
    counts = {initial: 1}
    for e in order:
        nxt = defaultdict(int)
        for state, multiplicity in counts.items():
            for coin in (0, 1):
                after, _, _ = sweep.native_event(state, EDGES[e], coin)
                nxt[after] += multiplicity
        counts = dict(nxt)
    assert sum(counts.values()) == 2 ** len(order)
    return {state: Fraction(n, 2 ** len(order)) for state, n in counts.items()}


def pointer_distribution(initial, commanded, retries):
    """Exact actual-native-policy quotient by stage, port and attempt count.

    Coin records are distinct fine leaves. Equal future states have their
    rational masses added; all timeout stages retain unconditional weight.
    """
    require(type(initial) is int and 0 <= initial < 12, 'I0: pointer column')
    require(type(commanded) is int and -2 <= commanded <= 2, 'I0: pointer command')
    require(type(retries) is int and 1 <= retries <= 128, 'I0: pointer retries')
    outputs = defaultdict(Fraction)
    discovered = defaultdict(Fraction)
    state = tuple(int(v == initial) for v in range(12))
    pending = Fraction(1)
    attempts = 0
    for edge in EDGES:
        for _ in range(retries):
            hold = Fraction()
            for coin in (0, 1):
                after, (_, j) = native(state, edge, coin)
                if j:
                    discovered[(after.index(1), attempts + 1)] += pending / 2
                else:
                    assert after == state
                    hold += pending / 2
            attempts += 1
            pending = hold
    outputs[('discovery_timeout', None, 0, attempts)] += pending
    for (port, used), weight in discovered.items():
        active = {used: weight}
        path = list(meter.route(port, 0))
        cursor = 0
        direction = 1 if commanded > 0 else -1
        commands = [(dest, 'routing_timeout', 0) for dest in path]
        for step in range(6 * abs(commanded)):
            cursor = (cursor + direction) % 6
            commands.append((meter.CYCLE[cursor], 'measurement_timeout', direction * step))
        known = port
        for dest, failure, signed_before in commands:
            edge = tuple(sorted((known, dest)))
            loads = tuple(int(v == known) for v in range(12))
            branches = [native(loads, edge, coin) for coin in (0, 1)]
            assert sum(bool(rec[1]) for _, rec in branches) == 1
            assert all(after.index(1) == (dest if rec[1] else known) for after, rec in branches)
            advanced = defaultdict(Fraction)
            for count, mass in active.items():
                for k in range(1, retries + 1):
                    advanced[count + k] += mass * Fraction(1, 2 ** k)
                outputs[(failure, None, signed_before, count + retries)] += mass * Fraction(1, 2 ** retries)
            active = dict(advanced)
            known = dest
        for count, mass in active.items():
            outputs[('complete_closed_pointer_word', commanded, 6 * commanded, count)] += mass
    outputs = {key: value for key, value in outputs.items() if value}
    assert sum(outputs.values()) == 1
    expected = strict.reset_success(initial, retries) * (1 - Fraction(1, 2 ** retries)) ** (6 * abs(commanded))
    assert sum(p for key, p in outputs.items() if key[0] == 'complete_closed_pointer_word') == expected
    return outputs
