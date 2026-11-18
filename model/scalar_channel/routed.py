"""Execute local ancilla transport, loop instruments and retained feed-forward.

The primitive domain is a star and a neighbor SWAP, with supplied classical
control. The simulation evolves all code and ancilla stabilizers. Only a
single local ancilla is measured; no multi-star Pauli measurement is primitive.
"""
from __future__ import annotations

from collections import Counter
import numpy as np

from .compiler import (Tableau, build, fundamental_cycles, local_cycles,
                       routed_measurement, bits, path)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def clifford(tab, operation, a, b=None):
    """Conjugate i**p X**x Z**z stabilizers by an explicit local gate."""
    require(0 <= a < tab.q, 'Gate address outside hardware')
    if operation in ('CX', 'SWAP'):
        require(b is not None and 0 <= b < tab.q and a != b,
                'Two-qubit gate needs distinct hardware addresses')
    if operation == 'SWAP':
        clifford(tab, 'CX', a, b)
        clifford(tab, 'CX', b, a)
        clifford(tab, 'CX', a, b)
        return
    rows = []
    for p, x, z in tab.rows:
        xa, za = (x >> a) & 1, (z >> a) & 1
        if operation == 'H':
            p += 2 * xa * za
            if xa != za:
                x ^= 1 << a
                z ^= 1 << a
        elif operation in ('S', 'SDG'):
            p += xa if operation == 'S' else -xa
            if xa:
                z ^= 1 << a
        elif operation == 'CX':
            if xa:
                x ^= 1 << b
            if (z >> b) & 1:
                z ^= 1 << a
        else:
            raise ValueError('Unsupported Clifford gate')
        rows.append((p % 4, x, z))
    tab.rows = rows


def controlled_factor(tab, ancilla, target, x, z):
    require((x, z) in ((1, 0), (0, 1), (1, 1)), 'Empty Pauli factor')
    if (x, z) == (1, 0):
        clifford(tab, 'CX', ancilla, target)
    elif (x, z) == (0, 1):
        clifford(tab, 'H', target)
        clifford(tab, 'CX', ancilla, target)
        clifford(tab, 'H', target)
    else:
        clifford(tab, 'SDG', target)
        clifford(tab, 'CX', ancilla, target)
        clifford(tab, 'S', target)


def tableau_record(tab):
    """Lossless pure stabilizer state, including unfinished apparatus bits."""
    return dict(qubits=tab.q,
                stabilizers=[[p, hex(x), hex(z)] for p, x, z in tab.rows])


def measure_loop(tab, graph, chord, pauli, sign, *, fault=None,
                 omit_correction=False, diagnostic_signs=None,
                 diagnostic_seed=0):
    """One complete branch; faults retain its unfinished quantum state.

    Sign is a sampled local pointer outcome, with its exact probability
    returned. 'missing_ack' and 'wrong_ack' are external fault controls at the
    first transport. A full controller uses the retained record to proceed
    or halt. Neither fault authorizes admission or resets the traveling bit.
    Scratch preparation is checked by sequential local Z measurements. A
    negative diagnostic is retained and halts without reset. Explicit
    diagnostic_signs select a branch; otherwise diagnostic_seed samples it.
    Deterministic diagnostics consume no random draws. The returned branch
    probability includes every diagnostic and the final pointer, if reached.
    """
    require(type(sign) is int and sign in (-1, 1), 'A pointer has exactly two outcomes')
    require(fault in (None, 'missing_ack', 'wrong_ack'), 'Unknown fault')
    M, N = graph['M'], graph['N']
    require(type(diagnostic_seed) is int and diagnostic_seed >= 0,
            'Diagnostic seed must be a nonnegative integer')
    require(diagnostic_signs is None or
            (len(diagnostic_signs) == N and
             all(type(s) is int and s in (-1, 1) for s in diagnostic_signs)),
            'Diagnostics require one sign per scratch bit')
    require(tab.q == M + N, 'Hardware must include one ancilla per star')
    require(0 <= chord < M, 'Chord outside graph')
    require(pauli[1] >> M == 0 and pauli[2] >> M == 0,
            'Data measurement cannot target an ancilla')
    p, x, z = pauli
    phase = (p - (x & z).bit_count()) % 4
    require(phase in (0, 2), 'Measured Pauli must be Hermitian')
    center = graph['owner'][chord]
    route = routed_measurement(graph, pauli, center)
    edges = {tuple(sorted(map(int, e))) for e in graph['data']['edges']}
    require(all(tuple(sorted((a, b))) in edges for a, b in zip(route, route[1:])),
            'Ancilla route crosses a non-neighbor')
    diagnostics, trace = [], []
    diagnostic_probability = 1.0
    rng = np.random.default_rng(diagnostic_seed)
    for star in range(N):
        observable = (0, 0, 1 << (M + star))
        expected = tab.expectation(observable)
        observed = (diagnostic_signs[star] if diagnostic_signs is not None else
                    expected if expected else int(rng.choice((-1, 1))))
        # A certain local read is the identity on this state. Avoid rebuilding
        # the stabilizer basis for the deterministic blank-preparation case.
        if expected:
            require(observed == expected, 'Impossible requested diagnostic')
            probability = 1.0
        else:
            probability = tab.measure(observable, observed)
        diagnostic_probability *= probability
        diagnostics.append(dict(star=star, outcome=observed, probability=probability))
        trace.append(['diagnostic_z', star, observed])
        if observed < 0:
            return dict(status='failed', reason='dirty_ancilla', outcome=None,
                        probability=diagnostic_probability, diagnostics=diagnostics,
                        route=[], gates=trace, reset=False, corrected=False)
    clifford(tab, 'H', M + center)
    trace.append(['h', center])
    if phase == 2:
        tab.pauli((0, 0, 1 << (M + center)))
        trace.append(['z', center])
    visited, applied = set(), set()
    current = center
    for step, star in enumerate(route):
        if step:
            if fault:
                return dict(status='failed', reason=fault, outcome=None,
                            probability=diagnostic_probability, diagnostics=diagnostics,
                            route=route[:step],
                            gates=trace, reset=False, corrected=False)
            clifford(tab, 'SWAP', M + current, M + star)
            trace.append(['swap', current, star])
            current = star
        if star not in visited:
            for k in bits(x | z):
                if graph['owner'][k] == star:
                    controlled_factor(tab, M + star, k, (x >> k) & 1, (z >> k) & 1)
                    trace.append(['controlled_pauli', star, k, (x >> k) & 1, (z >> k) & 1])
                    applied.add(k)
            visited.add(star)
    require(applied == set(bits(x | z)) and current == center,
            'Incomplete local measurement route')
    clifford(tab, 'H', M + center)
    trace.append(['h', center])
    probability = tab.measure((0, 0, 1 << (M + center)), sign)
    trace.append(['pointer_z', center, sign])
    # The sign is retained before feedback or ancilla reinitialization.
    record = dict(status='admitted', outcome=sign,
                  probability=diagnostic_probability * probability,
                  pointer_probability=probability, diagnostics=diagnostics,
                  route=route, gates=trace, reset=sign < 0,
                  corrected=sign < 0 and not omit_correction)
    if sign < 0:
        if not omit_correction:
            tab.pauli((0, 0, 1 << chord))
            trace.append(['feedback_z', chord])
        tab.pauli((0, 1 << (M + center), 0))
        trace.append(['reset_x', center])
    return record


def run(seed=193019):
    graph = build()
    tab = Tableau(graph['M'] + graph['N'])
    fundamental, _ = fundamental_cycles(graph)
    rng = np.random.default_rng(seed)
    records = []
    visits = Counter()
    for chord, _, S in fundamental:
        rec = measure_loop(tab, graph, chord, S, int(rng.choice((-1, 1))))
        require(rec['probability'] == .5, 'Loop branch probability differs from one half')
        require(tab.expectation(S) == 1, 'Feed-forward failed to prepare loop')
        visits.update(rec['route'])
        records.append(rec)
    require(all(tab.expectation(P) == 1 for P in graph['B']), 'Vacuum numbers changed')
    require(all(tab.expectation(r['pauli']) == 1 for r in local_cycles(graph)),
            'Local penalty loop is not prepared')
    require(all(tab.expectation((0, 0, 1 << (graph['M'] + b))) == 1
                for b in range(graph['N'])), 'Ancilla leaked or was not returned')
    # Occupy the physical transit reference by local even bilinears.
    occupied = list(range(graph['E'], graph['nsource']))
    pair_operations = 0
    for i, j in zip(occupied[::2], occupied[1::2]):
        route = path(graph, i, j)
        for a, b in zip(route, route[1:]):
            tab.pauli(graph['oriented'](a, b))
            pair_operations += 1
    require([tab.expectation(P) for P in graph['B']] ==
            [-1 if i in occupied else 1 for i in range(graph['n'])],
            'Transit preparation changed a complementary occupation')
    # A failed transport has a retained branch, no completed syndrome or reset.
    fault_results = []
    chord, _, S = fundamental[0]
    for fault in ('missing_ack', 'wrong_ack'):
        trial = Tableau(tab.q)
        rec = measure_loop(trial, graph, chord, S, 1, fault=fault)
        require(rec['status'] == 'failed' and not rec['reset'], 'Fault was admitted')
        rec['halted_state'] = tableau_record(trial)
        fault_results.append(rec)
    trial = Tableau(tab.q)
    trial.pauli((0, 1 << graph['M'], 0))
    dirty = measure_loop(trial, graph, chord, S, 1)
    require(dirty['reason'] == 'dirty_ancilla' and not dirty['reset'],
            'Dirty scratch was admitted or erased')
    dirty['halted_state'] = tableau_record(trial)
    # Competing control: identical hardware and measured sign without feedback.
    trial = Tableau(tab.q)
    bad = measure_loop(trial, graph, chord, S, -1, omit_correction=True)
    uncorrected = trial.expectation(S)
    require(uncorrected == -1, 'Syndrome detector missed omitted correction')
    return dict(schema='oph.routed_loop_preparation.v1', seed=seed,
                data_qubits=graph['M'], local_ancillas=graph['N'],
                retained_sign_records=records,
                retained_sign_bits=len(records),
                diagnostic_measurements=sum(len(r['diagnostics']) for r in records),
                pointer_measurements=len(records),
                transport_hops=sum(len(r['route']) - 1 for r in records),
                visits_by_star={str(k): v for k, v in sorted(visits.items())},
                pair_operations=pair_operations,
                all_local_loops_plus=True, all_transit_numbers_correct=True,
                ancillas_returned_zero=True, failure_controls=fault_results,
                dirty_ancilla_control=dirty, final_state=tableau_record(tab),
                no_feedback_control=dict(record=bad, detected_sign=uncorrected),
                quantum_measurements='local ancilla Z only',
                controller='supplied finite local program',
                dense_quantum_state_materialized=False)
