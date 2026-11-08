"""Complete branch, reference, current, erasure and source-attachment controls."""
from dataclasses import replace
from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest

from . import history as h


def trace(p, edges, coins, initial):
    state = initial
    events = []
    snapshots = [(0, sum(((initial >> v) & 1) << i for i, v in enumerate(p['ports'])))]
    for tick, (edge, coin) in enumerate(zip(edges, coins, strict=True), 1):
        a, b = edge
        rec = h.r.record(edge, coin, ((state >> a) & 1) - ((state >> b) & 1))
        state += rec.current * ((1 << b) - (1 << a))
        if edge in p['incident_seams']:
            events.append((tick, rec))
            snapshots.append((tick, sum(((state >> v) & 1) << i for i, v in enumerate(p['ports']))))
    return tuple(events), tuple(snapshots), state


@pytest.mark.parametrize('size', [2, 3, 4, 5])
def test_every_occupation_sector_all_coin_histories(size):
    ports = tuple(range(size))
    incident = tuple((i, i+1) for i in range(size-1))
    p = h.r.plan(ports, incident)
    for initial in range(1 << size):
        for coins in product((0, 1), repeat=size-1):
            events, expected, _ = trace(p, p['tree'], coins, initial)
            observed = h.observe(p, events, size-1)
            for weights in (tuple(F(1) for _ in ports), tuple(F(i-1, 3) for i in ports)):
                result = h.ward(p, observed, weights)
                if initial in (0, (1 << size)-1):
                    assert result['status'] == 'CONSTANT_ERASURE'
                    assert observed['snapshots'] == ()
                else:
                    assert observed['snapshots'] == expected
            assert h.ward(p, observed, tuple(F(1) for _ in ports))['net_change'] == 0


def test_boundary_flux_all_inputs_and_prefixes():
    p = h.r.plan((1, 2, 3), ((0, 1), (1, 2), (2, 3), (3, 4)))
    edges = ((0, 1), (3, 4), (1, 2), (0, 1), (2, 3))
    for initial in range(32):
        for coins in product((0, 1), repeat=len(edges)):
            events, expected, last = trace(p, edges, coins, initial)
            observed = h.observe(p, events, len(edges))
            for subset in range(8):
                weights = tuple(F((subset >> i) & 1) for i in range(3))
                result = h.ward(p, observed, weights)
                direct = sum((weights[i] * (((last >> v) & 1) - ((initial >> v) & 1))
                              for i, v in enumerate(p['ports'])), F())
                assert result['net_change'] == direct
            if observed['absolute_history_reconstructed']:
                assert observed['snapshots'] == expected


def test_operator_ward_identity_all_fine_branches_with_reference():
    p = h.r.plan((1, 2), ((0, 1), (1, 2)))
    edges = ((0, 1), (1, 2), (0, 1))
    blocks = {}
    for initial in range(8):
        for coins in product((0, 1), repeat=3):
            events, _, last = trace(p, edges, coins, initial)
            blocks.setdefault(events, []).append((initial, last))
    profile = (F(2), F(-3, 2))
    values = [sum((profile[i] * ((x >> v) & 1) for i, v in enumerate(p['ports'])), F())
              for x in range(8)]
    A = np.diag([float(v) for v in values])
    rng = np.random.default_rng(1141)
    psi = rng.normal(size=(8, 3)) + 1j*rng.normal(size=(8, 3))
    psi /= np.linalg.norm(psi)
    rho = np.outer(psi.ravel(), psi.ravel().conj())
    total_effect = np.zeros((8, 8))
    output_expectation = current_expectation = 0.
    for events, pairs in blocks.items():
        K = np.zeros((8, 8))
        delta = sum((x for _, x in h.current_increments(p, events, 3, profile)), F())
        for x, y in pairs:
            assert values[y] - values[x] == delta
            K[y, x] = 1/np.sqrt(8)
        E = K.T @ K
        assert np.allclose(K.T @ A @ K, E @ A + float(delta) * E, atol=1e-15)
        total_effect += E
        KR = np.kron(K, np.eye(3))
        block = KR @ rho @ KR.T
        output_expectation += np.trace(np.kron(A, np.eye(3)) @ block).real
        current_expectation += float(delta) * np.trace(block).real
        observed = h.observe(p, events, 3)
        if observed['absolute_history_reconstructed']:
            mask = observed['decoded']['initial_mask']
            P = np.diag([int(((x >> 1) & 3) == mask) for x in range(8)])
            assert np.array_equal(K, K @ P)
    assert np.allclose(total_effect, np.eye(8), atol=1e-15)
    assert abs(output_expectation - np.trace(np.kron(A, np.eye(3)) @ rho).real
               - current_expectation) < 1e-13


@pytest.mark.parametrize('total', range(13))
def test_original_source_all_local_number_sectors(total):
    p = h.r.source_plan(0, tuple(range(12)))
    value = h.r.fixed_program(0, total, p)
    initial = (1 << total)-1
    state = initial
    slots = []
    for t, edge in enumerate(p['tree']):
        a, b = edge
        rec = h.r.record(edge, t % 2, ((state >> a) & 1)-((state >> b) & 1), 'fixed')
        outcomes = h.r.source._outcomes(value, value['commands'][t], {})
        slots.append(next(i for i, out in enumerate(outcomes)
                          if h.r.inherited._same_record(out.record, rec)))
        state += rec.current*((1 << b)-(1 << a))
    out = h.query(value, p, tuple(slots), {(initial, 0, initial, 0): F(1)}, use_global_total=True)
    assert out['density'] == {(state, 0, state, 0): F(1, 2**11)}
    assert len(out['archive']) == 10 and len(out['memory']) == 1
    hist = out['charge_history']
    result = h.ward(p, hist, (F(1),)*12)
    assert result['net_change'] == 0
    if total < 12:
        assert hist['absolute_history_reconstructed']
        assert hist['snapshots'][0] == (0, initial)
        assert hist['snapshots'][-1] == (11, state)
    else:
        assert hist['decoded']['status'] == 'CONSTANT_ERASURE'
    zero = h.query(value, p, (0,)*11, {(initial, 0, initial, 0): F(1)}, use_global_total=True)
    assert zero['density'] == {}
    assert len(zero['archive']) == 10 and len(zero['memory']) == 1
    # The hold slot has zero weight at repair probability one. Its syntactic
    # label can be a coverage erasure; no conditional state is normalized.
    assert zero['charge_history']['decoded'] == zero['decoded']


def test_erasure_does_not_delete_observable_flux():
    p = h.r.plan((1, 2), ((0, 1), (1, 2)))
    events, _, _ = trace(p, ((0, 1),), (0,), 1)
    observed = h.observe(p, events, 1)
    assert observed['decoded']['status'] == 'COVERAGE_ERASURE'
    assert h.ward(p, observed, (F(1), F(1)))['net_change'] == 1
    assert observed['snapshots'] == ()


@pytest.mark.parametrize('mutation', ['current', 'order', 'weight', 'boolean'])
def test_invalid_record_and_profile_rejected(mutation):
    p = h.r.plan((0, 1, 2), ((0, 1), (1, 2)))
    events, _, _ = trace(p, p['tree'], (0, 1), 1)
    if mutation == 'current':
        events = ((1, replace(events[0][1], current=3)),) + events[1:]
    elif mutation == 'order':
        events = events[::-1]
    weights = (F(1),)*3 if mutation != 'weight' else (1, 1, 1)
    with pytest.raises(ValueError):
        h.current_increments(p, events, True if mutation == 'boolean' else 2, weights)


def test_missing_boundary_is_rejected_at_source_attachment():
    full = h.r.source_plan(0, tuple(range(12)))
    missing = h.r.plan(full['ports'], tuple(e for e in full['incident_seams'] if e[1] < 12))
    value = h.r.fixed_program(0, 1, full)
    with pytest.raises(ValueError, match='incident-edge'):
        h.query(value, missing, (0,)*11, {(1, 0, 1, 0): F(1)})
