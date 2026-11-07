"""Independent endpoint, complete instrument and backaction controls."""
from dataclasses import replace
from fractions import Fraction as F

import numpy as np
import pytest
import sympy as s

from . import response as r
from model.observer import native
from model.local_repair_family import local_bias as local


@pytest.fixture(scope='module')
def data():
    return r.data()


def independent_rows():
    states = tuple(native.states(binary=True)); index = {state: i for i, state in enumerate(states)}
    result = [('hold', None, None, None, None, F(), ())]
    theta = (25, 25)+(8,)*10
    for edge in local.SEAMS:
        a, b = edge
        for law in ('C', 'E'):
            for coin in (0, 1):
                slots = {raw: [] for raw in range(4 if law == 'C' else 8)}
                weight = F(theta[a if coin else b], 60*(theta[a]+theta[b]))
                for i, state in enumerate(states):
                    total = state[a]+state[b]; output = list(state)
                    output[a], output[b] = (total+coin)//2, (total+1-coin)//2
                    current = state[a]-output[a]
                    raw = {0: 0, 1: 1, -1: 2}[current]+(4*(total % 2) if law == 'E' else 0)
                    slots[raw].append((i, index[tuple(output)]))
                for raw, pairs in slots.items():
                    result.append(('repair', law, edge, coin, raw, weight, tuple(pairs)))
    return result


def test_all_native_outcomes_against_independent_endpoint_arithmetic(data):
    actual = [(x.kind, x.law, x.seam, x.coin, x.raw, x.weight, x.pairs) for x in data['branches']]
    assert actual == independent_rows()
    assert len(actual) == 721
    assert sum(not row[-1] for row in actual) == 421
    assert sum(x.kind == 'repair' and x.label is None for x in data['branches']) == 240
    assert r.validate(data)


def test_source_definition_cannot_be_reweighted(monkeypatch):
    original = r.interpreter.program
    def changed():
        value = original(); value['parameters']['port_weights'] = [24, 24]+[8]*10
        return value
    monkeypatch.setattr(r.interpreter, 'program', changed)
    with pytest.raises(ArithmeticError, match='authoritative native member'):
        r.source()


def test_zero_raw_slots_cannot_be_dropped(data):
    changed = dict(data); changed['branches'] = data['branches'][1:]
    with pytest.raises(ArithmeticError, match='slots were dropped'):
        r.validate(changed)


def test_native_response_includes_the_measuring_repair(data):
    changed = dict(data)
    changed['B'] = tuple(s.diag(*data['D'][e, :]) for e in range(30))
    with pytest.raises(ArithmeticError, match='backaction'):
        r.validate(changed)
    equal = [e for e, (a, b) in enumerate(local.SEAMS) if data['parameters'].port_weights[a] == data['parameters'].port_weights[b]]
    assert len(equal) == 22
    for lag in (1, 2, 5):
        C = r.covariance(data, lag)
        assert all(C[:, e] == s.zeros(30, 1) for e in equal)
    assert r.covariance(data, 0)[0, 0] > 0
    assert r.covariance(data, 1)[0, 1] == F(5440, 49533)
    passive = (data['D']*data['P']*s.diag(*data['pi'])*data['D'].T)[0, 0]
    assert passive > 0
    assert r.covariance(data, 1)[0, 0] == 0


def test_current_sign_is_checked_on_every_source_column(data):
    changed = dict(data); rows = list(data['branches'])
    at = next(i for i, row in enumerate(rows) if row.pairs and row.current)
    rows[at] = replace(rows[at], current=-rows[at].current)
    changed['branches'] = tuple(rows)
    with pytest.raises(ArithmeticError, match='continuity'):
        r.validate(changed)


def independent_three_record_moment(data, edges):
    total = F(); theta = (25, 25)+(8,)*10
    for initial, pi in zip(data['states'], data['pi'], strict=True):
        paths = [(initial, F(pi))]
        for edge in edges:
            a, b = local.SEAMS[edge]; expanded = []
            for state, weight in paths:
                for coin in (0, 1):
                    out = list(state); endpoint = state[a]+state[b]
                    out[a], out[b] = (endpoint+coin)//2, (endpoint+1-coin)//2
                    # 60 from the score cancels the tagged-seam probability.
                    score_weight = F(theta[a if coin else b], theta[a]+theta[b])*(state[a]-state[b])
                    expanded.append((tuple(out), weight*score_weight))
            paths = expanded
        total += sum(weight for _, weight in paths)
    return total


@pytest.mark.parametrize('edges', [(1, 0, 1), (0, 1, 0), (1, 3, 7), (8, 3, 2)])
def test_ordered_three_time_response_from_independent_paths(data, edges):
    events = tuple((t, edge) for t, edge in enumerate(edges))
    assert r.ordered_moment(data, events) == independent_three_record_moment(data, edges)


def test_quantum_reference_and_all_retained_blocks_match_interpreter(data):
    amplitudes = {(i, a): F((i+2*a) % 5-2, 7) for i in (0, 1, 6, 13, 20, 32, 65) for a in range(2)}
    rho = {(i, a, j, b): x*y for (i, a), x in amplitudes.items() for (j, b), y in amplitudes.items() if x*y}
    program = r.source()
    actual = r.interpreter._outcomes(program, program['commands'][0], rho)
    for row, out in zip(data['branches'], actual, strict=True):
        assert r.apply_branch(row, rho) == out.density
    assert sum(r.interpreter.trace(out.density) for out in actual) == r.interpreter.trace(rho)
    assert r.interpreter.reference_marginal(out.density for out in actual) == r.interpreter.reference_marginal((rho,))
    assert any(any(i != j and value for (i, a, j, b), value in out.density.items()) for out in actual)
    # Committing all flags erases no initial-reference bookkeeping from the full outputs.
    diagonal = {key: value for key, value in rho.items() if key[0] == key[2]}
    for row in data['branches']:
        assert r.interpreter.trace(r.apply_branch(row, rho)) == r.interpreter.trace(r.apply_branch(row, diagonal))


def test_full_response_rank_and_exact_nonclosure(data):
    ranks = r.response_ranks(data)
    assert ranks['observable_rank_by_depth'][-1] == 66
    assert ranks['reachable_rank_by_depth'][10] == 65
    h = r.hierarchy(data)
    assert h['mean_dimension_before'] == 12
    assert h['mean_dimension_after_one_attempt'] == 19
    assert F(h['drift_contrast']) == F(17, 990)
    for port in range(12):
        assert sum(sign*int(port in pair) for sign, pair in zip(h['signs'], h['pairs'], strict=True)) == 0


def test_characteristic_polynomial_against_independent_modular_power_traces(data):
    packet = r.spectral_certificate(data)
    # Newton identities in a finite field independently check every coefficient.
    prime = 65521
    C = np.array(1980*data['P'], dtype=np.int64) % prime
    power = np.eye(66, dtype=np.int64); traces = [0]
    for _ in range(66):
        power = power@C % prime; traces.append(int(np.trace(power)) % prime)
    coefficients = [1]
    for k in range(1, 67):
        coefficients.append((-sum(coefficients[k-j]*traces[j] for j in range(1, k+1))*pow(k, -1, prime)) % prime)
    assert coefficients == [int(v) % prime for v in packet['integer_transition_characteristic_coefficients']]
    assert packet['distinct_centered_record_response_poles'] == 65
    assert packet['physical_mass_interpretation'] is False


def test_exact_resolvent_solves_native_recurrence(data):
    z = F(2)
    response = r.resolvent(data, z)
    # Sum twelve coefficients, then compute the exact remaining native recurrence.
    partial = s.zeros(30)
    for k in range(1, 13):
        partial += r.covariance(data, k)/z**k
    A = s.diag(*data['pi'])*data['G'].T-data['pi']*data['mean'].T
    tail = data['D']*(data['P']**12)*(z*s.eye(66)-data['P']).inv()*A/z**12
    assert response == partial+tail


def test_preparation_bound_and_no_physical_premise_claim():
    packet = r.preparation_error()
    assert packet['attempts'] == 1600
    assert F(packet['half_diamond_upper']) < F(1, 10**6)
    assert packet['new_preparation_archive_in_replacement_target'] is False


@pytest.mark.parametrize('events', [(), ((0, 0), (0, 1)), ((True, 0),), ((0, -1),), [(0, 0)]])
def test_reject_invalid_event_domain(data, events):
    with pytest.raises(ValueError):
        r.ordered_moment(data, events)


@pytest.mark.parametrize('as_string', [False, True])
def test_receipt_checker_accepts_path_and_cli_string(tmp_path, monkeypatch, as_string):
    import json
    packet = {'source_sha256': {'source.py': 'abc'}, 'full_raw_slots': 721}
    path = tmp_path/'receipt.json'; path.write_text(json.dumps(packet))
    monkeypatch.setattr(r, 'report', lambda: packet)
    assert r.check(str(path) if as_string else path) == packet
    path.write_text(json.dumps({**packet, 'full_raw_slots': 720}))
    with pytest.raises(ValueError, match='receipt or source bindings differ'):
        r.check(str(path) if as_string else path)


def test_receipt_checker_rejects_duplicate_json_keys(tmp_path, monkeypatch):
    path = tmp_path/'receipt.json'; path.write_text('{"full_raw_slots":721,"full_raw_slots":720}')
    monkeypatch.setattr(r, 'report', lambda: {'full_raw_slots': 720})
    with pytest.raises(ValueError):
        r.check(path)
