from collections import Counter, defaultdict
from dataclasses import replace
from fractions import Fraction
from itertools import combinations

import numpy as np
import pytest
import sympy as sp

from .._source import ROOT
from . import native_record_dynamics as n
from .certificate import compare_arrays


@pytest.fixture(scope='module', params=('C', 'E'))
def native(request):
    return n.native_instrument(ROOT, request.param)


def test_complete_padded_instrument_and_full_carrier_escape(native):
    expected = 240 if native.metadata['candidate'] == 'C' else 480
    assert len(native.labels) == expected
    assert np.count_nonzero(native.amplitudes) == 3960
    np.testing.assert_allclose(n.effects_diagonal(native).sum(axis=0), 1, atol=2e-15, rtol=0)
    full = n.native_instrument(ROOT, native.metadata['candidate'], binary=False)
    double = np.array([max(x) == 2 for x in full.states])
    surviving = np.zeros(78, dtype=int)
    for row in full.destinations:
        columns = np.flatnonzero(row >= 0)
        surviving[columns] += double[row[columns]]
    np.testing.assert_array_equal(surviving, 50*double)


def test_all_first_record_outputs_match_exact_formula(native):
    packet = n.record_references(native, np.eye(66)/66, ROOT)
    P = packet['native_chart_projectors']
    # Derive the posterior directly from the admitted initial loads.
    for r, (edge, coin, raw) in enumerate(native.labels):
        supported = np.flatnonzero(native.destinations[r] >= 0)
        if not len(supported):
            assert not packet['complete_probabilities'][r].any()
            continue
        outputs = [native.states[native.destinations[r, x]] for x in supported]
        weights = np.array(outputs).sum(axis=0)/(2*len(outputs))
        expected = np.einsum('i,iab->ab', weights, P)
        actual = packet['subnormalized_family_outputs'][r]/packet['raw_record_probabilities'][r]
        np.testing.assert_allclose(actual, expected, atol=2e-15, rtol=0)
        if native.costs[r] == 1:
            a, b = native.seams[edge]
            start, end = (a, b) if raw % 4 == 1 else (b, a)
            np.testing.assert_allclose(actual, np.eye(3)/5+(9*P[end]-P[start])/20,
                                       atol=2e-15, rtol=0)
    np.testing.assert_allclose(packet['subnormalized_family_outputs'].sum(axis=0),
                               np.eye(3)/3, atol=2e-15, rtol=0)


def test_two_record_integer_law_by_independent_endpoint_enumeration(native):
    labels = {label: i for i, label in enumerate(native.labels)}
    expected = np.zeros((len(labels), len(labels)), dtype=np.int64)
    ports = np.zeros(expected.shape+(12,), dtype=np.int64)
    code = {0: 0, 1: 1, -1: 2}
    def steps(x):
        for e, (a, b) in enumerate(native.seams):
            for coin in (0, 1):
                # Direct endpoint balancing, independent of source native().
                out = list(x)
                total = x[a]+x[b]
                out[a] = (total+coin)//2
                out[b] = total-out[a]
                current = x[a]-out[a]
                raw = code[current]
                if native.metadata['candidate'] == 'E':
                    raw += 4*(abs(x[a]-x[b]) % 2)
                yield labels[e, coin, raw], tuple(out)
    first = {x: tuple(steps(x)) for x in native.states}
    for x in native.states:
        for r, y in first[x]:
            for s, z in first[y]:
                expected[r, s] += 1
                ports[r, s] += z
    packet = n.two_record_reference_law(native, ROOT)
    np.testing.assert_array_equal(packet['probability_numerators'], expected)
    np.testing.assert_array_equal(packet['reference_port_numerators'], ports)
    assert int(expected.sum()) == 237600
    # Pair antipodal rays exactly from their projectors, then count support.
    data = n.geometry_packet(ROOT)
    groups = []
    unused = set(range(12))
    while unused:
        i = min(unused)
        group = [j for j in unused if n.is_zero(data['projectors'][i]-data['projectors'][j])]
        assert len(group) == 2
        groups.append(group)
        unused.difference_update(group)
    frame = data['frame']
    for triple in combinations([group[0] for group in groups], 3):
        assert sp.simplify(sp.Matrix([list(frame[i, :]) for i in triple]).det()) != 0
    count = sum(ports[:, :, group].sum(axis=2) > 0 for group in groups)
    rank = np.minimum(count, 3)
    actual = [int(expected[rank == k].sum()) for k in (1, 2, 3)]
    assert actual == ([300, 2940, 234360] if native.metadata['candidate'] == 'C'
                      else [1200, 11760, 224640])


def test_clock_variance_from_full_66_state_poisson_equation():
    native = n.native_instrument(ROOT, 'C')
    T, A = sp.zeros(66), sp.zeros(66)
    for r, row in enumerate(native.destinations):
        for x in np.flatnonzero(row >= 0):
            T[int(row[x]), int(x)] += sp.Rational(1, 60)
            if native.costs[r]:
                A[int(row[x]), int(x)] += sp.Rational(1, 60)
    pi, one = sp.ones(66, 1)/66, sp.ones(1, 66)
    mean = (one*A*pi)[0]
    correction = (sp.eye(66)-T+pi*one).inv()*(A*pi-mean*pi)
    variance = sp.simplify(mean*(1-mean)+2*(one*A*correction)[0])
    assert mean == sp.Rational(5, 33)
    assert variance == sp.Rational(1597, 11979)
    assert variance != mean*(1-mean)


@pytest.mark.parametrize('z', [1, .5, .25])
def test_doob_complete_channel_and_all_filter_failures(native, z):
    driven = n.doob_instrument(native, z)
    channel = driven['instrument']
    np.testing.assert_allclose(n.effects_diagonal(channel).sum(axis=0), 1, atol=3e-15, rtol=0)
    stationary = driven['stationary_state']
    np.testing.assert_allclose(n.apply(channel, stationary), stationary, atol=3e-15, rtol=0)
    rng = np.random.default_rng(912)
    psi = rng.normal(size=66)+1j*rng.normal(size=66)
    psi /= np.linalg.norm(psi)
    state = np.outer(psi, psi.conj())
    filtered = n.filter_implementation(native, z, state)
    assert len(filtered['outcomes']) == 1+3*len(native.labels)
    assert filtered['total_probability'] == pytest.approx(1, abs=3e-15)
    assert filtered['success_probability'] == pytest.approx(driven['herald_success_probability'], abs=3e-15)
    np.testing.assert_allclose(filtered['successful_output'],
        driven['herald_success_probability']*n.apply(channel, state), atol=3e-15, rtol=0)
    if z < 1:
        assert filtered['success_probability'] < 1
        unconditional = sum(output for _, output in filtered['outcomes'])
        assert np.linalg.norm(unconditional-n.apply(channel, state)) > 1e-3


def test_equal_event_count_apis_reject_a_changed_law(native):
    weighted = n.doob_instrument(native, .5)['instrument']
    for call in (lambda: n.quotient_counts(weighted), lambda: n.native_blocks(weighted),
                 lambda: n.two_record_reference_law(weighted, ROOT),
                 lambda: n.doob_instrument(weighted, .5)):
        with pytest.raises(ValueError, match='unchanged native'):
            call()
    n.quotient_counts(n.doob_instrument(native, 1)['instrument'])


def test_array_replay_detects_count_and_probability_changes():
    original = {'counts': np.array([0, 1, 66]), 'values': np.array([.25, .5])}
    compare_arrays(original, {'counts': original['counts'].copy(),
                              'values': original['values']+1e-16})
    for changed in ({'counts': np.array([0, 2, 66]), 'values': original['values']},
                    {'counts': original['counts'], 'values': original['values']+1e-5}):
        with pytest.raises(ValueError):
            compare_arrays(original, changed)


def test_exact_chart_matches_executed_source_coordinates():
    from model.common_source.slow_geometry import raw_frame
    actual = n.geometry_packet(ROOT)['frame']
    expected = sp.Matrix([[sp.Rational(a)+sp.Rational(b)*sp.sqrt(5) for a, b in row]
                          for row in raw_frame()])
    assert n.is_zero(actual-expected)


def test_cached_chart_still_rejects_changed_source_bytes(tmp_path):
    from pathlib import Path
    import shutil
    packet = n.geometry_packet(ROOT)
    for name in packet['source_hashes']:
        target = tmp_path/name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT/name, target)
    n.geometry_packet(tmp_path)
    chart = tmp_path/'model/common_source/slow_geometry.py'
    chart.write_text(chart.read_text()+'\n# changed live source\n')
    with pytest.raises(ValueError, match='source bytes differ'):
        n.geometry_packet(tmp_path)
