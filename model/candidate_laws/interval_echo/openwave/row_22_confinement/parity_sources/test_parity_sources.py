import itertools
import json

import networkx as nx
import numpy as np
import pytest
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from model.executable_model import carrier
from . import measure as m
from . import released


def payload(name):
    return json.loads((m.HERE / name).read_text())


def test_native_physical_preparations_and_observable_dual():
    ex, co = candidate.native_law.local_term_parts()
    h = ex + candidate.KAPPA*co
    v = np.zeros(99)
    v[2] = 1
    out = expm(-1j*np.pi*h) @ v
    assert abs(out[6]) == pytest.approx(1)
    assert abs(abs(np.vdot(v, out))) < 1e-12
    for level in range(3):
        board = carrier.federation(level)
        rings, dual, gates, parked, partner = m.survey(board)
        assert all(len(r) in (5, 6) for r in rings)
        assert sum(len(r) == 5 for r in rings) == 12
        assert len(rings) == len(board.external) - board.cells + 2
        assert len(set(partner.values())) == len(partner)
        # Simultaneously toggle one side of every bridge. All six-load
        # guards still accept and every real bridge reads a mismatch.
        state = m.prepare(gates, gates, parked, partner)
        assert all(sum(12*a+i in state for i in range(12)) == 6 for a in range(board.cells))
        assert set(m.read_ports(state, gates)[1]) == set(gates)
        assert set(m.read_ports(state, gates)[0]) == {i for i, r in enumerate(rings) if len(r) == 5}
        for a, b in dual.edges:
            assert len(set(rings[a]) & set(rings[b])) == 2


def test_all_two_source_minima_without_using_ground_formula():
    # Independent exhaustive syndrome energies from XOR convolution of optional
    # physical edges. Each edge can be selected zero or one time.
    _, dual, gates, _, _ = m.survey(carrier.federation(0))
    best = np.full(1 << len(dual), np.inf)
    best[0] = 0
    state = np.arange(len(best))
    for f, g in gates:
        toggled = state ^ (1 << f) ^ (1 << g)
        best = np.minimum(best, best[toggled] + 1)
    assert np.isfinite(best).sum() == 2048
    for a, b in itertools.combinations(dual, 2):
        mask = (1 << a) + (1 << b)
        d = nx.shortest_path_length(dual, a, b)
        assert best[mask] == d
        assert best[state & mask == mask].min() == min(d, 2)


def test_records_source_sectors_and_exterior_ground_without_inventories():
    raw, unguarded = payload('records.json'), payload('released.json')
    # This independent read uses a single late appointment with branch fixed
    # by checking ALL time settings against each exact energy forecast.
    def compare(arm, prediction):
        angles = np.arctan2(*(2*np.array(arm['x_y'])[::-1] / raw['phase_shots'] - 1))
        residuals = np.angle(np.exp(1j*(angles - prediction * np.array(raw['times']))))
        assert abs(residuals).max() < raw['simultaneous_energy_half_width']
        assert abs(arm['measured_work'] - prediction) < raw['simultaneous_energy_half_width']
    for level in raw['levels']:
        _, dual, gates, _, _ = m.survey(carrier.federation(level['level']))
        other = [r for r in unguarded['rows'] if r['level'] == level['level']]
        for row, free in zip(level['rows'], other, strict=True):
            d = row['distance']
            assert nx.shortest_path_length(dual, *row['endpoints']) == d
            for name, n in [('least', d), ('detour', d+1), ('sources_only', min(d, 2))]:
                a = row['trials'][name]
                bits = 0
                for u, v in a['mismatch_edges']:
                    assert (u, v) in gates
                    bits ^= (1 << u) | (1 << v)
                assert a['odd_faces'] == [i for i in range(len(dual)) if bits & (1 << i)]
                assert set(row['endpoints']) <= set(a['odd_faces'])
                if name != 'sources_only':
                    assert a['odd_faces'] == row['endpoints']
                assert len(a['mismatch_edges']) == n
                compare(a, 32.5*n)
                compare(free['coherent_work'][name], 32*n)
    assert abs(unguarded['filter_accepted']/unguarded['filter_shots']-0.5) < 0.001


def test_release_records_and_full_board_dilute_bounds():
    raw, free = payload('records.json'), payload('released.json')
    for level in raw['levels']:
        result = level['release']
        p, q, b = result['ports_p_q_b']
        board = carrier.federation(level['level'])
        assert q in board.neighbours[p] and p//12 == q//12
        assert b in board.neighbours[p] and p//12 != b//12
        assert set(result['before']) < set(result['after'])
        assert len(set(result['after']) - set(result['before'])) == 2
        for name, c, hopping in [('registered', 32.5, 1), ('exchange_only', 0.5, 1),
                                 ('contrast_only', 32, 0), ('shuttered', 32.5, 0)]:
            trial = result['variants'][name]
            omega = np.hypot(c, hopping)
            exact = (hopping/omega * np.sin(np.array(trial['times'])*omega/2))**2
            np.testing.assert_allclose(exact, trial['extra_pair_probability'], atol=1e-13)
            counts = np.array(trial['extra_pair_counts'])
            shots = result['shots']
            assert np.all(abs(counts-shots*exact) < 1+6*np.sqrt(shots*exact*(1-exact)))
        for row in (r for r in free['rows'] if r['level'] == level['level']):
            i, j = row['dilute_double_ports']
            boundary = len(set(board.neighbours[i])) + len(set(board.neighbours[j])) - 2*(j in board.neighbours[i])
            private = 32 * (12 - board.degree[i] - board.degree[j])
            assert row['dilute_upper'] == 32.5 * boundary + private <= 390
            for ports_key, work_key in [('dilute_even_original_ports', 'dilute_even_original_trial'),
                                        ('dilute_even_adjacent_ports', 'dilute_even_upper')]:
                a, b = row[ports_key]
                internal = 2*int(b in board.neighbours[a])
                assert row[work_key] == (32.5*(board.degree[a]+board.degree[b]-internal)
                                         + 32*(12-board.degree[a]-board.degree[b]))
            assert row['dilute_even_original_trial'] == 389
            assert row['dilute_even_upper'] == 324
    x, z = m.source_block()
    h = x + candidate.KAPPA*z
    plus = (np.eye(9)[2] + np.eye(9)[6])/np.sqrt(2)
    np.testing.assert_allclose(h @ plus, 32*plus, atol=1e-12)
    assert np.linalg.eigvalsh(h).min() > -1e-13


def test_exact_acquisition_replay():
    fresh = json.loads(json.dumps(m.run()))
    frozen = payload('records.json')
    for new_level, old_level in zip(fresh['levels'], frozen['levels'], strict=True):
        for name, old in old_level['release']['variants'].items():
            new = new_level['release']['variants'][name]
            # expm produces binary64 approximations. Counts, contacts, clocks
            # and all other fields must replay exactly, including zero signals.
            actual = np.array(new.pop('extra_pair_probability'))
            expected = np.array(old.pop('extra_pair_probability'))
            np.testing.assert_array_equal(actual == 0, expected == 0)
            np.testing.assert_array_max_ulp(actual, expected, maxulp=16)
    assert fresh == frozen
    assert json.loads(json.dumps(released.run())) == payload('released.json')


def test_half_filled_permutations_on_a_coherent_ground():
    from .half_filled import both_odd, run

    n = 6
    masks = [sum(1 << i for i in sites) for sites in itertools.combinations(range(n), n//2)]
    positions = {m: i for i, m in enumerate(masks)}
    h = np.zeros((len(masks), len(masks)))
    for j, bits in enumerate(masks):
        for a in range(n):
            b = (a+1) % n
            if (bits >> a ^ bits >> b) & 1:
                h[j, j] += 32.5
                h[positions[bits ^ (1 << a) ^ (1 << b)], j] -= 0.5
    energies, frame = np.linalg.eigh(h)
    psi = frame[:, 0]
    assert max(abs(psi)) < 0.99
    a, b = 0b0011, 0b0110
    keep = np.array([bool((bits & a).bit_count() % 2 and (bits & b).bit_count() % 2)
                     for bits in masks])
    ground_pin = np.linalg.eigvalsh(h[np.ix_(keep, keep)])[0]
    probabilities, chosen = [], None
    certified = float(both_odd(n, a.bit_count(), b.bit_count(), (a ^ b).bit_count()))
    for perm in itertools.permutations(range(n)):
        rows = [positions[sum(1 << perm[i] for i in range(n) if bits >> i & 1)] for bits in masks]
        y = np.eye(len(masks))[rows].copy()
        y[~keep] = 0
        v = y @ psi
        p = v @ v
        probabilities.append(p)
        if p >= certified and chosen is None:
            chosen = y
    assert np.mean(probabilities) == pytest.approx(certified)
    assert np.mean(keep) == pytest.approx(certified)
    assert chosen is not None
    v = chosen @ psi
    rayleigh = v @ h @ v / (v @ v)
    bound = np.linalg.norm(h @ chosen - chosen @ h, 2) / np.linalg.norm(v)
    assert ground_pin - 1e-10 <= rayleigh <= energies[0] + bound + 1e-10
    report = run()
    assert report == payload('half_filled.json')
    assert [r['pairs'] for r in report['levels']] == [66, 861, 13041]
    assert all(r['source_work_upper'] < 19129 and r['private_service_width'] == 32
               and r['native_services_per_port'] == 6 for r in report['levels'])
