"""Independent geometry, native-column, pulse and public record audits."""

from functools import lru_cache
import gzip
from itertools import combinations
import json

import networkx as nx
import numpy as np
import pytest
from scipy.linalg import expm

from . import measure as m


@lru_cache(maxsize=1)
def selected():
    return m.source.source()


def archive():
    return json.loads(gzip.decompress((m.HERE / 'results.json.gz').read_bytes()))


def test_all_full_native_memory_columns_give_calibrated_spatial_contacts():
    hx, hc = m.law.local_term_parts()
    for strength in (8, 4, 7.5):
        full = expm(-1j * np.pi * (hx + strength * hc))[:, :9]
        independently = np.zeros((99, 9), complex)
        for j, (a, b) in enumerate(m.law.SRC):
            independently[m.law.SRC.index((b, a)), j] = np.exp(-1j * strength * np.pi * (a - b)**2)
        np.testing.assert_allclose(full, independently, atol=5e-14)
    # A no-contrast binary seam has the same positive swap. The native
    # contrast does not determine the representation of calibrated pi gates.
    swap = np.eye(9)[[m.law.SRC.index((b, a)) for a, b in m.law.SRC]]
    np.testing.assert_allclose(expm(-1j * np.pi * m.law.seam_term(0)), swap, atol=5e-14)


def test_selected_native_source_and_improper_character_from_geometry():
    sec, arms, _, _, _, ideal, finite, summary = selected()
    board, _, _, _, _, _, _, geometry = m.geometry(sec, arms)
    h = np.zeros((220, 220))
    for j, word in enumerate(sec.states):
        occupied = set(word)
        h[j, j] = 8 * sum(6 - board.degree[x] for x in word)
        for a, b in board.seams:
            if a >= 12 and b >= 12:
                continue
            left, right = a in occupied, b in occupied
            h[j, j] += 8 * (left != right)
            if left != right:
                h[j, j] += .5
                if b < 12:
                    changed = tuple(sorted(occupied ^ {a, b}))
                    h[sec.index[changed], j] -= .5
    np.testing.assert_allclose(h @ ideal, summary['energy'] * ideal, atol=2e-12)
    assert summary['gap_lower'] > .016
    points, edges = m.carrier.icosahedron_from_coordinates()
    xyz = points[m.carrier.isomorphism(m.carrier.cell().seams, edges)]
    for item in geometry['arm_stabilizer']:
        p = item['ports']
        transform = np.linalg.lstsq(xyz, xyz[p], rcond=None)[0]
        assert abs(np.linalg.det(transform) - item['spatial_determinant']) < 1e-12
        permuted = np.empty(220)
        for j, word in enumerate(sec.states):
            permuted[sec.index[tuple(sorted(p[x] for x in word))]] = ideal[j]
        np.testing.assert_allclose(permuted, item['spatial_determinant'] * ideal, atol=5e-11)
    np.testing.assert_allclose(np.linalg.norm(finite), 1, atol=1e-13)


def test_component_decoder_against_all_six_vertex_graphs():
    lookup = m.pair_lookup()
    independent = np.zeros_like(lookup)
    for code in range(2**15):
        graph = nx.Graph()
        graph.add_nodes_from(range(6))
        graph.add_edges_from(edge for i, edge in enumerate(m.EDGES6) if code >> i & 1)
        independent[code] = sorted(map(len, nx.connected_components(graph))) == [3, 3]
    np.testing.assert_array_equal(lookup, independent)


def test_real_contact_read_is_convex_through_entire_native_pulse():
    words = list(combinations(range(6), 3))
    swap = np.zeros((20, 20))
    contrast = np.zeros(20)
    for j, word in enumerate(words):
        occupied = set(word)
        change = (0 in occupied) != (1 in occupied)
        destination = tuple(sorted(occupied ^ {0, 1})) if change else word
        swap[words.index(destination), j] = 1
        contrast[j] = change
    h = (np.eye(20) - swap) / 2 + 8 * np.diag(contrast)
    rng = np.random.default_rng(701)
    state = rng.normal(size=20)
    state /= np.linalg.norm(state)
    read = rng.random(20) > .5
    endpoints = [float(np.sum(abs(x[read])**2)) for x in (state, swap @ state)]
    for duration in np.linspace(0, np.pi, 19):
        actual = expm(-1j * duration * h) @ state
        predicted = np.cos(duration / 2)**2 * endpoints[0] + np.sin(duration / 2)**2 * endpoints[1]
        assert abs(float(np.sum(abs(actual[read])**2)) - predicted) < 2e-14


def test_physical_words_frames_phases_and_complete_source_histories():
    data = archive()
    board = m.carrier.federation(0)
    seamset = set(board.seams)
    for panel in data['protocols']:
        occupied_labels = np.arange(board.n)
        for a, b in panel['physical_pi_contacts']:
            assert tuple(sorted((a, b))) in seamset
            occupied_labels[a], occupied_labels[b] = occupied_labels[b], occupied_labels[a]
        # Reconstruct the intended whole-pattern action independently of
        # the producer's map update and destination-array contraction.
        expected = np.arange(board.n)
        if panel['protocol'].startswith('exchange'):
            for x, y in enumerate(data['geometry']['second_frame']):
                expected[x] = 12 * data['geometry']['second_owner'] + y
                expected[12 * data['geometry']['second_owner'] + y] = x
            if panel['protocol'] == 'exchange_then_improper_frame':
                f = np.arange(board.n)
                f[:12] = data['geometry']['improper_reflection_ports']
                expected = f[expected]
        np.testing.assert_array_equal(occupied_labels[expected], np.arange(board.n))
        for phase in panel['phases']:
            bank = phase['complete_source_histories_and_phase']
            assert len(bank['counts']) == 26 and sum(bank['counts']) == m.SHOTS
            assert min(bank['probabilities']) >= 0 and abs(sum(bank['probabilities']) - 1) < 5e-11
        for stage in panel['stages']:
            bank = stage['public_component_read']
            assert len(bank['counts']) == 3 and sum(bank['counts']) == m.SHOTS
            assert abs(stage['finite_two_separate_connected_triples'] - stage['ideal_two_separate_connected_triples']) < data['finite_pair_vector_error_upper']
        if panel['protocol'] == 'exchange_then_improper_frame':
            assert panel['phases'][0]['conditional_minus'] > 1 - 1e-12
        else:
            assert panel['phases'][0]['conditional_plus'] > 1 - 1e-12
        assert abs(panel['phases'][0]['conditional_plus'] - panel['phases'][1]['conditional_plus']) < 1e-14
    first, second = data['protocols'][:2]
    assert data['initial_two_core_probability'] > .998
    assert first['minimum_public_two_core_probability'] < .077
    assert abs(first['phases'][2]['conditional_minus'] - second['phases'][2]['conditional_minus']) > .04


def test_sampler_and_optimized_comparator_reject_bad_inputs():
    rng = np.random.default_rng(0)
    with pytest.raises(ArithmeticError):
        m.draw([-.01, 1.01], rng)
    with pytest.raises(ArithmeticError):
        m.draw([.1, .2], rng)
    with pytest.raises(ArithmeticError):
        m.compare({'phase': 1}, {'phase': -1})


def test_counted_loss_on_both_routes_uses_all_stage_appointments():
    from scipy.stats import binom

    data = archive()
    banks = sum(len(panel['stages']) for panel in data['protocols'])
    assert banks == 554
    for panel, upper in zip(data['protocols'][:2], (.07846, .07808), strict=True):
        hits, other, refused = panel['stages'][36]['public_component_read']['counts']
        assert hits + other + refused == data['shots']
        assert binom.cdf(hits, hits + other, upper) < .01 / (2 * banks)


def test_every_proper_half_turn_changes_the_observed_source_and_polar_control():
    from . import half_turns

    data = json.loads(gzip.decompress(half_turns.RECEIPT.read_bytes()))
    points, edges = m.carrier.icosahedron_from_coordinates()
    xyz = points[m.carrier.isomorphism(m.carrier.cell().seams, edges)]
    half_perms = set()
    for perm in m.carrier.automorphisms(m.carrier.cell().seams):
        matrix = np.linalg.lstsq(xyz, xyz[list(perm)], rcond=None)[0]
        if np.linalg.det(matrix) > 0 and abs(np.trace(matrix) + 1) < 1e-12:
            half_perms.add(tuple(perm))
    assert len(half_perms) == 15
    assert {tuple(row['ports']) for row in data['panels'][1:]} == half_perms
    _, _, _, _, _, _, finite, summary = selected()
    sharp = np.zeros(len(finite))
    sharp[data['words'].index(summary['source_word'])] = 1
    words = np.array(data['words'])
    for row in data['panels']:
        perm = np.array(row['ports'])
        assert np.array_equal(perm[perm], np.arange(12))
        for name, state in (('selected', finite), ('polar', sharp)):
            counts = np.array(row[name]['address_counts'])
            assert counts.size == 224 and counts.sum() == data['starts']
            if name == 'polar':
                moved = sorted(perm[summary['source_word']])
                assert counts[4 + data['words'].index(moved)] == data['starts']
            sector = np.array([data['words'].index(sorted(t)) for t in perm[words]])
            # Separate simulated banks use the same source acceptance.
            for law, phase_bank in row[name]['phase_counts'].items():
                assert len(phase_bank) == 6 and sum(phase_bank) == data['starts']
                if law == 'native_8':
                    dark = np.linalg.norm((state - state[sector]) / 2)**2
                    hits = phase_bank[-1]
                    trials = sum(phase_bank[-2:])
                    assert abs(hits / trials - dark) < .00025
    bounds = half_turns.address_separation(data)
    assert min(bounds['selected']) > .0013 and min(bounds['polar']) > .9999
    data['panels'][1]['selected']['address_counts'] = data['panels'][0]['selected']['address_counts']
    assert half_turns.address_separation(data)['selected'][0] <= 0
