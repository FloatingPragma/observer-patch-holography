"""Synthetic source tests for delayed challenges, missed waves and tag witnesses."""

import json
import math

import numpy as np
import pytest
from scipy.linalg import expm

from . import radar, tagged_radar as tagged


def test_perfect_swap_echo_decodes_only_arrival_slots_and_silent_capture_has_no_reply():
    unitary = np.array([[0, 1], [1, 0]], complex)
    histories = tagged.pair_histories(unitary, 0, 1, 2)
    laws = tagged.decode_histories(histories, 2)
    expected = np.zeros((2, 2, 3, 3))
    expected[tagged.ACTIVE, 0, 0, 0] = expected[tagged.ACTIVE, 1, 1, 1] = .5
    expected[tagged.SILENT, 0, 0, 2] = expected[tagged.SILENT, 1, 1, 2] = .5
    np.testing.assert_allclose(laws, expected)
    assert histories[0, 0, 0, 0] == .5 and histories[0, 1, 2, 2] == .5
    active, silent = [tagged.information(x) for x in laws]
    assert active['source_bit_to_B'] == active['source_bit_to_A'] == 1
    assert silent['source_bit_to_A'] == 0
    # A known, deterministically echoed bit contains no fresh B tag beyond b.
    assert active['receiver_tag_to_A_given_bit_and_capture'] == 0


def test_missed_wave_survives_into_return_window_in_both_controls():
    angle = math.pi / 6
    unitary = expm(-1j * angle * np.array([[0, 1], [1, 0]]))
    histories = tagged.pair_histories(unitary, 0, 1, 1)
    c2, s2 = math.cos(angle)**2, math.sin(angle)**2
    expected = .5 * np.array([c2**3, c2**2 * s2**2, c2**3 * s2])
    np.testing.assert_allclose(histories[0, 0, -1], expected, atol=1e-15)
    np.testing.assert_array_equal(histories[0, :, -1], histories[1, :, -1])
    laws = tagged.decode_histories(histories, 1)
    assert laws[0, 0, 1, 1] > 0  # A late outbound click sets the late reply slot.
    assert laws[0, 0, 0, 1] > 0  # An early tag can arrive in the late A slot.
    assert laws[0, :, 1, 0].sum() == 0  # No late release arrives before release.
    assert tagged.event_values(laws[1])['false_echo'] > 0
    assert tagged.event_values(laws[1])['return_tag_correct'] == 0
    assert tagged.information(laws[1])['receiver_tag_to_A_given_bit_and_capture'] == 0
    assert tagged.information(laws[1])['source_bit_to_A'] > 0
    assert tagged.information(laws[0])['receiver_tag_to_A_given_bit_and_capture'] > 0


def test_captured_branch_matches_independent_first_hit_factorization():
    h = np.array([[.1, .7, .3], [.7, -.2, .4], [.3, .4, .5]])
    unitary = expm(-.4j * h)
    ticks = 3
    direct = tagged.pair_histories(unitary, 0, 2, ticks)
    long_legs = radar.first_hit(unitary, (0, 1, 2), 2 * ticks)
    short_legs = radar.first_hit(unitary, (0, 1, 2), ticks)
    independent = tagged.factorized_capture_histories(long_legs, short_legs, 0, 2, ticks)
    np.testing.assert_allclose(direct[0, :, :-1], independent[:, :-1], atol=1e-14)
    np.testing.assert_allclose(direct.sum(axis=(-3, -2, -1)), [1, 1], atol=1e-14)


def test_identity_source_has_false_echoes_without_any_receiver_capture():
    histories = tagged.pair_histories(np.eye(2), 0, 1, 2)
    laws = tagged.decode_histories(histories, 2)
    assert np.all(laws[:, :, 2, 0] == .5)
    np.testing.assert_array_equal(laws[0], laws[1])
    events = tagged.event_values(laws[0])
    assert events['false_echo'] == events['outbound_timeout'] == events['any_A_click'] == 1
    assert events['delivered_correct'] == events['receiver_capture'] == 0
    assert tagged.information(laws[0])['receiver_tag_to_A_given_bit_and_capture'] is None
    json.dumps(tagged.information(laws[0]), allow_nan=False)


def test_retained_event_partition_and_integer_record_replay():
    unitary = expm(-.6j * np.array([[0, 1], [1, 0]]))
    histories = tagged.pair_histories(unitary, 0, 1, 2)
    weights = radar.acquisition_weights(histories[0])
    counts = tagged._draw(weights, 40003, 12, 0, 0, 0, 1)
    np.testing.assert_array_equal(counts, tagged._draw(weights, 40003, 12, 0, 0, 0, 1))
    decoded = tagged.decode_histories(counts, 2).astype(np.int64)
    assert decoded.sum() == 40003 and counts.dtype == np.int64
    event = tagged.event_values(decoded)
    assert (event['delivered_correct'] + event['delivered_wrong'] + event['return_timeout_after_capture']
            + event['false_echo'] + event['both_timeouts']) == 40003
    assert event['receiver_capture'] + event['outbound_timeout'] == 40003
    assert np.all(counts[1, :2] == 0)


def test_conditional_information_bound_excludes_source_knowledge_and_false_branches():
    counts = np.zeros((2, 3, 3), dtype=np.int64)
    counts[0, 0, 0] = counts[0, 1, 1] = 250000
    counts[1, 2, 1] = 500000  # False echoes carry a source bit but are outside the witness.
    weights = radar.acquisition_weights(counts / counts.sum())
    bound = tagged.information_lower_bound(counts, weights, 18)
    assert 0 < bound['mutual_information_lower_bits'] <= 1
    assert 0 < bound['capture_weighted_information_lower_bits_per_trial'] <= .5
    assert bound['information_given_bit_and_capture_lower_bits'] > 0
    deterministic = np.zeros_like(counts)
    deterministic[0, 0, 0] = deterministic[1, 1, 1] = 500000
    weights = radar.acquisition_weights(deterministic / deterministic.sum())
    assert tagged.information_lower_bound(deterministic, weights, 18)['mutual_information_lower_bits'] == 0
    silent = counts.copy()
    silent[:, :2] = 0
    silent[0, 0, 2] = silent[0, 1, 2] = 250000
    weights = radar.acquisition_weights(silent / silent.sum())
    assert tagged.information_lower_bound(silent, weights, 18)['mutual_information_lower_bits'] == 0


def test_resource_budget_counts_vacuum_tests_and_source_on_intervals():
    histories = tagged.pair_histories(np.array([[0, 1], [1, 0]], complex), 0, 1, 2)
    counts = np.rint(histories * 100).astype(np.int64)
    resources = tagged._resources(histories, counts, 2, .5)
    assert resources[0]['released_echo_count'] == resources[0]['captured_scalar_count'] == 100
    assert resources[1]['released_echo_count'] == 0
    assert resources[0]['acquired_B_occupation_tests'] == 200
    assert resources[0]['acquired_A_occupation_tests'] == 200
    assert resources[1]['acquired_A_occupation_tests'] == 400
    assert resources[0]['acquired_propagation_action_until_A_record'] == 100
    assert resources[1]['acquired_propagation_action_until_A_record'] == 50
    assert resources[0]['elapsed_source_action_time'] == resources[1]['elapsed_source_action_time'] == 400


def test_synthetic_wrapper_retains_all_false_returns_and_no_json_nan(monkeypatch):
    monkeypatch.setattr(radar, '_frozen_unitary', lambda protocol, mark: np.eye(30, dtype=complex))

    def deterministic_counts(weights, shots, *unused):
        counts = (weights * shots // radar.LATTICE).astype(np.int64)
        counts.ravel()[int(np.argmax(weights))] += shots - int(counts.sum())
        return counts

    monkeypatch.setattr(tagged, '_draw', deterministic_counts)
    protocol = dict(endpoints=list(radar.ENDPOINTS), mark=2, dt=.5, ticks=12,
                    shots=tagged.SHOTS, seed=73)
    summary, arrays = tagged.run(protocol)
    json.dumps(summary, allow_nan=False)
    assert len(summary['dyads']) == 60
    assert arrays['tagged_time_counts'].dtype == np.int64
    assert arrays['tagged_time_counts'].sum() == 120 * tagged.SHOTS
    assert arrays['tagged_decoder_probabilities'].shape == (2, 2, 6, 6, 2, 3, 3)
    for dyad in summary['dyads']:
        assert dyad['active_events']['false_echo']['probability'] == 1
        assert dyad['active_events']['delivered_correct']['probability'] == 0
        assert dyad['active_minus_silent']['any_A_click']['probability_difference'] == 0
        assert dyad['active_information']['receiver_tag_to_A_given_bit_and_capture'] is None
        assert dyad['active_tag_information_lower_bound']['mutual_information_lower_bits'] == 0
