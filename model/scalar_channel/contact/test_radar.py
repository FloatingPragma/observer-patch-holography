"""Synthetic checks of projective histories, independent renewal and withheld fits."""

import json

import numpy as np
import pytest
from scipy.linalg import expm

from . import radar


def test_first_hit_is_a_measured_history_and_renewal_predicts_it():
    h = np.array([[.3, .8, 0], [.8, -.2, .6], [0, .6, .4]])
    u = expm(-.37j * h)
    actual = radar.first_hit(u, (0, 1, 2), 7)
    predicted = radar.renewal_prediction(u, (0, 1, 2), 7)
    np.testing.assert_allclose(actual, predicted, atol=2e-14)
    state = np.array([1., 0, 0], complex)
    hits = []
    for _ in range(7):
        state = u @ state
        hit = np.diag([0., 0, 1.]) @ state
        hits.append(float(np.vdot(hit, hit).real))
        state = np.diag([1., 1, 0.]) @ state
    np.testing.assert_allclose(actual[0, 2], hits + [float(np.vdot(state, state).real)], atol=2e-14)
    assert abs(actual[0, 2, 1] - abs((u @ u)[2, 0])**2) > 1e-3
    np.testing.assert_allclose(actual.sum(axis=-1), 1 - np.eye(3), atol=2e-14)


def test_outbound_timeout_has_no_return_and_return_timeout_is_retained():
    legs = np.zeros((2, 2, 3))
    legs[0, 1] = [.2, .3, .5]
    legs[1, 0] = [.4, .1, .5]
    joint = radar.round_trip_law(legs)
    np.testing.assert_allclose(joint[0, 1], [[.08, .02, .1], [.12, .03, .15], [0, 0, .5]])
    np.testing.assert_allclose(joint[1, 0], [[.08, .12, .2], [.02, .03, .05], [0, 0, .5]])
    distances, successes = radar.distance_statistics(joint, .5)
    assert distances[0, 1] == pytest.approx(.7)
    assert successes[0, 1] == pytest.approx(.5)
    counts = np.rint(joint * 100).astype(np.int64)
    pairs = radar._pair_summaries(joint, counts, .5)
    assert pairs[0]['acquired_projector_tests'] == 260
    assert pairs[0]['outbound_timeout_count'] == 50
    assert pairs[0]['return_timeout_count'] == 25
    assert pairs[0]['successful_round_trip_count'] == 25


def test_integer_acquisition_is_complete_and_replays_from_retained_weights():
    p = np.array([[.123456789012345, 0], [.333333333333333, .543209877654322]])
    weights = radar.acquisition_weights(p)
    assert weights.dtype == np.int64
    assert weights.sum() == radar.LATTICE
    assert weights[0, 1] == 0
    assert .5 * abs(weights / radar.LATTICE - p).sum() <= p.size / (2 * radar.LATTICE)
    first = radar.replay_counts(weights, 10003, 123, 1, 0, 5)
    second = radar.replay_counts(weights.copy(), 10003, 123, 1, 0, 5)
    np.testing.assert_array_equal(first, second)
    assert first.dtype == np.int64 and first.sum() == 10003 and first[0, 1] == 0
    with pytest.raises(ValueError, match='retained acquisition'):
        radar.replay_counts(np.zeros((2, 2), dtype=np.int64), 10, 1, 0, 0, 1)


def _synthetic_distances():
    points = np.array([[.2, .3, .7], [1.2, -.4, .5], [0, 0, 0], [2, 0, 0],
                       [0, 2, 0], [0, 0, 2]])
    return np.linalg.norm(points[:, None] - points[None, :], axis=-1)


def test_trilateration_withholds_pair_before_constructing_prediction():
    distances = _synthetic_distances()
    result = radar.predict_withheld_pair(distances, 0, 1)
    assert result['domain_valid']
    assert result['predicted_distance'] == pytest.approx(distances[0, 1])
    changed = distances.copy()
    changed[0, 1] = changed[1, 0] = 100
    assert radar.predict_withheld_pair(changed, 0, 1) == result
    evaluated = radar.geometry_test(changed)['pairs'][0]
    assert not evaluated['accepted'] and evaluated['relative_error'] > .9
    assert radar.geometry_test(distances)['all_pairs_accepted']


def test_geometry_rejects_singular_non_euclidean_and_radially_inconsistent_data():
    flat = np.arange(6.)
    collinear = abs(flat[:, None] - flat[None, :])
    singular = radar.predict_withheld_pair(collinear, 0, 1)
    assert singular['reason'] == 'singular_or_ill_conditioned_anchor_frame'
    distances = _synthetic_distances()
    distances[3, 4] = distances[4, 3] = 10
    assert radar.predict_withheld_pair(distances, 0, 1)['reason'] == 'non_euclidean_anchor_distances'
    distances = np.ones((6, 6)) - np.eye(6)
    assert radar.predict_withheld_pair(distances, 0, 1)['reason'] == 'endpoint_distances_incompatible_with_anchor_frame'


def test_zero_success_and_json_have_explicit_undefined_ranges():
    joint = np.zeros((6, 6, 3, 3))
    for i in range(6):
        for j in range(6):
            if i != j:
                joint[i, j, -1, -1] = 1
    distances, successes = radar.distance_statistics(joint, .5)
    assert successes.sum() == 0
    matrix = radar._json_matrix(distances)
    assert matrix[0][1] is None and matrix[0][0] == 0
    result = radar.geometry_test(distances)
    assert result['valid_domains'] == 0 and not result['all_pairs_accepted']
    json.dumps(dict(distances=matrix, geometry=result), allow_nan=False)


def test_frozen_menu_rejects_substituted_clock_and_endpoint():
    protocol = dict(endpoints=list(radar.ENDPOINTS), mark=2, dt=.5, ticks=12, shots=1000, seed=73)
    radar._settings(protocol)
    for changes in (dict(dt=.25), dict(ticks=24), dict(mark=3), dict(endpoints=[1, 5, 11, 16, 22, 29]),
                    dict(shots=True), dict(seed=-1)):
        with pytest.raises(ValueError):
            radar._settings(dict(protocol, **changes))


def test_wrapper_retains_complete_timeout_records_on_a_synthetic_source(monkeypatch):
    monkeypatch.setattr(radar, '_frozen_unitary', lambda protocol, coupling: np.eye(30, dtype=complex))
    protocol = dict(endpoints=list(radar.ENDPOINTS), mark=2, dt=.5, ticks=12, shots=7, seed=73)
    summary, arrays = radar.run(protocol)
    json.dumps(summary, allow_nan=False)
    assert arrays['radar_joint_counts'].dtype == np.int64
    assert arrays['radar_joint_counts'].shape == (2, 6, 6, 13, 13)
    np.testing.assert_array_equal(arrays['radar_joint_counts'].sum(axis=(-2, -1)),
                                  np.broadcast_to(7 * (1 - np.eye(6)), (2, 6, 6)))
    for setting in summary['settings']:
        assert setting['source_emissions'] == 210
        assert setting['return_turnarounds'] == 0
        assert setting['acquired_projector_tests'] == 2520
        assert setting['acquired_source_action_time'] == 1260
        assert setting['geometry_from_born_laws']['valid_domains'] == 0
        assert all(row['conditional_round_trip_mean'] is None for row in setting['pairs'])
        assert setting['distance_confidence']['lower'][0][1] == .5
        assert setting['distance_confidence']['upper'][0][1] == 6
        assert not setting['euclidean_error_from_born_laws']['defined']


def test_gram_certificate_excludes_low_dimensions_and_includes_negative_eigenvalues():
    simplex = np.ones((6, 6)) - np.eye(6)
    certificate = radar.euclidean_error_bound(simplex)
    expected = np.sqrt(1.2) - 1
    for entry in certificate['dimensions'][:4]:
        assert entry['minimum_possible_maximum_relative_error'] == pytest.approx(expected, abs=1e-9)
        assert entry['excludes_five_percent_fit']
    assert certificate['dimensions'][4]['minimum_possible_maximum_relative_error'] == 0
    invalid = simplex.copy()
    invalid[0, 1] = invalid[1, 0] = 4
    all_dimensions = radar.euclidean_error_bound(invalid)['dimensions']
    assert all(x['negative_eigenvalue_obstruction'] > 0 for x in all_dimensions)
    assert all(x['minimum_possible_maximum_relative_error'] > 0 for x in all_dimensions)
    euclidean = radar.euclidean_error_bound(_synthetic_distances())['dimensions']
    assert all(x['minimum_possible_maximum_relative_error'] == 0 for x in euclidean[2:])


def test_simultaneous_confidence_propagates_to_gram_and_common_scale():
    distances = np.ones((6, 6)) - np.eye(6)
    counts = np.full((6, 6), 10**9, dtype=np.int64)
    summary, lower, upper = radar.confidence_envelopes(distances, counts, .5, 12)
    radius = 5.5 * np.sqrt(np.log(60 / radar.CONFIDENCE_ALPHA) / (2 * 10**9))
    assert summary['statistical_half_width'][0][1] == pytest.approx(radius)
    assert lower[0, 1] == pytest.approx(1 - radius)
    exact = radar.euclidean_error_bound(distances)
    uncertain = radar.euclidean_error_bound(distances, lower, upper)
    assert uncertain['distance_uncertainty_operator_bound'] > 0
    assert (0 < uncertain['dimensions'][2]['minimum_possible_maximum_relative_error']
            < exact['dimensions'][2]['minimum_possible_maximum_relative_error'])
    assert uncertain['dimensions'][2]['excludes_five_percent_fit']
    uniform = radar.common_multiplier_intersection(lower, upper, 2 * lower, 2 * upper)
    assert not uniform['common_multiplier_excluded']
    assert uniform['interval'][0] <= 2 <= uniform['interval'][1]
    marked = distances.copy()
    marked[0, 1] = marked[1, 0] = 1.5
    _, marked_lower, marked_upper = radar.confidence_envelopes(marked, counts, .5, 12)
    nonuniform = radar.common_multiplier_intersection(lower, upper, marked_lower, marked_upper)
    assert nonuniform['common_multiplier_excluded'] and nonuniform['interval'] is None
    # Broad uncertainty removes a spectral exclusion; it never supplies an embedding.
    wide_lower, wide_upper = .5 * distances, 6 * distances
    uninformative = radar.euclidean_error_bound(distances, wide_lower, wide_upper)
    assert all(row['minimum_possible_maximum_relative_error'] == 0 for row in uninformative['dimensions'])


def test_lattice_bias_bounds_conditional_rounding_without_source_readings():
    legs = np.zeros((6, 6, 3))
    for left in range(6):
        for right in range(6):
            if left != right:
                legs[left, right] = [.2, .3, .5]
    joint = radar.round_trip_law(legs)
    weights = np.zeros(joint.shape, dtype=np.int64)
    for left in range(6):
        for right in range(6):
            if left != right:
                weights[left, right] = radar.acquisition_weights(joint[left, right])
    counts = np.rint(joint * 100).astype(np.int64)
    bound = radar._lattice_conditional_bias(weights, counts, .5, 2)
    exact, _ = radar.distance_statistics(joint, .5)
    rounded, _ = radar.distance_statistics(weights / radar.LATTICE, .5)
    assert np.all(abs(exact - rounded) <= bound)
    assert 0 < bound[0, 1] < 1e-7
