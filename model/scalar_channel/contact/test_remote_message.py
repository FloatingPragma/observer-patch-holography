"""Synthetic checks of restricted statistical decoding and large aggregate records."""

import inspect
import json

import numpy as np
import pytest

from . import radar, remote_message as message


def test_reader_abstains_under_null_and_too_small_contrast():
    trials = message.TRIALS
    payload = [trials // 3, trials // 2] * 4
    null = message.read_message([trials // 2, trials // 2], payload, trials)
    assert not null['confidence_separated']
    assert null['decoded_bits'] == [message.ABSTAIN] * 8
    small = message.read_message([trials // 2, trials // 2 + 10], payload, trials)
    assert not small['confidence_separated'] and small['threshold_click_numerator'] is None


def test_reader_uses_known_labels_with_either_signal_sign():
    trials = message.TRIALS
    secret = [0, 1, 1, 0, 1, 0, 0, 1]
    for calibration in ([trials // 4, 3 * trials // 4], [3 * trials // 4, trials // 4]):
        received = [calibration[bit] for bit in secret]
        result = message.read_message(calibration, received, trials)
        assert result['confidence_separated'] and result['decoded_bits'] == secret
        assert result['bit_one_has_higher_click_rate'] == (calibration[1] > calibration[0])
        assert result['threshold_click_numerator'] == trials
        assert result['threshold_click_denominator'] == 2


def test_reader_has_no_capture_payload_or_source_bit_argument():
    assert list(inspect.signature(message.read_message).parameters) == [
        'calibration_A_clicks', 'payload_A_clicks', 'trials']
    first = np.zeros((2, 25, 25), dtype=np.int64)
    second = np.zeros_like(first)
    first[0, 0, 0], first[1, 24, 24] = 123, 877
    second[1, 24, 12], second[0, 0, 24] = 123, 877
    # Change every hidden launch/capture/time label while preserving A's click total.
    assert message.A_click_total(first) == message.A_click_total(second) == 123
    result = message.read_message([0, 1000], [message.A_click_total(first)] * 8, 1000)
    changed = message.read_message([0, 1000], [message.A_click_total(second)] * 8, 1000)
    assert result == changed


def test_large_multinomial_counts_are_int64_exact_and_reproducible():
    weights = np.zeros((2, 25, 25), dtype=np.int64)
    weights[0, 0, 0] = radar.LATTICE // 4
    weights[1, 24, 24] = 3 * radar.LATTICE // 4
    first = message.multinomial_block(weights, message.TRIALS, 19, 1, 0, 0, 1, 3)
    second = message.multinomial_block(weights, message.TRIALS, 19, 1, 0, 0, 1, 3)
    np.testing.assert_array_equal(first, second)
    assert first.dtype == np.int64 and first.sum() == message.TRIALS
    assert first[0, 0, 0] > 2**32 and first[1, 24, 24] > 2**32
    assert np.count_nonzero(first) == 2
    changed = message.multinomial_block(weights, message.TRIALS, 19, 0, 0, 0, 1, 3)
    assert not np.array_equal(first, changed)


def test_hoeffding_allocation_and_sufficient_resource_bound():
    budget = message.allowances(message.TRIALS)
    assert budget['payload_sampling'] > budget['calibration_sampling'] > 0
    gap = 4 * budget['calibration'] + 2 * budget['payload']
    assert message.sufficient_block_trials(1.01 * gap) < message.TRIALS
    assert message.sufficient_block_trials(.99 * gap) > message.TRIALS
    assert message.sufficient_block_trials(0) is None
    assert message.sufficient_block_trials(6 * message.LATTICE_TV_BOUND) is None


def test_private_sign_stream_is_independent_of_calibration_and_payload_counts():
    secret = message._private_payload(73, 1, 0, 5)
    weights = np.array([radar.LATTICE // 2, radar.LATTICE // 2], dtype=np.int64)
    for block in range(10):
        message.multinomial_block(weights, message.TRIALS, 73, block % 2, 1, 0, 5, block)
    np.testing.assert_array_equal(secret, message._private_payload(73, 1, 0, 5))
    assert secret.dtype == np.int64 and np.all((secret == 0) | (secret == 1))


def test_null_wrapper_retains_all_timeouts_and_abstentions(monkeypatch):
    monkeypatch.setattr(message, 'sign_measurement_certificate', lambda: dict(synthetic=True))
    laws = np.zeros((2, 2, 6, 6, 2, 25, 25))
    for source in range(6):
        for receiver in range(6):
            if source != receiver:
                laws[:, :, source, receiver, :, 24, 24] = .5
    spec = dict(endpoints=list(radar.ENDPOINTS), mark=2, dt=.5, ticks=12, shots=1048576, seed=73)
    summary, arrays = message.run(spec, laws)
    json.dumps(summary, allow_nan=False)
    assert summary['total_correct_bits'] == summary['total_incorrect_bits'] == 0
    assert summary['total_abstentions'] == 480
    assert not summary['every_dyad_calibration_passed']
    assert not summary['all_payload_bits_correct']
    assert arrays['remote_message_time_counts'].dtype == np.int64
    assert arrays['remote_message_time_counts'].sum() == 600 * message.TRIALS
    assert arrays['remote_message_A_click_counts'].sum() == 0
    assert np.all(arrays['remote_message_decoded_bits'] == message.ABSTAIN)
    replay_summary, replay_arrays = message.run(spec, laws, arrays['remote_message_time_counts'])
    np.testing.assert_array_equal(replay_arrays['remote_message_private_bits'], arrays['remote_message_private_bits'])
    np.testing.assert_array_equal(replay_arrays['remote_message_decoded_bits'], arrays['remote_message_decoded_bits'])
    assert replay_summary['record_validation']['mode'] == 'retained histogram source compatibility'
    corrupted = arrays['remote_message_time_counts'].copy()
    corrupted[0, 0, 1, 0, 0, 24, 24] -= 1
    corrupted[0, 0, 1, 0, 0, 24, 0] += 1
    with pytest.raises(ValueError, match='forbidden acquisition bin'):
        message.run(spec, laws, corrupted)


def test_saved_histogram_shape_bounds_totals_and_category_tampering():
    weights = np.zeros((2, 25, 25), dtype=np.int64)
    weights[:, 24, 24] = radar.LATTICE // 2
    counts = np.zeros_like(weights)
    counts[:, 24, 24] = message.TRIALS // 2
    assert message.validate_block_record(counts, weights) == 0
    changed = counts.copy()
    changed[0, 24, 24] += 10**9
    changed[1, 24, 24] -= 10**9
    with pytest.raises(ValueError, match='Bernstein'):
        message.validate_block_record(changed, weights)
    with pytest.raises(ValueError, match='int64 histogram shape'):
        message._validate_saved_records(np.zeros((2, 6, 6, 10, 2, 25, 25), dtype=np.int32))
    wrong = np.zeros((2, 6, 6, 10, 2, 25, 25), dtype=np.int64)
    wrong[0, 0, 0, 0, 0, 0, 0] = 1
    with pytest.raises(ValueError, match='self-dyad'):
        message._validate_saved_records(wrong)
    wrong[0, 0, 0, 0, 0, 0, 0] = 0
    with pytest.raises(ValueError, match='wrong trial total'):
        message._validate_saved_records(wrong)


def test_full_synthetic_message_recovers_private_byte_from_A_counts_alone(monkeypatch):
    monkeypatch.setattr(message, 'sign_measurement_certificate', lambda: dict(synthetic=True))
    laws = np.zeros((2, 2, 6, 6, 2, 25, 25))
    for source in range(6):
        for receiver in range(6):
            if source != receiver:
                for bit in (0, 1):
                    laws[:, 0, source, receiver, bit, 12 * bit, 12 * bit] = .5
                    laws[:, 1, source, receiver, bit, 12 * bit, 24] = .5
    spec = dict(endpoints=list(radar.ENDPOINTS), mark=2, dt=.5, ticks=12, shots=1048576, seed=73)
    summary, arrays = message.run(spec, laws)
    assert summary['all_payload_bits_correct'] and summary['total_correct_bits'] == 480
    assert summary['total_incorrect_bits'] == summary['total_abstentions'] == 0
    replay, retained = message.run(spec, laws, arrays['remote_message_time_counts'])
    assert replay['all_payload_bits_correct']
    np.testing.assert_array_equal(retained['remote_message_A_click_counts'], arrays['remote_message_A_click_counts'])


def _boundary_laws():
    laws = np.zeros((2, 2, 6, 6, 2, 25, 25))
    for source in range(6):
        for receiver in range(6):
            if source != receiver:
                laws[:, :, source, receiver, :, 24, 0] = .25 + .49 / radar.LATTICE
                laws[:, :, source, receiver, :, 24, 1] = .125
                laws[:, :, source, receiver, :, 24, 24] = .125 - .49 / radar.LATTICE
    return laws


def test_frozen_original_weights_survive_tie_and_half_lattice_boundary(monkeypatch):
    monkeypatch.setattr(message, 'sign_measurement_certificate', lambda: dict(synthetic=True))
    original = _boundary_laws()
    changed = original.copy()
    for source in range(6):
        for receiver in range(6):
            if source != receiver:
                changed[:, :, source, receiver, 0, 24, 0] -= .02 / radar.LATTICE
                changed[:, :, source, receiver, 0, 24, 24] += .02 / radar.LATTICE
                changed[:, :, source, receiver, 1, 24, 0] += .02 / radar.LATTICE
                changed[:, :, source, receiver, 1, 24, 24] -= .02 / radar.LATTICE
    assert not np.array_equal(radar.acquisition_weights(original[0, 0, 0, 1]),
                              radar.acquisition_weights(changed[0, 0, 0, 1]))
    spec = dict(endpoints=list(radar.ENDPOINTS), mark=2, dt=.5, ticks=12, shots=1048576, seed=73)
    _, acquired = message.run(spec, original)
    frozen = acquired['remote_message_acquisition_weights']
    summary, replayed = message.run(spec, changed, frozen_weights=frozen)
    np.testing.assert_array_equal(replayed['remote_message_acquisition_weights'], frozen)
    np.testing.assert_array_equal(replayed['remote_message_time_counts'], acquired['remote_message_time_counts'])
    assert summary['acquisition_weight_validation']['mode'] == 'verified frozen acquisition masses'
    _, scored = message.run(spec, changed, acquired['remote_message_time_counts'], frozen_weights=frozen)
    np.testing.assert_array_equal(scored['remote_message_A_click_counts'], acquired['remote_message_A_click_counts'])


def test_frozen_weights_reject_wrong_shape_totals_zeros_and_bin_corruption():
    laws = _boundary_laws()
    weights = np.zeros(laws.shape, dtype=np.int64)
    for mark in (0, 1):
        for control in (0, 1):
            for source in range(6):
                for receiver in range(6):
                    if source != receiver:
                        weights[mark, control, source, receiver] = radar.acquisition_weights(
                            laws[mark, control, source, receiver])
    message._validate_frozen_weights(laws, weights)
    with pytest.raises(ValueError, match='complete integer shape'):
        message._validate_frozen_weights(laws, weights.astype(float))
    changed = weights.copy()
    changed[0, 0, 0, 0, 0, 24, 0] = 1
    with pytest.raises(ValueError, match='self-dyad'):
        message._validate_frozen_weights(laws, changed)
    changed = weights.copy()
    changed[0, 0, 0, 1, 0, 24, 0] += 1
    with pytest.raises(ValueError, match='wrong menu total'):
        message._validate_frozen_weights(laws, changed)
    changed[0, 0, 0, 1, 0, 24, 0] -= 2
    changed[0, 0, 0, 1, 0, 24, 2] += 1
    with pytest.raises(ValueError, match='exact physical zero'):
        message._validate_frozen_weights(laws, changed)
    for delta, error in ((200, 'bin bound'), (10000, 'total-variation')):
        changed = weights.copy()
        changed[0, 0, 0, 1, 0, 24, 1] += delta
        changed[0, 0, 0, 1, 0, 24, 24] -= delta
        with pytest.raises(ValueError, match=error):
            message._validate_frozen_weights(laws, changed)
    altered_source = laws.copy()
    altered_source[0, 0, 0, 1, 0, 24, 0] -= .01
    altered_source[0, 0, 0, 1, 0, 24, 1] += .01
    with pytest.raises(ValueError, match='total-variation'):
        message._validate_frozen_weights(altered_source, weights)
