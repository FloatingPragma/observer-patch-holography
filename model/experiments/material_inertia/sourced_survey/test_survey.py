"""Native source, separate meters, exact energy control and 30-clock prediction."""

import copy
import json

import numpy as np
import pytest
from scipy.sparse.linalg import expm_multiply

from . import measure, source


def test_native_write_passive_conservation_and_reversed_current():
    source.local_proof()
    _, route, (u, v), h = source.apparatus_map()
    assert not {u, v} & set(route[:2])
    pair = h[[u, v], :][:, [u, v]].toarray() / 2
    q = np.diag([0, 1])
    current = 1j * (pair @ q - q @ pair)
    for sign in (1, -1):
        wave = np.array([1, sign * 1j]) / np.sqrt(2)
        np.testing.assert_allclose(np.vdot(wave, current @ wave), sign / 2, atol=3e-15)
    assert h.dtype == np.int64 and (h - h.T).count_nonzero() == 0


def test_bounded_native_flight_matches_independent_propagation():
    record, err = source.flights()
    board, route, (u, v), h = source.apparatus_map()
    wave = np.zeros((board.n, 3), complex)
    wave[u] = [1, 1 / np.sqrt(2), 1 / np.sqrt(2)]
    wave[v] = [0, 1j / np.sqrt(2), -1j / np.sqrt(2)]
    for t in (1, 4):
        a, b = expm_multiply(-0.5j * t * h, wave)[list(route[:2])]
        for k, name in enumerate(source.NAMES):
            expected = [abs(a[k]) ** 2 + abs(b[k]) ** 2, abs(a[k] - b[k]) ** 2 / 2]
            np.testing.assert_allclose(
                list(record[f"{name}_{t}"].values()),
                expected,
                atol=err["capture_or_energy_probability_error_upper"],
                rtol=0,
            )


def test_counts_only_held_prediction_and_exact_spectrum_control():
    panel = json.loads(measure.RECORD.read_text())
    decoded = measure.decode(panel["counts"])
    assert decoded == panel["inference"]
    assert all(
        decoded[f"{k}_{t}"]["contractions"] == 10
        and decoded[f"{k}_{t}"]["expansions"] == 3
        for k in source.NAMES
        for t in (1, 4)
    )
    ref = decoded["equal_energy_time_control"]
    assert ref["exact_weight_volume_differences"][0][1] < -8e-5
    other = copy.deepcopy(panel["counts"])
    # Exchange and symmetric arrival are both represented; do not change total admissions.
    read = other["sharp_1"]["incoming_orientation"]
    read[2], read[1] = read[1], read[2]
    with pytest.raises(AssertionError):
        measure.decode(other)


def test_native_held_replay():
    result = measure.run()
    stable = json.loads(measure.RECORD.read_text())
    assert result["exact_equal_energy"] == stable["exact_equal_energy"]
    assert result["measured_forward_mixer"] == stable["measured_forward_mixer"]
    assert result["inference"]["sharp_1"]["volume_changes"][0][1] < -0.0049
    assert max(result["exact_equal_energy"]["volume_differences"][0]) < -0.00025
