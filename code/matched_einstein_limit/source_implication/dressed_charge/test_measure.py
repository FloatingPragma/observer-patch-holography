import json
from pathlib import Path

import numpy as np

from . import interferometer, measure


def test_native_compact_dressing_signed_current_and_compiled_meter():
    data = measure.run()
    stored = json.loads(Path(measure.__file__).with_name("results.json").read_text())
    assert (
        data["eigenvalues_of_integer_flight"]
        == stored["eigenvalues_of_integer_flight"]
        == {"-1": [-3, -3], "1": [29, 29]}
    )
    assert data["initial_current_rate_per_unit_head_charge"] == "-3/8"
    assert data["terminal_settings"] == 61
    assert data["active"] == [4, 5, 6, 7, 8, 11]
    assert np.allclose(
        data["twice_restricted_dressing"], stored["twice_restricted_dressing"]
    )
    np.testing.assert_allclose(
        [rec["signed_internal_current"] for rec in data["records"]],
        [0.01719640674658, -0.01155533027501],
        atol=1e-12,
        rtol=0,
    )
    for fresh, existing in zip(data["records"], stored["records"]):
        assert abs(fresh["current_into_six"]) < 1e-14
        assert fresh["finite_difference_error"] < 2e-9
        assert abs(fresh["terminal_read"] - data["parked_expectation"]) < 1e-12
        assert abs(fresh["bare_head_change"] - existing["bare_head_change"]) < 1e-12


def test_launch_disturbance_and_projector_variance():
    spatial = (
        np.array(
            json.loads(Path(measure.__file__).with_name("results.json").read_text())[
                "twice_restricted_dressing"
            ]
        )
        / 2
    )
    projector = spatial / 3
    np.testing.assert_allclose(projector @ projector, projector, atol=1e-14, rtol=0)
    for k in range(6):
        sharp = np.diag(np.eye(6)[k])
        assert np.linalg.norm(sharp @ spatial - spatial @ sharp, ord=2) > 1.4
        assert spatial[k, k] == 1 and (spatial @ spatial)[k, k] == 3


def test_held_pair_meter_with_unoccupied_and_occupied_spectator():
    board = measure.station.iso.pairing.formation.federation(1)
    adj = interferometer.register(board)
    target = np.zeros((7, 7))
    target[3, 6] = target[6, 3] = 1.0
    settings = interferometer.terminal_settings(adj, target)
    composed = sum(weight * effect for weight, effect in settings)
    np.testing.assert_allclose(composed, target, atol=2e-14, rtol=0)
    rng = np.random.default_rng(6430)
    for outside in (0, 0.4):
        amplitudes = rng.normal(size=7) + 1j * rng.normal(size=7)
        amplitudes *= np.sqrt(1 - outside) / np.linalg.norm(amplitudes)
        signal = (amplitudes.conj() @ composed @ amplitudes).real
        np.testing.assert_allclose(
            signal, 2 * (amplitudes[3].conj() * amplitudes[6]).real, atol=1e-14
        )
