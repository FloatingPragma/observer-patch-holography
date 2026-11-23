"""Check the matched marginals, constraint and independent held joint meters."""

import hashlib
import json
import math

import numpy as np

from . import correlation as trial


def test_matched_reduced_states_and_gravitational_constraint():
    # In each x sine/cosine oscillator, directly trace out the other field.
    X = np.array([[0.0, 1.0], [1.0, 0.0]]) / math.sqrt(2 * trial.REFERENCE)
    cross = []
    for sign in (1, -1):
        state = np.array([[0.0, 1.0], [float(sign), 0.0]]) / math.sqrt(2)
        for reduced in (state @ state.T, state.T @ state):
            np.testing.assert_allclose(reduced, np.eye(2) / 2)
        cross.append(
            2 * state.ravel() @ np.kron(X, X) @ state.ravel() / trial.scalar.VOL0
        )
    kappa = trial.observer_field.COUPLING
    for interaction in (0.0, 0.3):
        measured = []
        for sign in (1, -1):
            _, source = trial.field_source(sign, interaction)
            initial = source(trial.TIMES[0])
            np.testing.assert_array_equal(initial[:3], np.zeros(3))
            p = initial[3:6]
            h = kappa * (p - 0.5 * sum(p)) / trial.scalar.VOL0
            measured.append(h[0] * h[1] + h[1] * h[2] + h[2] * h[0])
        np.testing.assert_allclose(
            measured[0] - measured[1],
            kappa * interaction * (cross[0] - cross[1]),
            atol=5e-15,
        )


def test_preregistered_interaction_break_and_independent_geodesic_area():
    packet = json.loads((trial.HERE / "correlation.json").read_text())
    assert (
        packet["protocol_sha256"]
        == hashlib.sha256(
            (trial.HERE / "CORRELATION_PROTOCOL.md").read_bytes()
        ).hexdigest()
    )
    assert (
        packet["calibration_receipt_sha256"]
        == hashlib.sha256((trial.scalar.HERE / "results.json").read_bytes()).hexdigest()
    )
    exact_cross = 2 / (trial.REFERENCE * trial.scalar.VOL0)
    for name, arm in packet["arms"].items():
        result = trial.analyze(
            arm["records"],
            arm["apertures"],
            packet["calibrations"],
            packet["diagonal_mass_from_spectrometer"],
            arm["calibrated_interaction"],
        )
        held = result["radar"][-1]["times"]
        centre = held[0]
        assert centre["max_marginal_moment_difference"] < 0.002
        assert (
            abs(centre["measured_delta_cross_field_covariance"] - exact_cross)
            < 0.000003
        )
        assert max(row["max_complete_error"] for row in held) < 0.0002
        np.testing.assert_allclose(
            [row["geometry"] for row in held],
            [row["geometry"] for row in arm["result"]["radar"][-1]["times"]],
            atol=3e-10,
            rtol=1e-6,
        )
        area_contrast = np.array([row["delta_Guu"] for row in result["area"]])
        if name == "mixed":
            forecast = centre["complete_forecast"][0]
            assert centre["geometry"][0] > 0.001
            assert abs(centre["geometry"][0] - centre["marginal_forecast"][0]) > 0.001
            assert abs(forecast - area_contrast[-1]) < 0.0002
            assert all(np.diff(abs(area_contrast - forecast)) < 0)
            assert max(abs(v) for v in centre["delta_axial_null"]) < 0.00001
        else:
            for a, b in zip(
                arm["records"]["plus"], arm["records"]["minus"], strict=True
            ):
                assert a["radar"] == b["radar"]
            assert arm["apertures"]["plus"] == arm["apertures"]["minus"]
            assert abs(area_contrast[-1]) < 0.00001
