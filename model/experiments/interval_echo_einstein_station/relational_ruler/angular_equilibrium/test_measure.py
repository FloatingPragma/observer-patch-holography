"""Check the native sector, complete-board controls and the registered outcome."""

import hashlib
import json
from types import SimpleNamespace

import numpy as np

from model.executable_model import federation
from . import measure as m


def test_full_native_trial_sector():
    board = federation(0)
    routes = m.ActiveRoutes(
        SimpleNamespace(cells=board.cells, external=np.array(board.external))
    )
    blanks, energies = m.local_calibration()
    assert energies["coin_odd"]["source_energy"] == 32
    assert energies["coin_odd"]["measured_potential"] == [0.0, 0.0]
    native, parked, enabled, guards, _ = m.build_station(board, blanks)
    assert tuple(sorted((parked[0], routes.anchor()))) in guards
    operation, bins, multiplicities = m.stencil(routes, 2, 3, len(enabled))
    np.testing.assert_array_equal(multiplicities, np.bincount(bins))
    assert len(bins) == len(enabled)
    np.testing.assert_allclose(
        operation.toarray(), native.toarray() / 2.5 - np.eye(len(bins)), atol=1e-13
    )
    row = m.read_window(routes, 2, 3, 0)
    values, vectors = np.linalg.eigh(native.toarray())
    free = (vectors * np.exp(-1j * values)) @ vectors.T
    src = np.flatnonzero(bins == 3)
    for theta, p in zip(m.PHASES, row["phase_probabilities"], strict=True):
        applied = np.where(bins == 3, np.exp(1j * theta), 1)
        final = free[src, :] @ (applied[:, None] * free[:, src])
        np.testing.assert_allclose(
            np.square(abs(final)).sum() / len(src), p, rtol=0, atol=2e-12
        )
    row_probability = row["outcomes"]["uninterrupted"]["return_probability"]
    assert abs(row_probability - row["phase_probabilities"][0]) < 2e-12


def test_archive_and_registered_boundaries():
    packet = json.loads((m.HERE / "results.json").read_text())
    for fname, field in (
        ("PROTOCOL.md", "protocol_sha256"),
        ("PHASE_PROTOCOL.md", "phase_protocol_sha256"),
    ):
        assert (
            hashlib.sha256((m.HERE / fname).read_bytes()).hexdigest() == packet[field]
        )
    assert [p["cells"] for p in packet["panels"]] == [80, 1280, 20480, 327680]
    assert packet["complete_board_control"]["complete_mobile_ports"] == 880
    held_errors = []
    for panel in packet["panels"]:
        for row in panel["windows"]:
            assert row["single_flight_norm_tail"] < 1e-13
            assert (
                sum(row["one_leg_counts"])
                == sum(row["double_load_one_leg_counts"])
                == m.SHOTS
            )
            for arm in row["outcomes"].values():
                assert sum(arm["counts"]) == sum(arm["double_load_counts"]) == m.SHOTS
            outcomes = row["outcomes"]
            assert (
                abs(
                    row["refreshed_prediction_from_one_leg_counts"]
                    - outcomes["angular_refresh"]["counts"][0] / m.SHOTS
                )
                < 0.00020
            )
            assert row["phase_acquisitions"] == m.phase_receipts(
                row["phase_probabilities"], panel["level"], row["window"]
            )
            for acquisition in row["phase_acquisitions"]:
                for held in acquisition["holdouts"]:
                    assert held["absolute_error"] < held["simultaneous_bound"]
                    held_errors.append(held["absolute_error"])
    assert len(held_errors) == 32
    assert max(held_errors) < 0.000064
    for row in packet["panels"][-1]["windows"]:
        assert row["absolute_refresh_error"] > 0.01
        assert row["error_fraction_of_departure"] > 0.10
        assert max(row["phase_probabilities"]) - min(row["phase_probabilities"]) > 0.006
