"""Verify the held-pair boundary and the published acquisitions."""

from hashlib import sha256
import json

import numpy as np

from .station_chart import ANCHOR_PAIRS, HELD_PAIRS, HERE, PROTOCOL, TRAIN, chart


def test_held_records_cannot_move_the_prediction():
    points = np.array(
        [[0, 0, 0], [3, 0, 0], [0, 4, 0], [0, 0, 5], [1, 2, 1], [-2, 1, 3], [2, -1, 1]]
    )
    distances = np.linalg.norm(points[:, None] - points, axis=2)
    radar = [distances[i, j] for i, j in TRAIN]
    pulse = [1 + 0.8 * distances[i, j] for i, j in ANCHOR_PAIRS]
    held = np.array([distances[i, j] for i, j in HELD_PAIRS])
    spec = json.loads(PROTOCOL.read_text())
    correct = chart(radar, pulse, held, 1 + 0.8 * held, 0.0, spec)
    corrupted = chart(radar, pulse, 2 * held, 4 * held, 0.5, spec)
    assert correct["passed"] and not corrupted["passed"]
    for good, bad in zip(correct["held"], corrupted["held"], strict=True):
        assert good["echo_prediction"] == bad["echo_prediction"]
        assert good["pulse_prediction"] == bad["pulse_prediction"]


def test_archived_counts_cover_the_enrolled_menus():
    receipt = json.loads((HERE / "station_chart.json").read_text())
    spec = json.loads(PROTOCOL.read_text())
    assert receipt["protocol_sha256"] == sha256(PROTOCOL.read_bytes()).hexdigest()
    path = HERE / "station_chart_counts.npz"
    assert receipt["counts_sha256"] == sha256(path.read_bytes()).hexdigest()
    with np.load(path, allow_pickle=False) as archive:
        keys = set()
        for level in spec["levels"]:
            for rep in range(spec["replicas"]):
                prefix = f"l{level}_r{rep}_"
                for kind in ("calibration_hands", "independent_hands"):
                    key = prefix + kind
                    keys.add(key)
                    assert np.all(archive[key].sum(axis=-1) == spec["shots"])
                for pairs, kind in (
                    (TRAIN, "first"),
                    (ANCHOR_PAIRS + HELD_PAIRS, "free"),
                    (HELD_PAIRS, "echo_hands"),
                    (HELD_PAIRS, "echo_timeouts"),
                ):
                    for left, right in pairs:
                        for a, b in ((left, right), (right, left)):
                            key = f"{prefix}{kind}_{a}_{b}"
                            data = archive[key]
                            keys.add(key)
                            assert np.issubdtype(data.dtype, np.integer) and np.all(
                                data >= 0
                            )
                            if kind == "first":
                                assert data.shape == (spec["appointments"] + 1,)
                                assert data.sum() == spec["shots"]
                            if kind == "free":
                                assert data.shape == (spec["appointments"], 2, 4)
                                assert np.all(data.sum(axis=(1, 2)) == spec["shots"])
                            if kind == "echo_hands":
                                assert data.shape == (2 * spec["appointments"] - 1, 4)
                                assert data.sum() == spec["return_quota"]
                            if kind == "echo_timeouts":
                                assert data.shape == (2,)
        assert keys == set(archive.files)


def test_held_clock_and_flight_results_follow_from_the_saved_counts():
    from .shells import read_hand

    spec = json.loads(PROTOCOL.read_text())
    receipt = json.loads((HERE / "station_chart.json").read_text())
    assert receipt["flat_control"]["passed"]
    assert [b["cells"] for b in receipt["boards"]] == [20, 80, 320]
    with np.load(HERE / "station_chart_counts.npz", allow_pickle=False) as archive:
        for level, board in zip(spec["levels"], receipt["boards"], strict=True):
            assert not board["population"]["passed"]
            for rep, observation in enumerate(board["acquisitions"]):
                prefix = f"l{level}_r{rep}_"
                calibration = read_hand(archive[prefix + "calibration_hands"])
                fresh = archive[prefix + "independent_hands"]
                phase = read_hand(fresh)
                clock_error = np.max(np.abs(phase[4:] - calibration[:4].mean(axis=0)))
                np.testing.assert_allclose(clock_error, observation["stationary_clock_max_error"])
                assert clock_error <= spec["clock_tolerance_radians"]
                for (i, j), row in zip(HELD_PAIRS, observation["held"], strict=True):
                    assert row["stations"] == [i, j]
                    echo, pulse = [], []
                    for a, b in ((i, j), (j, i)):
                        returns = archive[f"{prefix}echo_hands_{a}_{b}"]
                        half_time = (read_hand(fresh[a, 2:] + returns) - phase[a, 0]) / 2
                        echo.append(np.average(half_time, weights=returns.sum(axis=1)))
                        panel = archive[f"{prefix}free_{a}_{b}"]
                        pulse_time = read_hand(panel.sum(axis=1)) - phase[a, 0]
                        pulse.append(np.average(pulse_time, weights=panel[:, 0].sum(axis=1)))
                    for kind, observed in (("echo", np.mean(echo)), ("pulse", np.mean(pulse))):
                        np.testing.assert_allclose(observed, row[kind], rtol=1e-12)
                        error = abs(row[f"{kind}_prediction"] / observed - 1)
                        np.testing.assert_allclose(error, row[f"{kind}_relative_error"], rtol=1e-12)
                        assert error > spec["relative_tolerance"]
                assert not observation["passed"]
    assert not receipt["source_free_pass"]
