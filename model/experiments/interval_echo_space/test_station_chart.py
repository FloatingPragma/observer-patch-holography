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
                np.testing.assert_allclose(
                    clock_error, observation["stationary_clock_max_error"]
                )
                assert clock_error <= spec["clock_tolerance_radians"]
                for (i, j), row in zip(HELD_PAIRS, observation["held"], strict=True):
                    assert row["stations"] == [i, j]
                    echo, pulse = [], []
                    for a, b in ((i, j), (j, i)):
                        returns = archive[f"{prefix}echo_hands_{a}_{b}"]
                        half_time = (
                            read_hand(fresh[a, 2:] + returns) - phase[a, 0]
                        ) / 2
                        echo.append(np.average(half_time, weights=returns.sum(axis=1)))
                        panel = archive[f"{prefix}free_{a}_{b}"]
                        pulse_time = read_hand(panel.sum(axis=1)) - phase[a, 0]
                        pulse.append(
                            np.average(pulse_time, weights=panel[:, 0].sum(axis=1))
                        )
                    for kind, observed in (
                        ("echo", np.mean(echo)),
                        ("pulse", np.mean(pulse)),
                    ):
                        np.testing.assert_allclose(observed, row[kind], rtol=1e-12)
                        error = abs(row[f"{kind}_prediction"] / observed - 1)
                        np.testing.assert_allclose(
                            error, row[f"{kind}_relative_error"], rtol=1e-12
                        )
                        assert error > spec["relative_tolerance"]
                assert not observation["passed"]
    assert not receipt["source_free_pass"]


def test_round_product_accepts_an_actual_radial_coordinate_and_rejects_inconsistent_ranges():
    from .radial_product import BOOK, control, triangle_certificate

    book = json.loads(BOOK.read_text())
    settings = json.loads(PROTOCOL.read_text())
    assert control(settings, book)["passed"]
    # Equal angular separations permit equal distances at a common radial level
    # when the angular scale is free. A fourth short distance restricts that scale.
    angles = np.ones((4, 4)) - np.eye(4)
    paths = [(0, 1), (1, 2), (0, 2)]
    assert not triangle_certificate(angles, paths, [1.0, 1.0, 1.0], 0.05, 256)[
        "excludes_every_radius"
    ]
    assert triangle_certificate(
        angles, [*paths, (0, 3)], [1.0, 1.0, 1.0, 0.1], 0.05, 256
    )["excludes_every_radius"]
    receipt = json.loads((HERE / "radial_product.json").read_text())
    assert (
        receipt["counts_sha256"]
        == sha256((HERE / "radial_product_counts.npz").read_bytes()).hexdigest()
    )
    for board in receipt["boards"]:
        for menu in board["menus"]:
            for row in [menu["population"], *menu["acquisitions"]]:
                assert (
                    row["training_certificate"]["excludes_every_radius"]
                    and not row["passed"]
                )


def test_collar_depth_uses_hidden_words_and_keeps_every_detector_outcome():
    from .collar_depth import BOOK, lifted_operator
    from .shell_source import Board
    from .shells import protocol

    spec = json.loads(BOOK.read_text())
    board = Board(0, protocol("shell_protocol.json")["matching"])
    h, sites = lifted_operator(board, True)
    n = board.car.n
    a, b = sites[0]
    assert h[3 * n + a, 2 * n + b] == -4
    assert np.max(abs(h[3 * n :, n + b].toarray())) == 0
    receipt = json.loads((HERE / "collar_depth.json").read_text())
    path = HERE / "collar_depth_counts.npz"
    assert receipt["counts_sha256"] == sha256(path.read_bytes()).hexdigest()
    with np.load(path, allow_pickle=False) as archive:
        for level, row in enumerate(receipt["boards"]):
            for panel in row["experiments"]:
                step = spec["appointments"].index(panel["hand_radians"])
                key = f"l{level}_t{step}_site{panel['launch_site']}_word{panel['word']}_transfer{int(panel['transfer'])}"
                table = archive[key]
                assert table.shape == (4, 12 * row["cells"], 4)
                assert (
                    table.sum() == spec["shots"]
                    and table[3].sum() == panel["depth_two_count"]
                )
