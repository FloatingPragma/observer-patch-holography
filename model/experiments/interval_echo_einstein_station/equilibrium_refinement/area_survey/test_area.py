"""Independent source, clock, interval, and geometric controls."""

import json
from fractions import Fraction
from math import log

import numpy as np

from . import acquire, certificate, measure


def load(name):
    return json.loads((measure.HERE / name).read_text())


def test_discounted_commute_clock_with_exact_two_state_controls():
    # Deterministic ping-pong has a two-step return; an immobile resident never returns.
    z = 1 - 1 / measure.MEAN_TIMEOUT
    occupied = np.eye(2, dtype=bool)
    for transition, expected in (
        (np.array([[0.0, 1.0], [1.0, 0.0]]), z * z),
        (np.eye(2), 0.0),
    ):
        for mask, _, success, _ in measure.discounted_pair(
            transition, *occupied, measure.MEAN_TIMEOUT
        ):
            np.testing.assert_allclose(success[mask], expected, atol=1e-12)


def test_same_source_multiplicities_and_transfer_catalogue():
    s = measure.survey(1)
    original = next(
        r
        for r in json.loads((measure.HERE.parent / "results.json").read_text())[
            "observations"
        ]
        if r["board_cells"] == 80 and len(r["cells"]) == 2
    )
    assert s.rim == original["cut_ports"] and s.cells == original["cells"]
    np.testing.assert_allclose(
        original["class_mean_native_energy"],
        [s.h.diagonal()[s.cut_count == j].mean() for j in range(3)],
    )
    np.testing.assert_array_equal(np.bincount(s.cut_count), [190, 80, 6])
    grad, _ = measure.differentials(s.uniform, s.starts, s.scores)
    dirs, ctrl, _ = measure.survey_directions(s, grad)
    assert abs(ctrl @ s.entropy_scores) < 1e-12 and abs(ctrl @ s.h.diagonal()) < 1e-12
    t = dirs["hidden_interior_frequency"]
    p = s.uniform + 1e-6 * t
    np.testing.assert_allclose(
        measure.distances(s.uniform, s.starts, s.scores),
        measure.distances(p, s.starts, s.scores),
    )
    assert abs(t @ s.h.diagonal()) < 1e-10 and t @ s.entropy_scores < -1.4
    uniform = load("results.json")["rows"][1]
    assert (
        uniform["timed"]["calibrated_eta"] > 0
        and uniform["all_attempts"]["calibrated_eta"] < 0
    )
    same = uniform["timed"]["directions"]["matched_energy_and_cut_load"]
    assert max(abs(same[k]) for k in ("energy", "cut", "volume")) < 1e-6
    np.testing.assert_allclose(
        measure.fifth(s, uniform["timed"]["edges"])["extra_squared_length"],
        uniform["held_fifth"]["extra_squared_length"],
        rtol=1e-6,
    )


def test_exact_native_all_records_and_independent_center_obstruction():
    receipt = certificate.run()
    assert receipt["trip_probability_error_upper"] < acquire.NATIVE_ERROR
    assert receipt["rows"] == 24 and receipt["right_inverse_defect_upper"] < 0.002
    assert receipt["positive_interval_halfwidth_lower"] > 3e-6
    assert receipt["held_geometry"] == load("certificate.json")["held_geometry"]
    assert receipt["held_geometry"]["flat_3d_extra_squared_length_lower"] >= 36000
    assert (4 * log(4) - 3 * log(6)) / 3 > 0.05663
    # With all launches fixed, 3*dP1+4*dP2=0. A shift in the two-rim class
    # preserves normalization and changes the independent first load moment by 2/3.
    d = (Fraction(1, 3), Fraction(-4, 3), Fraction(1))
    assert (
        sum(d) == 0 and 3 * d[1] + 4 * d[2] == 0 and d[1] + 2 * d[2] == Fraction(2, 3)
    )


def test_fresh_measurements_include_failed_starts_and_opposite_volume_arms():
    archive = load("counts.json")
    observed = acquire.run()
    for name, arms in archive["counts"].items():
        for sign, record in arms.items():
            assert sum(record["interior_counts"]) == acquire.SHOTS
            assert all(sum(row) == acquire.SHOTS for row in record["return_counts"])
            assert all(min(row) > 0 for row in record["return_counts"])
            assert abs(record["computed_volume_offset"]) < 1e-5
        for output in (archive, observed):
            assert (
                output["held"][name]["volume"][0]
                <= 0
                <= output["held"][name]["volume"][1]
            )
    for output in (archive, observed):
        assert output["calibrated_eta"][0] > 0
        assert output["held"]["zero_to_two"]["area"][1] < -39000
        assert output["held"]["zero_to_two"]["entropy"][1] < -1.78
        assert output["held"]["matched_energy_and_cut_load"]["residual"][1] < -6
