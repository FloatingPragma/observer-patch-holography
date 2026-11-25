"""Verify common reflection times, fixed cohorts and dimensional scaling."""

import json
import math
from fractions import Fraction as F

import pytest

from observer_zero_expansion import C, ECHO_FILE, HERE, read_bank, report, window


def test_windows_on_an_accelerating_ruler_and_an_incomplete_source():
    # R(t)=2+t^2/2 at t=1,2,3 for the first marker.
    tracks = {(1, 0): {1: F(5, 2), 2: F(4), 3: F(13, 2)}, (2, 0): {1: F(10), 3: F(10)}}
    result = window(tracks, 2, 1, F(1, 1000))
    assert result["sources"] == [(1, 0)]
    m = result["measurement"]
    assert (m["H_per_tick"], m["relative_acceleration_per_tick2"]) == (F(1, 2), F(1, 4))
    assert (m["speed_m_per_s"], m["range_acceleration_m_per_s2"]) == (2 * C, 1000 * C)
    assert (m["H_per_s"], m["relative_acceleration_per_s2"]) == (500, 250000)
    assert m["mean_range_m"][1] == F(C, 250)
    assert window(tracks, 5, 1, F(1))["measurement"] is None
    with pytest.raises(ValueError):
        window(tracks, 2, 0, F(1))
    with pytest.raises(ValueError):
        window(tracks, 2, 1, F(0))


def test_archived_early_and_complete_surveys():
    path = HERE / ECHO_FILE
    raw = json.loads(path.read_bytes())
    tracks = read_bank(path)
    result = report(tracks, F(1))
    # Check the early cohort directly from integer clock receipts. Twice the
    # assigned reflection time equals twice send time plus the round trip.
    rows = {t: {} for t in (12, 16, 20)}
    for send, waits in zip(
        raw["origin_emission_ticks"], raw["echo_waits"], strict=True
    ):
        for label, wait in zip(raw["reflectors"], waits, strict=True):
            twice_time = 2 * send + wait
            if twice_time % 2 == 0 and twice_time // 2 in rows:
                rows[twice_time // 2][tuple(label)] = wait
    labels = set.intersection(*(set(row) for row in rows.values()))
    assert len(labels) == 102
    early = result["early_common_time_window"]
    m = early["measurement"]
    assert set(early["sources"]) == labels
    assert m["mean_range_relay_units"] == [
        F(sum(row[p] for p in labels), 204) for row in rows.values()
    ]
    assert m["speed_over_c"] == F(1, 408)
    assert m["H_per_tick"] == F(1, 2844)
    assert m["range_acceleration_relay_units_per_tick2"] == -F(1, 816)
    assert m["relative_acceleration_per_tick2"] == -F(1, 5688)
    assert sorted(m["changing_sources"]) == [(67108865, 4), (67108865, 10)]
    assert result["integer_radar_ticks_with_full_coverage"] == list(range(20, 115))
    assert result["all_individual_ranges_constant_from_first_common_time"]
    for key in ("same_early_beacons_later", "full_common_time_window"):
        late = result[key]
        assert (
            late["measurement"]["speed_over_c"]
            == late["measurement"]["relative_acceleration_per_tick2"]
            == 0
        )
    assert result["full_common_time_window"]["beacons"] == 156
    assert result["full_common_time_window"]["remote_observers"] == 13
    assert len(result["full_cohort_send_history"]) == 113
    for row, total in zip(
        result["full_cohort_send_history"], raw["total_round_trip_ticks"], strict=True
    ):
        assert row["mean_range_relay_units"] == F(total, 312)
    launch = result["all_beacons_at_send_ticks_0_4_8"]["measurement"]
    assert launch["speed_over_c"] == F(5, 64)
    assert launch["relative_acceleration_per_tick2"] == -F(95, 43616)
    scaled = report(tracks, F(1, 1000), "radar", 16, 4)["window"]["measurement"]
    assert scaled["speed_m_per_s"] == m["speed_m_per_s"]
    assert scaled["H_per_s"] == 1000 * m["H_per_s"]
    assert (
        scaled["relative_acceleration_per_s2"]
        == 10**6 * m["relative_acceleration_per_s2"]
    )


def test_local_clock_and_planck_rescaling():
    import numpy as np
    from scipy.linalg import expm
    from observer_zero_information import LAW, candidate

    matrix = LAW.seam_term(candidate.KAPPA)
    ids = [LAW.SRC.index(s) for s in ((1, 0), (0, 1))]
    block = matrix[np.ix_(ids, ids)]
    assert np.all(np.delete(matrix[:, ids], ids, axis=0) == 0)
    eig = np.linalg.eigvalsh(block)
    np.testing.assert_allclose(eig, [8, 9])
    U = expm(-1j * math.pi * block)
    assert abs(U[1, 0]) ** 2 == pytest.approx(1)
    assert abs((U @ U)[0, 0]) ** 2 == pytest.approx(1)
    clocks = json.loads((HERE / "observer_zero_information.json").read_bytes())
    for record in clocks["records"]:
        for setting in ("immediate", "delayed"):
            flips = record[setting]["complete_flip_reads"]
            assert len(flips) >= 3 and flips == [i % 2 for i in range(len(flips))]
    tracks = read_bank(HERE / ECHO_FILE)
    tau = F("5.391247e-44")
    one = report(tracks, tau)
    milestones = one["elapsed_from_preparation"]
    assert [m["tick"] for m in milestones] == [
        max(r["archive_round"] for r in clocks["local_position_readout"]["records"]),
        96,
        144,
    ]
    assert milestones[1]["seconds_at_supplied_tick_duration"] == F("5.17559712e-42")
    hypothesis = one["one_planck_time_hypothesis"]
    assert hypothesis["clock_gap_radians_per_tick"] == pytest.approx(
        (eig[1] - eig[0]) * math.pi
    )
    m = one["early_common_time_window"]["measurement"]
    double = report(tracks, 2 * tau, "radar", 16, 4)["window"]["measurement"]
    assert float(m["H_per_s"]) == pytest.approx(hypothesis["early_window_H_per_s"])
    assert m["speed_m_per_s"] == double["speed_m_per_s"]
    assert m["mean_range_m"] == [r / 2 for r in double["mean_range_m"]]
    assert (
        m["relative_acceleration_per_s2"] == 4 * double["relative_acceleration_per_s2"]
    )
    assert hypothesis["measured_planck_times_per_tick"] is None
