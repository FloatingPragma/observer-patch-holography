"""Independent clock integrations and reconstruction of the stored observer records."""

from collections import Counter
from decimal import Decimal
from fractions import Fraction
import json
from math import exp, log
from pathlib import Path

import pytest

from level202_distance_redshift import clock_comparison
from level202_vacuum_density import readout

HERE = Path(__file__).resolve().parent


def test_late_clocks_admit_distinct_gravitating_vacuum_densities():
    raw = (HERE / "level202_distance_redshift.json").read_bytes()
    result = readout(raw)
    assert result["observer_vacuum_energy_density_J_m3"] is None
    assert result["beacons_with_range_and_period_at_calibration"] == 95
    assert result["range_in_light_ticks"] == [2, 19]
    assert result["static_geometries_matching_all_calibrated_beacons_at_last_query"] == [
        "Minkowski", "Einstein_static_dust"
    ]
    for model in result["static_geometries"].values():
        radius = model["curvature_radius_in_light_ticks"]
        k = Fraction(1, radius**2) if radius is not None else 0
        vacuum = Fraction(model["Lambda_times_light_tick_squared"])
        dust = Fraction(model["eight_pi_G_dust_mass_density_times_tick_squared"])
        assert (dust + vacuum) / 3 == k
        assert vacuum / 3 == dust / 6
        assert model["exact_timestamp_matches_by_query_tick"] == {
            "0": 60, "4": 62, "16": 95, "96": 95
        }
    assert result["static_geometries"]["Einstein_static_dust"]["Lambda_times_light_tick_squared"] == "1/10000"

    # Change only the second received pulse at one held-out beacon, keeping its radar range.
    survey = json.loads(raw)
    target = next(row for row in survey["panels"]["96"] if row["echo_at"] is not None)
    target["beacon_received"][1] += 1
    altered = readout(json.dumps(survey).encode())
    assert altered["static_geometries_matching_all_calibrated_beacons_at_last_query"] == []
    for model in altered["static_geometries"].values():
        assert model["exact_timestamp_matches_by_query_tick"]["96"] == 94


def test_integrated_null_rays():
    # Along either radial light ray the difference of exp(-Ht) equals H times
    # the beacon's comoving separation. Integrate from the original clock events.
    h = 0.02
    k = 0.12
    period = 4
    span = 20

    def receive(t):
        return -log(exp(-h * t) - k) / h

    reflection = receive(0)
    d = receive(reflection) / 2
    expected_z = (
        receive(reflection + period / 2) - receive(reflection - period / 2)
    ) / period - 1
    expected_next = (receive(receive(span)) - span) / 2
    expected_excess = receive(reflection + span) - receive(reflection) - span
    got = clock_comparison(d, Decimal(str(h)), span)
    assert got["z_centered"] == pytest.approx(expected_z, rel=1e-12)
    assert got["next_radar_distance"] == pytest.approx(expected_next, rel=1e-12)
    assert got["extra_received_ticks_for_long_source_period"] == pytest.approx(
        expected_excess, rel=1e-12
    )


def test_survey_clocks_and_matched_baselines():
    data = json.loads((HERE / "level202_distance_redshift.json").read_text())
    panels = data["panels"]
    actual_counts = {}
    for epoch, rows in panels.items():
        assert len({tuple(r["source"]) for r in rows}) == data["enrolled"] == 274
        counts = Counter()
        for row in rows:
            if row["radar_distance"] is not None:
                assert row["radar_distance"] == (row["echo_at"] - row["query_at"]) / 2
            if row["z_pulse"] is not None:
                a, b = row["beacon_received"]
                c, d = row["beacon_sent"]
                assert d - c == 4 and (c + d) / 2 == row["reflection_at"]
                assert (b - a) / (d - c) - 1 == row["z_pulse"]
                if row["radar_distance"] is not None:
                    counts[row["z_pulse"]] += 1
        actual_counts[epoch] = counts
    assert actual_counts == {
        "0": {0: 91, 0.25: 5},
        "4": {0: 94, 0.25: 1},
        "16": {0: 95},
        "96": {0: 95},
        "transparent": {0: 274},
    }
    assert 33 == data["last_closure_in_ball"] < 96
    early, middle, settled = [
        {
            tuple(r["source"]): r["radar_distance"]
            for r in panels[epoch]
            if r["radar_distance"] is not None and r["z_pulse"] is not None
        }
        for epoch in ("0", "16", "96")
    ]
    assert settled == middle and settled.keys() < early.keys()
    assert Counter(d - early[s] for s, d in settled.items()) == {
        0: 60,
        0.5: 28,
        1: 4,
        1.5: 3,
    }
    comparison = data["vacuum_comparison"]
    late = {
        tuple(r["source"]): r
        for r in panels["96"]
        if r["radar_distance"] is not None and r["z_pulse"] is not None
    }
    records = comparison["matched_long_baseline"]
    assert {tuple(r["source"]) for r in records} == set(late)
    assert len(records) == 95
    h = float(comparison["h_per_tick_if_tick_is_planck_time"])
    span = int(comparison["duration_ticks"])
    assert h * span == pytest.approx(log(2), rel=2e-15)
    for r in records:
        original = late[tuple(r["source"])]
        d = original["radar_distance"]
        assert r["first_source_pulse"] == original["reflection_at"]
        assert r["first_received"] == original["echo_at"]
        assert r["later_roundtrip_ticks"] == 2 * d
        assert r["late_source_pulse_offset"] == r["late_receipt_offset"] == 0
        assert r["predicted"]["z_centered"] / (h * d) == pytest.approx(1, abs=2e-15)
        assert r["predicted"]["next_radar_distance"] == 2 * d
        assert r["predicted"]["extra_received_ticks_for_long_source_period"] == d
