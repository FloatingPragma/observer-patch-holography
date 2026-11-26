import json
from fractions import Fraction

import pytest

from .exterior_reader import HERE, INNER, TICK, check_archive, run


@pytest.fixture(scope="module")
def survey():
    return run()


def test_independent_meters_predict_held_counts_and_resolve_blind_controls(survey):
    rows = survey["readouts"]
    held = rows["held"]
    reference = rows["ordinary"]
    matter, collar = map(Fraction, survey["ordinary_and_collar_counts"]["held"])
    assert (matter, collar) == (Fraction(1, 3), Fraction(8, 33))
    for meter, gain in survey["ordinary_gain_by_meter"].items():
        assert gain == pytest.approx(0.25, abs=1e-8)
        prediction = survey["held_forecasts_from_ordinary"][meter]
        assert prediction["mu"] == pytest.approx(held[meter], abs=3e-9)
        assert prediction["redshift_ratio"] == pytest.approx(held["redshift_ratio"], abs=3e-9)
        assert prediction["initial_radial_acceleration"] == pytest.approx(
            -held["mu_from_two_release_intervals"]/INNER**2, abs=3e-9
        )
        assert rows["blank"][meter] == 0
        assert rows["both_collar_blind"][meter] == reference[meter]
    split = rows["light_alone_collar_blind"]
    assert split["mu_from_clock_spacing"] == reference["mu_from_clock_spacing"]
    assert split["mu_from_two_release_intervals"] == held["mu_from_two_release_intervals"]
    assert split["mu_from_two_release_intervals"] - split["mu_from_clock_spacing"] > 0.060
    assert (held["echo_inner_ticks"] - reference["echo_inner_ticks"])*TICK > 0.048
    for row in (reference, held):
        assert abs(row["spectroscopic_ratio"] - row["redshift_ratio"]) < 2e-10
        assert row["mirror_recovered_frequency_relative_error"] < 2e-10
    assert rows["blank"]["echo_inner_ticks"]*TICK == 4


def test_flight_and_massive_trajectory_refinement(survey):
    for name in ("ordinary", "held"):
        row = survey["readouts"][name]
        rays = row["signed_rays"]
        assert [r["panels"] for r in rays] == [4, 8, 16, 32, 64]
        for r in rays:
            assert max(r["outbound_error_upper"], r["inbound_error_upper"]) <= r["exact_arithmetic_time_error_bound"] + 2e-14
        for key in ("outbound_error_upper", "inbound_error_upper"):
            for a, b in zip(rays, rays[1:]):
                assert 14 < a[key]/b[key] < 18
        assert rays[-1]["outbound_error_upper"] < 1.4e-10
        assert rays[-1]["inbound_error_upper"] < 1e-11
        releases = row["released_clocks"]
        for duration, m1 in releases["1"].items():
            m3 = releases["3"][duration]
            assert abs(m1["radius_ticks"] - m3["radius_ticks"]) <= 2
            assert abs(m1["action_time_ticks"] - m3["action_time_ticks"]) <= 2
            assert m1["probe_proper_ticks"] == m3["probe_proper_ticks"]
        levels = releases["1"]["1/2"]["radius_error_upper_by_refinement"]
        for a, b in zip(levels, levels[1:]):
            assert 14 < a/b < 20
        assert levels[-1] < 9e-13


def test_independent_reader_agrees_with_analytic_clock_archive(survey):
    previous = json.loads((HERE / "results.json").read_text())["exterior"]
    for name in ("ordinary", "held"):
        row = survey["readouts"][name]
        analytic = previous[name]
        assert analytic["tick"] == TICK
        for measured, target in zip(row["received_outer_ticks"], analytic["received"], strict=True):
            assert abs(measured - target)*TICK <= row["signed_rays"][-1]["outbound_error_upper"] + 2*TICK
        for duration, target in zip(("1/16", "1/32"), analytic["released_radius_ticks"], strict=True):
            for mass in ("1", "3"):
                assert abs(row["released_clocks"][mass][duration]["radius_ticks"] - target) <= 2
    check_archive(survey, json.loads((HERE / "exterior_reader.json").read_text()))
