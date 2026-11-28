"""Replay and verify the deterministic signs of the two matched material responses."""

from fractions import Fraction
import json

from . import time_four as t


def test_exact_connected_radius_counterexample_on_both_mirror_launches():
    data = t.run()
    assert data == json.loads((t.m.HERE / "time_four.json").read_text())
    assert [data["launch"], data["exact_reflected_launch"]] == [2, 13]
    assert [(x["dimension"], x["steps"]) for x in data["sectors"]] == [
        (96, 128),
        (1104, 256),
    ]
    known = {"target": [69104635, 69104665], "phase_feedback": [-134266663, -134266648]}
    ordinary = json.loads((t.m.HERE / "results.json").read_text())["standards"]
    for name, interval in known.items():
        receipt = data["contrasts"][name]
        lower, upper = [Fraction(i, 10**9) for i in interval]
        center = Fraction(receipt["fixed_point_center"])
        error = Fraction(receipt["absolute_error_bound"])
        assert receipt["radius_billionths"] == interval
        assert lower <= center - error < center + error <= upper
        for launch in ordinary[name]["probes"]:
            assert lower < launch["radial"][2] < upper


def test_billionth_rounding_handles_negative_endpoints_and_grid_points():
    q = Fraction(1, 2 * t.BILLION)
    assert t.bounds(0 * q, q) == [-1, 1]
    assert t.bounds(-3 * q, q) == [-2, -1]
