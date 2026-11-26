"""Held cross ranges on supplied smooth products and every admitted native angular triad."""
from fractions import Fraction
import json
import math
from pathlib import Path

import numpy as np
import pytest

from . import native, smooth


HERE = Path(__file__).parent


def test_all_independent_native_roundtrips_and_five_midpoint_registrations():
    actual = native.run()
    published = json.loads((HERE / "results.json").read_text())["native"]
    assert json.loads(json.dumps(actual["centers"])) == published["centers"]
    assert actual["read_gates"] == published["read_gates"]
    assert actual["native_exchange_hand"] == pytest.approx(math.pi)
    assert actual["distinct_roundtrips"] == 990
    assert actual["read_gates"]["echo_read_undo"] == 1980
    assert len(actual["centers"]) == 5
    for row in actual["centers"]:
        assert [(v["radius"], v["pairs"]) for v in row["antipodal_pairs_at_true_midpoint"]] == [(1, 4), (2, 7)]
        assert len(row["orthogonal_triads"]) == 2
        for triad in row["orthogonal_triads"]:
            assert triad["radii"] == [1, 2, 2]
            assert len(triad["fork_pairs"]) == 1
            assert triad["scalar"] == "-39/8"
            contractions = {Fraction(pair["contraction"]) for pair in triad["cross_pairs"]}
            assert contractions == {Fraction(-3), Fraction(3,4), Fraction(-3,16)}
            held = next(p for p in triad["cross_pairs"] if p["axes"] == triad["fork_pairs"][0])
            assert held["ranges"] == [3,3,3,3]


def test_signed_range_errors_and_unequal_radii_on_anisotropic_products():
    report = smooth.run()
    last = [r for r in report if r["source"] == "product" and r["short_radius"] == 1/32]
    for row in last:
        estimated = np.array([p["contraction"] for p in row["pairs"]])
        assert max(abs(estimated - row["plane_limits"])) < .00026
        assert abs(row["scalar"] - row["expected_scalar"]) < .00062
        lo, hi = row["scalar_clock_interval"]
        assert lo < row["scalar"] < hi
        assert hi - lo < .000005
    coarse = [r for r in report if r["source"] == "product" and r["short_radius"] == 1/4]
    assert max(abs(r["scalar"] - r["expected_scalar"]) for r in coarse) > .039


def test_radial_clocks_on_flat_space_detect_an_off_center_diameter():
    rows = [r for r in smooth.run() if r["source"] == "bent_flat"]
    for row in rows:
        ell = row["short_radius"]
        r, s = ell, 2*ell
        expected = 12 * math.sin(ell)**2 / (r*s)
        low, high = row["scalar_clock_interval"]
        assert low < expected < high
        assert 5.8 < low < high < 6.1
        assert row["axial_shortfall"][0] == pytest.approx(2*r*(1-math.cos(ell)), abs=row["clock_tick"])
