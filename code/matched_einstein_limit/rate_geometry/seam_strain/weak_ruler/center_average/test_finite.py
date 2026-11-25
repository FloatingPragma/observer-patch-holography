"""Positive finite preparations with nonlinear area and regional count held fixed."""

import json

import pytest

from . import finite


def test_finite_repairing_area_count_and_independent_energy_entropy():
    native = json.loads((finite.HERE / "native.json").read_text())
    predicted = native["equal_individual_area_coefficient_and_regional_number"]
    archived = json.loads((finite.HERE / "finite.json").read_text())
    replay = finite.run()
    for x in (replay, archived):
        assert x["enrolled_ports"] == native["regional"]["ports"]
        assert x["repair_probability"] == 1e-7
        assert x["rail_weights"] == pytest.approx([0.66666669170, 0.33333330830],
                                                  abs=2e-11, rel=0)
        assert x["input_trace"] == pytest.approx(1., abs=1e-14, rel=0)
        assert x["minimum_input_eigenvalue"] > 1.4e-8
        assert x["independently_assembled_input_discrepancy"] < 3e-16
        for obs in ("number", "energy"):
            ref_lo, ref_hi = native["regional"][obs]["reference"]
            assert ref_lo-1e-12 < x["reference_" + obs] < ref_hi+1e-12
        a_lo, a_hi = native["total_individually_fixed_volume_area_coefficient"]
        assert a_lo < x["reference_individual_area_coefficient"] < a_hi
        assert x["area_change_at_count_match_without_rotation"] > 1.8e-10
        assert abs(x["constrained_area_change"]) < 1e-17
        assert abs(x["constrained_number_change"]) < 1e-18
        assert abs(x["directly_subtracted_count"]) < 1e-15
        assert x["directly_subtracted_energy"] == pytest.approx(
            x["regional_energy_change"], rel=0, abs=2e-15
        )
        for key, field in [("ordinary_per_repairing", "rail_secant"),
                           ("right_angle_per_repairing", "angle_secant")]:
            lo, hi = predicted[key]
            assert lo < x[field] < hi
        for key, field in [("energy", "energy_secant"), ("regional_entropy", "entropy_secant")]:
            lo, hi = predicted["matched_slopes"][key]
            assert lo < x[field] < hi
    assert replay["right_angle"] == pytest.approx(archived["right_angle"], rel=1e-7)
    assert replay["regional_energy_change"] == pytest.approx(archived["regional_energy_change"],
                                                           rel=1e-6)
    assert replay["regional_entropy_change"] == pytest.approx(archived["regional_entropy_change"],
                                                            abs=5e-14)
