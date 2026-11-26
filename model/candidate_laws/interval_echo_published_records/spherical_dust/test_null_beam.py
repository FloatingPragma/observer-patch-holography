import copy
from fractions import Fraction as F

import pytest

from . import null_beam as n


def test_matched_screens_pointwise_work_and_fixed_affine_work():
    data = n.verify()
    for name in ("ordinary", "held"):
        readings = [data["stocks"][name][str(d)]["readings"] for d in n.SCALES]
        assert all(abs(row["signed_point_residual"]) < row["point_error_bound"] for row in readings)
        assert all(abs(row["signed_integral_residual"]) < row["integral_error_bound"] for row in readings)
        assert all(abs(b["signed_point_residual"]) < 0.26 * abs(a["signed_point_residual"])
                   for a, b in zip(readings, readings[1:]))
    signed = data["held_minus_ordinary"]["128"]
    assert signed["midpoint_focusing"] > 0.58
    assert abs(signed["twice_affine_deficit"] - signed["weighted_work_integral"]) < 1.1e-9


def test_controls_separate_geometry_and_matter_at_identical_gain():
    data = n.measure()["stocks"]
    h, weak, blind = [data[name]["128"] for name in ("held", "half_collar_force", "blind_collar_work")]
    assert h["acquisition"]["events"] == blind["acquisition"]["events"]
    assert h["acquisition"]["shell_mass_tickets"] == weak["acquisition"]["shell_mass_tickets"]
    assert h["readings"]["ordinary_release_gain"] == weak["readings"]["ordinary_release_gain"] == blind["readings"]["ordinary_release_gain"]
    assert weak["readings"]["signed_point_residual"] < -0.25
    assert blind["readings"]["signed_point_residual"] > 0.52
    assert weak["readings"]["signed_integral_residual"] < -0.00066
    assert blind["readings"]["signed_integral_residual"] > 0.00139


def test_affine_and_three_axis_optical_reader_reject_corruption():
    data = n.measure()
    packet = data["stocks"]["held"]["32"]["acquisition"]
    release = data["ordinary_release"]["32"]
    for index in (0, 3):
        changed = copy.deepcopy(packet)
        changed["events"][3][index] += 10**7
        with pytest.raises(ValueError, match="calibration"):
            n.decode(changed, 32, release, F(1, 3))


def test_raw_work_channel_can_change_without_affecting_the_screen():
    data = n.measure()
    arm = data["stocks"]["held"]["32"]
    changed = copy.deepcopy(arm["acquisition"])
    changed["shell_mass_tickets"] = [[t[0], 0, 0] for t in changed["shell_mass_tickets"]]
    read = n.decode(changed, 32, data["ordinary_release"]["32"], F(1, 3))
    assert read["midpoint_focusing"] == arm["readings"]["midpoint_focusing"]
    assert read["twice_affine_deficit"] == arm["readings"]["twice_affine_deficit"]
    assert read["midpoint_work_forecast"] == read["weighted_work_integral"] == 0
