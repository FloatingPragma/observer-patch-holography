import json
from math import log, pi

from . import held_calibration as held


def test_proper_area_coefficient_and_held_transverse_control():
    report = held.run()
    assert report == json.loads(held.RECORD.read_text())
    region = report["three_sectors"]
    bands = region["simultaneous_99_percent"]
    for n, (lo, hi) in zip(held.TILES, bands["target_area_bands"], strict=True):
        assert lo < n * pi < hi
    for name in ("eta_fit_3_to_4", "eta_held_4_to_5"):
        lo, hi = bands[name]
        assert lo < log(2) / pi < hi and lo < region[name] < hi
    lo, hi = bands["held_capacity_minus_calibrated_area_nats"]
    assert lo < 0 < hi and lo < region["held_capacity_minus_calibrated_area_nats"] < hi
    assert report["finest_radius"] == 1 / 32
    for name, by_meter in report["held_frame_point_estimates"].items():
        for meter in by_meter.values():
            statistic = meter["held_max_abs_capacity_minus_eta_area_over_ell4"]
            if name == "maximally_symmetric":
                assert statistic < 0.004
            else:
                assert name == "curved_plane_product" and statistic > 6.5
