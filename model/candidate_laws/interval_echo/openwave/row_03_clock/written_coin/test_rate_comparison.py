"""Check the two-species detector against a true common-rate change and the records."""

import json
import math

import numpy as np

from model.candidate_laws.interval_echo.openwave.row_03_clock.written_coin.rate_comparison import (
    HERE,
    compare,
)


def test_common_rate_control():
    starts, timer_starts = 2**20, 2**26
    panels = {}
    for name, frequency in {
        "8_even": 32,
        "8_odd": 33,
        "4_even": 16,
        "4_odd": 16.5,
    }.items():
        rows = []
        for t in np.linspace(0.02, 12.02, 193):
            timer = round(timer_starts * math.sin(t / 16) ** 2)
            flip = round(starts * math.sin(frequency * t / 2) ** 2)
            rows.append([timer, starts - flip, flip, 0])
        panels[name] = {"counts": rows}
    control = compare(
        {"panels": panels, "starts": starts, "timer_starts": timer_starts}
    )
    assert control["common_rate_compatible"]
    assert not control["equal_shift_compatible"]


def test_material_contrast_control():
    result = compare(json.loads((HERE / "records.json").read_text()))
    assert not result["common_rate_compatible"]
    assert result["equal_shift_compatible"]
    assert (
        result["predicted_4_odd_common_rate"][1] < result["frequency_bands"]["4_odd"][0]
    )
