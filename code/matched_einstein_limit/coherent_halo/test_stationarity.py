import json
from math import pi
from pathlib import Path
import numpy as np

from ..fixed_volume.chord_kernel.measure import flat_quadrature
from .stationarity import integrated_chords, run
from model.experiments.interval_echo_einstein_station.collar_pairing.measure import (
    compare,
)


def test_integer_sum_matches_independent_point_clock():
    for n in (3, 7, 16):
        np.testing.assert_allclose(
            integrated_chords(n), flat_quadrature(n)["integral_per_radius_fourth"][0]
        )


def test_fixed_clock_source_has_nonzero_stationarity_deficit():
    data = run()
    compare(data, json.loads(Path(__file__).with_name("stationarity.json").read_text()))
    rows = data["panels"]
    target = 48 * 8 * pi * pi / 15
    assert all(
        0 < r["maximum_switch_bound"] < 1 and r["normalized_area_derivative"] == 0
        for r in rows
    )
    assert abs(rows[-1]["normalized_entropy_derivative"] / target - 1) < 0.02
    assert all(r["normalized_entropy_derivative"] > 0.95 * target for r in rows)
