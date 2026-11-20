"""Independent checks of the native spherical atlas and measured radar read."""

import json
import math

import numpy as np

from .geometry import _pulse, run


def test_isolated_massless_link_matches_analytic_detector_populations():
    length, speed, wait = 0.8, 1.1, 0.4
    omega_wait = math.sqrt(2) * speed * wait / length
    cosine = math.cos(omega_wait)
    expected = [(1 + cosine)**2 / 4, (1 - cosine)**2 / 4, math.sin(omega_wait)**2 / 2]
    np.testing.assert_allclose(abs(_pulse(length, speed, wait)[:, 0])**2, expected,
                               atol=1e-13, rtol=0)


def test_independently_enrolled_length_and_clock_transfer_to_native_heldouts():
    reference_length, reference_frequency, speed = 1.75, 0.8, 1.3
    result = run(reference_length=reference_length, reference_frequency=reference_frequency, speed=speed)
    json.dumps(result, allow_nan=False)
    np.testing.assert_allclose(
        [result["calibrated_speed"], result["relay_speed"]],
        [speed, math.sqrt(2) * speed / math.pi],
        atol=1e-13, rtol=0,
    )
    assert result["calibration"]["reference_length"] == reference_length
    assert result["reference_frequency"] == reference_frequency
    assert [len(level["cell_centers"]) for level in result["levels"]] == [20, 80, 320]
    assert [len(level["seam_lengths"]) for level in result["levels"]] == [30, 120, 480]
    assert [level["held_out_seams"] for level in result["levels"]] == [29, 120, 480]
    beacon_a, beacon_b = result["base_beacons"]
    for level in result["levels"]:
        count = 4**level["level"]
        relay = level["relay"]
        assert relay["cells"][0] == count * (beacon_a + 1) - 1
        assert relay["cells"][-1] == count * (beacon_b + 1) - 1
        assert relay["great_circle_distance"] <= relay["atlas_route_distance"] + 1e-13
        np.testing.assert_allclose(
            [relay["receiver_probability"], relay["return_probability"]],
            1.0, atol=1e-13, rtol=0,
        )
        np.testing.assert_allclose(
            relay["relay_distance"], relay["atlas_route_distance"], atol=1e-13, rtol=0
        )
        assert level["maximum_swap_probability_error"] < 1e-13
        assert 0 < level["detuned_receiver_probability"] < 0.99
        np.testing.assert_allclose(
            np.asarray(level["cell_centers"])[count - 1::count],
            result["levels"][0]["cell_centers"], atol=1e-13, rtol=0
        )
    for obstruction in result["uniform_hop_obstructions"]:
        assert [(row["coarse_hops"], row["fine_hops"]) for row in obstruction] == [(1, 3), (2, 4)]
        assert all(row["checked_pairs"] > 0 for row in obstruction)
