from pathlib import Path
import json

import numpy as np

from . import transport


def test_native_bridge_interpolation_and_raw_pairing_enclosure():
    data = transport.run()
    stored = json.loads(
        Path(transport.__file__).with_name("transport.json").read_text()
    )
    assert data["unique_two_hop_path"] == [5, 1, 38]
    for sign, wanted in (("-1", [756, 886]), ("1", [637, 254])):
        assert [s["minor_mod_1009"] for s in data["sectors"][sign]] == wanted
        assert [s["minor_mod_1009"] for s in stored["sectors"][sign]] == wanted
    control = data["ordinary_calibrated_control"]
    assert_same_tree(control, stored["ordinary_calibrated_control"])
    ref, ordinary, held = control[
        "nominal_reference_and_ordinary_and_unpairing_changes"
    ]
    assert abs(ref) < 1e-12 and abs(ordinary) < 1e-12
    low, high = control["exact_source_unpairing_current_change_enclosure"]
    assert low < held < high < 0
    assert data["nominal_raw_pairing_current_rate"] < -0.0577
    assert data["nominal_fixed_volume_pairing_current_rate"] < -0.154
    np.testing.assert_allclose(
        data["nominal_initial_current_rate"],
        stored["nominal_initial_current_rate"],
        atol=1e-14,
        rtol=0,
    )


def assert_same_tree(actual, stored):
    """Equal structure and exact non-float entries; floats agree to the platform roundoff of the rate check."""
    if isinstance(actual, dict):
        assert actual.keys() == stored.keys()
        for key in actual:
            assert_same_tree(actual[key], stored[key])
    elif isinstance(actual, list):
        assert len(actual) == len(stored)
        for a, b in zip(actual, stored):
            assert_same_tree(a, b)
    elif isinstance(actual, float) or isinstance(stored, float):
        np.testing.assert_allclose(actual, stored, atol=1e-14, rtol=0)
    else:
        assert actual == stored
