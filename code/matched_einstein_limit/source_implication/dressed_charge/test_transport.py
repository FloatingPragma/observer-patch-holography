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
    assert control == stored["ordinary_calibrated_control"]
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
