import json
from pathlib import Path

from . import curved


def test_six_station_geometry_and_volume_interventions():
    fresh = curved.run()
    saved = json.loads(Path(curved.__file__).with_name("curved.json").read_text())
    assert fresh["directed_mean_range_bounds"] == saved["directed_mean_range_bounds"]
    assert fresh["positive_space_form"] == saved["positive_space_form"]
    assert fresh["negative_space_form"]["large_inverse_radius_min"] == 18
    assert fresh["negative_space_form"]["certified_subintervals"] > 0
    assert (
        fresh["energy_entropy_volume_jacobian"]
        == saved["energy_entropy_volume_jacobian"]
    )
    readings = fresh["volume_interventions"]["numeric_center"]
    assert readings[0]["volume_change"] < 0 < readings[2]["volume_change"]
    assert abs(readings[1]["volume_change"]) < 1e-14
    assert (
        max(row["entropy_change"] for row in readings)
        - min(row["entropy_change"] for row in readings)
        < 2e-14
    )
