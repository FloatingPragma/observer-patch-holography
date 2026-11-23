"""Check estimates and independent forecasts directly from the archived observer tapes."""

import json
from pathlib import Path

import numpy as np

from .observer_field import calibrate, forecast, held_errors


def test_probe_tapes_reproduce_frozen_equations_and_future_records():
    data = json.loads(Path(__file__).with_name("observer_field.json").read_text())
    for row in data["refinements"]:
        fitted = calibrate(row["calibration_tape"])
        np.testing.assert_allclose(list(fitted.values()), list(row["calibrated"].values()), atol=1e-12)
        for arm in row["arms"].values():
            errors = held_errors(arm["held_tapes"], fitted)
            for key, value in errors.items():
                np.testing.assert_allclose(value, arm[key], atol=1e-12)
            for start in (0, 3):
                later = forecast(arm["held_tapes"][start:start+3], fitted)
                for i, estimates in enumerate(later, start=2*(start//3)):
                    archived = arm["record_forecasts"][i]
                    np.testing.assert_allclose(estimates["arm_errors"], archived["arm_errors"], atol=1e-11)
                    np.testing.assert_allclose(estimates["scalar_error"], archived["scalar_error"], atol=1e-11)
    fine = data["refinements"][-1]["arms"]
    assert fine["coupled"]["max_einstein_component_error"] < 0.005
    assert fine["test_field_on_vacuum_geometry"]["max_einstein_component_error"] > 0.05
