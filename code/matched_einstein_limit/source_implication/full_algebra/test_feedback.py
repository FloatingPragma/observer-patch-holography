import json
from pathlib import Path

import numpy as np
import pytest

from . import feedback


def test_tracial_ensemble_selection_and_held_root():
    survey = feedback.station.Survey(feedback.station.radar_menus()[6])
    gain, volume = feedback.calibration()
    x, price, replay = feedback.select(survey, gain, volume)
    # The held continuation enters only after selecting the source and multiplier.
    held = json.loads(Path(feedback.__file__).with_name("results.json").read_text())
    y = np.array(held["x_real"])+1j*np.array(held["x_imag"])
    assert np.abs(np.linalg.eigvalsh(x-y)).sum() < 1e-9
    assert abs(held["price"]-price) < 3e-8
    assert replay["last_update_trace_norm"] < 3e-13
    assert abs(replay["volume_difference"]) < 2e-12
    with pytest.raises(ValueError, match="outside"):
        feedback.fixed_price(survey, gain, feedback.PRICE[1]+.01)


def test_physical_bracket_and_atom_reachability():
    result = feedback.run()
    lo, hi = result["physical_price_bracket"]
    assert feedback.PRICE[0] < lo < result["feedback"]["price"] < hi < feedback.PRICE[1]
    assert result["true_state_trace_distance_upper"] < 5e-5
    assert result["volume_parity_gradient_lower"] > .1
    prep = result["native_atoms"]
    assert len(prep["kept"])+len(prep["dephased"]) == prep["count"] == 72
    assert min(prep["corrected_minimum_weights_lower"].values()) > .003
    saved = json.loads(Path(feedback.__file__).with_name("feedback.json").read_text())
    assert result["physical_price_bracket"] == saved["physical_price_bracket"]
    assert result["native_odd_weight"] == pytest.approx(saved["native_odd_weight"], abs=2e-10, rel=0)
