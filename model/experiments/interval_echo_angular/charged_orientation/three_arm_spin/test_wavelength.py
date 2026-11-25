import json

import numpy as np
import pytest

from .arm_wavelength import HERE, read_counts


def test_local_counts_resolve_the_standing_wave_without_an_energy_input():
    receipt = json.loads((HERE / "wavelength.json").read_text())
    assert sum(receipt["original_start_counts"]) == receipt["launches"]
    data = read_counts(receipt["original_start_counts"][:4], receipt["launches"])
    assert data == receipt["inferred"]
    assert data["zero_phase_third_difference_99"][1] < -1e-5
    first, alias = data["phase_intervals_99"]
    theta = np.arccos(17 - receipt["selected_pole"])
    assert first[0] < theta < first[1]
    assert alias[0] < np.pi - theta < alias[1]
    with pytest.raises(ValueError, match="unresolved"):
        read_counts([100, 100, 100, 100], 10000)
