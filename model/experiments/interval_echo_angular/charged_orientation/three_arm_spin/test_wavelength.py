import json

import numpy as np
import pytest

from .arm_wavelength import HERE, read_counts


def test_local_counts_resolve_the_standing_wave_without_an_energy_input():
    receipt = json.loads((HERE / "wavelength.json").read_text())
    assert sum(receipt["original_start_counts"]) == receipt["launches"]
    data = read_counts(receipt["original_start_counts"][:4], receipt["launches"])
    for key, value in data.items():
        np.testing.assert_allclose(value, receipt["inferred"][key], rtol=1e-10)
    assert data["zero_phase_third_difference_99"][1] < -1e-5
    first, alias = data["phase_intervals_99"]
    theta = np.arccos(17 - receipt["selected_pole"])
    assert first[0] < theta < first[1]
    assert alias[0] < np.pi - theta < alias[1]
    with pytest.raises(ValueError, match="unresolved"):
        read_counts([100, 100, 100, 100], 10000)


def test_four_counters_forecast_eight_held_counters_with_both_spatial_aliases():
    from .smooth_profile import prediction, run

    data = json.loads((HERE / "smooth_profile.json").read_text())
    training, held = data["training"], data["held"]
    assert set(training["stations"]).isdisjoint(held["stations"])
    assert sum(held["counts_with_other_final"]) == held["original_starts"]
    fit = prediction(training["counts"], training["original_starts"])
    assert not fit["flat_compatible_with_training"]
    assert not fit["quadratic_compatible_with_training"]
    phase, alias = fit["phase_per_step_aliases"]
    np.testing.assert_allclose(alias, np.pi - np.array(phase[::-1]))
    bounds = np.array(fit["predicted"])
    observed = np.array(held["intervals"])
    assert np.all(bounds[:, 0] < observed[:, 0]) and np.all(
        observed[:, 1] < bounds[:, 1]
    )
    # The nearest node and antinode are resolved even at the outside limits.
    assert 30 * bounds[0, 1] < bounds[1, 0]
    replay = run()["held"]
    assert all(replay["forecast_overlaps"])
    assert replay["physical_port_groups"] == held["physical_port_groups"]
