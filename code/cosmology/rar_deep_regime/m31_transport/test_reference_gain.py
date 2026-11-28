"""The local reference change cannot produce the observed steepest shell gain."""

from pathlib import Path
import json

import numpy as np

from cosmology.rar_deep_regime.m31_transport import reference_gain


def test_all_allowed_reference_mixtures_and_paired_baryon_variants():
    result = reference_gain.run()
    archive = json.loads(
        Path(reference_gain.__file__).with_name("reference_gain.json").read_text()
    )
    assert result == archive
    bound = result["maximum_gain_ratio_with_local_ordinary_calibration"]
    assert 1.12 < bound < 1.13
    assert result["maximum_gain_ratio_with_one_overall_calibration"] < 1.01
    assert (
        result["surveys"]["Chemin"]["least_two_error_gain_over_tested_sources"] > 1.88
    )
    assert result["surveys"]["Corbelli"]["least_two_error_gain_over_tested_sources"] < 1
    for case in result["surveys"]["Chemin"]["cases"].values():
        assert all(np.diff(case["bounds_on_mean_outer_to_inner_gain"]) < 0)
        assert case["bounds_on_mean_outer_to_inner_gain"][-1] > bound
