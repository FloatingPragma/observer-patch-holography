import json

import numpy as np

from . import clock_scores, force_forecast


def test_frozen_force_forecasts_cover_the_whole_head_algebra():
    saved = json.loads((clock_scores.HERE / "clock_scores.json").read_text())
    rebuilt = clock_scores.run()
    assert saved["event_order"] == rebuilt["event_order"]
    for charge, score in rebuilt["whole_algebra_absolute_force"]["signed_force"].items():
        held = saved["whole_algebra_absolute_force"]["signed_force"][charge]
        assert score["weight_numerators"] == held["weight_numerators"]
        assert score["intercept_numerator"] == held["intercept_numerator"]
        assert score["absolute_error_all_normalized_heads_upper"] < 6e-14
    assert rebuilt["entropy_differential"]["target_error_per_unit_HS_change_upper"] < 2e-5


def test_original_denominator_clock_counts_and_independent_material_readings():
    saved = json.loads((force_forecast.HERE / "force_forecast.json").read_text())
    rebuilt = force_forecast.run()
    for name, arm in rebuilt["preparations"].items():
        counts = np.array(arm["kept_bin_counts"])
        for route, total in arm["all_other_counts"].items():
            pair = [int(v) for v in route.split("-")]
            mask = [record[:2] == pair for record in rebuilt["event_order"]]
            assert counts[mask].sum() + total == rebuilt["synthetic_draws_per_preparation_per_route"]
        for charge, force in arm["force"].items():
            observed = force["force_from_separate_material_operator"]
            assert abs(force["forecast_from_population_clock_probabilities"] - observed) < 1e-13
            assert abs(force["forecast_from_simulated_clock_counts"] - observed) < force[
                "simultaneous_99pct_half_width_including_calibration"
            ]
            archived = saved["preparations"][name]["force"][charge]["force_from_separate_material_operator"]
            assert abs(observed - archived) < 1e-13
