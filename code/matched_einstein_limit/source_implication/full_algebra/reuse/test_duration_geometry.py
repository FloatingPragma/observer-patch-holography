import json
from pathlib import Path

import numpy as np

from . import duration_geometry, metric_counts, metric_forecast


def test_all_deadlines_both_orientations_and_all_gram_corners():
    record = duration_geometry.run()
    # JSON canonicalizes tuples as arrays.
    assert json.loads(json.dumps(record)) == json.loads(
        Path(duration_geometry.__file__).with_suffix(".json").read_text()
    )
    counts = record["simulated_duration_counts"]
    for arm in counts["records"].values():
        assert all(
            sum(arm[key]) == counts["trials_per_arm_and_setting"]
            for key in ("minus_counts", "plus_counts")
        )
        assert (
            arm["smallest_signed_difference"] > counts["simultaneous_99pct_half_width"]
        )
    assert len(record["deadlines"]) == 11 and len(record["opposite_time_orders"]) == 4
    assert record["common_metric_box"]["volume_lower"] == "1.57"


def test_both_single_link_readers_predict_absolute_independent_lengths():
    record = metric_forecast.run()
    saved = json.loads(Path(metric_forecast.__file__).with_suffix(".json").read_text())
    assert record["joint_events"] == saved["joint_events"]
    assert record["absolute_readers"] == saved["absolute_readers"]
    held = record["held_population_comparisons"]
    volumes = []
    for name, data in held.items():
        np.testing.assert_allclose(
            data["separate_lengths"],
            saved["held_population_comparisons"][name]["separate_lengths"],
            rtol=0,
            atol=1e-10,
        )
        for p in data["single_edge_predictions"]:
            np.testing.assert_allclose(p, data["separate_lengths"], rtol=0, atol=2e-9)
        volumes.append(data["separate_volume"])
    assert max(volumes) - min(volumes) > 0.019


def test_single_link_count_intervals_separate_two_independent_volumes():
    result = metric_counts.run()
    held = result["records"]
    for q in (0, 1):
        low = held["dephased_Xplus_Yplus"][q]["volume_band"]
        high = held["kept_Zplus_Zminus"][q]["volume_band"]
        assert low[1] < high[0]
    for settings in held.values():
        assert all(sum(x["counts"]) == 2**72 for x in settings)
