import json

import numpy as np

from . import geometry, measure, meter


def test_native_compensated_targets_and_failed_calibrations():
    data = meter.two_target_calibration()
    assert len(data["inputs"]) == 7
    assert min(data["correct_probabilities"]) > 1 - 3e-13
    control = meter.two_target_calibration(omit_y=True)
    for i, q in enumerate(data["charges"]):
        assert abs(control["correct_probabilities"][i] - (0.5 if q else 1)) < 3e-13


def test_ternary_input_reaches_third_service_and_restores_requested_sign():
    r = json.loads(measure.RECEIPT.read_text())
    p = r["profiles"][0]
    # Contaminate the first and second meter, then independently put an excess
    # inside the region and vacancy on the first meter.
    a, b = p["sites"][:2]
    x = np.array([a["cross"], 0, a["park"]])
    y = np.array([b["park"], a["cross"], a["cross"]])
    bins, q = geometry.labels(x, y, p, r["pointer_decoder"])
    assert list(q) == [0, 1, 0]
    assert bins[0] // 17 >= 9 and bins[1] // 17 >= 1
    assert [r["pointer_decoder"][int(i % 17)] for i in bins] == list(q)


def test_complete_counts_and_operator_bounds():
    r = json.loads(measure.RECEIPT.read_text())
    assert r["total_output_bins"] == 1754 and len(r["profiles"]) == 60
    cmp, fallback = geometry.event_masks(r["pointer_decoder"])
    for source in r["launches"]:
        for panels in source["panels"].values():
            for row in panels:
                counts = dict(row["counts"])
                assert sum(counts.values()) == r["shots"]
                assert all(
                    i < geometry.EXHAUSTION or i == geometry.FILTER_MISS for i in counts
                )
                assert sum(v for i, v in counts.items() if cmp[i]) == row["charged"]
                assert (
                    sum(v for i, v in counts.items() if fallback[i]) == row["fallback"]
                )
                assert row["accepted"] + counts[geometry.FILTER_MISS] == r["shots"]
    assert r["tests"]["below_native_002"] and r["tests"]["controls_above_20"]
    universal = r["universal_distinct_owner_bound"]
    assert universal["upper"] < 0.001071 and all(
        p["charge_on_entire_cluster"]["lower"] <= universal["upper"]
        for p in r["profiles"]
    )
    # Small independent subspace whose known eigenvector saturates the maximum.
    B = np.array([[1, 0], [0, 0.6], [0, 0.8]])
    v = measure.operator_bound(B, np.array([False, True, False]), 0.0)
    assert v["lower"] <= 0.36 <= v["upper"]
