import json

import numpy as np

from . import certificate, measure


def test_capped_arrival_effect_counts_and_stored_direction_certificate():
    report = measure.run()
    saved = json.loads((measure.HERE / "results.json").read_text())
    assert len(report["panels"]) == len(saved["panels"])
    for panel, prior in zip(report["panels"], saved["panels"]):
        np.testing.assert_allclose(panel["readings"], prior["readings"], atol=2e-12)
        if panel["source"] == "unpair":
            assert not panel["gain_test"]["feasible_positive"]
        if panel["source"] == "matched_ordinary":
            assert panel["gain_test"]["best_minimum_gain"] > 0.007
        if panel["source"] == "matched_unpair":
            assert (
                not panel["gain_test"]["feasible_positive"]
                and not panel["opposite_rectangles"]
            )
        if panel["source"] == "ordinary":
            assert panel["gain_test"]["best_minimum_gain"] > 0.003
    ks = np.zeros((5, 2, 2, 4, 4), complex)
    for n in range(5):
        for i in range(2):
            ks[n, (i + n) % 2, i] = np.eye(4)
    # Exact alternating complete transfer. No branch misses the first appointment.
    out = measure.ranging.first_hits(ks, 0, 1)
    np.testing.assert_array_equal(out[1], np.eye(4))
    np.testing.assert_array_equal(out[2:], np.zeros((3, 4, 4)))


def test_rational_witnesses_and_early_transmission():
    bounded = certificate.certify()
    stored = json.loads((measure.HERE / "certificate.json").read_text())
    assert bounded == stored
    matched = bounded["matched_pairing"]
    assert matched["minimum_signed_lower"] > 7e-4
    assert matched["angle_interval"][0] < -294.13 < matched["angle_interval"][1]
    assert (
        matched["ordinary_angle_interval"][0]
        < -989.94
        < matched["ordinary_angle_interval"][1]
    )
    for label, expected in (
        ("paired_four_edge_splits", [1, 1, -1, 1]),
        ("ordinary_four_edge_splits", [1, -1, -1, 1]),
    ):
        signs = [np.sign(r["lower"]) for r in matched[label]]
        assert signs == expected
    assert len(bounded["shortest_head_projectors"]) == 8
    for row in bounded["readings"]:
        assert row["minimum_signed_lower"] > 0
        if row["source"] == "state":
            assert row["minimum_signed_lower"] > 0.198
    low, high = bounded["preprocessed_leading_unpair_slope"]
    assert -0.233 < low < high < -0.232
