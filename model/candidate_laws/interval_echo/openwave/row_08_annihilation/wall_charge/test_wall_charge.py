import json
from itertools import product

import numpy as np

from . import measure, pickup


def test_native_common_and_gradient_standards():
    calibration = json.loads((measure.HERE / "calibration.json").read_text())
    fresh = pickup.calibration()
    np.testing.assert_array_equal(fresh["counts"], calibration["standards"]["counts"])
    for gradient, reverse in product((0, 1), repeat=2):
        probs = np.array(fresh["probabilities"])[gradient, reverse]
        for j, (a, b) in enumerate(pickup.QWORDS):
            value = (b - a if gradient else a + b - 2) * (-1 if reverse else 1)
            if value in (-2, 0, 2):
                np.testing.assert_allclose(
                    probs[j], np.eye(17)[{-2: 10, 0: 5, 2: 14}[value]], atol=1e-13
                )
            else:
                np.testing.assert_allclose(probs[j, 5], 0.5)
                assert (probs[j, 14] > probs[j, 10]) == (value > 0)
    assert pickup.full_dimer_calibration() == calibration["full_dimer"]
    assert pickup.routing() == calibration["routes"]


def test_charge_result_covers_all_allowed_flight_states():
    _, anchors, _ = measure.wall_flight.extended_geometry()
    _, _, _, distal = measure.contact_masks(anchors)
    q, d = measure.pickups()
    for j, (a, b) in enumerate(measure.CONTACTS):
        # Read both pointer modes and reconstruct the two original ternary loads.
        np.testing.assert_array_equal(1 + (q[:, j] - d[:, j]) // 2, measure.LOADS[:, a])
        np.testing.assert_array_equal(1 + (q[:, j] + d[:, j]) // 2, measure.LOADS[:, b])
        assert np.all(q[d[:, j] != 0, j] == 0)
        complete = [2 * (a // 2), 2 * (a // 2) + 1, 2 * (b // 2), 2 * (b // 2) + 1]
        assert np.all(measure.LOADS[:, complete].sum(axis=1) == 4)
    # Independently invert the logical basis of all terminal records. Ring
    # bits plus the four remotely read contacts suffice for all ten dimers.
    words = np.column_stack([measure.BITS[:, :4], q[:, 2:], d[:, 2:]])
    assert len(np.unique(words, axis=0)) == 1024
    assert distal.sum() == 2
    np.testing.assert_array_equal(d[distal, 4], -d[distal, 5])


def test_portable_held_counts_conjugation_and_cut():
    with np.load(measure.HERE / "counts.npz") as saved:
        counts = saved["counts"]
        assert np.all(counts.sum(axis=-1) == measure.SHOTS)
        cdf = saved["cdf"]
        p = np.diff(cdf, prepend=0, axis=-1) / measure.QUANTA
        # Inverting every 0 and 2 permutes word w into 1023-w. Without
        # attachments the receiver law factorizes from both ring controls.
        np.testing.assert_allclose(
            p[:, 0, 0], p[:, 1, 1, :, :, ::-1], atol=2 / measure.QUANTA, rtol=0
        )
        cut = p[1].reshape(2, 2, 2, 41, 16, 64).sum(axis=-2)
        for r, w in product((0, 1), repeat=2):
            np.testing.assert_allclose(
                cut[r, :, w], cut[0, :, 0], atol=16 / measure.QUANTA, rtol=0
            )
        out = measure.decode(counts, saved["hand"], cdf)
    assert out["settings"]["k8_cut0_g0"]["far"]["earliest"] == 12
    assert out["settings"]["k8_cut1_g0"]["far"]["earliest"] is None
    assert out["at32"]["k8_cut0_g0_r0_wall1"]["far_gradient_counts"] == [
        0,
        646075,
        0,
        0,
    ]
    assert out["at32"]["k8_cut0_g0_r1_wall0"]["far_gradient_counts"] == [
        0,
        42894181,
        0,
        0,
    ]
