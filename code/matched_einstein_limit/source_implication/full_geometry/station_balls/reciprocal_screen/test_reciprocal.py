import json

import numpy as np
from sympy import isprime

from . import dwell, fixed, measure


def test_transpose_instrument_and_rejected_unmatched_state():
    # Arbitrary complex coin blocks, including states with a Y component.
    rng = np.random.default_rng(14897)
    raw = rng.normal(size=(2, 4, 4)) + 1j * rng.normal(size=(2, 4, 4))
    x = raw @ raw.conj().transpose(0, 2, 1)
    x = x / np.trace(x, axis1=1, axis2=2).sum()
    k = raw[:, None, None, :, :]
    p = measure.transmission(k, x, x.transpose(0, 2, 1))
    reverse = measure.transmission(k.swapaxes(-1, -2), x, x.transpose(0, 2, 1))
    np.testing.assert_allclose(p, reverse, atol=1e-14, rtol=0)
    plain = np.broadcast_to(np.eye(4), x.shape)
    assert (
        abs(
            measure.transmission(k, x, plain)
            - measure.transmission(k.swapaxes(-1, -2), x, plain)
        ).item()
        > 0.1
    )


def test_witnesses_and_series_control():
    assert isprime(fixed.PRIME)
    for name, fn in (("fixed", fixed.certificate), ("dwell", dwell.witness)):
        assert fn() == json.loads((measure.HERE / (name + ".json")).read_text())
    # A rectangular matrix with a pivot at degree six and a genuine delayed
    # pivot at degree nine after subtracting an identical leading row.
    jets = [np.zeros((3, 2), dtype="int64") for _ in range(4)]
    jets[0][:2] = [[1, 1], [1, 1]]
    jets[3][1, 1] = 1
    assert dwell.stages(jets) == ([[6, 1], [7, 0], [8, 0], [9, 1]], [0, 1])


def test_actual_both_direction_counts_and_nonzero_gain():
    saved = json.loads((measure.HERE / "results.json").read_text())
    fresh = measure.run()
    lo, hi = saved["zero_flight_sample_interval"]
    assert lo < fresh["zero_flight_hit"] < hi
    for q in ("-1", "1"):
        for method in ("matched", "held_recalibrated", "plain", "held_fixed"):
            a = saved["charges"][q][method]
            b = fresh["charges"][q][method]
            if method in ("matched", "held_recalibrated"):
                assert b["maximum_pair_asymmetry"] < 5e-14
                np.testing.assert_allclose(
                    3 * b["unequal_patch_forward"],
                    5 * b["unequal_patch_reverse"],
                    atol=1e-13,
                )
            else:
                assert (
                    b["maximum_pair_asymmetry"] >= a["sampled_pair_asymmetry_lower"] > 0
                )
            for side in ("forward", "reverse"):
                obs = np.array(a[side]["counts"])
                assert (
                    obs.shape == (5, 6)
                    and np.all(obs >= 0)
                    and np.all(obs.sum(1) == saved["starts_per_launch"])
                )
                lo, hi = a[side]["sampled_area_gain_interval"]
                assert lo < b[side]["determinant"] < hi
                if method == "matched":
                    assert lo > 0
                    wide_lo, wide_hi = a[side]["equilibrium_area_gain_interval"]
                    assert 0 < wide_lo < lo < hi < wide_hi


def test_fixed_lengths_affine_span_and_modular_work():
    from . import fixed_volume
    from ..geometry import ANCHORS, raw_bounds

    saved = json.loads((measure.HERE / "fixed_volume.json").read_text())
    assert fixed_volume.certificate() == saved
    assert saved["frozen_modular_coefficient"][1] < -236
    base = raw_bounds()[:4, :4]
    # Test the native four-clock convention independently on point-valued means.
    point = np.array([[float(v.mid) for v in row] for row in base])

    def calibrated(values):
        log = np.zeros((4, 4))
        for i in ANCHORS:
            for j in ANCHORS:
                if i != j:
                    log[i, j] = np.log(values[i, j] / values[j, i])
        gain = log.mean(axis=1)
        return gain, values * values.T * np.exp(-gain[:, None] - gain[None, :])

    w = np.zeros((4, 4))
    for (i, j), a in zip(
        fixed_volume.combinations(range(4), 2), fixed_volume.UPPER, strict=True
    ):
        w[i, j] = a
        w[j, i] = -a
    gain, sq = calibrated(point)
    for shift in (w, np.array([1, 1, -1, -1])[:, None]):
        g2, s2 = calibrated(point * np.exp(0.01 * shift))
        np.testing.assert_allclose(s2, sq, atol=1e-13)
        if shift is w:
            np.testing.assert_allclose(g2, gain, atol=1e-14)
