import numpy as np
import pytest

from . import measure


@pytest.fixture(scope="module")
def data():
    return measure.run()


def test_rest_and_flight_coincidence(data):
    phase = np.sin(measure.TIMES / 2) ** 2
    np.testing.assert_allclose(data["blank_dimer"], phase, atol=2e-14)
    np.testing.assert_allclose(data["parked_head"], phase, atol=2e-14)
    joint = np.array(data["mobile_joint_same_flip_by_home_away"])
    np.testing.assert_allclose(joint[:, 1, 1], phase / 8, atol=2e-14)
    np.testing.assert_allclose(joint[:, 1, 0], 5 * phase / 8, atol=2e-14)
    np.testing.assert_allclose(joint.sum(axis=(1, 2)), 1, atol=2e-14)


def test_certified_source_and_negative_controls(data):
    p = data["left_odd_center"]
    k = data["orbital"][8]
    assert np.isclose(k["left_defect"], 3 * p * (1 - p) / 8)
    assert np.isclose(k["left_arbitrary_phase"], 3 * p * (1 - p) / 16)
    assert (
        data["independent_simulated_endpoint_counts"]["arbitrary_phase_99pct_lower"]
        > 0.02
    )
    assert abs(k["right_arbitrary_phase"]) < 2e-14
    assert abs(data["orbital"][0]["left_arbitrary_phase"]) < 2e-14


def test_single_line_with_free_phase_offset_and_contrast():
    for rate in (0.02, 0.6, 1.0, 8.0, 9.0, 54.17):
        line = 0.45 + 0.4 * np.cos(rate * measure.TIMES + 0.341)
        assert abs(measure.arbitrary_phase_defect(line)) < 1e-13
    assert abs(measure.arbitrary_phase_defect(np.ones(5))) < 1e-13


def test_directed_rail_and_timing(data):
    j = 17 * data["left_coin_X_center"] / 8
    for r in data["directed"]:
        assert np.isclose(r["endpoint_rate"], r["orientation"] * j, atol=2e-9)
        np.testing.assert_allclose(
            r["all_flip"], 0.75 * np.array(data["blank_dimer"]), atol=2e-14
        )


def test_joint_times_and_conditional_arrival_near_release():
    source, *_ = measure.density_center()
    left = measure.marginal_blocks(source, "left")
    t = np.array([1e-4, 0.063, 0.413, 1.817])
    joint = measure.mixed_joint(left, t)
    far = (np.sin(4 * t) ** 2 + np.sin(4.5 * t) ** 2) / 4
    np.testing.assert_allclose(joint[:, :, 1].sum(axis=1), far, atol=2e-14, rtol=0)
    np.testing.assert_allclose(
        joint[:, 1, 1], np.sin(t / 2) ** 2 / 8, atol=2e-14, rtol=0
    )
    assert abs(joint[0, 1, 1] / far[0] - 1 / 290) < 1e-9


def test_clock_motion_spectral_control(data):
    levels = data["exact_level_multiplicities"]
    assert levels["native"] == {"0": 2, "8": 3, "9": 3}
    assert levels["two_independent_blank_dimers"] == {"16": 1, "17": 2, "18": 1}
