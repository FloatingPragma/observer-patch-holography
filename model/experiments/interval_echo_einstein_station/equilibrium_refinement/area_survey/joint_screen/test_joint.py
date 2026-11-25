"""Native wave, physical source, image calibration, row rank and complete counts."""

import json
from fractions import Fraction
import numpy as np

from . import acquire, certificate, measure, wave


def data(name):
    return json.loads((measure.HERE / name).read_text())


def test_transverse_shape_is_an_area_and_fixed_point_matrix_products():
    q = np.array([0.016, 0.007, 0.002])
    coords = np.array([[0.0, 0.0], [1.0, 0.8], [-0.6, 1.0]])

    def direct(coords):
        rays = np.array([coords[a] - coords[b] for a, b in ((0, 1), (0, 2), (1, 2))])
        m = (rays.T * q) @ rays
        return np.sqrt(3 * np.linalg.det(m)) / 2

    A = abs(np.linalg.det(coords[1:] - coords[0])) / 2
    np.testing.assert_allclose(direct(coords), A * measure.shape3(q[None])[0])
    frame = np.array([[0.4, -1.0], [0.7, 1.3]])
    np.testing.assert_allclose(
        direct(coords @ frame + 11), abs(np.linalg.det(frame)) * direct(coords)
    )
    np.testing.assert_allclose(measure.shape3(np.ones((1, 3)) / 3), [1])
    assert measure.shape3(np.array([[1.0, 0, 0]]))[0] == 0
    rng = np.random.default_rng(41)
    a = rng.integers(-(2**50), 2**50, (4, 5), dtype=np.int64).astype(object) * 2**15
    b = rng.integers(-(2**52), 2**52, (5, 4), dtype=np.int64).astype(object) * 2**16
    # The integer kernel multiplies square propagation matrices; this sample
    # exercises signs, carry chains and the complete digit range.
    a = a[:, :4]
    b = b[:4]
    np.testing.assert_array_equal(wave.product(wave.split(a), wave.split(b)), a @ b)


def test_joint_trajectory_recurrence_and_sham_are_distinct():
    s = measure.clock.survey(1)
    short = measure.instrument(s, 2)
    x = short.support.copy()
    expected = np.zeros_like(x)
    for t in range(1, 65):
        x = s.transition @ x
        expected += x * 2.0 ** (-t)
    np.testing.assert_allclose(short.terminal, expected, rtol=0, atol=2e-14)
    # At zero flight the apex and two distinct opposite pixels cannot all be occupied by two residents.
    sh = measure.Screens(
        s.occupied, short.support, short.support, short.rim_areas, s.starts, s.scores
    )
    h, true, sham = sh.face_probs(s.uniform)
    assert np.all(true == 0) and np.all(sham > 0) and np.all(h > 0)
    released = measure.run()
    r = released["movements"]
    assert released["ordinary_gain"][4] > 0
    assert r["zero_to_two"]["image_slopes"][4] < -34000
    assert r["matched_energy_and_cut_load"]["calibrated_residual"][4] < -8.8


def test_independent_native_integer_rank_and_finite_entropy_gap():
    cert = certificate.run()
    arch = data("certificate.json")
    for r in (cert, arch):
        assert r["rank"] == 36 and r["right_inverse_defect_upper"] < 0.00003
        assert r["reachable_halfwidth_lower"] > 1 / 65536
        assert r["correction_to_exact_joint_fiber_linf_upper"] < 1e-8
        assert r["certified_finite_entropy_loss_lower_nats"] > 0.0046
        assert r["numerical_event_probability_error_upper"] < acquire.TERMINAL_BIAS
        assert r["numerical_old_clock_error_upper"] < acquire.CLOCK_BIAS
        assert r["new_ticket_largest_joint_pixel_difference_upper"] < 1e-13
        assert (
            abs(
                float(Fraction(r["new_ticket_rim_two_probability_difference"]))
                + 1 / 65536
            )
            < 1e-12
        )
        assert r["old_clock_lottery_image_changes"][4][0] > 2.38
        assert r["old_clock_lottery_image_changes"][6][1] < -5.95


def test_independent_cohorts_holdouts_and_full_outcome_accounting():
    for r in (acquire.run(), data("counts.json")):
        for stock in r["counts"].values():
            assert sum(stock["bulk_word_counts"]) == acquire.SHOTS
            assert all(sum(c) == acquire.SHOTS for c in stock["complete_rim_clocks"])
            for face in stock["joint_faces"]:
                assert (
                    len(face["two_resident"]) == 8
                    and sum(face["two_resident"]) == acquire.SHOTS
                )
                assert (
                    len(face["two_independent"]) == 5
                    and sum(face["two_independent"]) == acquire.SHOTS
                )
        assert r["area_gains"][4][0] > 0
        assert r["held"]["zero_to_two"]["residuals"][4][1] < -19
        old = r["lotteries"]["old_clock_lottery"]
        new = r["lotteries"]["augmented_lottery"]
        assert old["image"][4][0] > 2 and old["image"][6][1] < -5
        assert (
            new["image"][4][0] < 0 < new["image"][4][1] and new["entropy"][1] < -0.00469
        )
        assert old["independent_capacity"][0] > 0 and new["independent_capacity"][0] > 0


def test_full_record_finite_entropy_separation_on_positive_native_tickets():
    s = measure.clock.survey(1)
    groups, tickets = measure.interior_reallocation(s)
    assert groups == {648: 30, 649: 10, 650: 1, 778: 123, 779: 26}
    assert sum(tickets) == 552 and min(tickets) == 1
    delta = tickets - 2
    assert not np.any(delta[s.occupied.any(axis=1)])
    cert = certificate.interior_control(s)
    assert cert["exact_bulk_entropy_loss_nats"][0] > 0.08910
    readout = data("counts.json")["lotteries"]["interior_reallocation"]
    assert readout["entropy"][1] < -0.08910
    assert readout["independent_capacity"][0] < 0 < readout["independent_capacity"][1]
    assert all(a < 0 < b for a, b in readout["image"])
