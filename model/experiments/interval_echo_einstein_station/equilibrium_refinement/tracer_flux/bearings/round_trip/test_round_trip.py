"""Check the free Hamiltonian, actual two-leg clock, radius holdouts and distortion bound."""

from functools import lru_cache
from math import sqrt
import json
import numpy as np
from scipy.linalg import expm
from . import dynamics, measure


@lru_cache(maxsize=1)
def fresh():
    return measure.run()


def saved():
    return json.loads(measure.HERE.joinpath("results.json").read_text())


def test_free_contact_wave_and_two_leg_geometric_deadline():
    h = np.array([[3, -1, 0], [-1, 4, -1], [0, -1, 5]], dtype=int)
    u, err = dynamics.wave(h, 1 / 128)
    expected = expm(-0.5j * h / 128)
    # The integer centering phase is common to every launch and detector.
    phase = np.exp(0.5j * 4 / 128)
    assert np.max(np.linalg.norm(u - phase * expected, axis=0)) < err
    tag = np.arange(12)
    p = np.ones((12, 12)) / 22
    np.fill_diagonal(p, 0.5)
    lam = 1 / 64
    fail, err = dynamics.returns(p, tag, lam, 0, 22)
    success_at_other = (1 - lam) / 22 / (1 - (1 - lam) * 21 / 22)
    cols = np.arange(12)
    assert (
        np.max(abs(fail[cols[:, None] != cols[None, :]] - (1 - success_at_other**2)))
        < err
        < 1e-8
    )


def test_all_released_clocks_recover_angular_subspace_and_radial_classes():
    for packet in (fresh(), saved()):
        assert (
            1.95
            < packet["all_resident_states_elapsed_time_ratio"][0]
            < 2
            < packet["all_resident_states_elapsed_time_ratio"][1]
            < 2.05
        )
        assert packet["zero_contrast_relative_drift_upper"] < 0.01
        for exp in packet["exposures"]:
            assert exp["native"]["conditional_drift_relative_upper"] < 0.015
            assert exp["two_leg_deadline_bias_fraction_upper"] < 0.0016
            assert exp["expired_probability_error_upper"] < 1e-8
            for bank in exp["banks"].values():
                for counts in bank["rows_expired_completed_blocked"]:
                    assert sum(c is None for c in counts) == 1
                    assert all(
                        sum(c) == measure.SHOTS and min(c) >= 0
                        for c in counts
                        if c is not None
                    )
                for bounds, h in zip(bank["classes"], [11, 14, 15], strict=True):
                    assert (
                        h - 0.11 < bounds[0] < h + 0.11
                        and h - 0.11 < bounds[1] < h + 0.11
                    )
                a, b = bank["slow_spectrum_99_enclosure"]
                assert bank["fourth_eigenvalue_upper"] < a < 12 / (1 - 1 / sqrt(5)) < b
                assert bank["population_slow_to_bearing_angle_sine_upper"] < 0.045


def test_raw_time_flat_rank_and_curved_small_ball_control():
    for packet in (fresh(), saved()):
        for exp in packet["exposures"]:
            for bank in exp["banks"].values():
                assert bank["squared_range_fourth_gram_eigenvalue_lower"] > 40
                assert bank["squared_range_smallest_centered_eigenvalue_lower"] > 21
                assert (
                    bank["minimum_relative_distance_error_for_flat_3d"] > 0.021 > 1 / 52
                )
