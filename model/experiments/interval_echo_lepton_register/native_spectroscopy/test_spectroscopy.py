"""Check the finite clock and independent controls on the registered source."""

from functools import cache

import numpy as np
from numpy.testing import assert_allclose

from .measure import SHOTS, clock_probabilities, run


@cache
def receipt():
    return run()


def test_fourier_pointer_matches_direct_clock_read():
    n = 64
    times = np.arange(n)
    roots = np.exp(2j * np.pi * np.outer(times, times) / n)
    frequencies = (5.0, 5.25, 23.84)
    for frequency in frequencies:
        response = np.exp(-2j * np.pi * times * frequency / n)
        clock_amplitude = roots.T @ response / n
        assert_allclose(
            clock_probabilities(response[:, None])[0],
            abs(clock_amplitude) ** 2,
            atol=2e-15,
        )


def test_held_spectroscopy_and_complete_support():
    data = receipt()
    for name in ("blank", "written_real", "face_band_real"):
        case = data["cases"][name]
        assert (
            abs(case["held_opposite_tag_contrast"])
            < case["held_opposite_tag_radius_99"]
        )
    head = data["cases"]["face_band_positive_phase"]
    assert (
        head["held_opposite_tag_contrast"] - head["held_opposite_tag_radius_99"] > 0.11
    )
    assert len(set(head["selected_bins"])) < 3
    counts = np.array(head["all_bins_at_held_receiver"]["guessed_tags_and_no_arrival"])
    assert counts.shape == (3, 4)
    np.testing.assert_array_equal(counts.sum(axis=1), [SHOTS] * 3)
    silent = data["cases"]["written_silent_reference"]["all_bins_at_held_receiver"]
    assert silent["excess_lower_99"] > 0.0004
    assert (
        silent["opposite_tag_contrast_given_arrival"]
        - silent["opposite_tag_radius_given_arrival_99"]
        > 0.02
    )
    for board in data["support"]:
        assert board["blank"]["formation_support_exact"] == ["32", "54"]
        for row in board["chiral"]:
            assert row["support"]["heaviest_to_lightest_ceiling"] <= 40
            assert row["support"]["koide_ceiling"] < 2 / 3


def test_two_blank_host_neutral_routes():
    data = receipt()
    floor = data["blank_neutral_formation"]
    assert floor["neutral_formation_floor"] == 72
    assert floor["neutral_to_lightest_negative_ratio_floor"] == "48/23"
    response = data["filled_memory_receiver"]
    assert response["filled_columns_checked"] == 11
    assert response["unfilled_matrix_elements"] == 0
    row = response["banks"][1]["probabilities"]
    assert_allclose(row[0], row[1], atol=2e-14, rtol=0)
    assert abs(row[0] - row[2]) > 0.1
