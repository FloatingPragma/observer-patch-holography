"""Checks for the ordinary thermometer and signed fixed-volume calibration."""

import numpy as np

from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture import (
    sector_ball_meter,
    sector_ball_source,
    sector_boundary_clock,
    sector_loaded_ball,
)


def test_three_native_readings_resolve_the_source_spectrum():
    receipt = sector_ball_meter.calibrate()
    assert (
        len(receipt["spectral_classes"]) == 9 and len(receipt["contact_triples"]) == 15
    )
    assert all(n > 0 for n in receipt["ambiguous_two_contact_readings"].values())
    assert sum(row["word_count"] for row in receipt["contact_triples"]) == 448


def test_occupied_home_rank_and_signed_double_transfer():
    words = sector_ball_source.thermal.apparatus()[0]
    board, homes, seams, *_ = sector_ball_source.balls.energy.returned(words)
    graph, _, _ = sector_ball_source.balls.radar.chart(board, seams)
    d, g, _ = sector_loaded_ball.routes(graph, words, homes, 13, 13)
    edge = graph.toarray()
    area = []
    for volume in (6, 11, 13):
        inside = np.argsort(np.argsort(d, kind="stable"), axis=1) < volume
        area.append(np.array([edge[z][:, ~z].sum() for z in inside])[g])
    empty = words[:, 1] == 0
    np.testing.assert_array_equal(area[0], np.where(empty, 10, 14))
    assert np.ptp(area[1]) == np.ptp(area[2]) == 0
    before, after, _ = sector_boundary_clock.move_double(graph, homes, words, words[18])
    assert area[0][after] - area[0][before] == 4


def test_terminal_contrast_normalization_and_positive_empty_gain():
    arms, *_ = sector_ball_source.sources((8,))
    words = sector_ball_source.thermal.apparatus()[0]
    d, err = arms[8]
    assert abs(sum(d)) < err < 1e-7
    assert abs(d @ (words[:, 1] == 0) - 0.01139636) < 1e-7
