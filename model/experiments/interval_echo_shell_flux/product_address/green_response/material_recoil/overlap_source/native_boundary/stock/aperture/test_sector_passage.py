"""Exterior source identifiability and physical slot-clock checks."""

import numpy as np

from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture import (
    sector_passage_chart as chart,
)
from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture import (
    sector_passage_inversion as inversion,
)
from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture import (
    sector_passage_tide as tide,
)


def test_relative_acceleration_includes_both_releases_on_observed_ruler():
    # At opposite ends of a three-station path either packet can hop inward.
    L = np.array([[1, -1, 0], [-1, 2, -1], [0, -1, 1]])
    R = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]], dtype=float)
    assert tide.project(L, R)[0, 2] == -1
    change = tide.response(R, 3 * R, root=True)
    np.testing.assert_allclose(change, np.sqrt(R))
    np.testing.assert_allclose(tide.project(L, change)[0, 2], 1 - np.sqrt(2))


def test_sparse_joint_frequency_and_connected_source_have_different_signs():
    source = chart.source
    words = source.thermal.apparatus()[0]
    x, y = (words[:, [0, 3]] > 0).astype(int).T
    (delta, eps) = source.sources((256,))[0][256]
    joint = delta @ (x * y)
    dx, dy = delta @ x, delta @ y
    connected = joint - x.mean() * dy - y.mean() * dx - dx * dy
    assert joint + eps < 0 < connected - 4 * eps


def test_serial_source_offsets_cancel_on_a_tour():
    adj = np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]])
    row = chart.serial(adj, np.array([6.0, 6.0, 6.0]), np.array([1.0, 2.0, 1.0]))
    assert row["relative_commute_change"] == 4 / 18
    assert row["maximum_offset_difference"] > 0


def test_joint_occupancy_is_an_independent_exterior_reading():
    _, _, ports, clearance, adj, uv, _, periods = chart.apparatus()
    ranges = []
    errors = []
    for ticks in periods.values():
        step, _, r, e, _ = chart.law(adj, uv, ticks)
        assert np.all(step.sum(1) == 12) and np.array_equal(step, step.T)
        ranges.append(r)
        errors.append(e)
    F, C, err, cert = inversion.basis(
        ports, adj, list(periods), np.array(ranges), errors
    )
    ix = np.flatnonzero(clearance >= 2)
    a, b = np.triu_indices(len(ix), 1)
    array = C[1:, ix[a], ix[b]].T
    picks = inversion.baselines(array)
    assert np.linalg.norm(np.linalg.inv(array[picks]), ord=np.inf) * 5 * err < 1
    assert C[-1, ix[a], ix[b]].min() < -err and C[-1, ix[a], ix[b]].max() > err
    assert cert["held_classes"] == 3 and F.shape == (9, 6)
