import json
from pathlib import Path

import numpy as np
import pytest

from . import live_sphere as sphere


def round_fixture():
    # Five regular four-simplex vertices and three independent points in S³.
    _, q = np.linalg.eigh(np.eye(5) - np.ones((5, 5)) / 5)
    vertices = q[:, 1:] * np.sqrt(5 / 4)
    held = np.array([[1., 2., -2., 1.], [2., -1., 3., 2.], [-1., 2., 0., 1.]])
    points = np.vstack([vertices, held / np.linalg.norm(held, axis=1)[:, None]])
    d = 1.6 * np.arccos(np.clip(points @ points.T, -1, 1))
    np.fill_diagonal(d, 0.)
    return d


def test_round_calibration_predicts_an_independent_held_triangle():
    d = round_fixture()
    out = sphere.sphere_fit(d)
    assert abs(out["curvature"] - 1 / 1.6**2) < 1e-13
    assert out["all_distances_on_shortest_branch"]
    assert out["maximum_held_diagonal_residual"] < 1e-13
    assert out["maximum_held_distance_error"] < 1e-13
    assert not out["all_radii"]["excluded"]
    change = np.zeros_like(d)
    change[1, 4] = change[4, 1] = .37
    h = 1e-5
    fd = (sphere.sphere_fit(d+h*change)["curvature"] - sphere.sphere_fit(d-h*change)["curvature"])/(2*h)
    assert abs(sphere.sphere_curvature_slope(d, change, out["curvature"]) - fd) < 1e-9
    assert abs(sphere.sphere_curvature_slope(d, d, out["curvature"]) + 2*out["curvature"]) < 1e-13


def test_six_nearly_equal_distances_exclude_every_radius():
    d = np.ones((8, 8)) - np.eye(8)
    out = sphere.sphere_fit(d)
    check = out["all_radii"]
    assert abs(out["curvature"] - np.arccos(-.25)**2) < 1e-13
    assert check["excluded"] and len(check["witness_indices"]) == 6
    assert .1010 < check["maximum_distance_change_lower"] < .1011
    assert out["maximum_held_distance_error"] > .99
    # All anchors prescribe the fifth vertex for every held point; each is
    # already unit length. Only comparisons outside the anchors expose this.
    assert out["maximum_held_diagonal_residual"] < 1e-13
    d[2, 3] = .8
    with pytest.raises(ValueError, match="symmetric"):
        sphere.all_radii_obstruction(d)


def test_area_response_integrates_the_independent_volume():
    from scipy.integrate import quad

    k, ell, step = .4, .5, 1e-5
    result = sphere.sphere_fixed_volume(k, ell)
    integral = quad(lambda r: 4*np.pi*np.sin(np.sqrt(k)*r)**2/k, 0, result["proper_radius"])[0]
    assert abs(integral - 4*np.pi*ell**3/3) < 1e-14
    fd = (sphere.sphere_fixed_volume(k+step, ell)["area"] - sphere.sphere_fixed_volume(k-step, ell)["area"])/(2*step)
    assert abs(result["area_per_curvature_at_fixed_volume"] - fd) < 1e-9
    assert abs(sphere.sphere_fixed_volume(k, 1e-3)["area_per_curvature_at_fixed_volume"] / (1e-3)**4
               + 4*np.pi/5) < 1e-6


def test_eight_live_port_ranges_fail_the_held_spherical_test():
    survey = json.loads(Path(__file__).with_name("live_survey.json").read_text())
    assert survey["stations"] == [0, 1, 38, 39, 5, 7, 40, 42]
    panel, = (row for row in survey["panels"]
              if row["coupled"] and row["deadline"] == 12 and row["read"] == "pooled")
    d = np.array(panel["lengths"])
    result = sphere.sphere_fit(d)
    assert .38 < result["curvature"] < .39 and result["all_distances_on_shortest_branch"]
    assert result["all_radii"]["maximum_distance_change_lower"] > .30
    assert result["maximum_held_distance_error"] > 2.99
    slopes = {name: sphere.sphere_curvature_slope(d, delta, result["curvature"])
              for name, delta in panel["length_slopes"].items()}
    assert slopes["ordinary"] < -1.2e-4 and slopes["decorrelated"] > 4.8e-5
