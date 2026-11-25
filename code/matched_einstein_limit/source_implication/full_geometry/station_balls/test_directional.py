import json

import numpy as np

from . import directional as d
from . import port_apertures
from .geometry import HERE


def test_complete_rotating_frame_and_missing_hemisphere_rank():
    n = d.frame()
    matrix = np.array([[1.0, 0.3, -0.2], [0.3, -0.5, 0.7], [-0.2, 0.7, 0.1]])
    rotation = np.array([[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]])
    observations = np.einsum("ij,jk,ik->i", n, matrix, n)
    np.testing.assert_allclose(d.tensor(n, observations), matrix, atol=1e-14)
    np.testing.assert_allclose(
        d.tensor(n @ rotation, observations), rotation.T @ matrix @ rotation, atol=1e-14
    )
    # A hemisphere cap (cos half angle = 0) erases the traceless quadratic.
    axis = n[0]
    hemi, weights = d.solid_directions(axis, 0.0, 8)
    np.testing.assert_allclose(
        hemi.T @ (weights[:, None] * hemi), np.eye(3) / 3, atol=1e-14
    )
    cap, weights = d.solid_directions(axis, 0.5, 8)
    np.testing.assert_allclose(weights @ cap, d.FIRST * axis, atol=1e-14)
    np.testing.assert_allclose(
        cap.T @ (weights[:, None] * cap),
        d.ALPHA * np.eye(3) + d.BETA * np.outer(axis, axis),
        atol=1e-14,
    )


def test_smooth_quantitative_convergence_and_quadrature():
    saved = json.loads((HERE / "directional.json").read_text())
    result = d.report()
    # Both quadratures and the affine inversion use floating-point math.
    for old, new in zip(
        saved["prescribed_smooth_comparisons"],
        result["prescribed_smooth_comparisons"],
        strict=True,
    ):
        for k in old:
            np.testing.assert_allclose(old[k], new[k], atol=3e-11, rtol=3e-12)
    for field in ("spatial_error", "null_even_error", "null_odd_error"):
        errors = [row[field] for row in result["prescribed_smooth_comparisons"]]
        assert all(3.8 < a / b < 4.8 for a, b in zip(errors, errors[1:]))
    n = d.frame()
    q = np.array(result["null_even_target"])
    p = np.array(result["null_odd_target"])
    fine = d.screen_areas(n, q, p, 0.25, order=12)
    np.testing.assert_allclose(
        fine,
        result["prescribed_smooth_comparisons"][1]["paired_screen_areas"],
        atol=1e-12,
        rtol=0.0,
    )
    # Flat spacetime gives twelve exact null readings, including the q=0 branch.
    flat = d.screen_areas(n, np.zeros((3, 3)), np.zeros(3), 0.25)
    np.testing.assert_allclose(flat, np.pi * 0.25**2, atol=1e-15, rtol=0.0)


def test_native_angular_panel_replay():
    saved = json.loads((HERE / "port_apertures.json").read_text())
    fresh = port_apertures.acquire()
    assert fresh == saved
    assert fresh["ordinary_rank"] == 4 and fresh["pole_weighted_rank"] == 9
    assert sum(fresh["counts_including_outside"]) == fresh["starts"]
    assert fresh["maximum_odd_amplitude"] > fresh["simultaneous_error"]
    assert fresh["even_panel_spread"] > 2 * fresh["simultaneous_error"]


def test_joint_signed_contraction_bound():
    axes = d.frame()
    rng = np.random.default_rng(96591)
    for _ in range(30):
        n = rng.standard_normal(3)
        n /= np.linalg.norm(n)
        z = axes @ n
        b = 10 * z**2 / 3 - 17 / 18
        coefficient = np.column_stack([b / 2 + z / 3, b / 2 - z / 3])
        np.testing.assert_allclose(np.sum(coefficient**2), 377 / 108, atol=1e-14)
        error = np.sign(coefficient)
        fake_screen = np.pi * (1 - error / 6)
        even, odd = d.null(axes, fake_screen, 1.0)
        contraction = n @ even @ n + 2 * n @ odd
        np.testing.assert_allclose(contraction, np.sum(abs(coefficient)), atol=2e-14)
        assert contraction <= np.sqrt(377) / 3
