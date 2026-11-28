"""Check the branch sign on a translated halo and the shear/control observations."""

import json

import numpy as np
from scipy.integrate import quad

from cosmology.rar_deep_regime.record_shear import measure


def test_off_center_harmonic_source_and_compact_core():
    # u=1/|r-d| gives surface density pi/|R-d|. For R>d the
    # cylinder integral can be evaluated in polar coordinates about the source.
    d = 0.7
    radii = np.array([1.3, 2.1, 3.0, 4.0, 6.0, 10.0])[::-1]

    def shear(r):
        enclosed = (
            quad(lambda t: np.sqrt(r * r - d * d * np.sin(t) ** 2), 0, 2 * np.pi)[0]
            / r**2
        )
        at_ring = (
            quad(
                lambda t: 1 / np.sqrt(r * r + d * d - 2 * d * r * np.cos(t)),
                0,
                2 * np.pi,
            )[0]
            / 2
        )
        return enclosed - at_ring

    x = radii**-2
    y = np.array([shear(r) for r in radii])
    cov = np.eye(len(x))
    result = measure.shape_test(x, x, x, y, cov)
    assert result["bin_centres"]["largest_chord_z"] < 0
    alternate = measure.shape_test(x, x, x, y + 8 + 100 * x, cov)
    np.testing.assert_allclose(
        [v["centre"]["chord_excess"] for v in result["triples"]],
        [v["centre"]["chord_excess"] for v in alternate["triples"]],
        atol=1e-13,
    )
    # A positive x^2 term is strictly convex and violates the proposed sign.
    assert all(
        v["centre"]["chord_excess"] > 0
        for v in measure.shape_test(x, x, x, x**2, cov)["triples"]
    )


def test_stationary_walk_and_all_published_contrasts():
    data = json.loads(measure.RECEIPT.read_text())
    measure.verify(data, measure.run(data["source_data"]))
    for row in data["cubic_walk"]:
        assert row["minimum_projected_laplacian_off_source"] > 0
        assert row["minimum_projected_laplacian_on_sources"] < 0
        assert row["minimum_after_reversing_modular_sign_off_source"] < 0
    for survey in data["measurements"].values():
        for read in survey.values():
            assert len(read["triples"]) == 20
            assert (
                read["bin_centres"]["largest_chord_z"]
                < read["simultaneous_95pct_z_threshold"]
            )
            assert (
                read["arbitrary_within_bin_weights"]["largest_chord_z"]
                <= read["bin_centres"]["largest_chord_z"]
            )
