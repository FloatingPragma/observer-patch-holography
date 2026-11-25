import json
from math import cos, pi, sin, sqrt
from pathlib import Path
import sys

from scipy.integrate import quad

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from cosmology.rar_deep_regime.double_ring_pressure import elliptical as e


def test_optical_apertures_in_mass_ellipse_coordinates():
    # Integrate the unaveraged convergence with the area-preserving change
    # x = s cos(phi)/sqrt(q), y = s sin(phi)*sqrt(q), out to the actual circle.
    c, eta, b = 3778.96, 1.027, 1.397
    for q in (0.4, 0.946, 1.0):

        def disk_mean(radius):
            def mass_per_angle(phi):
                edge = radius / sqrt(cos(phi) ** 2 / q + q * sin(phi) ** 2)
                return quad(lambda s: (2 - eta) / 2 * b**eta * s ** (1 - eta), 0, edge)[
                    0
                ]

            return c * 4 * quad(mass_per_angle, 0, pi / 2)[0] / (pi * radius**2)

        expected = disk_mean(1.43) - disk_mean(2.07)
        assert abs(e.circular_contrast(c, q, eta, b, 1.43, 2.07) - expected) < 2e-9


def test_unfavourable_radius_pairs_and_directed_bound():
    result = e.run()
    assert result == json.loads((e.HERE / "elliptical.json").read_text())
    c1 = e.geometry()[1][0]
    for name, panel in result["panels"].items():
        for part in panel["radius_partitions"]:
            inner, outer = part["ring_ranges_arcsec"]
            for r1 in inner:
                for r2 in outer:
                    for q in e.FIT["axis_ratio"]:
                        for eta in e.FIT["projected_density_slope"]:
                            for b in e.FIT["intermediate_axis_arcsec"]:
                                assert (
                                    e.circular_contrast(c1, q, eta, b, r1, r2)
                                    >= part["optical_contrast_lower"]
                                )
        assert panel["optical_minus_orbital_lower"] == (
            149 if name == "ring_errors_1" else 61
        )


def test_geometry_box_support_and_physical_sheet_cancellation():
    c1, c2 = e.geometry()[1]
    for row in e.run()["panels"].values():
        eps = row["aligned_geometry_threshold_percent_at_zero_secondary_bend"] / 100
        support = max(c2 * e2 - c1 * e1 for e1 in (-eps, eps) for e2 in (-eps, eps))
        assert abs(support - row["aligned_margin"]) < 4e-5
    for sheet in (-300, 500):
        assert abs(c2 * (sheet / c2) - c1 * (sheet / c1)) < 1e-10
