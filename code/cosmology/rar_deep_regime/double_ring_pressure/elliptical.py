"""Circular-aperture mass contrast of the two-source HST elliptical lens fit."""

import argparse
import json
from math import cos, floor, pi, sin

from mpmath import iv
from scipy.integrate import quad

from .measure import HERE, THETA, critical_density_bounds, decimal, geometry

# Median and marginal 68 percent ranges in Collett & Auger (2014), Table 1.
SOURCE = "https://doi.org/10.1093/mnras/stu1190"
FIT = {
    "axis_ratio": (0.946, 0.941, 0.955),
    "projected_density_slope": (1.027, 1.002, 1.050),
    "intermediate_axis_arcsec": (1.397, 1.396, 1.398),
}


def circular_contrast(critical, q, eta, b, theta1, theta2):
    def angular_factor(phi):
        return (q * cos(phi) ** 2 + sin(phi) ** 2 / q) ** (-eta / 2)

    average = 2 / pi * quad(angular_factor, 0, pi / 2)[0]
    return critical * average * ((b / theta1) ** eta - (b / theta2) ** eta)


def contrast_lower(critical, inner_max, outer_min):
    q_min = decimal(FIT["axis_ratio"][1])
    b = decimal(FIT["intermediate_axis_arcsec"][1])
    _, low, high = FIT["projected_density_slope"]
    eta_min, eta_max = decimal(low), decimal(high)
    # Jensen over the azimuth. On 0<q<=1, q+1/q decreases with q.
    average = ((q_min + 1 / q_min) / 2) ** (-eta_max / 2)
    inner, outer = decimal(inner_max), decimal(outer_min)
    # For these worst-case radii, t*exp(-eta*t) increases between the
    # two endpoints t=log(theta/b), so the difference increases with eta.
    assert b < inner < outer and eta_max * iv.ln(outer / b) < 1
    difference = (b / inner) ** eta_min - (b / outer) ** eta_min
    return floor(float((critical * average * difference).a))


def compare_box(critical, ring_ranges, orbital_box):
    inner, outer = ring_ranges
    certified_inner, certified_outer = orbital_box["ring_ranges_arcsec"]
    assert certified_inner[0] <= inner[0] <= inner[1] <= certified_inner[1]
    assert certified_outer[0] <= outer[0] <= outer[1] <= certified_outer[1]
    optical = contrast_lower(critical, inner[1], outer[0])
    orbital = orbital_box["orbital_contrast_upper"]
    return {
        "ring_ranges_arcsec": [inner, outer],
        "optical_contrast_lower": optical,
        "orbital_contrast_upper": orbital,
        "optical_minus_orbital_lower": optical - orbital,
    }


def run():
    canonical = json.loads((HERE / "results.json").read_text())[
        "all_positive_exteriors"
    ]
    _, critical = critical_density_bounds()
    _, (c1, c2), _ = geometry()
    params = [
        FIT[name][0]
        for name in (
            "axis_ratio",
            "projected_density_slope",
            "intermediate_axis_arcsec",
        )
    ]
    panels = {}
    for name, orbital in canonical.items():
        ranges = orbital["ring_ranges_arcsec"]
        if name == "extended_to_bolton_radii":
            doubled = canonical["ring_errors_2"]
            doubled_inner, doubled_outer = doubled["ring_ranges_arcsec"]
            assert ranges[0][1] == doubled_inner[1] and ranges[1] == doubled_outer
            small_inner = [ranges[0][0], doubled_inner[0]]
            pieces = [
                compare_box(critical[0], [small_inner, ranges[1]], orbital),
                compare_box(critical[0], doubled["ring_ranges_arcsec"], doubled),
            ]
        else:
            pieces = [compare_box(critical[0], ranges, orbital)]

        # Support of C2*e2-C1*e1 on the per-plane dimensionless error box.
        margin = -orbital["pressure_contrast_upper"]
        panels[name] = {
            "radius_partitions": pieces,
            "optical_minus_orbital_lower": min(
                piece["optical_minus_orbital_lower"] for piece in pieces
            ),
            "aligned_margin": margin,
            "aligned_geometry_threshold_percent_at_zero_secondary_bend": round(
                100 * margin / (c1 + c2), 6
            ),
        }

    return {
        "lens_fit": SOURCE,
        "fit_median_and_marginal_limits": {key: list(row) for key, row in FIT.items()},
        "surface_density_unit": "solar mass per square parsec",
        "nominal_primary_optical_contrast": round(
            circular_contrast(c1, *params, *THETA), 6
        ),
        "panels": panels,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    record = HERE / "elliptical.json"
    result = run()
    if args.write:
        record.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        assert result == json.loads(record.read_text())
        print("Replayed the circular-aperture contrast of the elliptical lens fit.")
