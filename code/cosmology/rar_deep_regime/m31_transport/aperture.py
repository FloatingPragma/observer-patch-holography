"""Bound outer M31 aperture growth for every survival and every positive gain."""

import argparse
import json
import numpy as np
from scipy.integrate import trapezoid

from cosmology.rar_deep_regime.m31_transport import measure as m


def aperture_bounds(src, weights, r1, r2, mu_max=1.0, panels=160):
    """Upper bound on M(r2)/M(r1) across the whole family of finite shell inputs."""
    # Split where a spherical shell injects. Between injections (r u)^2
    # is convex, so midpoint and trapezoid give certified volume-integral brackets.
    cuts = np.unique(np.r_[0.0, src.r[src.r < r2], r1, r2])
    index = np.searchsorted(cuts, r1)
    mid = (cuts[:index] + cuts[1 : index + 1]) / 2
    widths = np.diff(cuts[: index + 1])
    outer_edges = cuts[index:]
    mesh = mu_max * np.linspace(0, 1, panels + 1) ** 2
    outer_upper = []
    inner_lower = []
    for mu in mesh:
        inner_lower.append(widths @ m.ru(src, weights, mu, mid) ** 2)
        outer_upper.append(
            trapezoid(m.ru(src, weights, mu, outer_edges) ** 2, outer_edges)
        )
    # Both aperture integrals fall separately with mu (positive Yukawa kernel).
    finite_range_ratio = 1 + np.array(outer_upper[:-1]) / np.array(inner_lower[1:])
    worst = finite_range_ratio.argmax()

    # For mu >= mu_max use F=2 mu r u. Each shell i contributes
    # (p_i/s_i)[exp(-mu|r-s_i|)-exp(-mu(r+s_i))]. All terms are nonnegative.
    # Outer L2 norm via the triangle inequality and a one-sided exponential
    # integral for each shell below r1; use a full-line integral for the rest.
    amp = weights / src.r
    inner = src.r <= r1
    outer_norm = (
        amp[inner] @ np.exp(-mu_max * (r1 - src.r[inner])) / np.sqrt(2)
        + amp[~inner].sum()
    )
    # Lower bound uses just the positive diagonal terms in F^2. Integrate
    # each from s_i to s_i+1/mu, wholly inside r1 for the selected shells.
    selected = src.r <= r1 - 1 / mu_max
    inner_norm2 = (
        (1 - np.exp(-2))
        / 2
        * np.sum((amp[selected] * -np.expm1(-2 * mu_max * src.r[selected])) ** 2)
    )
    tail_ratio = 1 + outer_norm**2 / inner_norm2
    return {
        "upper_mass_ratio_all_lifetimes": max(finite_range_ratio[worst], tail_ratio),
        "upper_mass_ratio_high_loss": tail_ratio,
        "interval_for_largest_bound_per_kpc": mesh[worst : worst + 2].tolist(),
        "mu_split_per_kpc": mu_max,
        "panels": panels,
    }


def required(radii, observed, errors, vb2, e):
    v1, v2 = observed + errors * np.array([e, -e])
    if v1**2 <= vb2[0]:
        raise ValueError("Inner excess mass must be positive")
    return radii[1] * (v2**2 - vb2[1]) / (radii[0] * (v1**2 - vb2[0]))


def run(inputs, dr=0.05):
    src = m.source(inputs["HI_radius_kpc_and_surface_msun_pc2"], dr)
    # Only the separately reported absolute slope comparison uses SPARC.
    a0 = json.loads((m.HERE.parent / "finite_halo/results.json").read_text())[
        "sparc_a0_m_s2"
    ]
    surveys = {}
    for name, rows in inputs[
        "outer_radius_kpc_speed_error_kms_and_inclination_deg"
    ].items():
        radii, observed, errors, _ = np.array([rows[0], rows[-1]]).T
        results = {}
        variants = [(7.1e10, 6.4e10, "disk"), (7.1e10, 6.4e10, "spherical")]
        variants += [
            (7.1e10 * f, b, "disk")
            for f in (0.75, 1.0, 1.25)
            for b in (2.4e10, 0.75 * 6.4e10, 1.25 * 6.4e10)
        ]
        for dm, bm, geometry in variants:
            mass = dm * src.disk + bm * src.bulge + src.gas
            vb2 = m.ordinary_v2(src, radii, dm, bm, geometry)
            gain = m.G * m.calibrated_slope(mass.sum(), a0)
            bounds = aperture_bounds(src, mass / mass.sum(), *radii)
            ratios = [required(radii, observed, errors, vb2, e) for e in (0, 1, 2)]
            growth = []
            for e in (0, 1, 2):
                v = observed + errors * [e, -e]
                growth.append(
                    float(np.diff(radii * (v**2 - vb2))[0] / np.diff(radii)[0] / gain)
                )
            case = {
                **bounds,
                "required_mass_ratio_nominal_one_two_error_boxes": ratios,
                "required_slope_over_SPARC_ceiling_nominal_one_two_error_boxes": growth,
                "outer_speed_ceiling_conditioned_on_inner_measurement_kms": np.sqrt(
                    vb2[1]
                    + bounds["upper_mass_ratio_all_lifetimes"]
                    * radii[0]
                    / radii[1]
                    * (observed[0] ** 2 - vb2[0])
                ),
            }
            if bm != 2.4e10:
                # Identical angular apertures at the alternate distance.
                # Light and H I masses scale with D^2, sizes with D, and the
                # fixed central black-hole mass transforms oppositely in v^2.
                # The dimensionless all-mu shape bound remains the same.
                bh = m.G * 1.4e8 / radii
                case["two_error_mass_ratio_at_735_835_kpc"] = [
                    required(radii, observed, errors, d * (vb2 - bh) + bh / d, 2)
                    for d in (735 / 785, 835 / 785)
                ]
            results[f"disk_{dm / 1e10:g}_bulge_{bm / 1e10:g}_{geometry}"] = case
        surveys[name] = {"endpoint_radii_kpc": radii.tolist(), "source_cases": results}
    return {"dr_kpc": dr, "surveys": surveys}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    out = run(json.loads(m.INPUT.read_text()))
    path = m.HERE / "aperture.json"
    if args.write:
        with path.open("x") as handle:
            json.dump(out, handle, indent=2, allow_nan=False)
            handle.write("\n")
    else:
        m.stationary.check(json.loads(path.read_text()), out)
    for name, survey in out["surveys"].items():
        print(name, next(iter(survey["source_cases"].values())))
