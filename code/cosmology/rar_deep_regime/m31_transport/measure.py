"""Forward M31 orbits from a photometric source and stationary mark transport."""

import argparse
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path
import re

import numpy as np
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.optimize import brentq
from scipy.special import gammaincinv, ive, kve, roots_legendre

from cosmology.record_abundance import clustering, stationary

HERE = Path(__file__).parent
INPUT = HERE / "inputs.json"
RESULT = HERE / "results.json"
PRIMARY = [
    (
        "0909.3846v1",
        "lc1.tex",
        "d0f9414303f906e55820e5dc2fec2b032cf188a945b33d529b6a6d752b3f0cd1",
    ),
    (
        "0912.4133v1",
        "ec_m31.tex",
        "c3eff6565a80820c27d6110c466485f72a39b7c5703bb78057a1ef2d582d60eb",
    ),
]
G = 4.30091e-6  # (km/s)^2 kpc / solar mass
KPC_M = 3.085677581491367e19
DISK_SCALE, BULGE_RE, BULGE_INDEX = 5.6, 1.3, 1.1
R_MIN = 5 * DISK_SCALE
ETA = 0.11933 / 0.02242
MOLECULAR = 3.5e8
HELIUM = 1.3


def published_inputs(texts):
    chemin, corbelli = texts
    ctable = (
        chemin.split(r"\tablecaption{\hi\ rotation curve, surface density", 1)[1]
        .split(r"\startdata", 1)[1]
        .split(r"\enddata", 1)[0]
    )
    gas = []
    orbits = {"Chemin": [], "Corbelli": []}
    for line in ctable.strip().splitlines():
        cells = [s.strip() for s in line.split(r"\\", 1)[0].split("&")]
        if not len(cells) == 11:
            raise ValueError("Unexpected Chemin row")
        r = float(cells[1])
        gas.append([r, float(cells[10])])
        if r >= R_MIN:
            orbits["Chemin"].append(
                [r, float(cells[8]), float(cells[9]), float(cells[7])]
            )
    table = corbelli.split(
        r"\caption{ The HI rotation curve of M31 and the parameters", 1
    )[1].split(r"\end{tabular}", 1)[0]
    for line in table.splitlines():
        if re.match(r"^\s*\d+\.\d+\s*&", line):
            cells = [float(v) for v in line.split(r"\\", 1)[0].split("&")]
            if cells[0] >= R_MIN:
                orbits["Corbelli"].append([cells[0], cells[-2], cells[-1], cells[1]])
    assert (
        len(gas) == 100 and len(orbits["Corbelli"]) == 9 and len(orbits["Chemin"]) > 20
    )
    assert all(
        value in chemin
        for value in (
            "5.6  \\pm 0.1",
            "1.3 \\pm 0.1",
            "1.1 \\pm 0.1",
            r"6.4 \times 10^{10}",
            r"7.1 \times 10^{10}",
            r"3.5\times 10^8",
        )
    )
    return {
        "primary_tex_sha256": {ref: digest for ref, _, digest in PRIMARY},
        "HI_radius_kpc_and_surface_msun_pc2": gas,
        "outer_radius_kpc_speed_error_kms_and_inclination_deg": orbits,
    }


@dataclass
class Source:
    r: np.ndarray  # Shell midpoints in kpc
    edge: np.ndarray
    disk: np.ndarray  # Normalized thin exponential-disk mass in shells
    bulge: np.ndarray  # Deprojection of the adopted spherical Sersic profile
    gas: np.ndarray  # Measured radial weights, including H2 and helium, in Msun


def source(gas_profile, dr):
    edge = np.arange(round(140 / dr) + 1) * dr
    r = (edge[:-1] + edge[1:]) / 2
    disk = np.diff(1 - (1 + edge / DISK_SCALE) * np.exp(-edge / DISK_SCALE))
    # Abel inversion with projected radius R = r cosh(t) removes the
    # square-root endpoint singularity in the integral of -dI/dR.
    x, w = roots_legendre(80)
    proj = r[:, None] * np.cosh(8 * (1 + x))[None, :]
    n = BULGE_INDEX
    b = gammaincinv(2 * n, 0.5)
    y = (proj / BULGE_RE) ** (1 / n)
    bulge = r**2 * ((np.exp(-b * y) * y / proj) @ w)
    # Treat a quoted annulus as constant up to the midpoints to its neighbors.
    # The last annulus extends half a bin beyond the last quoted center.
    radius, sigma = np.array(gas_profile).T
    bounds = np.r_[
        0, (radius[1:] + radius[:-1]) / 2, radius[-1] + np.diff(radius)[-1] / 2
    ]
    cumulative = np.r_[0, np.cumsum(np.pi * 1e6 * np.diff(bounds**2) * sigma)]
    ix = np.searchsorted(bounds, edge, side="right").clip(1, len(bounds) - 1) - 1
    enclosed = cumulative[ix] + np.pi * 1e6 * sigma[ix] * (
        np.minimum(edge, bounds[-1]) ** 2 - bounds[ix] ** 2
    )
    hi = np.diff(enclosed)
    assert hi.min() >= 0
    gas = HELIUM * hi * (1 + MOLECULAR / hi.sum())
    return Source(r, edge, disk / disk.sum(), bulge / bulge.sum(), gas)


def ru(src, weights, mu, radii=None):
    """r times the spherical Yukawa potential of a normalized shell distribution."""
    if radii is None:
        radii = src.edge
        ix = slice(None)
    else:
        ix = np.searchsorted(src.r, radii)
    if mu == 0:
        inside = np.r_[0, np.cumsum(weights)]
        outside = np.r_[np.cumsum((weights / src.r)[::-1])[::-1], 0]
        return inside[ix] + radii * outside[ix]
    x = mu * src.r
    inside = np.r_[0, np.cumsum(weights * np.sinh(x) / x)]
    outside = np.r_[np.cumsum((weights * np.exp(-x) / x)[::-1])[::-1], 0]
    return inside[ix] * np.exp(-mu * radii) + outside[ix] * np.sinh(mu * radii)


def total_length(src, weights, mu):
    radial_field = ru(src, weights, mu)
    # Beyond the last shell r*u is an exact decaying exponential.
    return trapezoid(radial_field**2, src.edge) + radial_field[-1] ** 2 / (2 * mu)


def gas_force(src, radii, height=0.5):
    angles, w = roots_legendre(80)
    cos_phi = np.cos(np.pi / 2 * (angles + 1))
    keep = src.gas > 0
    s, mass = src.r[keep, None], src.gas[keep]
    v2 = []
    for r in radii:
        separation = (r - s) ** 2 + 2 * r * s * (1 - cos_phi) + height**2
        mean = ((r - s * cos_phi) / separation**1.5) @ (w / 2)
        v2.append(G * r * (mass @ mean))
    return np.array(v2)


def ordinary_v2(src, radii, disk_mass, bulge_mass, geometry):
    # The 1.4e8-Msun central black hole is kept in the ordinary force only.
    if geometry == "disk":
        y = radii / (2 * DISK_SCALE)
        disk = (
            (2 * G * disk_mass / DISK_SCALE)
            * y**2
            * (ive(0, y) * kve(0, y) - ive(1, y) * kve(1, y))
        )
        enclosed_bulge = np.interp(radii, src.edge, np.r_[0, np.cumsum(src.bulge)])
        return (
            disk
            + G * (1.4e8 + bulge_mass * enclosed_bulge) / radii
            + gas_force(src, radii)
        )
    if geometry == "spherical":
        mass = disk_mass * src.disk + bulge_mass * src.bulge + src.gas
        return (
            G * (1.4e8 + np.interp(radii, src.edge, np.r_[0, np.cumsum(mass)])) / radii
        )
    raise ValueError(f"Unknown baryon geometry: {geometry}")


def calibrated_slope(mass, a0_si):
    return np.sqrt(mass * a0_si * KPC_M / (1e6 * G))


def predict(src, radii, disk_mass, bulge_mass, a0_si, geometry="disk"):
    mass = disk_mass * src.disk + bulge_mass * src.bulge + src.gas
    weights = mass / mass.sum()
    slope = calibrated_slope(mass.sum(), a0_si)
    target = ETA * mass.sum() / slope
    mu = brentq(
        lambda value: total_length(src, weights, value) - target, 1e-6, 0.5, xtol=1e-14
    )
    vb2 = ordinary_v2(src, radii, disk_mass, bulge_mass, geometry)
    output = {}
    for condition, inverse_length in (("cosmic_share", mu), ("lossless", 0)):
        integrand = ru(src, weights, inverse_length) ** 2
        cumulative = np.interp(
            radii, src.edge, cumulative_trapezoid(integrand, src.edge, initial=0)
        )
        output[condition] = np.sqrt(vb2 + G * slope * cumulative / radii)
    return {
        "mu_per_kpc": mu,
        "point_source_density_efold_kpc": 1 / (2 * mu),
        "total_baryonic_injector_mass_msun": mass.sum(),
        "gas_mass_msun": src.gas.sum(),
        "lossless_outer_mass_slope_msun_kpc": slope,
        "baryon_enclosure_fraction": np.interp(
            radii, src.edge, np.r_[0, np.cumsum(weights)]
        ),
        "ordinary_speed": np.sqrt(vb2),
        "all_radial_injectors_speed_ceiling": np.sqrt(vb2 + G * slope),
        **output,
    }


def run(inputs, dr=0.05):
    src = source(inputs["HI_radius_kpc_and_surface_msun_pc2"], dr)
    calibration = json.loads(
        (HERE.parent / "receipts/lensing_profile.json").read_text()
    )
    a0 = json.loads((HERE.parent / "finite_halo/results.json").read_text())[
        "sparc_a0_m_s2"
    ]
    limits = calibration["sparc_cluster_bootstrap_95pct_m_s2"]
    banks = inputs["outer_radius_kpc_speed_error_kms_and_inclination_deg"]
    radii = np.unique([row[0] for rows in banks.values() for row in rows])
    central = predict(src, radii, 7.1e10, 6.4e10, a0)
    hybrid = predict(src, radii, 7.1e10, 2.4e10, a0)
    spherical = predict(src, radii, 7.1e10, 6.4e10, a0, "spherical")
    # These scenarios are shared across all rings, frozen before examining velocities.
    variants = [
        predict(src, radii, 7.1e10 * factor, bulge, accel)
        for factor in (0.75, 1.25)
        for bulge in (2.4e10, 6.4e10 * 0.75, 6.4e10 * 1.25)
        for accel in limits
    ]
    surveys = {}
    for name, rows in banks.items():
        r, v, e, inc = np.array(rows).T
        index = np.searchsorted(radii, r)
        forecasts = {
            s: central[s][index]
            for s in (
                "cosmic_share",
                "lossless",
                "all_radial_injectors_speed_ceiling",
                "ordinary_speed",
            )
        }
        lo = np.min([x["cosmic_share"][index] for x in variants], axis=0)
        hi = np.max([x["cosmic_share"][index] for x in variants], axis=0)
        upper = np.max(
            [x["all_radial_injectors_speed_ceiling"][index] for x in variants], axis=0
        )
        surveys[name] = {
            "rows": [
                {
                    "r_kpc": rr,
                    "observed_kms": vv,
                    "quoted_error_kms": ee,
                    "inclination_deg": ii,
                    **{s: vals[j] for s, vals in forecasts.items()},
                    "hybrid_cosmic_share_kms": hybrid["cosmic_share"][index[j]],
                    "spherical_cosmic_share_kms": spherical["cosmic_share"][index[j]],
                    "enclosed_baryons": central["baryon_enclosure_fraction"][index[j]],
                    "finite_sensitivity_kms": [lo[j], hi[j]],
                    "upper_sensitivity_kms": upper[j],
                }
                for j, (rr, vv, ee, ii) in enumerate(rows)
            ],
            "central_residual_rms_kms": {
                s: np.sqrt(np.mean((v - vals) ** 2)) for s, vals in forecasts.items()
            },
            "central_mean_observed_minus_predicted_kms": {
                s: np.mean(v - vals) for s, vals in forecasts.items()
            },
            "inside_one_quoted_error": {
                s: int(np.sum(abs(v - vals) <= e)) for s, vals in forecasts.items()
            },
            "above_entire_finite_envelope_plus_one_error": int(np.sum(v - e > hi)),
            "below_entire_finite_envelope_minus_one_error": int(np.sum(v + e < lo)),
            "above_all_injector_ceiling_plus_one_error_favorable_corner": int(
                np.sum(v - e > upper)
            ),
        }
    return {
        "protocol_sha256": sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest(),
        "source_sha256": inputs["primary_tex_sha256"],
        "dr_kpc": dr,
        "SPARC_scale_m_per_s2": a0,
        "SPARC_scale_95_bootstrap_m_per_s2": limits,
        "external_halo_dark_to_baryon_mass": ETA,
        **{
            k: central[k]
            for k in (
                "mu_per_kpc",
                "point_source_density_efold_kpc",
                "total_baryonic_injector_mass_msun",
                "gas_mass_msun",
                "lossless_outer_mass_slope_msun_kpc",
            )
        },
        "surveys": surveys,
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--fetch",
        action="store_true",
        help="verify against both original primary TeX sources",
    )
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    if args.write and not args.fetch:
        p.error("writing requires original sources")
    if args.fetch:
        with ThreadPoolExecutor(max_workers=2) as pool:
            inputs = published_inputs(pool.map(clustering.source, PRIMARY))
        if args.write:
            with INPUT.open("x") as out:
                json.dump(inputs, out, indent=2)
                out.write("\n")
        else:
            stationary.check(json.loads(INPUT.read_text()), inputs)
    result = run(json.loads(INPUT.read_text()))
    if args.write:
        with RESULT.open("x") as out:
            json.dump(result, out, indent=2, allow_nan=False)
            out.write("\n")
    else:
        stationary.check(json.loads(RESULT.read_text()), result)
    for name, survey in result["surveys"].items():
        print(
            name,
            len(survey["rows"]),
            survey["central_residual_rms_kms"],
            survey["rows"][-1],
        )
