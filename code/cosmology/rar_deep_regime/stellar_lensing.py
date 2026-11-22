"""Predict inner SLACS Einstein angles from the same galaxies' stellar spectra."""

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad, simpson
from scipy.special import chndtr, gammaincinv, roots_legendre

HERE = Path(__file__).resolve().parent
RECEIPT = HERE / "receipts/stellar_lensing.json"
CATALOGUE = (
    "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?"
    "-source=J%2FApJ%2F682%2F964&-out.all&-out.max=unlimited"
)
TABLES = {
    4: (
        "Name zFG zBG Re sigma e_sigma Mph Mul Lens",
        "597131ad3336f5e363ff6812953862b4fd4274ea453437e6a3b7580ef0fc6c43",
    ),
    5: (
        "Name bSIE bLTM Good?",
        "f2db13ad32978d9d0ea25167ee0999a4ef4192e3d03cc7e4ea0925136cb85e7b",
    ),
}
C_KM_S = 299792.458
ARCSEC_PER_RADIAN = 180.0 * 3600.0 / np.pi
FIBRE_RADIUS = 1.5
PSF_SIGMA = 1.4 / np.sqrt(8 * np.log(2))


def catalogue_table(source, number):
    fields, digest = TABLES[number]
    block = source.split(f"#Table\tJ_ApJ_682_964_table{number}:")[1].split("\n#Table")[
        0
    ]
    header, *lines = block[block.index("\nrecno\t") + 1 :].splitlines()
    raw = csv.DictReader([header, *lines[2:]], delimiter="\t")
    rows = [{k: row[k].strip() for k in fields.split()} for row in raw]
    canonical = "\n".join(sorted("\t".join(row.values()) for row in rows))
    if hashlib.sha256(canonical.encode()).hexdigest() != digest:
        raise ValueError(f"Changed scientific fields in Bolton et al. Table {number}")
    return rows


def stellar_moments(order=96, size=801):
    """I and I*sigma_los²/v_c² from exact Abel deprojection and isotropic Jeans."""
    radius = np.geomspace(1e-5, 1e4, size)  # in effective radii
    b4 = gammaincinv(8, 0.5)
    light = np.exp(-b4 * radius**0.25)
    nodes, weights = roots_legendre(order)
    angle, weights = (nodes + 1) * np.pi / 4, weights * np.pi / 4
    sec = 1 / np.cos(angle)
    pressure = np.empty_like(radius)
    for i, projected in enumerate(radius):
        # r = R sec(theta) on the sight line. For each r, Abel-invert
        # I(s) = exp[-b4 (s/Re)^(1/4)] with s = r sec(phi).
        s = projected * sec[:, None] * sec[None, :]
        minus_derivative = (b4 / 4) * s ** (-0.75) * np.exp(-b4 * s**0.25)
        nu = (minus_derivative @ (weights * sec)) / np.pi
        # Interchanging the Jeans and projection integrals gives
        # I sigma_los²/v_c² = 2 integral_R^inf nu(r) sqrt(r²-R²) dr/r.
        pressure[i] = 2 * projected * (nu @ (weights * np.tan(angle) ** 2))
    # The all-light dispersion equals v_c/sqrt(3) by the scalar virial theorem.
    area = radius**2
    log_radius = np.log(radius)
    global_k = simpson(pressure * area, x=log_radius) / simpson(
        light * area, x=log_radius
    )
    np.testing.assert_allclose(global_k, 1 / 3, rtol=2e-6)
    return radius, light, pressure


def aperture_factor(re, moments):
    radius, light, pressure = moments
    projected = radius * re
    # Gaussian PSF encircled by an offset circular fibre, common spectroscopic
    # seeing from Bolton et al. §7. The truncated tail is below machine accuracy.
    in_range = projected < FIBRE_RADIUS + 12 * PSF_SIGMA
    probability = np.zeros_like(radius)
    probability[in_range] = chndtr(
        (FIBRE_RADIUS / PSF_SIGMA) ** 2, 2, (projected[in_range] / PSF_SIGMA) ** 2
    )
    area = probability * radius**2
    x = np.log(radius)
    return float(simpson(area * pressure, x=x) / simpson(area * light, x=x))


def comoving(z):
    return quad(lambda t: 1 / np.sqrt(0.3 * (1 + t) ** 3 + 0.7), 0, z)[0]


def summary(rows, lens_key, aperture_model):
    ratios, pulls = [], []
    for row in rows:
        k = row["k"] if aperture_model else 0.5
        measured = row[lens_key]
        # Score in the observed dispersion coordinate. Its uncertainty is
        # reported directly and avoids squaring a Gaussian speed error.
        expected = C_KM_S * np.sqrt(measured * k / row["geometry"])
        dlog_k = row["re_error_dlog_k"] if aperture_model else 0.0
        se = np.hypot(row["se_sigma"], 0.5 * expected * np.hypot(0.02, dlog_k))
        pulls.append((row["sigma"] - expected) / se)
        ratios.append(row["geometry"] * (row["sigma"] / C_KM_S) ** 2 / (k * measured))
    ratios, pulls = np.array(ratios), np.array(pulls)
    worst = int(np.argmax(abs(pulls)))
    return {
        "galaxies": len(rows),
        "predicted_to_measured_angle_median": round(float(np.median(ratios)), 6),
        "median_absolute_fractional_angle_error": round(
            float(np.median(abs(ratios - 1))), 6
        ),
        "conditional_chisq": round(float(pulls @ pulls), 5),
        "within_two_assigned_errors": int(np.sum(abs(pulls) <= 2)),
        "worst_object": rows[worst]["name"],
        "worst_signed_speed_residual_in_assigned_errors": round(float(pulls[worst]), 5),
    }


def outer_shear(moments):
    """Transfer the same 22 lenses' independent stellar masses and spectra to weak shear."""
    # Gavazzi et al. Table 2: centre in h^-1 kpc, (tangential, cross, se) in h Msun/pc².
    bins = np.array(
        [
            [3.3, 2307, -1015, 1570],
            [5.8, 918, -526, 614],
            [10.1, 115, 141, 264],
            [17.6, 81, -180, 153],
            [30.8, 232, -36, 80],
            [53.9, 100, 36, 46],
            [94.1, 90, -45, 27],
            [164.5, 52, 52, 17],
            [287.5, 60, 34, 17],
        ]
    )
    data = HERE / "data/slacs_weak_lenses.dat"
    lenses = np.loadtxt(data, dtype=str, skiprows=2)
    if len(set(lenses[:, 0])) != 22:
        raise ValueError("The entire 22-lens shear sample must be matched.")
    z, sigma, re, efficiency, chab, _, salp, e_salp = lenses[:, 1:].astype(float).T
    g = 4.300917270e-6  # kpc (km/s)^2/Msun
    h = 0.7
    kpc_to_km = 3.08567758149e16
    rotation = json.loads(
        (HERE / "receipts/sparc_deep_regime_diagnostic.json").read_text()
    )
    scale = rotation["deep_radial_acceleration"]
    a0 = scale["a0_fixed_exponent_m_s2"]
    a_lo, a_hi = scale["a0_galaxy_cluster_bootstrap_95pct_m_s2"]
    # SI acceleration to (km/s)^2/kpc.
    a0, a_lo, a_hi = np.array([a0, a_lo, a_hi]) * kpc_to_km / 1000
    distances = C_KM_S / (100 * h) * 1000 * np.array([comoving(v) for v in z]) / (1 + z)
    max_mass = 10 ** (salp + e_salp)
    domain = max(
        10 * max(re * distances / ARCSEC_PER_RADIAN),
        np.sqrt(g * max(max_mass) / (0.1 * a_lo)),
    )
    speed2 = sigma**2 / np.array([aperture_factor(r, moments) for r in re])
    # For complete physical annuli of uniform sky-pair density and shape variance,
    # source count scales as D_l^-2, cancelling D_l² from Sigma_crit^-2 in Eq. 10.
    weight = efficiency**2 / np.sum(efficiency**2)

    def projected(radius, mass, acceleration):
        anomaly = np.sqrt(acceleration * g * mass) / (4 * g * radius)
        return mass / (np.pi * radius**2) + anomaly, anomaly

    outputs = []
    for i, (centre, measured, cross, se) in enumerate(bins):
        radius = centre / h
        if radius < domain:
            continue
        chab_profile, _ = projected(radius, 10**chab, a0)
        salp_profile, _ = projected(radius, 10**salp, a0)
        div = 1e6 * h  # Msun/kpc² to h Msun/pc²
        flat_profile = speed2 / (4 * g * radius)
        # Non-overlapping ordered bins imply all contributing radii are above
        # the preceding bin's representative radius. Only use this envelope
        # where even that conservative radius is in the declared exterior.
        previous = bins[i - 1, 0] / h
        bound = None
        if previous >= domain:
            profile, surface = projected(previous, max_mass, a_hi)
            kappa = 4 * np.pi * g * distances * surface / C_KM_S**2  # D_ls/D_s <= 1
            if np.max(kappa) >= 1:
                raise ValueError(
                    "The reduced-shear bound requires a subcritical exterior."
                )
            # Equation 10 converts each shape using the lens's catalogue
            # mean efficiency. If the bin's shape-weighted efficiency agrees,
            # E[beta / mean_beta] = 1 for that lens. Without this assumption
            # beta / mean_beta can reach 1 / mean_beta.
            calibrated = profile / ((1 - kappa) * div)
            bound = {}
            for name, limits in (
                ("per_lens_bin_mean_efficiency_calibrated", calibrated),
                ("only_source_efficiency_at_most_one", calibrated / efficiency),
            ):
                largest = int(np.argmax(limits))
                upper = float(limits[largest])
                bound[name] = dict(
                    bounding_lens=lenses[largest, 0],
                    upper_esd=round(upper, 4),
                    measured_minus_upper_in_quoted_shear_errors=round(
                        (measured - upper) / se, 4
                    ),
                )
        outputs.append(
            dict(
                centre_hinv_kpc=centre,
                tangential=measured,
                cross=cross,
                quoted_error=se,
                chabrier=round(float(weight @ chab_profile / div), 4),
                salpeter=round(float(weight @ salp_profile / div), 4),
                inner_speed_extended_flat=round(float(weight @ flat_profile / div), 4),
                ordered_bin_mass_box_bound=bound,
            )
        )
    return {
        "sources": [
            "https://doi.org/10.1086/519237",
            "https://doi.org/10.1088/0004-637X/705/2/1099",
        ],
        "source_pdf_sha256": {
            "Gavazzi published": "9d2fbe8b2a5ea0de7734113d1645222c2ff43b3122cf84cb5033cc32cd05012e",
            "astro-ph/0701589v2": "562896a54435d24f1beae72266b593ba899e2fbcf4de400a87dbe486da9c4a12",
            "0911.2471": "31f2d6b38c8143a9d8216a0cdfda3ca95e7f5536d3f784f479f47c01228e1668",
        },
        "spectra_and_stellar_population_table_sha256": hashlib.sha256(
            data.read_bytes()
        ).hexdigest(),
        "a0_m_s2": scale["a0_fixed_exponent_m_s2"],
        "outer_centre_floor_physical_kpc": round(float(domain), 4),
        "selection": "All 22 lenses from Gavazzi Table 1, their unchanged radial shear stack in Table 2 and their exact SDSS matches in Auger Table 4. Use bin centres beyond 10 times the largest measured effective radius and g_bar < 0.1*a_lo even for the largest Salpeter + 1 quoted log-mass-error. a_lo is the lower SPARC scale interval endpoint.",
        "channels": "Compact stellar baryons plus frozen SPARC square-root exterior, Chabrier and Salpeter stellar population models. A third control extends the measured central stellar field flat to the outer stack. All use zero slip, linear shear and complete uniform sky coverage at representative radii. Actual shape weights, per-lens fields of view, pair radii, perturbers and extended gas are unavailable here.",
        "bound": "Separate feasibility bound for compact stars plus anomaly. Give each galaxy up to its Salpeter + 1 quoted error stellar mass and the SPARC upper interval scale. Allow any positive lens weights and radii in the ordered disjoint bin. Use the previous centre as a lower radius and D_ls/D_s <= 1 to bound reduced shear. The tighter upper bound assumes each lens's shape-weighted mean source efficiency in this bin equals its Table 1 value. The wider bound allows arbitrary positive source weights and any efficiency in [0,1]. These parameter ranges describe sensitivities.",
        "bins_in_h_Msun_pc2": outputs,
        "centre_conditional_chisq": {
            key: round(
                sum(
                    ((r["tangential"] - r[key]) / r["quoted_error"]) ** 2
                    for r in outputs
                ),
                4,
            )
            for key in ("chabrier", "salpeter", "inner_speed_extended_flat")
        },
        "cross_zero_chisq": round(
            sum((r["cross"] / r["quoted_error"]) ** 2 for r in outputs), 4
        ),
    }


def build(path):
    source = path.read_text()
    spectra = catalogue_table(source, 4)
    images = catalogue_table(source, 5)
    eligible = {
        r["Name"]: r
        for r in spectra
        if r["Mph"] == "E"
        and r["Mul"] == "S"
        and r["Lens"] == "A"
        and r["Re"]
        and r["sigma"]
        and r["e_sigma"]
    }
    moments = stellar_moments()
    rows = []
    for image in images:
        if image["Good?"] != "Yes" or image["Name"] not in eligible:
            continue
        spec = eligible[image["Name"]]
        re, sigma = float(spec["Re"]), float(spec["sigma"])
        k = aperture_factor(re, moments)
        lens, background = comoving(float(spec["zFG"])), comoving(float(spec["zBG"]))
        # Deflection of a spherical logarithmic potential with Phi=Psi:
        # theta_E = 2 pi (v_c/c)^2 D_ls/D_s, k = sigma_ap²/v_c².
        geometry = 2 * np.pi * (1 - lens / background) * ARCSEC_PER_RADIAN
        rows.append(
            dict(
                name=image["Name"],
                k=k,
                sigma=sigma,
                se_sigma=float(spec["e_sigma"]),
                geometry=geometry,
                bSIE=float(image["bSIE"]),
                bLTM=float(image["bLTM"]),
                re_error_dlog_k=0.5
                * np.log(
                    aperture_factor(re * 1.035, moments)
                    / aperture_factor(re * 0.965, moments)
                ),
                physical_radius_kpc=(
                    C_KM_S
                    / 70
                    * 1000
                    * lens
                    / (1 + float(spec["zFG"]))
                    * float(image["bSIE"])
                    / ARCSEC_PER_RADIAN
                ),
            )
        )
    if len(rows) != len({row["name"] for row in rows}):
        raise ValueError("Matched systems must be unique")
    return {
        "primary_paper": "https://doi.org/10.1086/589327",
        "catalogue": CATALOGUE,
        "outer_shear_same_22_lenses": outer_shear(moments),
        "sha256_scientific_fields_sorted": {
            str(i): digest for i, (_, digest) in TABLES.items()
        },
        "selection": "All Table 5 modeled A-grade lenses with good SDSS dispersion and a Table 4 single early-type galaxy with measured effective radius, dispersion and error. No selection on residuals.",
        "inputs": "Flat Omega_m=0.3, Omega_Lambda=0.7. SDSS fibre radius 1.5 arcsec, survey median seeing FWHM 1.4 arcsec. Condition on a spherical isothermal total field, zero slip and isotropic spherical de Vaucouleurs starlight with each galaxy's measured intermediate-axis half-light radius. The spectral prediction precedes comparison with the arcs. The SLACS target selection itself favours strong lenses.",
        "error_scope": "Independent quoted SDSS speed error and first-order propagation of the survey 2% Einstein-angle and 3.5% light-radius empirical errors. Shared photometric covariance, sight-line structure, anisotropy, rotation, individual seeing, profile variation and selection scatter are unmodeled. Conditional scores have no intrinsic-scatter fit or population significance assigned.",
        "sie_einstein_radius_range_kpc": [
            round(float(f(r["physical_radius_kpc"] for r in rows)), 4)
            for f in (min, max)
        ],
        "primary_isothermal_ellipsoid_imaging": summary(rows, "bSIE", True),
        "same_arcs_light_traces_mass_imaging": summary(rows, "bLTM", True),
        "inverse_square_starlight_control": summary(rows, "bSIE", False),
        "object_columns": [
            "name",
            "predicted_angle_arcsec",
            "SIE_measured_arcsec",
            "LTM_measured_arcsec",
        ],
        "objects": [
            [
                r["name"],
                round(r["geometry"] * (r["sigma"] / C_KM_S) ** 2 / r["k"], 5),
                r["bSIE"],
                r["bLTM"],
            ]
            for r in rows
        ],
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("catalogue", type=Path)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = build(args.catalogue)
    if args.write:
        RECEIPT.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    elif result != json.loads(RECEIPT.read_text()):
        raise ValueError("Stellar-lensing receipt differs")
    print(json.dumps(result["primary_isothermal_ellipsoid_imaging"], indent=2))
