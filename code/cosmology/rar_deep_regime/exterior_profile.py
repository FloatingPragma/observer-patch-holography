"""Bound the common isothermal acceleration scale in the observed SPARC outskirts."""

from collections import defaultdict
from statistics import NormalDist

import numpy as np
from scipy.optimize import brentq

import rar_deep_regime_diagnostic as base


def exterior_rings(table1, table2, include_bulges=False):
    """Select by baryons and catalogue geometry, retaining all measured speeds."""
    meta = {galaxy["name"]: galaxy for galaxy in table1}
    bulges = {ring["name"] for ring in table2 if ring["vbul_kms"] != 0}
    groups = defaultdict(list)
    for ring in table2:
        name = ring["name"]
        galaxy = meta[name]
        disk_radius, hi_radius = galaxy["disk_scale_kpc"], galaxy["hi_radius_kpc"]
        r = ring["rad_kpc"]
        if (
            (name in bulges and not include_bulges)
            or galaxy["quality"] > 2
            or galaxy["inclination_deg"] < 30
            or min(disk_radius, hi_radius) <= 0
            or r < max(4 * disk_radius, hi_radius)
        ):
            continue
        gas_v, disk_v = ring["vgas_kms"], ring["vdisk_kms"]
        fixed_baryons, disk = gas_v * abs(gas_v), disk_v * abs(disk_v)
        bulge_v = ring["vbul_kms"]
        fixed_baryons += (
            base.UPSILON_BULGE * bulge_v * abs(bulge_v)
        )  # Fixed nondisk force.
        gbar = (
            (fixed_baryons + base.UPSILON_DISK * disk) * base.KM**2 / (r * base.KPC_M)
        )
        if 0 < gbar < 0.1 * base.DEEP_CUT_REFERENCE_A0_M_S2:
            groups[name].append(
                (r, ring["vobs_kms"], ring["e_vobs_kms"], fixed_baryons, disk)
            )
    return [(meta[name], np.asarray(sorted(groups[name])).T) for name in sorted(groups)]


def sin_deg(angle):
    return np.sin(np.deg2rad(angle))


def nuisance_bounds(galaxy, width):
    """Bound t = sin(i0)^2 / [(D/D0) sin(i)^2] and the disk mass-to-light ratio."""
    distance = galaxy["dist_mpc"]
    fraction = width * galaxy["e_dist_mpc"] / distance
    inclination = galaxy["inclination_deg"]
    tilt = width * galaxy["e_inclination_deg"]
    i_lo, i_hi = inclination - tilt, min(90.0, inclination + tilt)
    if fraction >= 1 or not 0 < inclination <= 90:
        raise ValueError("Unbounded distance or deprojection in the selected sample")
    sin0_sq = sin_deg(inclination) ** 2
    t_lo = sin0_sq / (sin_deg(i_hi) ** 2 * (1 + fraction))
    t_hi = sin0_sq / (sin_deg(i_lo) ** 2 * (1 - fraction)) if i_lo > 0 else np.inf
    offset = 0.1 if width else 0.0
    return t_lo, t_hi, 0.5 - offset, 0.5 + offset


def mass_coefficient(galaxy, mass_to_light):
    mass = mass_to_light * galaxy["L36_gsun"] + base.GAS_HELIUM * galaxy["MHI_gsun"]
    return np.sqrt(base.G_SI * base.M_SUN * 1e9 * mass) / base.KM**2


def predicted_speed_squared(galaxy, radius, baryons, ml, scale, model):
    if model == "compact_mass":
        extra = mass_coefficient(galaxy, ml) * np.sqrt(scale)
    elif model == "local_force":
        extra = np.sqrt(radius * base.KPC_M / base.KM**2 * baryons * scale)
    else:
        raise ValueError(f"Unknown baryonic prescription {model}")
    return baryons + extra


def speed_envelopes(galaxy, rings, width, scale_interval):
    """Predict catalogue-frame speed intervals from baryons, geometry and a scale box."""
    radius, _, _, fixed, disk = rings
    t_lo, t_hi, ml_lo, ml_hi = nuisance_bounds(galaxy, width)
    a_lo, a_hi = scale_interval
    if a_lo <= 0 or a_lo > a_hi:
        raise ValueError("Require an ordered positive acceleration interval")
    baryons = [fixed + f(ml_lo * disk, ml_hi * disk) for f in (np.minimum, np.maximum)]
    if np.any(baryons[0] <= 0):
        raise ValueError("Local-force square root crosses zero on this nuisance box")
    # Force-equivalent mass at a ring is R Vbar² / G in a disk force model.
    # Both squared-speed prescriptions scale as D at fixed angular size and flux.
    return {
        model: np.array(
            [
                np.sqrt(
                    predicted_speed_squared(
                        galaxy, radius, baryons[0], ml_lo, a_lo, model
                    )
                    / t_hi
                ),
                np.sqrt(
                    predicted_speed_squared(
                        galaxy, radius, baryons[1], ml_hi, a_hi, model
                    )
                    / t_lo
                ),
            ]
        )
        for model in ("compact_mass", "local_force")
    }


def compare_baryon_sources(groups, scale, training_interval, tolerance):
    """One frozen scale, identical radial samples and catalogue-frame speed bounds."""
    cases = {
        "catalogue": (0, (scale, scale)),
        "one_error_box": (1, (scale, scale)),
        "two_error_box": (2, (scale, scale)),
        "two_error_box_and_training_interval": (2, training_interval),
    }
    models = ("compact_mass", "local_force")
    limits = {case: {name: {} for name in models} for case in cases}
    names, sizes, radii = [], [], []
    predictions = {name: [] for name in models}
    for galaxy, rings in groups:
        name = galaxy["name"]
        radius, observed, error, _, _ = rings
        names.append(name)
        sizes.append(len(radius))
        radii.extend(radius)
        for case, (width, scales) in cases.items():
            for model, (lower, upper) in speed_envelopes(
                galaxy, rings, width, scales
            ).items():
                if case == "catalogue":
                    np.testing.assert_array_equal(lower, upper)
                    predictions[model].extend(lower)
                # Nearest speed in the enclosure gives a necessary error for each ring,
                # even if that ring could choose its own independent nuisance parameters.
                discrepancy = (
                    np.maximum.reduce(
                        (lower - observed, observed - upper, np.zeros_like(observed))
                    )
                    / error
                )
                j = int(np.argmax(discrepancy))
                summary = limits[case][model]
                count = int(np.count_nonzero(discrepancy > tolerance))
                summary["rings_above_threshold"] = (
                    summary.get("rings_above_threshold", 0) + count
                )
                summary["galaxies_with_ring_above_threshold"] = summary.get(
                    "galaxies_with_ring_above_threshold", 0
                ) + int(count > 0)
                if (
                    not summary.get("limiting_ring")
                    or discrepancy[j] > summary["worst_ring_lower_bound"]
                ):
                    summary["worst_ring_lower_bound"] = float(discrepancy[j])
                    summary["limiting_ring"] = {
                        "galaxy": name,
                        "radius_kpc": radius[j],
                        "observed_kms": observed[j],
                        "quoted_error_kms": error[j],
                        "predicted_speed_interval_kms": [lower[j], upper[j]],
                    }
    for model, summary in limits["two_error_box"].items():
        # Diagnose coexistence on all the other rings after ranking on the full sample.
        limiting = summary["limiting_ring"]
        retained = []
        for galaxy, rings in groups:
            include = np.ones(len(rings[0]), dtype=bool)
            if galaxy["name"] == limiting["galaxy"]:
                include = rings[0] != limiting["radius_kpc"]
            if include.any():
                retained.append((galaxy, rings[:, include]))
        witness = joint_witness(retained, 2, tolerance, scale, model)
        if witness is not None:
            summary["shared_witness_after_omitting_limiting_ring"] = {
                "galaxies": len(witness),
                "rings": sum(row["rings"] for row in witness),
                "maximum_residual": max(
                    row["maximum_absolute_velocity_residual_in_quoted_errors"]
                    for row in witness
                ),
                "columns": [
                    "galaxy",
                    "distance_mpc",
                    "inclination_deg",
                    "disk_mass_to_light",
                ],
                "values": [
                    [
                        row["galaxy"],
                        row["distance_mpc"],
                        row["inclination_deg"],
                        row["disk_mass_to_light"],
                    ]
                    for row in witness
                ],
            }
        else:
            summary["shared_witness_after_omitting_limiting_ring"] = None
    return {
        "comparison": "Fixed training scale and identical exterior selection. Compact extra squared speed sqrt(G Mphot a0); local-force extra sqrt(R Vbar² a0) in consistent units. Both add to the same Vbar². The force-equivalent mass R Vbar²/G does not measure the physical mass enclosed by a disk or a modular charge.",
        "bounds": "Necessary per-ring residuals in catalogue speed units. The wide case allows every nuisance and every calibration inside the training bootstrap interval independently at each ring. Error fractions condition on the catalogue model. Shared galaxy constraints can only increase the required errors.",
        "scale_m_s2": scale,
        "training_interval_m_s2": list(training_interval),
        "quoted_velocity_error_threshold": tolerance,
        "geometry_and_calibration_sensitivity": limits,
        "nominal_predictions_kms": {
            "row_order": "Use the galaxy list consecutively, taking the listed number of radius and speed entries for each galaxy.",
            "galaxies": names,
            "ring_counts": sizes,
            "radius_kpc": radii,
            **predictions,
        },
    }


def ngc0024_projection_control(groups, scale, training_high, tolerance):
    """Check the bottleneck against the two morphology inclinations in Chemin et al."""
    galaxy, rings = next((g, r) for g, r in groups if g["name"] == "NGC0024")
    radius, observed, error, gas, disk = rings[:, np.argmax(rings[0])]
    distance = galaxy["dist_mpc"]
    projection = np.sqrt(1 + 2 * galaxy["e_dist_mpc"] / distance) / sin_deg(
        galaxy["inclination_deg"]
    )
    baryons = gas + 0.6 * disk
    trials = []
    for value in (scale, training_high):
        intrinsic = np.sqrt(
            predicted_speed_squared(galaxy, radius, baryons, 0.6, value, "local_force")
        )
        trials.append(
            {
                "scale_m_s2": value,
                "highest_catalogue_frame_speed_by_inclination_deg": {
                    str(angle): float(intrinsic * projection * sin_deg(angle))
                    for angle in (70, 76, 78)
                },
                "minimum_inclination_to_enter_velocity_band_deg": float(
                    np.rad2deg(
                        np.arcsin(
                            (observed - tolerance * error) / (intrinsic * projection)
                        )
                    )
                ),
            }
        )
    return {
        "source": "Chemin et al. (2006), section 5.1 and table 7, https://arxiv.org/abs/astro-ph/0609148",
        "HI_table_preprint_pdf_sha256": "56562396fbd9930ac3df95c0ef05114206b30992ed558e6ec13fe854be52100f",
        "antecedent_HI_record": {
            "angular_radius_arcsec": 320,
            "rotational_speed_kms": 109.7,
            "quoted_error_kms": 2.0,
            "kinematic_inclination_deg": 64,
            "kinematic_inclination_error_deg": 3,
            "HI_isophotal_inclination_deg": 76,
            "optical_inclination_deg": 78,
        },
        "catalogue_ring_radius_kpc": radius,
        "catalogue_angular_radius_arcsec": float(
            np.rad2deg(radius / (distance * 1000)) * 3600
        ),
        "lower_edge_of_speed_acceptance_band_kms": float(observed - tolerance * error),
        "model": "Projection sensitivity on the limiting ring only. Hold the tabulated angular baryonic-force profile fixed; maximize over its catalogue distance box and disk M/L. Use the published morphologies as independent angle estimates on the same target. A new joint velocity-field and baryonic-image reduction would have to redetermine the radii and deprojected material distributions.",
        "readings": trials,
    }


def alternate_geometry(groups, scale, tolerance):
    """Use the full table intersections with the enrolled galaxies; retain every ring."""
    # Intersections of Schombert et al. Table 2 (TRGB entries only) and Dicaire
    # et al. Table 4 (kinematic inclinations) with the 41 fixed catalogue names.
    distances = {
        "DDO154": (4.04, 0.15),
        "NGC0024": (7.67, 0.32),
        "NGC2915": (4.29, 0.20),
        "NGC3741": (3.23, 0.12),
        "NGC6503": (6.25, 0.29),
        "UGCA442": (4.37, 0.20),
    }
    inclinations = {"DDO154": (59.0, 32.0), "NGC0024": (75.0, 3.0)}
    prescribed = []
    summary = {}
    for mode in ("distance_only", "inclination_only", "both"):
        worst, violating_rings, solutions, failed_galaxies = 0.0, 0, [], []
        for galaxy, rings in groups:
            name = galaxy["name"]
            prior = dict(galaxy)
            if mode != "inclination_only" and name in distances:
                prior["dist_mpc"], prior["e_dist_mpc"] = distances[name]
            if mode != "distance_only" and name in inclinations:
                prior["inclination_deg"], prior["e_inclination_deg"] = inclinations[
                    name
                ]
            # Re-express the same fixed angular ring, flux and projected velocity
            # in the alternate central frame. Its velocity residual/error is invariant.
            ratio = prior["dist_mpc"] / galaxy["dist_mpc"]
            velocity = sin_deg(galaxy["inclination_deg"]) / sin_deg(
                prior["inclination_deg"]
            )
            shifted = rings.copy()
            shifted[0] *= ratio
            shifted[1:3] *= velocity
            shifted[3:] *= ratio
            prior["L36_gsun"] *= ratio**2
            prior["MHI_gsun"] *= ratio**2
            lower, upper = speed_envelopes(prior, shifted, 2, (scale, scale))[
                "local_force"
            ]
            miss = (
                np.maximum.reduce(
                    (shifted[1] - upper, lower - shifted[1], np.zeros(len(lower)))
                )
                / shifted[2]
            )
            worst = max(worst, max(miss))
            violating_rings += int(np.count_nonzero(miss > tolerance))
            witness = joint_witness(
                [(prior, shifted)], 2, tolerance, scale, "local_force"
            )
            if witness is None:
                failed_galaxies.append(name)
            else:
                solutions.extend(witness)
            if mode == "both" and name == "NGC0024":
                for d, inc, ml in (
                    (7.67, 75.0, 0.5),
                    (7.99, 78.0, 0.6),
                    (8.31, 78.0, 0.6),
                ):
                    vsq = predicted_speed_squared(
                        galaxy,
                        rings[0],
                        rings[3] + ml * rings[4],
                        ml,
                        scale,
                        "local_force",
                    )
                    pred = (
                        np.sqrt(vsq * d / galaxy["dist_mpc"])
                        * sin_deg(inc)
                        / sin_deg(galaxy["inclination_deg"])
                    )
                    prescribed.append(
                        {
                            "distance_mpc": d,
                            "inclination_deg": inc,
                            "disk_mass_to_light": ml,
                            "predictions_in_SPARC_frame_kms": pred.tolist(),
                            "worst_quoted_speed_errors": float(
                                np.max(abs(pred - rings[1]) / rings[2])
                            ),
                        }
                    )
        summary[mode] = {
            "rings_with_unavoidable_miss": violating_rings,
            "worst_ring_lower_bound": worst,
            "galaxies_without_shared_witness": failed_galaxies,
            "rings_in_solved_galaxies": sum(v["rings"] for v in solutions),
            "maximum_residual_over_solved_galaxies": max(
                v["maximum_absolute_velocity_residual_in_quoted_errors"]
                for v in solutions
            ),
        }
        if mode == "both":
            summary[mode]["matched_galaxy_witnesses"] = [
                v
                for v in solutions
                if v["galaxy"] in distances or v["galaxy"] in inclinations
            ]
    return {
        "distance_source": "Schombert, McGaugh & Lelli (2020), Table 2, TRGB entries; section 2.3 attributes distances and errors to Tully et al. (2019). https://doi.org/10.3847/1538-3881/ab9d88",
        "distance_table_journal_pdf": "https://astroweb.case.edu/ssm/papers/AJv160n71.pdf",
        "distance_table_pdf_sha256": "adada9659ecdf2625a2ba21ce729e46258a4d89994263070ece245875803d1d5",
        "inclination_source": "Dicaire et al. (2008), Table 4 and Appendix A, inner H-alpha tilted-ring inclinations. https://doi.org/10.1111/j.1365-2966.2008.12868.x",
        "inclination_table_preprint_pdf_sha256": "fb717555180ec414ea49f1afdeeb0dc908ef571a2cfa1174910cd5d3792f00b9",
        "application": "Retain the original outer sample, force profiles and angular bins. Replace a geometry input for every available match within the stated table and method. Otherwise retain the catalogue input. Use twice the tabulated distance/inclination errors and disk M/L between 0.4 and 0.6. These conditional projection tests use a disjoint fixed scale. Stellar distances are independent of the ring velocities; independent imaging between the distance compilations is unverified. The angle of an inner H-alpha field may differ from that of its outer H I disk. Surface-density profiles would have to be redetermined in a complete image and velocity-field reduction.",
        "TRGB_distance_Mpc_and_quoted_error": distances,
        "Halpha_inclination_deg_and_quoted_error": inclinations,
        "NGC0024_three_test_rings_in_ascending_radius": prescribed,
        "matched_catalogue_controls": summary,
    }


def scale_constraints(groups, width, k):
    """Necessary lower and upper bounds on sqrt(a0), relaxing ties between rings."""
    lower, upper = (0.0, ""), (np.inf, "")
    for galaxy, rings in groups:
        _, speed, error, fixed_baryons, disk = rings
        t_lo, t_hi, ml_lo, ml_hi = nuisance_bounds(galaxy, width)
        baryon_lo = fixed_baryons + np.minimum(ml_lo * disk, ml_hi * disk)
        baryon_hi = fixed_baryons + np.maximum(ml_lo * disk, ml_hi * disk)
        c_lo, c_hi = (mass_coefficient(galaxy, ml) for ml in (ml_lo, ml_hi))
        lo = (np.maximum(0.0, speed - k * error) ** 2 * t_lo - baryon_hi).max()
        hi = ((speed + k * error) ** 2 * t_hi - baryon_lo).min()
        lower = max(lower, (lo / c_hi, galaxy["name"]))
        upper = min(upper, (hi / (c_lo if hi >= 0 else c_hi), galaxy["name"]))
    return lower, upper


def scale_gap(groups, width, k):
    lo, hi = scale_constraints(groups, width, k)
    return lo[0] - hi[0]


def fixed_scale_gap(groups, width, k, a0):
    lo, hi = scale_constraints(groups, width, k)
    return max(lo[0] - np.sqrt(a0), np.sqrt(a0) - hi[0])


def minimum_error(gap):
    if gap(0) <= 0:
        return 0.0
    hi = 1.0
    while gap(hi) > 0:
        hi *= 2
        if hi > 1024:
            raise ValueError("No enclosing velocity-error bracket found")
    return brentq(gap, 0.0, hi)


def joint_witness(groups, width, k, a0, model="compact_mass"):
    """One physical D, inclination and disk M/L for all the rings of each galaxy."""
    results = []
    for galaxy, rings in groups:
        radius, speed, error, fixed_baryons, disk = rings
        t_min, t_max, ml_lo, ml_hi = nuisance_bounds(galaxy, width)
        floor = np.maximum(0.0, speed - k * error)
        for ml in np.linspace(ml_lo, ml_hi, 201):
            squared = predicted_speed_squared(
                galaxy, radius, fixed_baryons + ml * disk, ml, a0, model
            )
            lo = max(t_min, np.max(squared / (speed + k * error) ** 2))
            hi = min(
                t_max,
                np.min(squared[floor > 0] / floor[floor > 0] ** 2, initial=np.inf),
            )
            if lo <= hi:
                t = (lo + hi) / 2
                break
        else:
            return None
        i0 = galaxy["inclination_deg"]
        delta_i = width * galaxy["e_inclination_deg"]
        sin_sq = sin_deg(i0) ** 2
        d0, delta_d = galaxy["dist_mpc"], width * galaxy["e_dist_mpc"]
        d_lo = max(
            d0 - delta_d, d0 * sin_sq / (t * sin_deg(min(90.0, i0 + delta_i)) ** 2)
        )
        d_hi = d0 + delta_d
        if i0 > delta_i:
            d_hi = min(d_hi, d0 * sin_sq / (t * sin_deg(i0 - delta_i) ** 2))
        if d_lo > d_hi + 1e-12:
            raise ValueError("Projection has no common distance and inclination")
        distance = (d_lo + d_hi) / 2
        inclination = np.rad2deg(
            np.arcsin(np.sqrt(min(1.0, d0 * sin_sq / (distance * t))))
        )
        # Recover the measured line-of-sight velocity from the physical witness.
        measured_frame = (
            np.sqrt(squared * distance / d0) * sin_deg(inclination) / sin_deg(i0)
        )
        discrepancy = np.max(np.abs(measured_frame - speed) / error)
        if (
            discrepancy > k + 1e-9
            or abs(distance - d0) > delta_d + 1e-9
            or abs(inclination - i0) > delta_i + 1e-7
        ):
            raise ValueError("Shared galaxy parameters fail direct reprojection")
        results.append(
            {
                "galaxy": galaxy["name"],
                "rings": len(speed),
                "distance_mpc": distance,
                "inclination_deg": inclination,
                "disk_mass_to_light": ml,
                "maximum_absolute_velocity_residual_in_quoted_errors": discrepancy,
            }
        )
    return results


def monotone_mass_gap(groups, k):
    """Is a nonnegative nondecreasing M_extra(r) possible, at catalogue geometry?"""
    gaps = []
    for _, rings in groups:
        radius, speed, error, fixed_baryons, disk = rings
        baryons = fixed_baryons + base.UPSILON_DISK * disk
        lower = np.maximum(
            0.0, radius * (np.maximum(0.0, speed - k * error) ** 2 - baryons)
        )
        upper = radius * ((speed + k * error) ** 2 - baryons)
        # Equal radii share one mass. At successive radii the previous lower bound persists.
        _, first = np.unique(radius, return_index=True)
        prefix = np.maximum.accumulate(np.maximum.reduceat(lower, first))
        gaps.append(np.max(prefix - np.minimum.reduceat(upper, first)))
    return max(gaps)


def measure(groups, table1, table2):
    n = sum(len(rings[0]) for _, rings in groups)
    if not n or any(np.any(rings[2] <= 0) for _, rings in groups):
        raise ValueError("A nonempty sample and positive quoted errors are required")
    k = NormalDist().inv_cdf(1 - 0.05 / (2 * n))
    sensitivity = {}
    for width in (0, 1, 2):
        low, high = scale_constraints(groups, width, k)
        sensitivity[str(width)] = {
            "minimum_maximum_velocity_error_lower_bound": minimum_error(
                lambda z: scale_gap(groups, width, z)
            ),
            "required_min_a0_m_s2": low[0] ** 2,
            "largest_lower_bound_galaxy": low[1],
            "allowed_max_a0_m_s2": np.sign(high[0]) * high[0] ** 2,
            "smallest_upper_bound_galaxy": high[1],
        }
    low, high = scale_constraints(groups, 2, k)
    if low[0] > high[0]:
        raise ValueError("The archived two-error-box interval has closed")
    a0 = ((low[0] + high[0]) / 2) ** 2
    witness = joint_witness(groups, 2, k, a0)
    if witness is None:
        raise ValueError("The relaxed common scale lacks a shared-parameter witness")
    all_exterior = exterior_rings(table1, table2, include_bulges=True)
    return {
        "selection": "Bulge absent at every catalogue ring, quality <= 2, inclination >= 30 deg, measured positive disk scale length and HI isodensity radius, R >= max(4 Rdisk, RHI), nominal 0 < gbar < 1.2e-11 m/s^2. Includes every qualifying speed regardless of error or excess sign.",
        "model": "V^2 = Vgas|Vgas| + Upsilon Vdisk|Vdisk| + sqrt(G Mphot a0), in consistent units. Mphot = Upsilon L_3.6 + 1.33 M_HI. Angular geometry and photometric flux held fixed as D varies. One common a0. RHI is the 1 solar mass/pc^2 crossing. An exterior compact source and settled spherical cuts are additional hypotheses.",
        "geometry_boxes": "Widths multiply quoted distance and inclination errors. Width 0 uses Upsilon=0.5. Nonzero widths use [0.4,0.6]. Necessary intervals release shared parameters by ring. Witness uses one distance, inclination and mass-to-light ratio per galaxy.",
        "simultaneous_velocity_threshold": "Bonferroni two-sided 95% threshold conditional on marginal Gaussian catalogue speed errors, selection and nuisance boxes. Includes no joint probability model for geometry, baryonic systematics or noncircular motions.",
        "galaxies": len(groups),
        "rings": n,
        "quoted_velocity_error_threshold": k,
        "sensitivity_by_nuisance_box_width": sensitivity,
        "two_error_box_constructive_scale_m_s2": a0,
        "two_error_box_shared_parameter_witness": witness,
        "nonnegative_nondecreasing_enclosed_mass_control": {
            "method": "At nominal geometry permit any nonnegative nondecreasing spherical mass per galaxy. At each ordered radius R, intersect the interval R[(V +/- k e)^2 - Vbar^2] with all preceding lower bounds and zero. Truncate the lower velocity endpoint at zero. Duplicate radii must share a mass. Fold fixed bulge M/L 0.7 into Vbar for the extended sample.",
            "compact_sample_minimum_maximum_velocity_error": minimum_error(
                lambda z: monotone_mass_gap(groups, z)
            ),
            "including_bulges_galaxies": len(all_exterior),
            "including_bulges_rings": sum(len(rings[0]) for _, rings in all_exterior),
            "including_bulges_minimum_maximum_velocity_error": minimum_error(
                lambda z: monotone_mass_gap(all_exterior, z)
            ),
        },
    }
