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
    if fraction >= 1 or i_lo <= 0 or inclination > 90:
        raise ValueError("Unbounded distance or deprojection in the selected sample")
    sin0_sq = sin_deg(inclination) ** 2
    t_lo = sin0_sq / (sin_deg(i_hi) ** 2 * (1 + fraction))
    t_hi = sin0_sq / (sin_deg(i_lo) ** 2 * (1 - fraction))
    offset = 0.1 if width else 0.0
    return t_lo, t_hi, 0.5 - offset, 0.5 + offset


def mass_coefficient(galaxy, mass_to_light):
    mass = mass_to_light * galaxy["L36_gsun"] + base.GAS_HELIUM * galaxy["MHI_gsun"]
    return np.sqrt(base.G_SI * base.M_SUN * 1e9 * mass) / base.KM**2


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


def joint_witness(groups, width, k, a0):
    """One physical D, inclination and disk M/L for all the rings of each galaxy."""
    results = []
    for galaxy, rings in groups:
        _, speed, error, fixed_baryons, disk = rings
        t_min, t_max, ml_lo, ml_hi = nuisance_bounds(galaxy, width)
        floor = np.maximum(0.0, speed - k * error)
        for ml in np.linspace(ml_lo, ml_hi, 201):
            squared = (
                fixed_baryons + ml * disk + mass_coefficient(galaxy, ml) * np.sqrt(a0)
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
        d_hi = min(d0 + delta_d, d0 * sin_sq / (t * sin_deg(i0 - delta_i) ** 2))
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
