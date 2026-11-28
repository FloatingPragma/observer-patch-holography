"""Direct stellar distance indicators for the frozen pair and analytic rebin bounds."""

import json

import numpy as np

from .measure import run as primary
from .source import AREA_UNITS, HERE, edges, enroll, filled, profile

PAIR = ("NGC3521", "NGC3198")


def coarse_tv_envelope(shells, sigma, width, dmin, dmax):
    """Conservative radial-TV interval for any distance multiplier in the range."""
    if min(sigma) < 0 or max(sigma) > 10_000:
        raise ValueError("These bounds require no threshold saturation")
    jumps = -np.diff(np.r_[sigma, 0])
    # Projecting a filled disk of radius r onto a coarse annular grid replaces its
    # perimeter radius by r - (r-a)(b-r)/(a+b) for the shell [a,b] containing r.
    min_r = shells * dmin
    max_loss = width**2 / (4 * np.maximum(width, 2 * min_r - width))
    signed_lower = dmin * (shells @ jumps) - np.maximum(jumps, 0) @ max_loss
    upper = dmax * (shells @ np.abs(jumps))
    return 2 * np.pi * AREA_UNITS * np.array([max(0, signed_lower), upper])


def evaluate(n, kept, target, distance, dist_bounds, error_multiplier):
    meta, lrow, rows, r, _ = kept[n]
    d0 = meta["dist_mpc"]
    drel = distance / d0
    dlo, dhi = np.array(dist_bounds) / d0
    alt = profile(meta, lrow, rows, r * drel, distance=distance)
    radii = np.array([float(x["r"]) for x in rows]) * d0 / float(lrow["Dist"])
    sigma = sum(filled(radii, rows, c)[0] for c in ("Sigma*", "SigmaHI", "SigmaH2"))
    slo, shi = coarse_tv_envelope(edges(radii)[1:], sigma, 0.8, dlo, dhi)
    v, ev = target["observed_kms"], error_multiplier * target["e_kms"]
    i = np.deg2rad(meta["inclination_deg"])
    ei = np.deg2rad(error_multiplier * meta["e_inclination_deg"])
    ys = (
        (v + np.array([-ev, ev])) * np.sin(i) / np.sin([i + ei, i - ei])
    ) ** 2 - np.array([dhi, dlo]) * target["baryon_v2"]
    if ys[0] <= 0 or slo <= 0 or i + ei >= np.pi / 2:
        raise ValueError("Positive-excess sub-edge-on control required")
    y = v * v - drel * target["baryon_v2"]
    total_hi = 1.36e9 * meta["MHI_gsun"] * drel**2
    return {
        "distance_mpc": distance,
        "distance_range": dist_bounds,
        "orbit_kpc": r * drel,
        "edge_kpc": alt["edge_kpc"],
        "mass_msun": alt["mass_msun"],
        "S": alt["S_msun_kpc"],
        "F": alt["F"],
        "Y": y,
        "S_envelope": [slo, shi],
        "Y_box": ys.tolist(),
        "required_gain": y / alt["S_msun_kpc"],
        "gain_envelope": [ys[0] / shi, ys[1] / slo],
        "catalogue_minus_profile_atomic_He_msun": total_hi - alt["atomic_with_he_msun"],
        "all_catalogue_stars_at_0p7_plus_gas_minus_measured_stars_atomic_msun": 0.7e9
        * meta["L36_gsun"]
        * drel**2
        + total_hi
        - alt["stellar_msun"]
        - alt["atomic_with_he_msun"],
    }


def run():
    kept, _ = enroll()
    targets = primary()["targets"]
    result = {}
    # Table 4 of Freedman et al. 2001 lists 0.08 mag random error. Add the 0.10
    # mag LMC-anchor error linearly in the primary box. Double both in the stress box.
    for mult in (1, 2):
        d = 13.8
        dist = {
            "NGC3521": (11.3, [11.3 - mult * 1.1, 11.3 + mult * 1.1]),
            "NGC3198": (d, (d * 10 ** (np.array([-1, 1]) * 0.18 * mult / 5)).tolist()),
        }
        a, b = [evaluate(n, kept, targets[n], *dist[n], mult) for n in PAIR]
        added = max(0, a["S"] * b["Y"] / a["Y"] - b["S"])
        result[str(mult)] = {
            "galaxies": dict(zip(PAIR, (a, b), strict=True)),
            "S_ratio": a["S"] / b["S"],
            "sqrt_mass_ratio": np.sqrt(a["mass_msun"] / b["mass_msun"]),
            "SF_ratio": a["S"] * a["F"] / (b["S"] * b["F"]),
            "Y_ratio": a["Y"] / b["Y"],
            "required_gain_ratio_envelope": [
                a["gain_envelope"][0] / b["gain_envelope"][1],
                a["gain_envelope"][1] / b["gain_envelope"][0],
            ],
            "central_missing_S_required_for_NGC3198_if_no_added_NGC3521": added,
            "central_minimum_monotone_tail_mass_msun": added * b["edge_kpc"],
            "central_minimum_arbitrary_annuli_mass_at_0p8kpc_msun": added * 0.8 / 2,
        }
    return result


if __name__ == "__main__":
    (HERE / "distance.json").write_text(
        json.dumps(run(), indent=2, allow_nan=False) + "\n"
    )
