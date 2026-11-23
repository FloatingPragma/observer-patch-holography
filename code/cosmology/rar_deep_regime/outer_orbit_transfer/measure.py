"""Test a shared exterior mass profile with H I orbit intervals and weak shear."""

import argparse
import hashlib
import json
import math
from pathlib import Path
from statistics import NormalDist

import numpy as np

from .. import lensing_profile
from ..finite_halo.transfer import profile_surface_bound
from ..udg_transfer.measure import GM_SUN, KPC

HERE = Path(__file__).resolve().parent


def enclosed_fraction(row):
    star = row["stellar_unit_msun"]
    return (
        0.6
        * (row["stellar_inside_unit_msun"] - 0.01 * star)
        / (0.6 * star + 10**0.10 * row["helium_hi_msun"])
    )


def orbit(row, speed_errors=0.0, generous=False, subtract_stars=True):
    """Upper corner of the nuisance box; a negative speed multiplier gives a lower central interval."""
    m_to_l = 0.6 if generous else 0.8
    gas_scale = 10**-0.10 if generous else 1.0
    mass = m_to_l * row["stellar_unit_msun"] + gas_scale * row["helium_hi_msun"]
    rad = KPC * row["radius_kpc"]
    g = GM_SUN * mass / rad**2
    v = max(0.0, row["speed_kms"] + speed_errors * row["speed_error_kms"])
    scale = 1.0
    if generous:
        inc = math.radians(row["inc_deg"])
        low_inc = math.radians(row["inc_deg"] - 2 * row["inc_error_deg"])
        scale -= 2 * row["distance_error_mpc"] / row["distance_mpc"]
        if low_inc <= 0 or scale <= 0:
            raise ValueError("Unbounded distance or inclination box")
        v *= math.sin(inc) / math.sin(low_inc)
    vsq = (
        m_to_l * (row["vdisk_unit_kms"] ** 2 + row["vbul_unit_kms"] ** 2)
        if subtract_stars
        else 0.0
    )
    ratio = (v**2 / scale - vsq) * 1e6 / (rad * g)
    return g, ratio


def shape_feasibility(rows, k):
    bounds = {}
    for row in rows:
        g, upper = orbit(row, k)
        _, lower = orbit(row, -k)
        bounds[row["name"]] = dict(
            x=1 / math.sqrt(g), lower=max(0.0, lower), upper=upper
        )
    gap, supplier, recipient = max(
        (a["lower"] * min(1.0, b["x"] / a["x"]) - b["upper"], i, j)
        for i, a in bounds.items()
        for j, b in bounds.items()
    )
    return dict(
        worst_excess=float(gap),
        lower_source=supplier,
        upper_source=recipient,
        knots=bounds,
    )


def min_errors(rows):
    hi = 1.0
    while shape_feasibility(rows, hi)["worst_excess"] > 0:
        hi *= 2
    lo = 0.0
    if shape_feasibility(rows, 0.0)["worst_excess"] <= 0:
        hi = 0.0
    for _ in range(60):
        mid = (hi + lo) / 2
        if shape_feasibility(rows, mid)["worst_excess"] > 0:
            lo = mid
        else:
            hi = mid
    return dict(
        min_speed_error_multiplier=hi,
        witness=shape_feasibility(rows, hi),
        one_speed_error=shape_feasibility(rows, 1.0),
    )


def compare(profile, weights, caps, mass_factor):
    low = np.array(profile["lower_bin_edge_m_s2"])
    select = low >= 1e-13
    upper = mass_factor * (2 * np.array(profile["g_bar_m_s2"]) - low)[select]
    signal = np.asarray(profile["tangential_esd"])[select]
    cov = np.asarray(profile["covariance"])[np.ix_(select, select)]
    names = list(caps)
    bounds = np.array(
        [profile_surface_bound(upper, g, max(ratio, 0.0)) for g, ratio in caps.values()]
    )

    def reduce(mask):
        if not mask.any():
            return dict(bins=0)
        w = weights[mask] / weights[mask].sum()
        b = bounds[:, mask]
        sigma = float(np.sqrt(w @ cov[np.ix_(mask, mask)] @ w))
        value = float(w @ signal[mask])
        cap = float(w @ b.min(axis=0))
        return dict(
            bins=int(mask.sum()),
            largest_bin_acceleration_m_s2=upper[mask].tolist(),
            upper_by_bin=b.min(axis=0).tolist(),
            tightest_by_bin=[names[i] for i in b.argmin(axis=0)],
            weighted_upper=cap,
            observed=value,
            statistical_se=sigma,
            gap_in_statistical_errors=(value - cap) / sigma,
            leave_one_out={
                n: float(w @ np.delete(b, i, axis=0).min(axis=0))
                for i, n in enumerate(names)
            }
            if len(names) > 1
            else {},
        )

    return dict(
        all_bins=reduce(np.ones(len(upper), dtype=bool)),
        exterior_scaled_bins=reduce(upper <= min(g for g, _ in caps.values())),
    )


def run():
    ip = HERE / "inputs.json"
    wp = HERE.parent / "finite_halo/results.json"
    rows = json.loads(ip.read_text())["galaxies"]
    primary = json.loads(lensing_profile.RECEIPT.read_text())
    surveys = json.loads(wp.read_text())["surveys"]
    enroll = {}
    for row in rows:
        name = row["name"]
        f = enclosed_fraction(row)
        re = row["radius_kpc"] / row["re_kpc"]
        g = orbit(row)[0]
        enroll[name] = dict(
            enclosed_fraction_lower=f,
            radius_in_re=re,
            central_gbar_m_s2=g,
            radial_pass=bool(re >= 5),
            mass90_pass=bool(f >= 0.9),
            mass80_pass=bool(f >= 0.8),
            deep=bool(g <= 1.2e-11),
        )
    panels = {}
    for label, cut, deep in (
        ("primary", 0.90, False),
        ("secondary", 0.80, False),
        ("primary_deep", 0.90, True),
        ("secondary_deep", 0.80, True),
    ):
        selected = [
            r
            for r in rows
            if enroll[r["name"]]["radial_pass"]
            and enclosed_fraction(r) >= cut
            and (not deep or enroll[r["name"]]["deep"])
        ]
        cases = {}
        if selected:
            bonf = NormalDist().inv_cdf(1 - 0.05 / len(selected))
            for name, k, wide, subtract in (
                ("central", 0, False, True),
                ("central_all_speed", 0, False, False),
                ("two_error_box", 2, True, True),
                ("two_error_all_speed", 2, True, False),
                ("bonferroni_box", bonf, True, True),
                ("bonferroni_all_speed", bonf, True, False),
            ):
                caps = {r["name"]: orbit(r, k, wide, subtract) for r in selected}
                mass_factor = 10**0.2 if wide else 1.0
                cases[name] = dict(
                    speed_error_multiplier=k,
                    lens_mass_factor=mass_factor,
                    negative_force_residuals=[n for n, (_, u) in caps.items() if u < 0],
                    surveys={
                        s: compare(
                            primary["profiles"][s],
                            np.array(w["positive_weights"]),
                            caps,
                            mass_factor,
                        )
                        for s, w in surveys.items()
                    },
                )
        panels[label] = dict(
            names=[r["name"] for r in selected],
            cases=cases,
            orbital_shape=min_errors(selected) if selected else None,
        )
    return dict(
        input_sha256=hashlib.sha256(ip.read_bytes()).hexdigest(),
        weight_receipt_sha256=hashlib.sha256(wp.read_bytes()).hexdigest(),
        lens_archive_sha256=primary["archive_sha256"],
        enrollment=enroll,
        panels=panels,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run()
    with args.output.open("x") as out:
        json.dump(result, out, indent=2, allow_nan=False)
        out.write("\n")
    for key, panel in result["panels"].items():
        print(key, panel["names"])
        if panel["names"]:
            print(
                "orbit multiplier",
                panel["orbital_shape"]["min_speed_error_multiplier"],
                panel["orbital_shape"]["witness"]["lower_source"],
                panel["orbital_shape"]["witness"]["upper_source"],
            )
            for label, case in panel["cases"].items():
                print(
                    label,
                    [
                        round(s["all_bins"]["weighted_upper"], 2)
                        for s in case["surveys"].values()
                    ],
                )
