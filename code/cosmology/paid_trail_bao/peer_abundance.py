"""Propagate a finite peer-agreement source interval through published galaxy images."""

import argparse
from fractions import Fraction as Q
from hashlib import sha256
from itertools import product
import json

import numpy as np

from . import pooled_fluctuations as fluctuations
from . import pooled_ordinary as pooled


def run():
    lab = pooled.scales.budget.LAB
    cert = lab / "peer_concentration.json"
    bound = json.loads(cert.read_bytes())
    a, d = (
        [float(Q(v)) for v in bound[key]]
        for key in ("ordinary_per_port_bounds", "single_count_change_bounds")
    )
    corners = [(u, v, -v) for u, v in product(a, d)]
    assert d[1] < 0 < a[0] and bound["joint_acceptance_upper_power_of_two"] < -1000
    peaks_path = pooled.HERE / "record_profile_maps.json"
    peaks = json.loads(peaks_path.read_bytes())["peak_msun_pc2"]
    lib = pooled.sky.neighbor.inventory.rp.galaxy.transfer.budget
    G = lib.local.profiles.galaxies.G
    _, mu = pooled.scales.source_data()
    mu = float(mu)
    inputs, *_ = lib.load_inputs()
    scores = {}
    profile = {}
    cosmic = {}
    passive = {}
    passive_upper = 0.0
    # Compute fields at all four corners. Mass, radial field and squared
    # speed are affine in the two calibration coefficients at fixed occupancy.
    with np.load(pooled.sky.MAPS) as raw:
        for name, tab in inputs.items():
            r, v, e, b = (
                np.array(tab[k])
                for k in ("radius_kpc", "observed_kms", "error_kms", "baryons_kms")
            )
            low = np.maximum(0, v - 2 * e) ** 2
            high = (v + 2 * e) ** 2
            passive[name] = int(((b * b + 1e-8 >= low) & (b * b - 1e-8 <= high)).sum())
            for tag, peak in peaks.items():
                radii = raw[f"{name}_{tag}_radius_kpc"]
                unique, rr = np.unique(radii, return_inverse=True)
                rings = np.searchsorted(r, radii)
                for cap in pooled.sky.neighbor.inventory.rp.CAPACITY:
                    p = raw[f"{name}_{tag}"] / cap
                    w = 6 * p**2 * (1 - p) ** 2 * (448 / 2110)
                    L = len(p)
                    conversion = pooled.WIDTH**2 * 1e6 * peak * cap * mu / 4
                    # At most L eligible sites. Each successful site can
                    # produce at most this many donors even in the largest pool.
                    mass_limit = (
                        2.0**-1000 * L * (48 * a[1] * L + abs(d[0])) * conversion
                    )
                    gain = G * max(1 / r.min(), r.max() / (pooled.WIDTH / 2) ** 2)
                    passive_upper = max(passive_upper, mass_limit * gain)
                    for pool, group in (
                        ("tile", np.arange(L)),
                        ("annulus", rings),
                        ("disk", np.zeros(L, int)),
                    ):
                        fields = [
                            pooled.field_from_tiles(
                                r,
                                rings,
                                rr,
                                unique,
                                conversion * pooled.expected_per_tile(w, group, coef),
                                G,
                            )[1]
                            for coef in corners
                        ]
                        for mode in fields[0]:
                            values = np.array([row[mode] for row in fields])
                            # Outward allowances for floating disk quadrature
                            # and arithmetic at the cached image resolution.
                            pad = 1e-8 + 1e-12 * np.abs(values).max(0)
                            lo = b * b + values.min(0) - pad
                            hi = b * b + values.max(0) + pad
                            possible = (hi >= low) & (lo <= high)
                            key = f"{pool}/{tag}/{cap}/{mode}"
                            score = scores.setdefault(
                                key, dict(possible_points=0, possible_complete=[])
                            )
                            score["possible_points"] += int(possible.sum())
                            if possible.all():
                                score["possible_complete"].append(name)
                            if tag == "mapped_baryons" and cap == 1:
                                profile.setdefault(f"{pool}/{mode}", {})[name] = dict(
                                    radius_kpc=r.tolist(),
                                    extra_v2_low=(lo - b * b).tolist(),
                                    extra_v2_high=(hi - b * b).tolist(),
                                )
    # Use the same BBN baryon density and homogeneous pool volumes as the
    # ordinary calibration. Each trial releases one member per eligible site.
    rho = 3 * 0.1**2 / (8 * np.pi * G) * 0.02235
    for tag, peak in peaks.items():
        for cap in pooled.sky.neighbor.inventory.rp.CAPACITY:
            p = rho * pooled.WIDTH / (1e6 * peak * cap)
            for n in (1, 5, 20, 80):
                ratios = [pooled.homogeneous(n, p, c, mu) for c in corners]
                cosmic[f"{tag}/{cap}/{n}"] = [min(ratios), max(ratios)]
    assert passive_upper < 1e-8
    material = lab / "peer_register_lab.json"
    measured = json.loads(material.read_bytes())
    coeff = [float(Q(x)) for x in measured["ordinary_coefficients"]]
    finite = {}
    for label, rate in (
        ("one_register", Q(measured["passive_probability"])),
        ("successful_preparation", Q(1)),
    ):
        sky = pooled.run(material, float(rate))
        row = dict(
            acceptance=str(rate),
            scores=sky["scores"],
            homogeneous=sky["homogeneous"],
            mapped_baryons=sky["mapped_baryons_first_capacity"],
        )
        if rate == 1:
            row["matter_fluctuations"] = fluctuations.forecasts(
                sky["homogeneous"], coeff, mu
            )
        finite[label] = row
    return dict(
        sha256={
            f.name: sha256(f.read_bytes()).hexdigest()
            for f in (
                cert,
                material,
                peaks_path,
                pooled.sky.MAPS,
                pooled.HERE / "source_profile.json",
            )
        },
        register_896=finite,
        renewed_scores=scores,
        renewed_homogeneous_ratios=cosmic,
        renewed_mapped_baryons=profile,
        passive_max_abs_extra_v2=passive_upper,
        passive_possible_points_by_galaxy=passive,
        point_count=sum(len(t["radius_kpc"]) for t in inputs.values()),
        planck_base_lcdm_ratio=0.12 / 0.02237,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    out = run()
    path = pooled.HERE / "peer_abundance.json"
    if args.write:
        path.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
    else:
        pooled.sky.neighbor.inventory.rp.galaxy.check(
            json.loads(path.read_bytes()), out
        )
    print(
        max(v["possible_points"] for v in out["renewed_scores"].values()),
        sum(out["passive_possible_points_by_galaxy"].values()),
        max(v[-1] for v in out["renewed_homogeneous_ratios"].values()),
    )
