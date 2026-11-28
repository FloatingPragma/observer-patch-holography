"""Predict orbital speeds from the ordinary relay at fixed return reporting rates."""

import argparse
from hashlib import sha256
import json

import numpy as np

from . import pooled_ordinary as pooled

HERE = pooled.HERE


def run():
    lab = pooled.scales.budget.LAB / "sector_balls.json"
    src = json.loads(lab.read_bytes())["first_return_rates"]
    libs = pooled.sky.neighbor.inventory.rp
    data = libs.galaxy.transfer.budget
    G = data.local.profiles.galaxies.G
    galaxies, *_ = data.load_inputs()
    _, mass = pooled.scales.source_data()
    peak_path = HERE / "record_profile_maps.json"
    peaks = json.loads(peak_path.read_bytes())["peak_msun_pc2"]
    coefficients = {
        rate: np.array(t["ordinary_relay"]["coefficient_bands"])
        for rate, t in src.items()
    }
    scores = {k: {} for k in coefficients}
    profiles = {k: {} for k in coefficients}
    complete_curves = {k: {} for k in coefficients}
    cosmic = {k: [] for k in coefficients}
    with np.load(pooled.sky.MAPS) as maps:
        for name, tab in galaxies.items():
            radii, vel, err, baryon = [
                np.array(tab[k])
                for k in ("radius_kpc", "observed_kms", "error_kms", "baryons_kms")
            ]
            low = np.maximum(0, vel - 2 * err) ** 2
            high = (vel + 2 * err) ** 2
            for tracer, peak in peaks.items():
                radial = maps[f"{name}_{tracer}_radius_kpc"]
                rings = np.searchsorted(radii, radial)
                unique, ring_ids = np.unique(radial, return_inverse=True)
                image = maps[f"{name}_{tracer}"]
                for cap in libs.CAPACITY:
                    p = image / cap
                    w = 6 * p**2 * (1 - p) ** 2 * 448 / 2110
                    scale = pooled.WIDTH**2 * 1e6 * peak * cap / 4 * float(mass)
                    for pool, groups in (
                        ("tile", np.arange(len(p))),
                        ("annulus", rings),
                        ("disk", np.zeros(len(p), int)),
                    ):
                        # All three mass and force kernels are linear in the
                        # independently calibrated donor coefficients.
                        basis = [
                            pooled.field_from_tiles(
                                radii,
                                rings,
                                ring_ids,
                                unique,
                                scale * pooled.expected_per_tile(w, groups, unit),
                                G,
                            )[1]
                            for unit in np.eye(3)
                        ]
                        for mode in basis[0]:
                            kernel = np.array([f[mode] for f in basis])
                            for rate, band in coefficients.items():
                                mid = band.mean(1) @ kernel
                                width = (
                                    np.ptp(band, axis=1) / 2 @ abs(kernel)
                                    + 1e-8
                                    + 1e-12 * abs(kernel).sum(0)
                                )
                                predicted = baryon**2 + mid
                                hit = (predicted + width >= low) & (
                                    predicted - width <= high
                                )
                                tag = f"{pool}/{tracer}/{cap}/{mode}"
                                s = scores[rate].setdefault(
                                    tag, dict(possible=0, complete=[])
                                )
                                s["possible"] += int(hit.sum())
                                if hit.all():
                                    s["complete"].append(name)
                                    complete_curves[rate].setdefault(tag, {})[name] = (
                                        dict(
                                            radius_kpc=radii.tolist(),
                                            observed_speed_kms=vel.tolist(),
                                            observed_error_kms=err.tolist(),
                                            predicted_speed_kms=np.sqrt(
                                                predicted
                                            ).tolist(),
                                            model_v2_error=width.tolist(),
                                        )
                                    )
                                if (
                                    pool == "tile"
                                    and tracer == "mapped_baryons"
                                    and cap == 1
                                    and mode == "spherical"
                                ):
                                    profiles[rate][name] = dict(
                                        radius_kpc=radii.tolist(),
                                        predicted_speed_kms=np.sqrt(
                                            np.maximum(0, predicted)
                                        ).tolist(),
                                        model_v2_error=width.tolist(),
                                    )
    rho = 3 * 0.1**2 / (8 * np.pi * G) * 0.02235
    for peak in peaks.values():
        for cap in libs.CAPACITY:
            filling = rho * pooled.WIDTH / (1e6 * peak * cap)
            for tiles in (1, 5, 20, 80):
                kernel = np.array(
                    [
                        pooled.homogeneous(tiles, filling, u, float(mass))
                        for u in np.eye(3)
                    ]
                )
                for rate, band in coefficients.items():
                    cosmic[rate].append((kernel @ band).tolist())
    sign_only = {}
    for galaxy, t in galaxies.items():
        obs = np.array(t["observed_kms"])
        b = np.array(t["baryons_kms"])
        errs = np.array(t["error_kms"])
        sign_only[galaxy] = dict(
            total=len(obs), requires_positive=int(np.count_nonzero(obs - 2 * errs > b))
        )
    return dict(
        six_marker_signed_double=dict(
            rule="a local six-marker area decrease has the sign opposite to moving a double into the home",
            spherical_speed_needs=sign_only,
            requires_positive=sum(z["requires_positive"] for z in sign_only.values()),
        ),
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (lab, pooled.sky.MAPS, peak_path, HERE / "source_profile.json")
        },
        assumptions="one successful preparation at every eligible site, one common reporting probability per forecast",
        speeds=sum(len(v["radius_kpc"]) for v in galaxies.values()),
        galaxies=len(galaxies),
        base_lcdm_density_ratio=0.1200 / 0.02237,
        forecasts={
            k: dict(
                scores=v,
                complete_profiles=complete_curves[k],
                mapped_baryons_tile=profiles[k],
                cosmic_ratio_1_to_80_tiles=[
                    float(np.min(cosmic[k])),
                    float(np.max(cosmic[k])),
                ],
            )
            for k, v in scores.items()
        },
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    out = run()
    path = HERE / "sector_geometry.json"
    if args.write:
        path.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
    else:
        pooled.sky.neighbor.inventory.rp.galaxy.check(json.loads(path.read_text()), out)
    print(
        {
            k: (
                max(v["possible"] for v in r["scores"].values()),
                r["cosmic_ratio_1_to_80_tiles"],
            )
            for k, r in out["forecasts"].items()
        }
    )
