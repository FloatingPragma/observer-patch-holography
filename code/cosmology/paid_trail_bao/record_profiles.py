"""Transfer a measured two-head repair law to held galactic apertures."""

import argparse
from hashlib import sha256
import json
from math import comb
from pathlib import Path

import numpy as np

from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture import (
    record_density as source,
)
from . import pair_galaxies as galaxy

HERE = galaxy.HERE
maps = galaxy.maps
CAPACITY = (1, 2, 4, 8, 16, 64, 256)
MENUS = ("none", "third_head")
TRACERS = ("atomic", "mapped_baryons")


def source_laws(filling, bank, *, per_joint_change=None):
    """Per-attempt connected logarithm or an independently calibrated linear receiver."""
    p = np.asarray(filling)
    assert np.all((0 <= p) & (p <= 1))
    interior = (0 < p) & (p < 1)
    out = {m: np.zeros_like(p, dtype=float) for m in MENUS}
    if not np.any(interior):
        return out
    z = p[interior]
    weights = np.array(
        [
            comb(4, k) * z**k * (1 - z) ** (4 - k) / row["words"]
            for k, row in enumerate(bank["classes"])
        ]
    )
    a, b, c, d = np.array([v["magnitude_counts"] for v in bank["classes"]]).T @ weights
    log_odds = np.log(a) + np.log(d) - np.log(b) - np.log(c)
    delta = (b * c - a * d) / (a + b + c + d)
    response = -log_odds if per_joint_change is None else per_joint_change
    out["none"][interior] = delta * response
    witness = bank["third_head_counts"]
    total = np.zeros_like(z)
    for label in sorted(set().union(*(r.keys() for r in witness))):
        a, b, c, d = np.array([r.get(label, [0] * 4) for r in witness]).T @ weights
        total += (b * c - a * d) / (a + b + c + d)
    out["third_head"][interior] = total * response
    return out


def observed_cells(cache, inputs):
    enrolled, _ = maps.enroll()
    prior = json.loads((HERE / "pair_galaxy_maps.json").read_bytes())["galaxies"]
    for name, row in inputs.items():
        meta, survey, stellar, _, _ = enrolled[name]
        grid, info = maps.raw_map(name, meta["dist_mpc"], cache)
        assert info["raw_fits_sha256"] == prior[name]["raw_sha256"]
        rings, _ = maps.ring_surface(grid, survey, stellar, meta["dist_mpc"])
        hi = grid / (maps.FINE**2 * maps.AREA_UNITS)
        for tag, tracer in zip(TRACERS, (hi, hi + rings), strict=True):
            img = maps.coarse(tracer, round(galaxy.WIDTH / maps.FINE), 0)
            axis = (np.arange(len(img)) + 0.5 - len(img) // 2) * galaxy.WIDTH
            shell = np.searchsorted(
                row["radius_kpc"], np.hypot(axis[:, None], axis[None, :])
            )
            selected = (shell < len(row["radius_kpc"])) & (img > 0)
            yield name, tag, img, shell, selected


def acquire(cache, bank):
    inputs, _, _, _ = galaxy.transfer.budget.load_inputs()
    peak = dict.fromkeys(TRACERS, 0.0)
    measured = {name: {} for name in inputs}
    for name, tag, img, shell, selected in observed_cells(cache, inputs):
        measured[name][tag] = (shell[selected], img[selected])
        peak[tag] = max(peak[tag], float(img[selected].max()))
    profiles = {}
    access = {}
    for tag in TRACERS:
        profiles[tag] = {}
        access[tag] = {}
        for cap in CAPACITY:
            curves = {m: {} for m in MENUS}
            for name, row in inputs.items():
                shell, density = measured[name][tag]
                laws = source_laws(density / (cap * peak[tag]), bank)
                if name == "NGC2841":
                    endpoints = np.array(list(laws.values()))
                    access[tag][str(cap)] = {
                        key: (
                            np.bincount(shell, value, minlength=len(row["radius_kpc"]))
                            * galaxy.WIDTH**2
                            * maps.AREA_UNITS
                        ).tolist()
                        for key, value in (
                            ("lower", endpoints.min(axis=0)),
                            ("upper", endpoints.max(axis=0)),
                        )
                    }
                for menu, values in laws.items():
                    sums = np.bincount(shell, values, minlength=len(row["radius_kpc"]))
                    curves[menu][name] = (
                        sums.cumsum() * galaxy.WIDTH**2 * maps.AREA_UNITS
                    ).tolist()
            profiles[tag][str(cap)] = curves
    return dict(
        peak_msun_pc2=peak,
        unit="nat pc^2 per attempted patch",
        enclosed=profiles,
        ngc2841_shell_access=access,
    )


def scale_scores(inputs, fields):
    target = {
        n: np.asarray(v["observed_kms"]) ** 2 - np.asarray(v["baryons_kms"]) ** 2
        for n, v in inputs.items()
    }
    folds = {}
    errors = []
    for held, row in inputs.items():
        train = sorted(set(inputs) - {held})
        gain = max(
            0.0,
            sum(np.mean(fields[n] * target[n]) for n in train)
            / sum(np.mean(fields[n] ** 2) for n in train),
        )
        f = fields[held]
        v2 = np.asarray(row["baryons_kms"]) ** 2 + gain * f
        valid = v2 >= 0
        diff = np.full(len(v2), np.inf)
        diff[valid] = abs(np.sqrt(v2[valid]) - np.asarray(row["observed_kms"])[valid])
        errors.extend(diff)
        lo, hi = galaxy.transfer.budget.bounds(row)
        lower, upper = 0.0, np.inf
        for value, a, b in zip(f, lo, hi, strict=True):
            if value != 0:
                x, y = sorted((a / value, b / value))
                lower, upper = max(lower, x), min(upper, y)
            elif a > 0 or b < 0:
                upper = -np.inf
        folds[held] = dict(
            gain_msun_per_nat_pc2=float(gain),
            within_two=int(np.sum(diff <= 2 * np.asarray(row["error_kms"]))),
            points=len(v2),
            undefined=int(np.sum(~valid)),
            any_gain=bool(lower <= upper),
        )
    median = np.median(errors)
    return dict(
        galaxies=folds,
        rings=len(errors),
        within_two=sum(x["within_two"] for x in folds.values()),
        undefined=sum(x["undefined"] for x in folds.values()),
        held_full_curves=sum(x["within_two"] == x["points"] for x in folds.values()),
        possible_full_curves=sum(x["any_gain"] for x in folds.values()),
        median_abs_kms=float(median) if np.isfinite(median) else None,
    )


def run(profiles):
    inputs, _, _, _ = galaxy.transfer.budget.load_inputs()
    G = galaxy.transfer.budget.local.profiles.galaxies.G
    tests = {}
    for tag, scales in profiles["enclosed"].items():
        tests[tag] = {}
        for cap, laws in scales.items():
            tests[tag][cap] = {}
            for menu, enclosed in laws.items():
                fields = {
                    n: G * np.asarray(e) / inputs[n]["radius_kpc"]
                    for n, e in enclosed.items()
                }
                tests[tag][cap][menu] = scale_scores(inputs, fields)
    refs = [
        HERE / "RECORD_DENSITY.md",
        HERE / "record_profile_maps.json",
        HERE / "source_profile.json",
        HERE / "pair_galaxy_maps.json",
        source.HERE / "record_density.json",
        maps.HERE / "leroy_table7.dat",
        maps.HERE / "leroy_sample.dat",
    ]
    residual = np.concatenate(
        [abs(np.asarray(r["observed_kms"]) - r["baryons_kms"]) for r in inputs.values()]
    )
    sigma = np.concatenate([r["error_kms"] for r in inputs.values()])
    return dict(
        sha256={p.name: sha256(p.read_bytes()).hexdigest() for p in refs},
        tests=tests,
        baryons_only=dict(
            within_two=int(sum(residual <= 2 * sigma)),
            median_abs_kms=float(np.median(residual)),
        ),
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--image-cache", type=Path)
    cli.add_argument("--write", action="store_true")
    options = cli.parse_args()
    path = HERE / "record_profile_maps.json"
    bank = json.loads((source.HERE / "record_density.json").read_bytes())
    for label, p in (("1/16", 1 / 16), ("1/4", 0.25), ("1/2", 0.5)):
        for menu, v in source_laws(np.array([p]), bank).items():
            np.testing.assert_allclose(
                v[0],
                bank["cases"][label][menu]["receiver_nats_per_attempt"],
                atol=1e-15,
            )
    if options.image_cache:
        value = acquire(options.image_cache, bank)
        if options.write:
            path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
        else:
            galaxy.check(json.loads(path.read_bytes()), value)
    result = run(json.loads(path.read_bytes()))
    out = HERE / "record_profiles.json"
    if options.write:
        out.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        galaxy.check(json.loads(out.read_bytes()), result)
    for name, scales in result["tests"].items():
        print(
            name,
            {
                m: [
                    (c, t[m]["within_two"], t[m]["possible_full_curves"])
                    for c, t in scales.items()
                ]
                for m in MENUS
            },
        )
