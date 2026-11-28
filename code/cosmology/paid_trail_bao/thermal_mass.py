"""Convert the measured overlap-induced internal energy on ordinary mass standards."""

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from . import record_profiles as rp

HERE = rp.HERE
LAB = rp.source.HERE


def calibration():
    energy = json.loads((LAB / "thermal_energy.json").read_bytes())
    ordinary = json.loads((LAB / "thermal_standard.json").read_bytes())
    # From the signed source and the native spectrum, independent of any orbit.
    low, high = energy["energy_mixed"]
    a, b = ordinary["double"]["energy_band"]
    return low / b, high / a


def acquire(cache):
    inputs, *_ = rp.galaxy.transfer.budget.load_inputs()
    bank = json.loads((LAB / "record_density.json").read_bytes())
    peaks = json.loads((HERE / "record_profile_maps.json").read_bytes())[
        "peak_msun_pc2"
    ]
    assert bank["classes"][0]["magnitude_counts"][1:] == [0, 0, 0]
    out = {
        tag: {str(k): {m: {} for m in rp.MENUS} for k in rp.CAPACITY}
        for tag in rp.TRACERS
    }
    for n, tag, img, ring, selected in rp.observed_cells(cache, inputs):
        for cap in rp.CAPACITY:
            full = cap * peaks[tag]
            p = img[selected] / full
            delta = rp.source_laws(p, bank, per_joint_change=1.0)
            for menu, d in delta.items():
                assert np.all(abs(d) <= 4 * p + 1e-15)
                by_shell = np.bincount(
                    ring[selected], weights=d, minlength=len(inputs[n]["radius_kpc"])
                )
                # At full filling, four double energies correspond to full*area in solar masses.
                out[tag][str(cap)][menu][n] = (
                    by_shell.cumsum()
                    * rp.galaxy.WIDTH**2
                    * rp.maps.AREA_UNITS
                    * full
                    / 4
                ).tolist()
    return out


def forecast(maps, scale, inputs):
    low, high = scale
    G = rp.galaxy.transfer.budget.local.profiles.galaxies.G
    tests = {}
    for tag, scales in maps.items():
        tests[tag] = {}
        for c, menu_maps in scales.items():
            tests[tag][c] = {}
            for menu, rows in menu_maps.items():
                total = 0
                residual = []
                worst = 0.0
                ring_scores = {}
                for name, v in inputs.items():
                    source = G * np.asarray(rows[name]) / v["radius_kpc"]
                    bary = np.asarray(v["baryons_kms"])
                    obs = np.asarray(v["observed_kms"])
                    endpoints = bary**2 + np.array([low, high])[:, None] * source
                    assert np.min(endpoints) >= 0
                    lower, upper = np.sort(np.sqrt(endpoints), axis=0)
                    errors = np.asarray(v["error_kms"])
                    allowed = (upper >= obs - 2 * errors) & (lower <= obs + 2 * errors)
                    total += int(allowed.sum())
                    mid = (lower + upper) / 2
                    residual.extend(abs(mid - obs))
                    worst = max(worst, max(abs(lower - bary)), max(abs(upper - bary)))
                    ring_scores[name] = dict(
                        within_two=int(allowed.sum()),
                        signed_residual_kms=(mid - obs).tolist(),
                    )
                tests[tag][c][menu] = dict(
                    within_two=total,
                    median_abs_kms=float(np.median(residual)),
                    maximum_baryon_speed_change_kms=float(worst),
                    galaxies=ring_scores,
                )
    return tests


def run(maps):
    inputs, _, meta, _ = rp.galaxy.transfer.budget.load_inputs()
    scale = calibration()
    high = scale[1]
    tests = forecast(maps, scale, inputs)
    prior = json.loads((HERE / "outer_disk.json").read_bytes())
    all_gas = prior["all_atomic_He_msun"] + prior["all_molecular_He_msun"]
    stars = (
        1e9 * meta["NGC2841"]["L36_gsun"]
    )  # Grant one solar mass per solar 3.6-micron luminosity.
    need = prior["strongest"]["required_extra_mass_floor_msun"]
    ratio = 0.1200 / 0.02237
    single = json.loads((LAB / "thermal_standard.json").read_bytes())["single"][
        "energy_band"
    ][0]
    thermal = json.loads((LAB / "thermal_energy.json").read_bytes())
    bath_floor = 2 * single + thermal["energy_detector"]["energy_range"][0]
    assert bath_floor > 0
    # A conditional product repair changes the joint cell by at most one
    # quarter per attempt. Count every extra pair's ordinary energy as well.
    repeated = thermal["energy_mixed"][1] / (4 * bath_floor)
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                HERE / "THERMAL_MASS.md",
                HERE / "thermal_mass_maps.json",
                HERE / "source_profile.json",
                HERE / "record_profile_maps.json",
                HERE / "outer_disk.json",
                LAB / "thermal_standard.json",
                LAB / "thermal_energy.json",
                LAB / "record_density.json",
                rp.maps.HERE / "leroy_table7.dat",
                rp.maps.HERE / "leroy_sample.dat",
            )
        },
        energy_per_double=scale,
        absolute_mass_ratio_ceiling=high,
        planck_cold_per_baryon=ratio,
        cosmic_minimum_factor=ratio / high,
        independent_baths=dict(
            minimum_pair_energy=bath_floor,
            mass_ratio_ceiling=repeated,
            cosmic_minimum_factor=ratio / repeated,
        ),
        ngc2841=dict(
            all_gas_msun=all_gas,
            granted_stars_msun=stars,
            all_extra_ceiling_msun=(all_gas + stars) * high,
            required_shell_msun=need,
            minimum_factor=need / ((all_gas + stars) * high),
        ),
        tests=tests,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    cli.add_argument("--image-cache", type=Path)
    args = cli.parse_args()
    path = HERE / "thermal_mass_maps.json"
    if args.image_cache:
        data = acquire(args.image_cache)
        if args.write:
            path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
        else:
            rp.galaxy.check(json.loads(path.read_bytes()), data)
    result = run(json.loads(path.read_bytes()))
    dst = HERE / "thermal_mass.json"
    if args.write:
        dst.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        rp.galaxy.check(json.loads(dst.read_bytes()), result)
    print({k: v for k, v in result.items() if k not in ("tests", "sha256")})
    print(
        {
            tag: [
                (
                    min(v[m]["within_two"] for v in caps.values()),
                    max(v[m]["within_two"] for v in caps.values()),
                    max(v[m]["maximum_baryon_speed_change_kms"] for v in caps.values()),
                )
                for m in rp.MENUS
            ]
            for tag, caps in result["tests"].items()
        }
    )
