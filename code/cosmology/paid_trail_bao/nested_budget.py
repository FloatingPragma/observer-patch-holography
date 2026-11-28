"""Test the complete closed inventory's spectral span against resolved disk masses."""

import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from . import closed_inventory as inventory

rp = inventory.rp
HERE = rp.HERE
LAB = inventory.LAB


def acquire(cache):
    inputs, *_ = rp.galaxy.transfer.budget.load_inputs()
    prior = json.loads((HERE / "pair_galaxy_maps.json").read_bytes())["galaxies"]
    out = {}
    for name, tracer, img, _, _ in rp.observed_cells(cache, inputs):
        r = np.asarray(inputs[name]["radius_kpc"])
        low = r[int(np.argmin(abs(r - r[-1] / 2)))]
        high = r[-1]
        width = rp.galaxy.WIDTH
        axis = (np.arange(len(img)) + 0.5 - len(img) // 2) * width
        squared = axis[:, None] ** 2 + axis[None, :] ** 2
        mass = img * width**2 * rp.maps.AREA_UNITS
        assert mass.min() >= 0
        row = out.setdefault(
            name,
            dict(
                raw_fits_sha256=prior[name]["raw_sha256"],
                radii_kpc=r.tolist(),
                start_kpc=float(low),
                tracers={},
            ),
        )
        row["tracers"][tracer] = dict(
            enclosed_msun=[float(mass[squared <= x * x].sum()) for x in r],
            imaged_msun=float(mass.sum()),
            integrated_rM_msun_kpc2=float(
                np.sum(mass * np.maximum(0, high**2 - np.maximum(low**2, squared))) / 2
            ),
        )
    return out


def evaluate(data):
    lab = json.loads((LAB / "nested_source.json").read_bytes())
    lines = {
        k: Q(v["span_per_double_upper"])
        for k, v in lab["inventory_spans"].items()
        if int(k)
    }
    assert Q(lab["inventory_spans"]["0"]["span_upper"]) == 0
    cal = json.loads((LAB / "thermal_standard.json").read_bytes())["double"][
        "energy_band"
    ]
    ratio = float(np.nextafter(float(max(lines.values()) / Q(cal[0])), np.inf))
    per_fuel = [
        [Q(x) / int(k) for x in row["energy_band"]]
        for k, row in lab["inventory_spans"].items()
        if int(k)
    ]
    low_energy = min(lo for lo, _ in per_fuel)
    high_energy = max(hi for _, hi in per_fuel)
    assert high_energy - low_energy == max(lines.values())
    absolute_ratio = float(np.nextafter(float(high_energy / Q(cal[0])), np.inf))
    inputs, *_ = rp.galaxy.transfer.budget.load_inputs()
    G = rp.galaxy.transfer.budget.local.profiles.galaxies.G
    out = {}
    for name, tab in inputs.items():
        img = data[name]
        r, v, b, e = [
            np.asarray(tab[k])
            for k in ("radius_kpc", "observed_kms", "baryons_kms", "error_kms")
        ]
        np.testing.assert_array_equal(img["radii_kpc"], r)
        j = list(r).index(img["start_kpc"])
        low = r / G * (np.maximum(0, v - 2 * e) ** 2 - b * b)
        high = r / G * ((v + 2 * e) ** 2 - b * b)
        # A nonnegative density gives a nondecreasing enclosed source mass.
        # Every earlier measured lower mass therefore constrains the next interval.
        floor = np.maximum.accumulate(np.maximum(low, 0))
        area = np.diff(r * r) / 2
        monotone = float(floor[j:-1] @ area[j:])
        # A separate signed comparison interpolates the low mass linearly in r.
        linear = float(
            np.sum(
                np.diff(r[j:])
                / 6
                * (
                    (2 * r[j:-1] + r[j + 1 :]) * low[j:-1]
                    + (r[j:-1] + 2 * r[j + 1 :]) * low[j + 1 :]
                )
            )
        )
        rows = {}
        for tag, m in img["tracers"].items():
            enclosed = np.array(m["enclosed_msun"])
            limit = ratio * enclosed
            allowed = (low <= limit) & (high >= -limit)
            local = m["integrated_rM_msun_kpc2"] * ratio
            # Even put the *entire* imaged fuel reservoir at the center.
            moved = m["imaged_msun"] * ratio * (r[-1] ** 2 - r[j] ** 2) / 2
            rows[tag] = dict(
                rings_compatible_signed=int(allowed.sum()),
                whole_curve_compatible=bool(allowed.all()),
                integrated_local_ceiling=local,
                integrated_relocated_ceiling=moved,
                required_to_local=monotone / local,
                required_to_relocated=monotone / moved,
                linear_to_local=linear / local,
                # The entire final stock can be moved to any single aperture.
                # Dropping every positive initial energy only increases its change.
                required_to_all_final_energy=float(
                    low.max() / (absolute_ratio * m["imaged_msun"])
                ),
            )
        out[name] = dict(
            radii_kpc=[float(r[j]), float(r[-1])],
            rings=len(r),
            nonnegative_mass_integral_lower=monotone,
            linear_mass_integral_lower=linear,
            tracers=rows,
        )
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                LAB / "NESTED_SOURCE.md",
                LAB / "nested_source.json",
                LAB / "thermal_standard.json",
                HERE / "nested_budget_maps.json",
                HERE / "source_profile.json",
                rp.maps.HERE / "leroy_table7.dat",
                rp.maps.HERE / "leroy_sample.dat",
            )
        },
        span_per_double_max=str(max(lines.values())),
        ordinary_energy=cal,
        mass_ratio_upper=ratio,
        absolute_energy_per_double=[str(low_energy), str(high_energy)],
        absolute_final_mass_ratio_upper=absolute_ratio,
        planck_base_lcdm_central_ratio=0.1200 / 0.02237,
        planck_central_to_span=(0.1200 / 0.02237) / ratio,
        galaxies=out,
        summary={
            t: dict(
                ring_count=sum(g["rings"] for g in out.values()),
                rings_compatible_signed=sum(
                    g["tracers"][t]["rings_compatible_signed"] for g in out.values()
                ),
                whole_curves_compatible=sum(
                    g["tracers"][t]["whole_curve_compatible"] for g in out.values()
                ),
                whole_curves_excluded_all_final_energy=sum(
                    g["tracers"][t]["required_to_all_final_energy"] > 1
                    for g in out.values()
                ),
                local_positive_integral_excluded=sum(
                    g["tracers"][t]["required_to_local"] > 1 for g in out.values()
                ),
                relocated_positive_integral_excluded=sum(
                    bool(g["tracers"][t]["required_to_relocated"] > 1)
                    for g in out.values()
                ),
            )
            for t in rp.TRACERS
        },
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--image-cache", type=Path)
    cli.add_argument("--write", action="store_true")
    opt = cli.parse_args()
    path = HERE / "nested_budget_maps.json"
    if opt.image_cache:
        raw = acquire(opt.image_cache)
        if opt.write:
            path.write_text(json.dumps(raw, indent=2, allow_nan=False) + "\n")
        else:
            rp.galaxy.check(json.loads(path.read_bytes()), raw)
    out = evaluate(json.loads(path.read_bytes()))
    result = HERE / "nested_budget.json"
    if opt.write:
        result.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
    else:
        rp.galaxy.check(json.loads(result.read_bytes()), out)
    print(out["mass_ratio_upper"], out["summary"])
