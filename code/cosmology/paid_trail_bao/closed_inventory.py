"""Use the closed source's eligible inventories and measured ordinary energy on galaxy images."""

import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from . import thermal_mass as mass

rp = mass.rp
HERE = mass.HERE
LAB = mass.LAB


def acquire(cache):
    inputs, *_ = rp.galaxy.transfer.budget.load_inputs()
    peaks = json.loads((HERE / "record_profile_maps.json").read_bytes())[
        "peak_msun_pc2"
    ]
    sources = {
        tracer: {str(c): {"closed": {}} for c in rp.CAPACITY} for tracer in rp.TRACERS
    }
    for name, tag, img, rings, mask in rp.observed_cells(cache, inputs):
        for cap in rp.CAPACITY:
            full = peaks[tag] * cap
            p = img[mask] / full
            assert np.all((0 <= p) & (p <= 1))
            eligible = 6 * p**2 * (1 - p) ** 2 * (448 / 2110)
            shell = np.bincount(
                rings[mask], weights=eligible, minlength=len(inputs[name]["radius_kpc"])
            )
            sources[tag][str(cap)]["closed"][name] = (
                shell.cumsum() * rp.galaxy.WIDTH**2 * rp.maps.AREA_UNITS * full / 4
            ).tolist()
    return sources


def run(sources):
    inputs, *_ = rp.galaxy.transfer.budget.load_inputs()
    data = json.loads((LAB / "closed_energy.json").read_bytes())
    interval = json.loads((LAB / "thermal_standard.json").read_bytes())["double"][
        "energy_band"
    ]
    priced = {
        name: [Q(row["energy_increment"])] * 2
        for name, row in data["matched_energy_extremes"].items()
    }
    priced.update(
        whole=[Q(data["canonical_total_increment"])] * 2,
        one_bond=[Q(data["mandatory_bond_increment"])] * 2,
    )
    bound = max(
        abs(priced["minimum"][0]), abs(priced["maximum"][0]), abs(priced["whole"][0])
    )
    priced["matched_envelope"] = [-bound, bound]
    supported = data["all_supported_prices"]
    bound = max(
        abs(Q(supported["coherent_lower"])), abs(Q(supported["coherent_upper"]))
    )
    priced["all_supported_envelope"] = [-bound, bound]
    result = {}
    for label, (lo, hi) in priced.items():
        limits = [e / Q(s) for e in (lo, hi) for s in interval]
        scaled = [
            float(np.nextafter(float(min(limits)), -np.inf)),
            float(np.nextafter(float(max(limits)), np.inf)),
        ]
        panels = mass.forecast(sources, scaled, inputs)
        result[label] = {}
        for tracer, cs in panels.items():
            for cap, k in cs.items():
                row = k["closed"]
                names = {n: r["within_two"] for n, r in row["galaxies"].items()}
                result[label][f"{tracer}/{cap}"] = dict(
                    within_two=row["within_two"],
                    complete=sum(
                        k == len(inputs[n]["radius_kpc"]) for n, k in names.items()
                    ),
                    max_speed_change=row["maximum_baryon_speed_change_kms"],
                    galaxies=names,
                )
    max_ratio = Q(2, 9) * Q(448, 2110) * bound / Q(interval[0])
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                HERE / "CLOSED_INVENTORY.md",
                HERE / "closed_inventory_maps.json",
                LAB / "closed_energy.json",
                LAB / "thermal_standard.json",
                HERE / "record_profile_maps.json",
                HERE / "source_profile.json",
            )
        },
        response_energy={k: [str(x) for x in v] for k, v in priced.items()},
        maximum_extra_per_ordinary=float(np.nextafter(float(max_ratio), np.inf)),
        planck_minimum_shortfall=float(
            np.nextafter(float(Q(12000, 2237) / max_ratio), -np.inf)
        ),
        all_accepted_stock_ratio_upper=float(
            np.nextafter(float(bound / (2 * Q(interval[0]))), np.inf)
        ),
        ordinary_double_energy=interval,
        galaxies=len(inputs),
        rings=sum(len(v["radius_kpc"]) for v in inputs.values()),
        forecasts=result,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--image-cache", type=Path)
    cli.add_argument("--write", action="store_true")
    opt = cli.parse_args()
    path = HERE / "closed_inventory_maps.json"
    if opt.image_cache:
        values = acquire(opt.image_cache)
        if opt.write:
            path.write_text(json.dumps(values, indent=2, allow_nan=False) + "\n")
        else:
            rp.galaxy.check(json.loads(path.read_bytes()), values)
    out = run(json.loads(path.read_bytes()))
    dst = HERE / "closed_inventory.json"
    if opt.write:
        dst.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
    else:
        rp.galaxy.check(json.loads(dst.read_bytes()), out)
    print(out["maximum_extra_per_ordinary"], out["planck_minimum_shortfall"])
    for k, rows in out["forecasts"].items():
        print(
            k,
            sorted({(p["within_two"], p["complete"]) for p in rows.values()}),
            max(p["max_speed_change"] for p in rows.values()),
        )
