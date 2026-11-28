"""Transfer an imaged pair-abundance profile to each withheld SPARC galaxy."""

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.closed_region.measure import (
    check,
)
from . import pair_occupancy as pairs
from . import poisson_transfer as transfer

HERE = transfer.HERE
maps = pairs.profile.maps
WIDTH = 0.8


def acquire(cache):
    inputs, _, _, _ = transfer.budget.load_inputs()
    enrolled, _ = maps.enroll()
    known = {}
    for filename in ("gas_image_inputs.json", "map_inputs.json"):
        for name, record in json.loads((maps.HERE / filename).read_bytes()).items():
            known[name] = record.get("radio", record)["raw_fits_sha256"]
    profiles = {}
    for name, rows in inputs.items():
        meta, survey, stellar, _, _ = enrolled[name]
        grid, provenance = maps.raw_map(name, meta["dist_mpc"], cache)
        assert provenance["raw_fits_sha256"] == known[name]
        hi = grid / (maps.FINE**2 * maps.AREA_UNITS)
        rings, _ = maps.ring_surface(grid, survey, stellar, meta["dist_mpc"])
        output = {}
        for tag, tracer in (("atomic", hi), ("mapped_baryons", hi + rings)):
            surface = maps.coarse(tracer, round(WIDTH / maps.FINE), 0)
            axis = (np.arange(len(surface)) + 0.5 - len(surface) // 2) * WIDTH
            radius = np.hypot(axis[:, None], axis[None, :])
            frame = np.pad(surface, 1)
            near = (
                frame[:-2, 1:-1] + frame[2:, 1:-1] + frame[1:-1, :-2] + frame[1:-1, 2:]
            ) / 4
            values = WIDTH**2 * maps.AREA_UNITS * surface * near
            output[tag] = [float(values[radius <= r].sum()) for r in rows["radius_kpc"]]
        profiles[name] = dict(raw_sha256=known[name], enclosed=output)
    return dict(units="Msun^2/pc^2", width_kpc=WIDTH, galaxies=profiles)


def run(source):
    inputs, _, _, _ = transfer.budget.load_inputs()
    assert source["galaxies"].keys() == inputs.keys()
    predictions = {}
    for tag in ("atomic", "mapped_baryons"):
        fields = {
            n: transfer.budget.local.profiles.galaxies.G
            * np.array(row["enclosed"][tag])
            / np.array(inputs[n]["radius_kpc"])
            for n, row in source["galaxies"].items()
        }
        predictions[tag] = transfer.scale_test(inputs, fields)
    residual = np.concatenate(
        [abs(np.asarray(r["observed_kms"]) - r["baryons_kms"]) for r in inputs.values()]
    )
    errors = np.concatenate([r["error_kms"] for r in inputs.values()])
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                HERE / "PAIR_GALAXIES.md",
                HERE / "pair_galaxy_maps.json",
                HERE / "source_profile.json",
                maps.HERE / "leroy_table7.dat",
                maps.HERE / "leroy_sample.dat",
            )
        },
        tests=predictions,
        baryons_only=dict(
            rings=len(errors),
            within_two=int(np.sum(residual <= 2 * errors)),
            median_abs_kms=float(np.median(residual)),
        ),
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--image-cache", type=Path)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    src = HERE / "pair_galaxy_maps.json"
    if args.image_cache:
        data = acquire(args.image_cache)
        if args.write:
            src.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
        else:
            check(json.loads(src.read_bytes()), data)
    result = run(json.loads(src.read_bytes()))
    path = HERE / "pair_galaxies.json"
    if args.write:
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        check(json.loads(path.read_bytes()), result)
    print(
        result["baryons_only"],
        {
            k: (v["within_two"], v["median_abs_kms"], v["possible_full_curves"])
            for k, v in result["tests"].items()
        },
    )
