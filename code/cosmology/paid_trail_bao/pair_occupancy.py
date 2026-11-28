"""Test a local ordinary-pair source against the measured NGC 2841 gas image."""

import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from . import dragonfly_profile as profile

HERE = profile.HERE


def acquire(cache):
    grid, info = profile.maps.raw_map(profile.outer.GALAXY, 14.1, cache)
    prior = json.loads((HERE / "dragonfly_gas.json").read_bytes())
    assert prior["sha256"] == info["raw_fits_sha256"]
    out = {}
    for width in (0.4, 0.8, 1.6):
        surf = profile.maps.coarse(
            grid / profile.maps.FINE**2, round(width / profile.maps.FINE), 0
        )
        axis = (np.arange(len(surf)) + 0.5 - len(surf) // 2) * width
        rad = np.hypot(axis[:, None], axis[None, :])
        frame = np.pad(surf, 1)
        neighbors = (
            frame[:-2, 1:-1] + frame[2:, 1:-1] + frame[1:-1, :-2] + frame[1:-1, 2:]
        ) / 4
        assert abs(surf.sum() * width**2 / grid.sum() - 1) < 1e-10
        for label, field in (("square", surf * surf), ("neighbor", surf * neighbors)):
            rings = []
            for r, _ in prior["enclosed_atomic_He_msun"]:
                # Units of solar masses squared per square kpc. Store in units
                # of 10^6 to normalize subsequent exact-rational computations.
                rings.append([r, round(width**2 * field[rad <= r].sum() / 10**6)])
            out[f"{width}_{label}"] = rings
    return dict(
        url=prior["url"], sha256=prior["sha256"], units="1e6 Msun^2/kpc^2", panels=out
    )


def evaluate(data):
    mass = profile.rings()
    out = {}
    for tag, apertures in data["panels"].items():
        x = {Q(str(r)): Q(n) for r, n in apertures}
        rows = [(r, x[r] - x[20], lo, hi) for r, lo, hi in mass]
        assert all(a[1] < b[1] for a, b in zip(rows, rows[1:]))
        train = [v for v in rows if v[0] <= 35]
        held = [v for v in rows if v[0] > 35]
        parameters = profile.vertices(train)
        bands = []
        for r, y, lo, hi in held:
            pred = [a + g * y for a, g in parameters]
            bands.append(
                dict(
                    radius_kpc=float(r),
                    overlap=bool(pred and max(pred) >= lo and min(pred) <= hi),
                    predicted_msun=[float(min(pred)), float(max(pred))]
                    if pred
                    else None,
                    orbital_msun=[float(lo), float(hi)],
                )
            )
        lookup = {str(float(r)): (v, lo, hi) for r, v, lo, hi in rows}
        a, b, c, d = [lookup[k] for k in ("20.57", "34.88", "43.08", "49.19")]
        inner = (b[2] - a[1]) / (b[0] - a[0])
        outer = (d[1] - c[2]) / (d[0] - c[0])
        assert min(inner, outer) > 0
        out[tag] = dict(
            inside_has_fit=bool(parameters),
            held=bands,
            joint_has_fit=bool(profile.vertices(rows)),
            inner_shell_count=float(b[0] - a[0]),
            outer_shell_count=float(d[0] - c[0]),
            max_inner_mass_per_pair_unit=float(inner),
            min_outer_mass_per_pair_unit=float(outer),
            minimum_outer_to_inner_conversion=float(outer / inner),
        )
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                HERE / "PAIR_OCCUPANCY.md",
                HERE / "pair_occupancy_image.json",
                HERE / "source_profile.json",
            )
        },
        panels=out,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--image-cache", type=Path)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    src = HERE / "pair_occupancy_image.json"
    if args.image_cache:
        result = acquire(args.image_cache)
        if args.write:
            src.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
        else:
            assert json.loads(src.read_bytes()) == result
    result = evaluate(json.loads(src.read_bytes()))
    path = HERE / "pair_occupancy.json"
    if args.write:
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        assert json.loads(path.read_bytes()) == result
    print(
        {
            k: (
                v["inside_has_fit"],
                sum(x["overlap"] for x in v["held"]),
                v["joint_has_fit"],
                v["minimum_outer_to_inner_conversion"],
            )
            for k, v in result["panels"].items()
        }
    )
