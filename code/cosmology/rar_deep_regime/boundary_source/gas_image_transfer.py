"""Full measured two-dimensional gas sources versus eight outside orbits."""

import argparse
import hashlib
import json
from pathlib import Path
from urllib.request import urlretrieve

import numpy as np

from .gas_transfer import CALIBRATOR, OTHER_PAIR_MEMBER
from .maps import (
    FINE,
    GEOMETRY,
    coarse,
    mapped_shell_response,
    raw_map,
    ring_surface,
    square_tv,
)
from .measure import run as primary
from .source import AREA_UNITS, HERE, enroll
from .stellar_maps import SUPERSAMPLE, sample_image, sky_grid


def molecular(n, grid, edge, distance, cache):
    info = {}
    if n == "NGC7793":  # absent from the 48-object HERACLES public release
        return np.zeros_like(grid), {"coverage": "absent"}
    extent = int(np.ceil(edge / FINE)) + 1
    fine = SUPERSAMPLE
    idx = (np.arange(2 * extent * fine) + 0.5 - extent * fine) * FINE / fine
    x, y = np.meshgrid(idx, idx, indexing="ij")
    use = np.hypot(x, y) < edge
    ra, dec = sky_grid(n, distance, x[use], y[use])
    values = {}
    for k in ("mom0", "emom0"):
        filename = f"{n.lower()}_heracles_{k}.fits.gz"
        path = cache / filename
        if not path.exists():
            urlretrieve("https://www.iram.fr/ILPA/LP001/" + filename, path)
        values[k] = sample_image(path, ra, dec, "K KM/S")
        info[k + "_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    co, err = values["mom0"], values["emom0"]
    buf = np.zeros_like(x)
    buf[use] = (
        np.where((co >= 3 * err) & (err > 0), co, 0)
        * 4.35
        * np.cos(np.deg2rad(GEOMETRY[n]["inc"]))
        / 0.8
    )
    result = buf.reshape(2 * extent, fine, 2 * extent, fine).mean(axis=(1, 3))
    info["sampled_fraction_of_aperture"] = float(np.isfinite(co + err).mean())
    info["coverage"] = "measured"
    return np.pad(result, len(grid) // 2 - extent), info


def acquire(cache):
    kept, _ = enroll()
    sources = {}
    for n, (meta, _lrow, _rows, orbit, source) in kept.items():
        if n == OTHER_PAIR_MEMBER:
            continue
        grid, provenance = raw_map(n, meta["dist_mpc"], cache)
        if grid.sum() <= 0:
            raise ValueError(n)
        h2, co_info = molecular(n, grid, source["edge_kpc"], meta["dist_mpc"], cache)
        cat_mass = 1.36e9 * meta["MHI_gsun"]
        hi = grid / (FINE**2 * AREA_UNITS)
        stars, _ = ring_surface(
            grid, _lrow, _rows, meta["dist_mpc"], columns=("Sigma*",)
        )
        scale = cat_mass / grid.sum()
        axes = (np.arange(len(grid)) + 0.5 - len(grid) // 2) * FINE
        rr = np.hypot(axes[:, None], axes[None, :])
        values = {}
        for width in (0.8, 1.6):
            for label, fields in [
                ("atomic", hi),
                ("both_gases", hi + h2),
                ("radial_stars_plus_mapped_gases", stars + hi + h2),
                ("global_HI_both_gases", hi * scale + h2),
            ]:
                field = coarse(fields, round(width / FINE), 0)
                _, finite = mapped_shell_response(field, width, 0.0, orbit)
                values[f"{width:g}_{label}_compact"] = square_tv(field, width)
                values[f"{width:g}_{label}_finite"] = finite
        sources[n] = {
            "radio": provenance,
            "CO": co_info,
            "atomic_He_msun": grid.sum(),
            "fraction_catalogue_HI": grid.sum() / cat_mass,
            "H2_He_msun": h2.sum() * FINE**2 * AREA_UNITS,
            "atomic_mass_fraction_beyond_orbit": grid[rr > orbit].sum() / grid.sum(),
            "source_msun_kpc": values,
        }
    return sources


def run(sources):
    target = primary()["targets"]
    held = sorted(set(sources) - {CALIBRATOR})
    result = {}
    for law, cal in sources[CALIBRATOR]["source_msun_kpc"].items():
        gain = target[CALIBRATOR]["excess_v2"] / cal
        rows = {}
        for n in held:
            excess = gain * sources[n]["source_msun_kpc"][law]
            speed = np.sqrt(excess + target[n]["baryon_v2"])
            boxes = {}
            for mult in (1, 2):
                required = {}
                for key in (CALIBRATOR, n):
                    t = target[key]
                    v = t["observed_kms"] + mult * t["e_kms"] * np.array([-1, 1])
                    required[key] = (v**2 - t["baryon_v2"]) / sources[key][
                        "source_msun_kpc"
                    ][law]
                boxes[str(mult)] = bool(
                    max(required[n][0], required[CALIBRATOR][0], 0)
                    <= min(required[n][1], required[CALIBRATOR][1])
                )
            rows[n] = {
                "common_gain_possible_with_quoted_speed_errors": boxes,
                "predicted_excess_v2": excess,
                "predicted_over_observed_excess": excess / target[n]["excess_v2"],
                "observed_kms": target[n]["observed_kms"],
                "total_residual_kms": speed - target[n]["observed_kms"],
            }
        errors = np.array([r["total_residual_kms"] for r in rows.values()])
        result[law] = {
            "calibration": CALIBRATOR,
            "gain": gain,
            "held": rows,
            "median_abs_kms": np.median(abs(errors)),
            "rms_kms": np.sqrt(np.mean(errors**2)),
        }
    return {
        "protocol_sha256": hashlib.sha256(
            (HERE / "GAS_IMAGE_TRANSFER_PROTOCOL.md").read_bytes()
        ).hexdigest(),
        "predictions": result,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rebuild-from-fits", type=Path)
    args = parser.parse_args()
    inputs = HERE / "gas_image_inputs.json"
    if args.rebuild_from_fits:
        inputs.write_text(
            json.dumps(acquire(args.rebuild_from_fits), indent=2, allow_nan=False)
            + "\n"
        )
    (HERE / "gas_image_transfer.json").write_text(
        json.dumps(run(json.loads(inputs.read_text())), indent=2, allow_nan=False)
        + "\n"
    )
