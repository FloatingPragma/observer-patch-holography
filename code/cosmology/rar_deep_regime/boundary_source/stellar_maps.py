"""Assess directly mapped baryonic phases in the locked orbit pair."""

import argparse
import hashlib
import json
from pathlib import Path
from urllib.request import urlretrieve

import numpy as np
from scipy.special import ellipe, ellipk

from .measure import G, run as primary_results

from .distance import PAIR, run as distance_run
from .maps import FINE, GEOMETRY, coarse, mapped_shell_response, square_tv
from .source import AREA_UNITS, HERE, edges, enroll

SUPERSAMPLE = 8


def files(name, cache):
    paths = {}
    for k in ("stellar", "ICAmask", "mom0", "emom0"):
        if k in ("stellar", "ICAmask"):
            filename = f"{name}.{k}.fits"
            url = f"https://irsa.ipac.caltech.edu/data/SPITZER/S4G/galaxies/{name}/P5/{filename}"
        else:
            filename = f"{name.lower()}_heracles_{k}.fits.gz"
            url = f"https://www.iram.fr/ILPA/LP001/{filename}"
        path = cache / filename
        if not path.exists():
            cache.mkdir(parents=True, exist_ok=True)
            urlretrieve(url, path)
        paths[k] = path
    return paths


def sky_grid(name, distance, x, y):
    """Rotate a disk tangent plane onto the celestial sphere (gnomonic)."""
    geom = GEOMETRY[name]
    pa, inc, ra, dec = np.deg2rad([geom[k] for k in ("pa", "inc", "ra", "dec")])
    east = (x * np.sin(pa) + y * np.cos(inc) * np.cos(pa)) / (distance * 1e3)
    north = (x * np.cos(pa) - y * np.cos(inc) * np.sin(pa)) / (distance * 1e3)
    cord = np.cos(dec) - north * np.sin(dec)
    sky_ra = ra + np.arctan2(east, cord)
    sky_dec = np.arctan2(np.sin(dec) + north * np.cos(dec), np.hypot(east, cord))
    return np.rad2deg(sky_ra), np.rad2deg(sky_dec)


def sample_image(path, ra, dec, unit, apply_detector_sip=False):
    from astropy.io import fits
    from astropy.wcs import WCS

    image, hdr = fits.getdata(path, header=True)
    assert hdr["BUNIT"] == unit
    # North-up mosaics advertise TAN. Some copied detector headers contain a
    # SIP block without a SIP projection; retain it only for the explicit control.
    if unit == "MJy/sr":
        if apply_detector_sip:
            hdr["CTYPE1"] += "-SIP"
            hdr["CTYPE2"] += "-SIP"
        else:
            for k in list(hdr):
                if k.startswith(("A_", "B_", "AP_", "BP_")):
                    del hdr[k]
    x, y = WCS(hdr, naxis=2).all_world2pix(ra, dec, 0)
    ix, iy = np.floor(x + 0.5).astype(int), np.floor(y + 0.5).astype(int)
    valid = (ix >= 0) & (iy >= 0) & (ix < image.shape[1]) & (iy < image.shape[0])
    samples = np.full(len(ra), np.nan)
    samples[valid] = image[iy[valid], ix[valid]]
    return samples


def acquire(name, occupied_grid, direct, cache, body):
    d = direct["distance_mpc"]
    _, lrow, rows, _, _ = body
    paths = files(name, cache)
    shell_edges = edges(
        np.array([float(r["r"]) for r in rows]) * d / float(lrow["Dist"])
    )
    ring_density = np.array([float(r["Sigma*"]) for r in rows])
    extent = int(np.ceil(shell_edges[-1] / FINE)) + 1
    n = SUPERSAMPLE
    axis = (np.arange(2 * extent * n) + 0.5 - extent * n) * FINE / n
    x, y = np.meshgrid(axis, axis, indexing="ij")
    r = np.hypot(x, y)
    use = r < shell_edges[-1]
    shell = np.searchsorted(shell_edges[1:], r[use])
    ra, dec = sky_grid(name, d, x[use], y[use])
    light = sample_image(paths["stellar"], ra, dec, "MJy/sr")
    alt_light = sample_image(paths["stellar"], ra, dec, "MJy/sr", True)
    alt_mask = sample_image(paths["ICAmask"], ra, dec, "MJy/sr", True)
    mask = sample_image(paths["ICAmask"], ra, dec, "MJy/sr")
    co = sample_image(paths["mom0"], ra, dec, "K KM/S")
    err = sample_image(paths["emom0"], ra, dec, "K KM/S")
    valid = np.isfinite(light) & (mask == 0)
    # Querejeta et al. eq. 6 per 0.75 arcsec pixel, divided by its face-on area.
    cos_i = np.cos(np.deg2rad(GEOMETRY[name]["inc"]))
    conversion = 9308.23 * 0.6 * cos_i / (1e6 * np.deg2rad(0.75 / 3600)) ** 2
    surface = np.maximum(light, 0) * conversion
    count = np.bincount(shell[valid], minlength=len(ring_density))
    if np.any(count == 0):
        raise ValueError("An entire source ring is missing old-stellar light")
    mean = (
        np.bincount(shell[valid], weights=surface[valid], minlength=len(ring_density))
        / count
    )
    imputed = np.where(valid, surface, mean[shell])
    av = np.isfinite(alt_light) & (alt_mask == 0)
    alt_density = np.maximum(alt_light, 0) * conversion
    alt_mean = np.bincount(shell[av], weights=alt_density[av]) / np.bincount(shell[av])
    alt_star = np.where(av, alt_density, alt_mean[shell])
    info = {
        "sha256": {
            k: hashlib.sha256(p.read_bytes()).hexdigest() for k, p in paths.items()
        },
        "masked_stellar_area_fraction": 1 - valid.mean(),
        "missing_stellar_light_fraction_after_infill": imputed[~valid].sum()
        / imputed.sum(),
        "CO_valid_area_fraction": float((np.isfinite(co) & np.isfinite(err)).mean()),
        "max_absolute_fractional_stellar_ring_imputation": float(
            max(1 - count / np.bincount(shell))
        ),
        "resolved_stellar_msun": float(imputed.sum() * AREA_UNITS * (FINE / n) ** 2),
    }
    fields = {
        "star": imputed,
        "star_detector_sip": alt_star,
        "star_ring_mass": imputed * (ring_density / mean)[shell],
        "CO3": np.where((co >= 3 * err) & (err > 0), co, 0) * 4.35 * cos_i / 0.8,
        "CO5": np.where((co >= 5 * err) & (err > 0), co, 0) * 4.35 * cos_i / 0.8,
    }
    packed = {}
    for key, val in fields.items():
        buf = np.zeros_like(r)
        buf[use] = val
        packed[key] = buf.reshape(2 * extent, n, 2 * extent, n).mean(axis=(1, 3))
    info["crop_extent_cells"] = extent
    assert len(occupied_grid) // 2 >= extent
    return packed, info


def build(cache):
    distance = distance_run()["1"]["galaxies"]
    kept, _ = enroll()
    arrays, info = {}, {}
    with np.load(HERE / "atomic_grids.npz") as atomic:
        for name in PAIR:
            fields, info[name] = acquire(
                name, atomic[name], distance[name], cache, kept[name]
            )
            arrays.update({f"{name}_{k}": v for k, v in fields.items()})
    np.savez_compressed(HERE / "stellar_co_grids.npz", **arrays)
    (HERE / "stellar_co_inputs.json").write_text(json.dumps(info, indent=2) + "\n")


def ordinary_disk_v2(radius, mass, rings, height):
    """Azimuthally averaged inward force times R for softened source rings."""
    ap = (radius + rings) ** 2 + height**2
    am = (radius - rings) ** 2 + height**2
    m = 4 * radius * rings / ap
    kernel = (ellipk(m) + (radius**2 - rings**2 - height**2) * ellipe(m) / am) / (
        np.pi * np.sqrt(ap)
    )
    return float(G * np.sum(mass * kernel))


def run():
    direct = distance_run()["1"]
    result = {}
    primary = primary_results()["targets"]
    enrolled, _ = enroll()
    inputs = json.loads((HERE / "stellar_co_inputs.json").read_text())
    with (
        np.load(HERE / "atomic_grids.npz") as atomic,
        np.load(HERE / "stellar_co_grids.npz") as fields,
    ):
        for name in PAIR:
            d = direct["galaxies"][name]
            hi = atomic[name] / (FINE**2 * AREA_UNITS)
            total = len(hi) // 2
            extent = inputs[name]["crop_extent_cells"]
            star, shaped, co3, co5, detector = [
                np.pad(fields[f"{name}_{k}"], total - extent)
                for k in ("star", "star_ring_mass", "CO3", "CO5", "star_detector_sip")
            ]
            meta = enrolled[name][0]
            inventory = (
                1.36e9 * meta["MHI_gsun"] * (d["distance_mpc"] / meta["dist_mpc"]) ** 2
            )
            hi_scale = inventory / atomic[name].sum()
            recipes = {
                "mapped_all": star + hi + co3,
                "mapped_all_detector_SIP": detector + hi + co3,
                "mapped_all_global_HI": star + hi * hi_scale + co3,
                "leroy_ring_stellar_mass": shaped + hi + co3,
                "mapped_all_5sigma_CO": star + hi + co5,
                "mapped_stars": star,
                "mapped_both_gases": hi + co3,
            }
            grids = {}
            for width in (0.8, 1.6):
                factor = round(width / FINE)
                for offset in (0, factor // 2):
                    cell = {k: coarse(v, factor, offset) for k, v in recipes.items()}
                    _, finite = mapped_shell_response(
                        cell["mapped_all"], width, offset * FINE, d["orbit_kpc"]
                    )
                    grids[f"{width:g}kpc_{offset * FINE:g}offset"] = {
                        k: square_tv(v, width) for k, v in cell.items()
                    } | {"mapped_all_finite": finite}
            speed2 = primary[name]["observed_kms"] ** 2
            axis = (np.arange(len(star)) + 0.5 - total) * FINE
            radii = np.hypot(axis[:, None], axis[None, :])
            ordinary = {}
            for recipe in ("mapped_all", "mapped_all_global_HI"):
                density = recipes[recipe]
                use = density > 0
                ordinary[recipe] = {}
                for h in (0.2, 0.4, 0.8):
                    v2 = ordinary_disk_v2(
                        d["orbit_kpc"],
                        density[use] * FINE**2 * AREA_UNITS,
                        radii[use],
                        h,
                    )
                    ordinary[recipe][str(h)] = {
                        "baryonic_v2": v2,
                        "excess_v2": speed2 - v2,
                    }
            result[name] = {
                "ordinary_disk_force": ordinary,
                "grids": grids,
                "stellar_msun": star.sum() * FINE**2 * AREA_UNITS,
                "ring_mass_scaled_stellar_msun": shaped.sum() * FINE**2 * AREA_UNITS,
                "CO_3sigma_msun_with_He": co3.sum() * FINE**2 * AREA_UNITS,
            }
    a, b = PAIR
    ratios = {
        grid: {k: v / result[b]["grids"][grid][k] for k, v in recipe.items()}
        for grid, recipe in result[a]["grids"].items()
    }
    ratio_v2 = {
        recipe: {
            height: row["excess_v2"]
            / result[b]["ordinary_disk_force"][recipe][height]["excess_v2"]
            for height, row in heights.items()
        }
        for recipe, heights in result[a]["ordinary_disk_force"].items()
    }
    return {
        "profiles": result,
        "ratios": ratios,
        "observed_excess_ratio": direct["Y_ratio"],
        "excess_ratio_using_mapped_baryons": ratio_v2,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rebuild-from-fits", type=Path)
    args = parser.parse_args()
    if args.rebuild_from_fits:
        build(args.rebuild_from_fits)
    (HERE / "stellar_co.json").write_text(
        json.dumps(run(), indent=2, allow_nan=False) + "\n"
    )
