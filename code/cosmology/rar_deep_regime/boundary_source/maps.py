"""Two full THINGS atomic-gas maps plus independently published stellar and CO rings."""

import argparse
import hashlib
import json
from pathlib import Path
from urllib.request import urlretrieve

import numpy as np

from .distance import PAIR, run as distance_run
from .source import AREA_UNITS, HERE, THRESHOLD_MAX, edges, enroll, filled

FINE = 0.4
# Walter et al. survey geometry, beam and global H I flux tables; exact centres from de Blok et al. Table 2.
GEOMETRY = {
    "NGC3198": dict(
        ra=15 * (10 + 19 / 60 + 55 / 3600),
        dec=45 + 32 / 60 + 58.9 / 3600,
        beam=[13.01, 11.56],
        pa=215,
        inc=72,
        flux=227.0,
    ),
    "NGC3521": dict(
        ra=15 * (11 + 5 / 60 + 48.6 / 3600),
        dec=-(2 / 60 + 9.2 / 3600),
        beam=[14.14, 11.15],
        pa=340,
        inc=73,
        flux=297.0,
    ),
}


def sexagesimal(d, m, s):
    return np.sign(d) * (abs(d) + m / 60 + s / 3600)


# Walter et al. survey centres, beams (major/minor arcsec), position angles
# and global gas fluxes, also adopting the same Leroy sample inclinations.
GEOMETRY.update(
    {
        name: dict(
            ra=sexagesimal(*ra) * 15,
            dec=sexagesimal(*dec),
            beam=beam,
            pa=pa,
            inc=inc,
            flux=flux,
        )
        for name, ra, dec, beam, pa, inc, flux in [
            ("DDO154", (12, 54, 5.9), (27, 9, 10), (14.09, 12.62), 230, 66, 82.1),
            ("IC2574", (10, 28, 27.7), (68, 24, 59), (12.81, 11.90), 56, 53, 387),
            ("NGC2403", (7, 36, 51.1), (65, 36, 3), (8.75, 7.65), 124, 63, 1054),
            ("NGC2841", (9, 22, 2.6), (50, 58, 35), (11.06, 9.37), 153, 74, 183),
            ("NGC5055", (13, 15, 49.2), (42, 1, 45), (10.06, 8.66), 102, 59, 379),
            ("NGC6946", (20, 34, 52.2), (60, 9, 14), (6.04, 5.61), 243, 33, 506),
            ("NGC7331", (22, 37, 4.1), (34, 24, 57), (6.13, 5.61), 168, 76, 179),
            ("NGC7793", (23, 57, 49.7), (-32, 35, 28), (15.60, 10.85), 290, 50, 246),
        ]
    }
)


def raw_map(name, distance, cache):
    # Astropy is needed only when reacquiring original FITS, not for offline replay.
    from astropy.io import fits

    props = GEOMETRY[name]
    stem = name.replace("NGC", "NGC_").replace("IC", "IC_")
    path = cache / f"{stem}.fits"
    if not path.exists():
        cache.mkdir(parents=True, exist_ok=True)
        urlretrieve(
            f"https://www2.mpia-hd.mpg.de/THINGS/Data_files/{stem}_NA_MOM0_THINGS.FITS",
            path,
        )
    image, hdr = fits.getdata(path, header=True)
    image = image.squeeze().astype(float)
    assert (
        hdr["BUNIT"] == "JY/B*M/S"
        and hdr["CTYPE1"] == "RA---SIN"
        and hdr["CTYPE2"] == "DEC--SIN"
    )
    assert hdr["CROTA1"] == 0 and hdr["EPOCH"] == 2000
    theta = np.deg2rad(hdr["CROTA2"])
    iy, ix = np.nonzero(np.isfinite(image) & (image > 0))
    delta = np.deg2rad(np.array([hdr["CDELT1"], hdr["CDELT2"]]))
    # Orthographic (SIN) projection of adopted centre, relative to pointing.
    ra0, dec0, ra, dec = np.deg2rad(
        [hdr["CRVAL1"], hdr["CRVAL2"], props["ra"], props["dec"]]
    )
    east0 = np.cos(dec) * np.sin(ra - ra0)
    north0 = np.sin(dec) * np.cos(dec0) - np.cos(dec) * np.sin(dec0) * np.cos(ra - ra0)
    px = (ix + 1 - hdr["CRPIX1"]) * delta[0]
    py = (iy + 1 - hdr["CRPIX2"]) * delta[1]
    east = (px * np.cos(theta) - py * np.sin(theta) - east0) * distance * 1000
    north = (px * np.sin(theta) + py * np.cos(theta) - north0) * distance * 1000
    pa, inc = np.deg2rad([props["pa"], props["inc"]])
    major = east * np.sin(pa) + north * np.cos(pa)
    minor = (east * np.cos(pa) - north * np.sin(pa)) / np.cos(inc)
    pix = abs(hdr["CDELT1"] * hdr["CDELT2"]) * 3600**2
    beam = np.prod(props["beam"]) * np.pi / (4 * np.log(2))
    mass_per_jy = 1.36 * 2.356e5 * distance**2
    flux = image[iy, ix] * (pix / beam / 1000)
    # Half-step quadrature inside a pixel prevents checkerboards when a deprojected
    # sky pixel approaches one fine cell along the projected minor axis.
    size = (
        int(np.ceil((max(abs(major).max(), abs(minor).max()) + 4.0) / (FINE * 4))) * 4
    )
    grid = np.zeros((2 * size, 2 * size))
    offsets = (np.arange(4) + 0.5) / 4 - 0.5
    for dx in offsets:
        for dy in offsets:
            px, py = np.array([dx, dy]) * delta * distance * 1000
            de = px * np.cos(theta) - py * np.sin(theta)
            dn = px * np.sin(theta) + py * np.cos(theta)
            x = major + de * np.sin(pa) + dn * np.cos(pa)
            y = minor + (de * np.cos(pa) - dn * np.sin(pa)) / np.cos(inc)
            bins = (
                np.floor(x / FINE).astype(int) + size,
                np.floor(y / FINE).astype(int) + size,
            )
            np.add.at(grid, bins, flux * (mass_per_jy / 16))
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    return grid, dict(
        raw_fits_sha256=digest,
        signed_flux_jy_kms=np.nansum(image) * pix / beam / 1000,
        positive_flux_jy_kms=flux.sum(),
        reference_things_flux_jy_kms=props["flux"],
        distance_mpc=distance,
    )


def make_grids(cache):
    direct = distance_run()["1"]["galaxies"]
    provenance = {}
    arrays = {}
    for name in PAIR:
        grid, info = raw_map(name, direct[name]["distance_mpc"], cache)
        provenance[name] = info
        arrays[name] = grid
    np.savez_compressed(HERE / "atomic_grids.npz", **arrays)
    (HERE / "map_inputs.json").write_text(json.dumps(provenance, indent=2) + "\n")


def ring_surface(grid, lrow, rows, distance, columns=("Sigma*", "SigmaH2")):
    coords = (np.arange(len(grid)) - len(grid) // 2 + 0.5) * FINE
    radius = np.array([float(b["r"]) for b in rows]) * distance / float(lrow["Dist"])
    shell_edges = edges(radius)
    non_atomic = np.r_[sum(filled(radius, rows, c)[0] for c in columns), 0]
    surface = np.zeros_like(grid)
    offsets = ((np.arange(8) + 0.5) / 8 - 0.5) * FINE
    for dx in offsets:
        for dy in offsets:
            r = np.hypot(coords[:, None] + dx, coords[None, :] + dy)
            surface += non_atomic[np.searchsorted(shell_edges[1:], r)] / 64
    return surface, np.hypot(coords[:, None], coords[None, :])


def coarse(field, factor, offset):
    # The measured grids have a zero margin everywhere, so roll preserves exterior.
    data = np.roll(field, (offset, offset), axis=(0, 1))
    return data.reshape(
        data.shape[0] // factor, factor, data.shape[1] // factor, factor
    ).mean(axis=(1, 3))


def square_tv(surface, width):
    frame = np.pad(np.clip(surface, 0, THRESHOLD_MAX), 1)
    return (
        width * AREA_UNITS * sum(np.abs(np.diff(frame, axis=a)).sum() for a in (0, 1))
    )


def mapped_shell_response(surface, width, offset, orbit):
    """Assign each bond to its midpoint radius and integrate the spherical field."""
    n = len(surface)
    bounds = (np.arange(n + 1) - n / 2) * width - offset
    centres = (bounds[:-1] + bounds[1:]) / 2
    frame = np.pad(np.clip(surface, 0, THRESHOLD_MAX), 1)
    xjump = abs(np.diff(frame[:, 1:-1], axis=0))
    yjump = abs(np.diff(frame[1:-1, :], axis=1))
    ex, cy = np.nonzero(xjump)
    cx, ey = np.nonzero(yjump)
    radii = np.r_[np.hypot(bounds[ex], centres[cy]), np.hypot(centres[cx], bounds[ey])]
    raw = np.r_[xjump[ex, cy], yjump[cx, ey]] * width * AREA_UNITS
    rank = np.argsort(radii)
    r = radii[rank]
    w = raw[rank] / raw.sum()
    if min(r) <= 0:
        raise ValueError("Positive grid-edge radius required")
    knots = np.r_[0, np.unique(r[r < orbit]), orbit]
    inside = np.searchsorted(r, knots, side="right")
    inner = np.r_[0, np.cumsum(w)]
    outer = np.r_[np.cumsum((w / r)[::-1])[::-1], 0]
    field = inner[inside] + knots * outer[inside]
    factor = (
        np.diff(knots)
        @ (field[:-1] ** 2 + field[:-1] * field[1:] + field[1:] ** 2)
        / (3 * orbit)
    )
    return float(raw[radii <= orbit].sum()), float(raw.sum() * factor)


def run():
    provenance = json.loads((HERE / "map_inputs.json").read_text())
    direct = distance_run()["1"]
    profiles, _ = enroll()
    reading = {}
    with np.load(HERE / "atomic_grids.npz") as data:
        for n in PAIR:
            raw = data[n]
            info = provenance[n]
            d = direct["galaxies"][n]
            meta, lrow, rows, _, _ = profiles[n]
            stars, rr = ring_surface(raw, lrow, rows, d["distance_mpc"])
            hi = raw / (FINE**2 * AREA_UNITS)
            total_hi = (
                1.36e9 * meta["MHI_gsun"] * (d["distance_mpc"] / meta["dist_mpc"]) ** 2
            )
            result = {
                "atomic_He_msun": raw.sum(),
                "fraction_total_sparc_atomic": raw.sum() / total_hi,
                "atomic_mass_fraction_beyond_orbit": raw[rr > d["orbit_kpc"]].sum()
                / raw.sum(),
                "atomic_mass_fraction_beyond_optical_aperture": raw[
                    rr > d["edge_kpc"]
                ].sum()
                / raw.sum(),
                "angular_field_stellar_molecular_msun": stars.sum()
                * FINE**2
                * AREA_UNITS,
                "grids": {},
            }
            for width in (0.8, 1.6):
                factor = round(width / FINE)
                for offset in (0, factor // 2):
                    b = coarse(stars, factor, offset)
                    h = coarse(hi, factor, offset)
                    gas_norm = (
                        info["reference_things_flux_jy_kms"]
                        / info["positive_flux_jy_kms"]
                    )
                    combined = b + h
                    enclosed, transported = mapped_shell_response(
                        combined, width, offset * FINE, d["orbit_kpc"]
                    )
                    result["grids"][f"{width:g}kpc_{offset * FINE:g}offset"] = {
                        "baryon_TV": square_tv(combined, width),
                        "baryon_TV_interior_seams": enclosed,
                        "baryon_TV_finite_shell": transported,
                        "atomic_TV": square_tv(h, width),
                        "stellar_molecular_TV": square_tv(b, width),
                        "max_surface_msun_pc2": combined.max(),
                        "baryon_TV_with_global_THINGS_flux": square_tv(
                            b + h * gas_norm, width
                        ),
                        "baryon_TV_with_global_SPARC_HI": square_tv(
                            b + h * total_hi / raw.sum(), width
                        ),
                    }
            reading[n] = result
    a, b = PAIR
    ratios = {
        grid: {
            key: entry[key] / reading[b]["grids"][grid][key]
            for key in entry
            if "TV" in key
        }
        for grid, entry in reading[a]["grids"].items()
    }
    return {
        "profiles": reading,
        "source_ratios": ratios,
        "observed_excess_ratio": direct["Y_ratio"],
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rebuild-from-fits", type=Path)
    args = parser.parse_args()
    if args.rebuild_from_fits:
        make_grids(args.rebuild_from_fits)
    (HERE / "maps.json").write_text(json.dumps(run(), indent=2, allow_nan=False) + "\n")
