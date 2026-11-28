"""Match rotation-curve clock deficits to detected H I line widths in the same disks."""

import argparse
from hashlib import sha256
import json
from pathlib import Path
from urllib.request import urlretrieve

import numpy as np

from . import pair_galaxies as galaxies

HERE = galaxies.HERE
C = 299792.458
CATALOG = "https://www2.mpia-hd.mpg.de/THINGS/Data_files/"


def maps(cache):
    from astropy.io import fits
    from astropy.wcs import WCS

    inputs, _, meta, _ = galaxies.transfer.budget.load_inputs()
    result = {}
    for name, tab in inputs.items():
        stem = name.replace("NGC", "NGC_").replace("IC", "IC_")
        fluxfile = cache / (stem + ".fits")
        widthfile = cache / (stem + "_NA_MOM2.fits")
        cache.mkdir(parents=True, exist_ok=True)
        for order, path in ((0, fluxfile), (2, widthfile)):
            if not path.exists():
                urlretrieve(CATALOG + f"{stem}_NA_MOM{order}_THINGS.FITS", path)
        a, ha = fits.getdata(fluxfile, header=True)
        b, hb = fits.getdata(widthfile, header=True)
        for k in (
            "CTYPE1",
            "CTYPE2",
            "CRVAL1",
            "CRVAL2",
            "CRPIX1",
            "CRPIX2",
            "CDELT1",
            "CDELT2",
            "CROTA2",
            "EPOCH",
        ):
            assert ha[k] == hb[k]
        assert (
            a.shape == b.shape
            and ha["BUNIT"].strip() == "JY/B*M/S"
            and hb["BUNIT"].strip() == "METR/SEC"
        )
        image = a.squeeze()
        y, x = np.nonzero(np.isfinite(image) & (image > 0))
        flux = image[y, x].astype(float)
        width = b.squeeze()[y, x].astype(float) / 1000
        ra, dec = np.deg2rad(WCS(ha, naxis=2, fix=False).all_pix2world(x, y, 0))
        geom = galaxies.maps.GEOMETRY[name]
        ra0, dec0, pa, inc = np.deg2rad([geom[k] for k in ("ra", "dec", "pa", "inc")])
        distance = 1000 * meta[name]["dist_mpc"]
        east = distance * np.cos(dec) * np.sin(ra - ra0)
        north = distance * (
            np.sin(dec) * np.cos(dec0) - np.cos(dec) * np.sin(dec0) * np.cos(ra - ra0)
        )
        radius = np.hypot(
            east * np.sin(pa) + north * np.cos(pa),
            (east * np.cos(pa) - north * np.sin(pa)) / np.cos(inc),
        )
        rings = np.asarray(tab["radius_kpc"])
        edges = np.r_[
            max(0, rings[0] - (rings[1] - rings[0]) / 2),
            (rings[:-1] + rings[1:]) / 2,
            rings[-1] + (rings[-1] - rings[-2]) / 2,
        ]
        rows = []
        idx = np.searchsorted(edges, radius, side="right") - 1
        for j, r in enumerate(rings):
            where = idx == j
            valid = where & np.isfinite(width) & (width >= 0)
            row = dict(
                radius_kpc=r,
                positive_flux_pixels=int(where.sum()),
                width_pixels=int(valid.sum()),
            )
            if np.any(valid):
                row.update(
                    covered_positive_flux=float(sum(flux[valid]) / sum(flux[where])),
                    width_median_kms=float(np.median(width[valid])),
                    width_rms_kms=float(
                        np.sqrt(np.average(width[valid] ** 2, weights=flux[valid]))
                    ),
                    width_max_kms=float(np.max(width[valid])),
                )
            rows.append(row)
        result[name] = dict(
            sha256={
                p.name: sha256(p.read_bytes()).hexdigest()
                for p in (fluxfile, widthfile)
            },
            width_url=CATALOG + f"{stem}_NA_MOM2_THINGS.FITS",
            distance_Mpc=distance / 1000,
            rings=rows,
        )
    return result


def inward_integrals(radius, residual):
    """Piecewise linear v^2-v_b^2 in log(r), retaining negative segments."""
    area = np.diff(np.log(radius)) * (residual[1:] + residual[:-1]) / 2
    return np.r_[np.cumsum(area[::-1])[::-1], 0.0]


def run(data):
    from .closed_inventory import LAB

    inputs, *_ = galaxies.transfer.budget.load_inputs()
    clock = json.loads((LAB / "renewal_clock.json").read_bytes())
    multipliers = {}
    for t, reader in clock["readers"].items():
        flags = reader["exact"]["flags"]
        ratio = flags[0][1] / flags[2][0]
        multipliers[t] = np.nextafter((ratio + 2 + 1 / ratio) / 4, np.inf)
    result = {}
    for name, tab in inputs.items():
        r = np.array(tab["radius_kpc"])
        v, e, b = [
            np.asarray(tab[k]) for k in ("observed_kms", "error_kms", "baryons_kms")
        ]
        low = np.maximum(0, v - 2 * e) ** 2 - b * b
        central = v * v - b * b
        central, lower = [inward_integrals(r, values) for values in (central, low)]
        start = int(np.argmin(abs(r - r[-1] / 2)))
        rows = []
        for j in range(start, len(r)):
            image = data[name]["rings"][j]
            assert image["radius_kpc"] == r[j]
            row = dict(
                radius_kpc=r[j],
                integral_central_kms2=central[j],
                integral_lower_kms2=lower[j],
            )
            for power, tag in ((1, "renewal"), (2, "rms_lapse")):
                need = C * np.sqrt(np.expm1(power * max(0, lower[j]) / C**2))
                out = dict(required_rms_kms=need)
                if "width_rms_kms" in image:
                    sigma = image["width_rms_kms"]
                    measured = (sigma / C) ** 2
                    outer = data[name]["rings"][-1].get("width_rms_kms")
                    if outer is not None:
                        out["potential_from_both_detected_widths_kms2"] = (
                            C**2
                            / power
                            * (np.log1p(measured) - np.log1p((outer / C) ** 2))
                        )
                    out.update(
                        observed_rms_kms=sigma,
                        observed_max_kms=image["width_max_kms"],
                        effective_samples_needed=(need / sigma) ** 2 if sigma else None,
                        allowed_potential_fraction=C**2
                        * np.log1p(measured)
                        / (power * lower[j])
                        if lower[j] > 0
                        else None,
                    )
                row[tag] = out
            # Tick-weighted emission samples the event law. For this ternary
            # clock the ratio of excess arithmetic/harmonic rate to event
            # CV^2 increases with f and peaks at equal endpoint weights.
            row["tick_weighted_minimum_rms_kms"] = {
                t: row["renewal"]["required_rms_kms"] / np.sqrt(m)
                for t, m in multipliers.items()
            }
            rows.append(row)
        result[name] = dict(
            last_radius_kpc=r[-1], selected_inner_index=start, outer_rings=rows
        )
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                HERE / "CLOCK_WIDTH.md",
                HERE / "clock_width_maps.json",
                HERE / "source_profile.json",
                LAB / "renewal_clock.json",
            )
        },
        native_tick_weighted_variance_multipliers=multipliers,
        galaxies=result,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--image-cache", type=Path)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    src = HERE / "clock_width_maps.json"
    if args.image_cache:
        measurements = maps(args.image_cache)
        if args.write:
            src.write_text(json.dumps(measurements, indent=2, allow_nan=False) + "\n")
        else:
            galaxies.check(json.loads(src.read_bytes()), measurements)
    result = run(json.loads(src.read_bytes()))
    dst = HERE / "clock_width.json"
    if args.write:
        dst.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        galaxies.check(json.loads(dst.read_bytes()), result)
    for name, row in result["galaxies"].items():
        print(name, row["outer_rings"][0]["renewal"])
