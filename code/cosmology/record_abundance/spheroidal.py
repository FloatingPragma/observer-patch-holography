"""Translate ordinary Draco apertures into bounds on site capacity and recovery."""

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import gzip
from hashlib import sha256
import json
from math import ceil, comb, cos, log, pi, radians, sqrt, tan
from urllib.request import urlopen

import numpy as np
from scipy.optimize import brentq
from scipy.spatial.distance import cdist

from . import clustering, mass_ratio, stationary

HERE = mass_ratio.HERE
RECEIPT = HERE / "spheroidal.json"
LAB = HERE / "draco_lab.npz"
BASE = "https://vizier.cfa.harvard.edu/ftp/cats/VIII/76/bvmaps/"
MAPS = {
    85.5: "57cd8c0f2abcb597b06ee39a8e7dbab9e7250da7df8fedf6ebdc49e310c147d4",
    86.0: "b6c22ad2cf577952ae118e53b4b221b77b52cb8b30c01ef74300d544fe60b0ca",
    86.5: "32a80c3a9437de5babef9e7bc41538c1fbf3aedc83cdc5073fb87ef523f26e4b",
    87.0: "4a9931d8df8e51c9dbea594db60062dec48a700075ddca643bb5ccec655920b0",
    87.5: "68d4d235466df5e042bcfab3341355fe1903277c7d797e1ab6b3371da93dcc65",
}
LAT = np.arange(33.5, 36.1, 0.5)
SOURCE = (
    (
        "0908.2995v3",
        "ms-mn2e.tex",
        "429edae10d5d9fd2673573ed72fbb2a803f28b5e016b7ca47ef4b0525f6f1a21",
    ),
    (
        "astro-ph/0301228v1",
        "ms.tex",
        "6965673a5960b9ec35906cc55a88ac52624a57700ae0b60b00b6a796f9d42868",
    ),
    (
        "1410.0028v1",
        "ms.tex",
        "dc6249c083df4ed9b73cbb78ea2a3ded5f5a073a6595ae74e309419d36ed24f4",
    ),
    (
        "2407.07769v1",
        "main.tex",
        "bb7cb4b6cc468fd748b2a3bb73d5f663efc145115e4b8d9438457743f0f49c03",
    ),
)
# Rounded Galactic position of Draco from ICRS 17h20m12.4s, +57d54m55s.
CENTER = np.array([86.367829, 34.721708])
DIST = 76000.0  # pc
PC = 3.08567758e18  # cm
MP = 1.6726219e-24
MSUN = 1.98847e33
FWHM = 40 / 60  # degrees, the regridded northern cube's effective FWHM
WINDOWS = (-274, -400, -335, -210)


def fetch_maps():
    def retrieve(lon):
        name = f"L{int(10 * lon):04}.fits"
        with urlopen(BASE + name + ".gz", timeout=45) as src:
            raw = gzip.decompress(src.read())
        if sha256(raw).hexdigest() != MAPS[lon]:
            raise ValueError(f"Unexpected LAB map {name}")
        return fits_slice(raw)

    with ThreadPoolExecutor(max_workers=5) as pool:
        maps = list(pool.map(retrieve, MAPS))
    velocity = maps[0][0]
    for v, _ in maps:
        np.testing.assert_array_equal(v, velocity)
    return velocity, np.stack([info[1] for info in maps])


def fits_slice(raw):
    """Read just the signed 16-bit primary image of the hash-pinned LAB product."""
    cards = {}
    for i in range(0, len(raw), 80):
        card = raw[i : i + 80].decode("ascii")
        if card[:8].strip() == "END":
            break
        if card[8:10] == "= ":
            cards[card[:8].strip()] = card[10:30].strip().strip("'").strip()
    assert (
        cards["CTYPE1"] == "VELO-LSR"
        and cards["BITPIX"] == "16"
        and cards["BUNIT"] == "K"
    )
    assert [int(cards[k]) for k in ("NAXIS1", "NAXIS2", "NAXIS3")] == [893, 361, 1]
    start = 2880 * (i // 2880 + 1)
    raw_image = np.frombuffer(raw, dtype=">i2", offset=start, count=893 * 361).reshape(
        361, 893
    )
    row = (
        (LAT - float(cards["CRVAL2"])) / float(cards["CDELT2"])
        + float(cards["CRPIX2"])
        - 1
    ).astype(int)
    v = (
        (np.arange(893) + 1 - float(cards["CRPIX1"])) * float(cards["CDELT1"])
        + float(cards["CRVAL1"])
    ) / 1000
    select = (v > -440) & (v < -175)
    image = raw_image[row, :][:, select]
    assert not (image == int(cards["BLANK"])).any()
    return v[select], image * float(cards["BSCALE"]) + float(cards["BZERO"])


def verify_papers():
    with ThreadPoolExecutor(max_workers=4) as pool:
        wolf, wham, gbt, pm = pool.map(clustering.source, SOURCE)
    row = next(line for line in wolf.splitlines() if line.startswith("Draco  (206)"))
    assert all(
        s in row for s in (r"2.2^{+0.7}_{-0.6}", r"2.11^{+0.31}_{-0.31}", "291", "76")
    )
    assert all(
        s in wham
        for s in ("0.024", "0.003 and 0.002 R", "less than 10", "2.5~", "8000")
    )
    row = next(line for line in gbt.splitlines() if line.startswith("Draco "))
    assert [s.strip() for s in row.split("&")][2:5] == ["-274", "2.8E5", "10"]
    assert "133" in row and r"24.19^{+ 6.31}_{- 2.97}" in pm


def warm_mass(radius, distance):
    # The tracer's velocity-resolved limit is beam-averaged. Add maximum
    # stellar absorption, correct 10% foreground attenuation, and fill all volume.
    beam_radius = distance * tan(radians(0.5))
    projected_min = (
        1 if radius <= beam_radius else 1 - sqrt(1 - (beam_radius / radius) ** 2)
    )
    em = 2.5 * (0.024 + 0.003) / 0.9
    volume = 4 * pi * radius**3 / 3
    upper = (
        1.33
        * MP
        * PC**3
        / MSUN
        * sqrt(em * pi * beam_radius**2 * volume / projected_min)
    )
    return {
        "warm_mass_ceiling_msun": upper,
        "minimum_spherical_shell_light_fraction_in_beam": projected_min,
    }


def neutral_budget(radius, distance, velocity, brightness):
    """Use an explicit lower bound on the beam sum everywhere in the aperture."""
    theta = np.degrees(np.arcsin(radius / distance))
    # Bound the longitude range of a great-circle disk. A clipped .01-degree
    # mesh covers it. Each patch gets a conservative upper distance to each beam.
    lon_extent = theta / cos(radians(CENTER[1] + theta))
    step = 0.01
    offsets = [
        np.arange(-ceil(limit / step), ceil(limit / step) + 1) * step
        for limit in (lon_extent, theta)
    ]
    dl, db = np.meshgrid(*offsets)
    # Retain all patches that could meet the disk by a lower distance bound.
    bmax = radians(CENTER[1] + theta + step)
    # Use the spherical chord as a distance lower bound for patch centers.
    lat = np.radians(CENTER[1] + db)
    hav = (
        np.sin(np.radians(db) / 2) ** 2
        + np.cos(lat) * cos(radians(CENTER[1])) * np.sin(np.radians(dl) / 2) ** 2
    )
    keep = 2 * np.degrees(np.arcsin(np.sqrt(hav))) <= theta + step
    centers = np.column_stack((dl[keep] + CENTER[0], db[keep] + CENTER[1]))
    lon_grid, lat_grid = np.meshgrid(list(MAPS), LAT, indexing="ij")
    positions = np.column_stack((lon_grid.ravel(), lat_grid.ravel()))
    # Include nearest grid samples needed by patches; no beam is used twice.
    selected = np.unique(
        cdist(centers * [cos(bmax), 1], positions * [cos(bmax), 1]).argmin(axis=1)
    )
    beams = positions[selected]
    # Straight interpolation in lon/lat bounds great-circle distance above.
    lon_scale = cos(radians(min(centers[:, 1].min(), beams[:, 1].min()) - step))
    distance_upper = cdist(centers * [lon_scale, 1], beams * [lon_scale, 1]) + step
    gain = float(np.exp(-4 * log(2) * (distance_upper / FWHM) ** 2).sum(axis=1).min())
    beam_area_pc2 = pi / (4 * log(2)) * (distance * radians(FWHM)) ** 2
    # Allow extreme positive covariance across positions. Off-line windows are
    # evaluated on the same beams without subtracting foreground or fitting a baseline.
    spectra = brightness.reshape(-1, len(velocity))[selected]
    dv = float(np.median(np.diff(velocity)))
    windows = {}
    for center in WINDOWS:
        mask = abs(velocity - center) <= 30
        width = int(mask.sum()) * dv
        line = float(dv * spectra[:, mask].sum())
        random = 5 * 0.09 * sqrt(2 * width)
        baseline = 0.04 * width
        scale = 1.823e18 * beam_area_pc2 * PC**2 * MP / (gain * MSUN)
        windows[str(center)] = {
            "beam_sum_K_km_per_s": line,
            "width_km_per_s": width,
            "hi_upper_msun": scale * (max(0, line) + len(beams) * (random + baseline)),
            "hi_upper_fully_correlated_msun": scale
            * (max(0, line) + len(beams) * (5 * 0.09 + 0.04) * width),
        }
    return {
        "effective_beam_fwhm_deg": FWHM,
        "beams_lon_lat_deg": beams.tolist(),
        "minimum_beam_sum": gain,
        "windows": windows,
    }


def native_witnesses(n, edges):
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import shortest_path

    graph = csr_matrix((np.ones(len(edges)), edges.T), shape=(n, n))
    distance = shortest_path(graph, directed=False, unweighted=True)
    return [np.minimum(distance[i], distance[j]) for i, j in edges]


def native_census(level, k):
    board = stationary.measure.federation(level)
    edges = np.asarray(board.external) // 12
    m = [
        int(np.count_nonzero((row > 0) & (row <= 2)))
        for row in native_witnesses(board.cells, edges)
    ]
    counts = Counter(m)
    # Given one specified endpoint active, the opposite endpoint and m unique
    # witnesses must be empty. Either endpoint may own the cut.
    per_active = (
        3
        * sum(count * comb(board.cells - j - 2, k - 1) for j, count in counts.items())
        / (len(m) * comb(board.cells - 1, k - 1))
    )
    return {
        "cells": board.cells,
        "fixed_ordinary_count": k,
        "p": k / board.cells,
        "witness_counts": {str(j): counts[j] for j in sorted(counts)},
        "expected_failed_mismatches_per_ordinary": per_active,
    }


def budget(
    radius, dyn, luminous, distance, velocity, brightness, *, guaranteed_light=None
):
    if guaranteed_light is None:
        guaranteed_light = luminous
    gas = warm_mass(radius, distance)
    atomic = neutral_budget(radius, distance, velocity, brightness)
    hi = atomic["windows"]["-274"]["hi_upper_msun"]
    gas_proxy = gas["warm_mass_ceiling_msun"] + 1.33 * hi
    p0 = stationary.screen.DENSITY
    f2 = (1 - p0) ** 11 * (2 - p0) / 2
    c = 0.11933 / 0.02242 / (1 - p0)
    baryon_ref = 0.02242 * 2.77536627e11 * 4 * pi / 3 * (radius / 1e6) ** 3
    limits = {}
    for name, f in (("constant", 1), ("one_step", (1 - p0) ** 4), ("two_step", f2)):
        threshold = dyn / (1 + c / f)
        limits[name] = {
            "maximum_dark_per_baryon": c / f,
            "minimum_baryon_msun": threshold,
            "required_stellar_ml_after_allowed_gas": max(0, threshold - gas_proxy)
            / luminous,
            "extra_baryons_at_stellar_ml_3_msun": max(
                0, threshold - gas_proxy - 3 * luminous
            ),
            "required_stellar_ml_if_no_gas": threshold / luminous,
        }
    # This inversion uses the dynamical target only to derive a necessary
    # inequality. Occupancy is still unmeasured. The mass-to-light choice is explicit.
    stars = 3 * luminous
    minimum_dark_ratio = dyn / (stars + gas_proxy) - 1

    def dark_ratio(p):
        return c * ((1 - p) ** 12 + (1 - p) ** 13) / (2 * f2)

    if minimum_dark_ratio < dark_ratio(0):
        maximum_p = brentq(lambda p: dark_ratio(p) - minimum_dark_ratio, 0, 1)
        limits["two_step"]["maximum_ordinary_occupation_at_ml_3"] = maximum_p
        limits["two_step"]["required_site_density_ratio_from_stars_at_ml_3"] = (
            p0 * (3 * guaranteed_light) / (maximum_p * baryon_ref)
        )
    return {
        "radius_pc": radius,
        "mass_proxy_msun": dyn,
        "light_used_for_inventory_lsun": luminous,
        "light_known_inside_for_density_lsun": guaranteed_light,
        "warm": gas,
        "neutral": atomic,
        "allowed_atomic_plus_warm_msun": gas_proxy,
        "proper_site_density_enhancement_minimum_from_stars_at_ml_1": p0
        * guaranteed_light
        / baryon_ref,
        "homogeneous_cosmic_baryon_mass_in_volume_msun": baryon_ref,
        "limits": limits,
    }


def run(v, data):
    outer_m = (
        900 * 24.19**2 / 4.30091e-3
    )  # Wolf uses a 291 pc sphere; Vitral uses eq. (14).
    return {
        "protocol_sha256": sha256(
            (HERE / "SPHEROIDAL_PROTOCOL.md").read_bytes()
        ).hexdigest(),
        "primary_tex": {d[0]: d[2] for d in SOURCE},
        "lab_sha256": {str(k): val for k, val in MAPS.items()},
        "native_fixed_census": [native_census(level, 8) for level in (1, 2)],
        "half_light": budget(291, 2.11e7, 1.1e5, DIST, v, data),
        "favorable_half_light_input_corner": budget(305, 1.8e7, 1.45e5, 81000, v, data),
        "outer": budget(900, outer_m, 2.2e5, DIST, v, data, guaranteed_light=1.1e5),
        "outer_lower_velocity": budget(
            900,
            900 * (24.19 - 2.97) ** 2 / 4.30091e-3,
            2.2e5,
            DIST,
            v,
            data,
            guaranteed_light=1.1e5,
        ),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--fetch",
        action="store_true",
        help="retrieve and compare the published survey maps and manuscript entries",
    )
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    if args.write and not args.fetch:
        parser.error("--write requires --fetch")
    if args.fetch:
        verify_papers()
        velocity, brightness = fetch_maps()
        if args.write:
            with LAB.open("xb") as out:
                np.savez_compressed(out, velocity=velocity, brightness=brightness)
    with np.load(LAB) as bank:
        if args.fetch:
            np.testing.assert_array_equal(velocity, bank["velocity"])
            np.testing.assert_array_equal(brightness, bank["brightness"])
        result = run(bank["velocity"], bank["brightness"])
    if args.write:
        with RECEIPT.open("x") as dest:
            json.dump(result, dest, indent=2, allow_nan=False)
            dest.write("\n")
    else:
        stationary.check(json.loads(RECEIPT.read_text()), result)
    print(
        json.dumps(
            {key: result[key]["limits"] for key in ("half_light", "outer")}, indent=2
        )
    )
