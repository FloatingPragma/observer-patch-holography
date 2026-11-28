"""Predict patch scales from the first occupied search aperture on each image."""

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path

import numpy as np
from scipy.signal import fftconvolve

from . import closed_inventory as inventory
from . import tidal_scales as scales
from .tidal_scale_kernel import annular_force

HERE = scales.HERE
WIDTH = inventory.rp.galaxy.WIDTH
MAPS = HERE / "tidal_neighbor_maps.npz"


def search_bins(shape):
    # The instrument reports the first occupied aperture, rounded up to
    # a quarter octave in distance on the square lattice.
    end = sum((s - 1) ** 2 for s in shape)
    if end == 0:
        return np.array([], dtype=int)
    return np.unique(
        [min(end, int(2 ** (j / 2))) for j in range(math.ceil(2 * math.log2(end)) + 1)]
        + [end]
    )


def peer_mass(w, rings, n, cell_mass, bins):
    """Expected annular mass in each reported bin conditional on the emitting cell.
    The last row retains emitters with no second return anywhere in the field.
    """
    assert np.all((w >= 0) & (w < 1))
    log_miss = np.log1p(-w)
    surviving = np.ones_like(w)
    stock = w * cell_mass
    out = []

    def collect(p):
        return np.bincount(rings.ravel(), (stock * p).ravel(), minlength=n)[:n]

    for k in bins:
        if k == bins[-1]:
            log_survival = log_miss.sum() - log_miss
        else:
            half = math.isqrt(k)
            a = np.arange(-half, half + 1)
            r2 = a[:, None] ** 2 + a[None, :] ** 2
            contacts = (0 < r2) & (r2 <= k)
            log_survival = fftconvolve(log_miss, contacts, mode="same")
        next_survival = np.exp(np.minimum(0, log_survival))
        assert np.max(next_survival - surviving) < 3e-12
        next_survival = np.minimum(surviving, next_survival)
        out.append(collect(surviving - next_survival))
        surviving = next_survival
    out.append(collect(surviving))
    result = np.array(out)
    np.testing.assert_allclose(result.sum(0), collect(np.ones_like(w)), rtol=1e-12)
    return result


def draw_check(w, bins, expected, seed):
    # Size-bias the emitting cell by its return probability, independently
    # draw all other cells, and read the first occupied aperture on that board.
    coords = np.argwhere(w > 0)
    p = w[w > 0]
    rng = np.random.default_rng(seed)
    counts = np.zeros(len(bins) + 1, int)
    draws = 65536
    for _ in range(draws // 256):
        root = rng.choice(len(p), 256, p=p / p.sum())
        returned = rng.random((256, len(p))) < p
        returned[np.arange(256), root] = False
        distance = ((coords[root, None] - coords) ** 2).sum(2)
        first = np.where(returned, distance, bins[-1] + 1).min(1)
        counts += np.bincount(np.searchsorted(bins, first), minlength=len(counts))
    cdf = np.cumsum(expected) / expected.sum()
    difference = float(np.max(abs(counts.cumsum() / draws - cdf)))
    allowance = math.sqrt(math.log(2 * (len(bins) + 1) / 0.0001) / (2 * draws))
    assert difference < allowance
    return dict(
        seed=seed,
        draws=draws,
        histogram=counts.tolist(),
        max_cdf_error=difference,
        simultaneous_allowance=allowance,
    )


def acquire(cache, inputs):
    rp = inventory.rp
    peak = json.loads((HERE / "record_profile_maps.json").read_bytes())["peak_msun_pc2"]
    maps = {}
    checks = {}
    for name, tag, img, rings, mask in rp.observed_cells(cache, inputs):
        y, x = np.nonzero(mask)
        box = np.s_[y.min() : y.max() + 1, x.min() : x.max() + 1]
        img = np.where(mask, img, 0)[box]
        rings = rings[box]
        bins = search_bins(img.shape)
        maps[f"{name}_{tag}_bins"] = bins
        for cap in rp.CAPACITY:
            full = cap * peak[tag]
            p = img / full
            w = 6 * p**2 * (1 - p) ** 2 * (448 / 2110)
            raw = peer_mass(
                w,
                rings,
                len(inputs[name]["radius_kpc"]),
                WIDTH**2 * rp.maps.AREA_UNITS * full / 4,
                bins,
            )
            key = f"{name}_{tag}_{cap}"
            maps[key] = raw
            if name == "DDO154" and tag == "mapped_baryons" and cap in (1, 8):
                checks[key] = draw_check(w, bins, raw.sum(1), 7192026 + cap)
    return maps, checks


def score(maps, inputs):
    profiles, mu = scales.source_data()
    stock_maps = scales.pixel_maps()
    G = inventory.rp.galaxy.transfer.budget.local.profiles.galaxies.G
    scores, details, free_gains = {}, {}, {}
    max_delta = 0.0
    for name, tab in inputs.items():
        r, v, e, baryon = [
            np.array(tab[k])
            for k in ("radius_kpc", "observed_kms", "error_kms", "baryons_kms")
        ]
        # Tracers can have distinct image bounds, so use the union of bin edges.
        bins = np.unique(
            np.concatenate([maps[f"{name}_{tag}_bins"] for tag in inventory.rp.TRACERS])
        )
        for label, (knots, K, refs) in profiles.items():
            for native in (1, 6):
                fine, coarse = (
                    np.array(
                        [
                            annular_force(
                                r, WIDTH * math.sqrt(int(k)) / native, knots, K, order
                            )
                            for k in bins
                        ]
                    )
                    for order in (12, 6)
                )
                for tag in inventory.rp.TRACERS:
                    b = maps[f"{name}_{tag}_bins"]
                    select = np.searchsorted(bins, b)
                    for cap in inventory.rp.CAPACITY:
                        stock = maps[f"{name}_{tag}_{cap}"]
                        expected = stock_maps[f"{tag}/{cap}"][name]
                        np.testing.assert_allclose(
                            stock.sum(0).cumsum(), expected, rtol=2e-11, atol=1e-8
                        )
                        conv = np.einsum("bij,bj->i", fine[select], stock[:-1])
                        coarse_conv = np.einsum("bij,bj->i", coarse[select], stock[:-1])
                        # Diagnostic only: does the imaged radial shape allow
                        # any one nonnegative normalization across a galaxy?
                        bands = scales.budget.passage_images.intervals(
                            G * r * conv,
                            np.maximum(v - 2 * e, 0) ** 2 - baryon**2,
                            (v + 2 * e) ** 2 - baryon**2,
                        )
                        shape_key = f"{label}/{native}/{tag}/{cap}"
                        free_gains.setdefault(shape_key, {})[name] = [
                            (max(0.0, lo), hi) for lo, hi in bands
                        ]
                        for calibration, D in refs.items():
                            factor = r * (G * float(mu / D))
                            v2 = baryon**2 + factor * conv
                            prediction = np.sqrt(np.maximum(0, v2))
                            hit = (v2 >= 0) & (abs(prediction - v) <= 2 * e)
                            low_order = baryon**2 + factor * coarse_conv
                            other_hit = (low_order >= 0) & (
                                abs(np.sqrt(np.maximum(0, low_order)) - v) <= 2 * e
                            )
                            assert np.array_equal(hit, other_hit)
                            max_delta = max(
                                max_delta,
                                float(np.max(abs(factor * (conv - coarse_conv)))),
                            )
                            key = f"{label}/{native}/{calibration}/{tag}/{cap}"
                            row = scores.setdefault(
                                key, dict(within_two=0, complete=0, undefined=0)
                            )
                            row["within_two"] += int(hit.sum())
                            row["complete"] += int(hit.all())
                            row["undefined"] += int((v2 < 0).sum())
                            if (label, native, calibration, tag, cap) == (
                                "source",
                                6,
                                6,
                                "mapped_baryons",
                                1,
                            ):
                                mass = stock.sum(1)
                                q = np.searchsorted(
                                    mass[:-1].cumsum(), mass[:-1].sum() * 0.5
                                )
                                details[name] = dict(
                                    registered_peer_fraction=float(
                                        1 - mass[-1] / mass.sum()
                                    ),
                                    median_registered_hop_kpc=float(
                                        WIDTH * math.sqrt(int(b[q])) / native
                                    ),
                                    radius_kpc=r.tolist(),
                                    predicted_v2_kms2=v2.tolist(),
                                    observed_kms=v.tolist(),
                                    error_kms=e.tolist(),
                                    within_two=int(hit.sum()),
                                )
    shape = {}
    for key, galaxies in free_gains.items():
        all_bands = [pair for row in galaxies.values() for pair in row]
        best, _ = scales.budget.passage_images.best_gain(all_bands)
        whole = [
            name
            for name, bands in galaxies.items()
            if max(lo for lo, hi in bands) <= min(hi for lo, hi in bands)
        ]
        largest_ordinary = max(
            float(mu / d) for d in profiles[key.split("/")[0]][2].values()
        )
        shape[key] = dict(
            best_common_gain_rings=int(best),
            whole_with_individual_gain=whole,
            gain_relative_to_largest_ordinary={
                name: [
                    max(lo for lo, hi in galaxies[name]) / largest_ordinary,
                    min(hi for lo, hi in galaxies[name]) / largest_ordinary,
                ]
                for name in whole
            },
        )
    return dict(
        summary=scores,
        shape_only=shape,
        nominated=details,
        largest_quadrature_shift_in_v2=max_delta,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--image-cache", type=Path)
    args = parser.parse_args()
    inputs, *_ = inventory.rp.galaxy.transfer.budget.load_inputs()
    path = HERE / "tidal_neighbors.json"
    if args.image_cache:
        maps, checks = acquire(args.image_cache, inputs)
        if args.write:
            np.savez_compressed(MAPS, **maps)
        else:
            with np.load(MAPS) as old:
                assert set(old.files) == set(maps)
                for key, v in maps.items():
                    np.testing.assert_allclose(v, old[key], rtol=3e-11, atol=1e-5)
    else:
        checks = json.loads(path.read_bytes())["independent_draws"]
    with np.load(MAPS) as archive:
        result = score(dict(archive), inputs)
    result.update(
        independent_draws=checks,
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                MAPS,
                HERE / "closed_inventory_maps.json",
                HERE / "TIDAL_SCALES.md",
                HERE / "pair_galaxy_maps.json",
                HERE / "record_profile_maps.json",
                HERE / "source_profile.json",
                scales.budget.LAB / "exterior_tide.json",
                scales.budget.LAB / "moving_records.json",
            )
        },
    )
    if args.write:
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        inventory.rp.galaxy.check(json.loads(path.read_bytes()), result)
    print(
        {
            k: max(v[k] for v in result["summary"].values())
            for k in ("within_two", "complete")
        }
    )


if __name__ == "__main__":
    main()
