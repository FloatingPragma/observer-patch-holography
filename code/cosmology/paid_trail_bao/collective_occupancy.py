"""Do the independently imaged enrollment opportunities supply simultaneous sources?"""

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import norm

from . import tidal_neighbors as neighbor

HERE = neighbor.HERE
MAPS = HERE / "collective_fillings.npz"
MULTIPLICITIES = (1, 2, 4, 8)


def acquire(cache, inputs):
    rp = neighbor.inventory.rp
    peaks = json.loads((HERE / "record_profile_maps.json").read_bytes())[
        "peak_msun_pc2"
    ]
    out = {}
    for name, tag, img, rings, mask in rp.observed_cells(cache, inputs):
        key = f"{name}_{tag}"
        axis = (np.arange(len(img)) + 0.5 - len(img) // 2) * neighbor.WIDTH
        radii = np.hypot(axis[:, None], axis[None, :])[mask]
        assert np.array_equal(
            rings[mask], np.searchsorted(inputs[name]["radius_kpc"], radii)
        )
        out[key] = img[mask] / peaks[tag]
        out[key + "_radius_kpc"] = radii
    return out


def count_law(prob, limit=8):
    """The last bin accumulates the tail without subtracting two nearly equal values."""
    law = np.zeros(limit + 1)
    law[0] = 1.0
    for p in prob:
        law[-1] += p * law[-2]
        law[1:-1] = law[1:-1] * (1 - p) + law[:-2] * p
        law[0] *= 1 - p
    return law


def upper_dependent(prob, k):
    # Remove the j most probable cells. At least k-j of the remaining
    # indicators must be one whenever the full field contains k sources.
    descending = np.sort(prob)[::-1]
    if len(descending) < k:
        return 0.0
    return min(1.0, *(descending[j:].sum() / (k - j) for j in range(k)))


def retained_pair(prob, t):
    # P(zero or one success) for independent cells with retained returns.
    # Evaluate the singleton terms in log space, including near-certain cells.
    miss = t * np.log1p(-prob)
    total = miss.sum()
    return -np.expm1(total) - np.sum(-np.expm1(miss) * np.exp(total - miss))


def retention_time(prob):
    high = 1
    while retained_pair(prob, high) < 0.95:
        high *= 2
    low = 0
    while high - low > 1:
        mid = (high + low) // 2
        if retained_pair(prob, mid) >= 0.95:
            high = mid
        else:
            low = mid
    return dict(
        opportunities=high,
        below=float(retained_pair(prob, low)),
        at=float(retained_pair(prob, high)),
    )


def run(maps):
    rp = neighbor.inventory.rp
    inputs, *_ = rp.galaxy.transfer.budget.load_inputs()
    peaks = json.loads((HERE / "record_profile_maps.json").read_bytes())[
        "peak_msun_pc2"
    ]
    stock = neighbor.scales.pixel_maps()
    rays = sum(len(d["radius_kpc"]) for d in inputs.values())
    sigma = float(norm.isf(0.01 / rays))
    result = {}
    with np.load(neighbor.MAPS) as peers:
        for name, tab in inputs.items():
            z = (np.array(tab["observed_kms"]) - tab["baryons_kms"]) / tab["error_kms"]
            j = int(z.argmax())
            tests = {}
            for tag in rp.TRACERS:
                for cap in rp.CAPACITY:
                    fill = maps[f"{name}_{tag}"] / cap
                    p = 6 * fill**2 * (1 - fill) ** 2 * (448 / 2110)
                    law = count_law(p)
                    assert np.all(law >= 0) and abs(law.sum() - 1) < 1e-12
                    per_cell = (
                        neighbor.WIDTH**2 * rp.maps.AREA_UNITS * peaks[tag] * cap / 4
                    )
                    np.testing.assert_allclose(
                        p.sum() * per_cell, stock[f"{tag}/{cap}"][name][-1], rtol=2e-11
                    )
                    # A singleton is exactly a returned emitter whose entire
                    # remaining image is empty. Verify against the spatial search.
                    no_peer = peers[f"{name}_{tag}_{cap}"][-1].sum() / per_cell
                    np.testing.assert_allclose(law[1], no_peer, rtol=2e-11, atol=1e-13)
                    row = dict(
                        cells=len(p),
                        mean=float(p.sum()),
                        independent={
                            str(k): float(sum(law[k:])) for k in MULTIPLICITIES
                        },
                        arbitrary_correlations={
                            str(k): float(upper_dependent(p, k)) for k in MULTIPLICITIES
                        },
                    )
                    if name == "DDO154":
                        row["retained_pair_95"] = retention_time(p)
                    tests[f"{tag}/{cap}"] = row
            result[name] = dict(
                peak_excess_in_quoted_errors=float(z[j]),
                radius_kpc=tab["radius_kpc"][j],
                above_two_errors=int((z > 2).sum()),
                above_adjusted_errors=int((z > sigma).sum()),
                menus=tests,
            )
    all_disks = {}
    for k in map(str, MULTIPLICITIES):
        # Favor each galaxy separately when picking a map and capacity.
        each = {
            n: {
                key: max(v[key][k] for v in panel["menus"].values())
                for key in ("independent", "arbitrary_correlations")
            }
            for n, panel in result.items()
            if panel["above_adjusted_errors"] > 0
        }
        all_disks[k] = dict(
            per_galaxy_upper=each,
            independent_cells_and_galaxies=math.prod(
                v["independent"] for v in each.values()
            ),
            any_cross_galaxy_dependence=min(
                v["arbitrary_correlations"] for v in each.values()
            ),
        )
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                MAPS,
                neighbor.MAPS,
                HERE / "source_profile.json",
                HERE / "record_profile_maps.json",
                neighbor.scales.budget.LAB / "COLLECTIVE_TIDE.md",
            )
        },
        one_sided_bonferroni_z=sigma,
        total_rings=rays,
        galaxies=result,
        all_excess_disks=all_disks,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    cli.add_argument("--image-cache", type=Path)
    args = cli.parse_args()
    if args.image_cache:
        inputs, *_ = neighbor.inventory.rp.galaxy.transfer.budget.load_inputs()
        maps = acquire(args.image_cache, inputs)
        if args.write:
            np.savez_compressed(MAPS, **maps)
        else:
            with np.load(MAPS) as saved:
                assert set(saved) == set(maps)
                for key, v in maps.items():
                    np.testing.assert_allclose(saved[key], v, rtol=2e-12, atol=1e-15)
    with np.load(MAPS) as saved:
        result = run(saved)
    dst = HERE / "collective_occupancy.json"
    if args.write:
        dst.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        neighbor.inventory.rp.galaxy.check(json.loads(dst.read_bytes()), result)
    print(
        {
            k: [v["independent_cells_and_galaxies"], v["any_cross_galaxy_dependence"]]
            for k, v in result["all_excess_disks"].items()
        }
    )
