"""Place the measured compensated passage-time profile on the enrolled galaxy images."""

import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from . import closed_inventory as inventory
from .source_profile import galaxies

rp = inventory.rp
HERE = inventory.HERE
SCALES = (0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 10.25, 12.8, 25.6, 51.2)


def enclosed(radius, seed, h):
    # One isotropic shell about each off-center seed. The image centers
    # straddle the origin, so every enrolled galactocentric offset is nonzero.
    d = seed[None, :]
    R = radius[:, None]
    assert np.all(d > 0)
    if h == 0:
        return (d <= R).astype(float)
    partial = np.clip((R**2 - (d - h) ** 2) / (4 * h * d), 0, 1)
    return np.where(R >= d + h, 1, np.where(R <= abs(d - h), 0, partial))


def intervals(field, lo, hi):
    rows = []
    for x, a, b in zip(field, lo, hi, strict=True):
        if x:
            rows.append(sorted((float(a / x), float(b / x))))
        elif a <= 0 <= b:
            rows.append([-np.inf, np.inf])
        else:
            rows.append([np.inf, -np.inf])
    return rows


def best_gain(rows):
    base = sum(a == -np.inf and b == np.inf for a, b in rows)
    # Count ends after starts at a tied endpoint, since intervals are closed.
    events = sorted(
        (value, order)
        for a, b in rows
        if np.isfinite(a) and a <= b
        for value, order in ((a, 0), (b, 1))
    )
    best, count, gain = base, base, 0.0
    for x, order in events:
        count += 1 if order == 0 else -1
        if count > best:
            best, gain = count, x
    assert sum(a <= gain <= b for a, b in rows) == best
    return best, gain


def run(cache):
    inputs, *_ = rp.galaxy.transfer.budget.load_inputs()
    clock = inventory.LAB / "closed_passage.json"
    M = list(map(Q, json.loads(clock.read_bytes())["enclosed_flux"]))
    shell = [M[0] / M[0], *[(y - x) / M[0] for x, y in zip(M[:-1], M[1:], strict=True)]]
    assert shell[0] == 1 and max(shell[1:]) <= 0 and sum(shell) == 0
    peak = json.loads((HERE / "record_profile_maps.json").read_bytes())["peak_msun_pc2"]
    previous = json.loads((HERE / "closed_inventory_maps.json").read_bytes())
    panels = {}
    for name, tracer, img, _, mask in rp.observed_cells(cache, inputs):
        tab = inputs[name]
        radius = np.array(tab["radius_kpc"])
        lo, hi = rp.galaxy.transfer.budget.bounds(tab)
        axis = (np.arange(len(img)) + 0.5 - len(img) // 2) * rp.galaxy.WIDTH
        seed, indices = np.unique(
            np.hypot(axis[:, None], axis[None, :])[mask], return_inverse=True
        )
        weights = []
        for cap in rp.CAPACITY:
            full = peak[tracer] * cap
            p = img[mask] / full
            assert np.all((0 <= p) & (p <= 1))
            mass = (
                6
                * p**2
                * (1 - p) ** 2
                * (448 / 2110)
                * full
                / 4
                * rp.galaxy.WIDTH**2
                * rp.maps.AREA_UNITS
            )
            weights.append(np.bincount(indices, weights=mass, minlength=len(seed)))
        weights = np.array(weights).T
        core = enclosed(radius, seed, 0)
        isolated = core @ weights
        for j, cap in enumerate(rp.CAPACITY):
            np.testing.assert_allclose(
                isolated[:, j],
                previous[tracer][str(cap)]["closed"][name],
                rtol=2e-12,
                atol=1e-6,
            )
        assert not enclosed(radius, seed, 2 * (radius[-1] + seed[-1])).any()
        for length in SCALES:
            # Exact compensation makes a fully enclosed seed contribute zero.
            # Subtract the two aperture fractions before applying each weight.
            fraction = np.zeros_like(core)
            for hop, w in enumerate(shell[1:], 1):
                if hop * length > radius[-1] + seed[-1]:
                    fraction += float(-sum(shell[hop:])) * core
                    break
                fraction += float(-w) * (core - enclosed(radius, seed, hop * length))
            fields = (fraction @ weights) * (galaxies.G / radius[:, None])
            for j, cap in enumerate(rp.CAPACITY):
                key = f"{tracer}/{cap}/{length}"
                panels.setdefault(key, {})[name] = intervals(fields[:, j], lo, hi)
    report = {}
    for key, curves in panels.items():
        rows = [v for gal in curves.values() for v in gal]
        best, gain = best_gain(rows)
        report[key] = dict(
            best_common_gain=gain,
            within_two=best,
            total=len(rows),
            complete_with_own_gain=[
                name for name, gal in curves.items() if best_gain(gal)[0] == len(gal)
            ],
            complete_at_best_gain=[
                name
                for name, gal in curves.items()
                if all(a <= gain <= b for a, b in gal)
            ],
        )
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                clock,
                HERE / "closed_inventory_maps.json",
                HERE / "record_profile_maps.json",
                HERE / "source_profile.json",
                HERE / "pair_galaxy_maps.json",
                HERE / "PASSAGE_MASS.md",
            )
        },
        models=report,
        best_points=max(v["within_two"] for v in report.values()),
        most_complete_with_own_gains=max(
            len(v["complete_with_own_gain"]) for v in report.values()
        ),
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--image-cache", type=Path, required=True)
    cli.add_argument("--write", action="store_true")
    opt = cli.parse_args()
    result = run(opt.image_cache)
    path = HERE / "passage_images.json"
    if opt.write:
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        rp.galaxy.check(json.loads(path.read_bytes()), result)
    print(result["best_points"], result["most_complete_with_own_gains"])
