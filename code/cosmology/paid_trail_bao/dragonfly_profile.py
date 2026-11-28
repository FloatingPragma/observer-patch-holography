"""Forward acceleration bands from NGC 2841's measured outer disk."""

import argparse
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

import numpy as np

from cosmology.rar_deep_regime.boundary_source import maps
from . import outer_disk as outer

HERE = outer.HERE


def rings():
    inputs, _, meta, masses = outer.budget.load_inputs()
    n = outer.GALAXY
    assert meta[n]["dist_mpc"] == 14.1
    scale = Q(str(masses[n])) * 10**9
    return [
        (r, lo * scale, hi * scale)
        for r, lo, hi in outer.transport.mass_intervals(inputs[n], masses[n])
        if 20 < r <= 50
    ]


def acquire(cache):
    grid, info = maps.raw_map(outer.GALAXY, 14.1, cache)
    prior = json.loads((maps.HERE / "gas_image_inputs.json").read_bytes())[outer.GALAXY]
    assert info["raw_fits_sha256"] == prior["radio"]["raw_fits_sha256"]
    np.testing.assert_allclose(grid.sum(), prior["atomic_He_msun"], rtol=1e-10, atol=0)
    axis = (np.arange(len(grid)) + 0.5 - len(grid) // 2) * maps.FINE
    rr = np.hypot(axis[:, None], axis[None, :])
    cuts = [20.0, *(float(r) for r, _, _ in rings())]
    return dict(
        url="https://www2.mpia-hd.mpg.de/THINGS/Data_files/NGC_2841_NA_MOM0_THINGS.FITS",
        sha256=info["raw_fits_sha256"],
        distance_mpc=info["distance_mpc"],
        signed_flux_jy_kms=info["signed_flux_jy_kms"],
        total_atomic_He_msun=round(grid.sum()),
        enclosed_atomic_He_msun=[[r, round(grid[rr <= r].sum())] for r in cuts],
    )


def vertices(rows):
    # Each exact inequality is a*A + b*gain <= c. Two distinct measured
    # radii with positive exterior gas and finite upper masses bound the polygon.
    limits = [(Q(-1), Q(0), Q(0)), (Q(0), Q(-1), Q(0))]
    for _, x, lo, hi in rows:
        limits.extend([(Q(1), x, hi), (Q(-1), -x, -lo)])
    points = set()
    for (a, b, c), (d, e, f) in combinations(limits, 2):
        det = a * e - b * d
        if det:
            A, g = (c * e - b * f) / det, (a * f - c * d) / det
            if all(u * A + v * g <= t for u, v, t in limits):
                points.add((A, g))
    return sorted(points)


def run(source):
    gas = {Q(str(r)): Q(m) for r, m in source["enclosed_atomic_He_msun"]}
    rows = [(r, 4 * (gas[r] - gas[20]), lo, hi) for r, lo, hi in rings()]
    assert all(a[1] < b[1] for a, b in zip(rows, rows[1:]))
    train = [a for a in rows if a[0] <= 35]
    held = [a for a in rows if a[0] > 35]
    points = vertices(train)

    def slope(pair, lower):
        a, b = pair
        return ((b[2] - a[3]) if lower else (b[3] - a[2])) / (b[1] - a[1])

    inner_pair = min(combinations(train, 2), key=lambda p: slope(p, False))
    outer_pair = max(combinations(held, 2), key=lambda p: slope(p, True))
    ceiling, floor = slope(inner_pair, False), slope(outer_pair, True)
    # Cross-check the *transplanted* stellar normalization with the separate
    # published integrated mass. Compare only radii actually in both apertures.
    outer_radii = [r for r, *_ in rows if r > 30]
    stellar_shell = 3 * (gas[outer_radii[-1]] - gas[outer_radii[0]])
    published_upper = outer.run()["outer_stars_upper_msun"]
    forecasts = []
    for r, x, lo, hi in held:
        band = [A + g * x for A, g in points]
        forecasts.append(
            dict(
                radius_kpc=float(r),
                orbit_msun=[float(lo), float(hi)],
                forecast_msun=[float(min(band)), float(max(band))] if points else None,
                overlap=bool(points and min(band) <= hi and max(band) >= lo),
            )
        )
    conflict = None
    if not vertices(rows):
        # A convex planar intersection is empty on at most three enrolled
        # strips (including the fixed nonnegative quadrant in each strip).
        for size in (2, 3):
            conflict = next(
                (
                    [float(e[0]) for e in part]
                    for part in combinations(rows, size)
                    if not vertices(part)
                ),
                None,
            )
            if conflict:
                break
        assert conflict
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                HERE / "DRAGONFLY_PROFILE.md",
                HERE / "dragonfly_gas.json",
                HERE / "source_profile.json",
            )
        },
        train_kpc=[float(e[0]) for e in train],
        shell_gain_comparison=dict(
            inner_kpc=[float(row[0]) for row in inner_pair],
            inner_upper=float(ceiling),
            outer_kpc=[float(row[0]) for row in outer_pair],
            outer_lower=float(floor),
            required_ratio=float(floor / ceiling),
        ),
        stellar_normalization_check=dict(
            radii_kpc=[float(outer_radii[0]), float(outer_radii[-1])],
            transplanted_stars_msun=int(stellar_shell),
            published_outer_stars_plus_two_errors_msun=published_upper,
            ratio=float(stellar_shell) / published_upper,
        ),
        interior_mass_and_gain=[[float(A), float(g)] for A, g in points],
        fit_gain_range=[float(f(g for _, g in points)) for f in (min, max)]
        if points
        else None,
        held=forecasts,
        all_rings_share_one_fit=bool(vertices(rows)),
        conflict_kpc=conflict,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument(
        "--image-cache",
        type=Path,
        help="Reacquire the enclosed radial gas masses from FITS",
    )
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    source_path = HERE / "dragonfly_gas.json"
    if args.image_cache:
        source = acquire(args.image_cache)
        if args.write:
            source_path.write_text(json.dumps(source, indent=2, allow_nan=False) + "\n")
        else:
            assert json.loads(source_path.read_bytes()) == source
    result = run(json.loads(source_path.read_bytes()))
    path = HERE / "dragonfly_profile.json"
    if args.write:
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        assert json.loads(path.read_bytes()) == result
    print(
        result["fit_gain_range"],
        sum(p["overlap"] for p in result["held"]),
        len(result["held"]),
        result["all_rings_share_one_fit"],
        result["conflict_kpc"],
    )
