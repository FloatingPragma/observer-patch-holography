"""Conditional cut-amplitude profiles from independent images and the native return receiver."""

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

from . import closed_inventory as inventory
from . import passage_images as images
from . import source_profile as profiles

rp = inventory.rp
HERE = inventory.HERE


def convex_test(tab, outer):
    r = np.array(tab["radius_kpc"], float)
    lo, hi = rp.galaxy.transfer.budget.bounds(tab)
    lo, hi = r * lo, r * hi
    if outer:
        keep = 2 * r >= r[-1]
        r, lo, hi = r[keep], lo[keep], hi[keep]
    if len(r) < (3 if outer else 1):
        return dict(included=False, rings=len(r))
    rings = len(r)
    if not outer:
        r, lo, hi = (np.r_[0, x] for x in (r, lo, hi))
    lo = np.maximum(0, lo)
    # A falling prefix of the lower convex hull of all upper endpoints is
    # flattened at their minimum. Every boundary segment is a two-point chord.
    n = len(r)
    gap = 0.0
    witness = {}

    def compare(g, kind, at):
        nonlocal gap, witness
        if g > gap:
            gap = float(g)
            witness = dict(
                kind=kind,
                radii_kpc=[float(r[i]) for i in at],
                mass_gap_msun=gap / profiles.galaxies.G,
            )

    for j in range(n):
        i = j + int(np.argmin(hi[j:]))
        compare(lo[j] - hi[i], "nonnegative" if i == j else "nondecreasing", (j, i))
        if 0 < j < n - 1:
            weight = (r[j] - r[:j, None]) / (r[None, j + 1 :] - r[:j, None])
            chord = (1 - weight) * hi[:j, None] + weight * hi[None, j + 1 :]
            left, right = np.unravel_index(np.argmin(chord), chord.shape)
            compare(lo[j] - chord[left, right], "convex", (left, j, j + 1 + right))
    # Independent feasibility test using adjacent slope inequalities only.
    dr = np.diff(r)
    slope = np.diff(np.eye(n), axis=0) / dr[:, None]
    A = np.vstack((-slope, -np.diff(slope, axis=0)))
    scale = max(1, abs(hi).max(), abs(lo).max())
    opt = linprog(
        np.zeros(n),
        A_ub=A,
        b_ub=np.zeros(len(A)),
        bounds=list(zip(lo / scale, hi / scale, strict=True)),
    )
    if opt.success:
        assert gap <= 1e-8 * scale
    else:
        assert opt.status == 2 and gap > 1e-9 * scale
    return dict(included=True, rings=rings, feasible=bool(opt.success), witness=witness)


def shape_catalog(inputs):
    meta, groups = profiles.galaxies.catalogue()
    kept, _ = profiles.galaxies.select(meta, groups)
    bulgeless = {
        name: dict(
            radius_kpc=[p["rad_kpc"] for p in rows],
            observed_kms=[p["vobs_kms"] for p in rows],
            error_kms=[p["e_vobs_kms"] for p in rows],
            baryons_kms=np.sqrt(profiles.bare_speed2(rows)),
        )
        for name, rows in kept.items()
    }
    result = {}
    for label, tabs in (("bulgeless", bulgeless), ("imaged", inputs)):
        result[label] = {}
        for domain in ("all", "outer"):
            out = {k: convex_test(v, domain == "outer") for k, v in tabs.items()}
            good = [v for v in out.values() if v["included"]]
            result[label][domain] = dict(
                enrolled=len(good),
                failures=sum(not v["feasible"] for v in good),
                galaxies=out,
            )
    return result


def rate(h, n):
    root = np.sqrt(h)
    return root if not n else root * np.tanh(n * np.arctanh(root))


def integrated(stations, rate_at_station, r):
    # Right-continuous cut amplitude jumps at the measured pixel center.
    knots = np.r_[0, stations]
    value = np.r_[0, rate_at_station]
    area = np.r_[0, np.cumsum(value[:-1] * np.diff(knots))]
    i = np.searchsorted(knots, r, side="right") - 1
    return (area[i] + (r - knots[i]) * value[i]) / r


def predict(curves, inputs):
    rows = []
    xx = {}
    xy = {}
    targets = {}
    for name, field in curves.items():
        tab = inputs[name]
        lo, hi = rp.galaxy.transfer.budget.bounds(tab)
        v = np.array(tab["observed_kms"])
        y = v * v - np.array(tab["baryons_kms"]) ** 2
        xx[name] = np.mean(field**2)
        xy[name] = np.mean(field * y)
        targets[name] = (lo, hi)
        for a, b in images.intervals(field, lo, hi):
            if b < 0:
                rows.append([np.inf, -np.inf])
            else:
                rows.append([max(0.0, a), b])
    # Include infinite upper endpoints from vanishing fields explicitly.
    constant = sum(a == 0 and b == np.inf for a, b in rows)
    finite = [(a, b) for a, b in rows if np.isfinite(a) and np.isfinite(b)]
    best, gain = images.best_gain(finite)
    count = complete = held_full = 0
    held = {}
    for name, field in curves.items():
        gain9 = max(
            0.0,
            sum(v for k, v in xy.items() if k != name)
            / sum(v for k, v in xx.items() if k != name),
        )
        lo, hi = targets[name]
        ok = (gain9 * field >= lo) & (gain9 * field <= hi)
        count += int(ok.sum())
        held_full += int(ok.all())
        held[name] = dict(gain=float(gain9), within_two=int(ok.sum()), rings=len(field))
        complete += int(
            np.all((gain * field >= lo - 1e-8) & (gain * field <= hi + 1e-8))
        )
    return dict(
        best_possible_count=best + constant,
        fitted_gain=gain,
        all_disk_least_squares_gain=float(max(0, sum(xy.values()) / sum(xx.values()))),
        complete_at_best_gain=complete,
        held_count=count,
        held_complete=held_full,
        held=held,
    )


def run(cache):
    inputs, *_ = rp.galaxy.transfer.budget.load_inputs()
    source = inventory.LAB / "fluctuation_return.json"
    lab = json.loads(source.read_bytes())
    sizes = [v["half_length"] for v in lab["boards"].values()] + [0]
    peak = json.loads((HERE / "record_profile_maps.json").read_bytes())["peak_msun_pc2"]
    width = rp.galaxy.WIDTH
    max_r = max(tab["radius_kpc"][-1] for tab in inputs.values())
    m = int(np.ceil(max_r / width))
    axis = (np.arange(-m, m) + 0.5) * width
    slots = int(np.count_nonzero(np.hypot(axis[:, None], axis[None, :]) <= max_r))
    predictions = {}
    for name, tracer, img, _, mask in rp.observed_cells(cache, inputs):
        r = np.array(inputs[name]["radius_kpc"])
        axis = (np.arange(len(img)) + 0.5 - len(img) // 2) * width
        stations, idx = np.unique(
            np.hypot(axis[:, None], axis[None, :])[mask], return_inverse=True
        )
        for cap in rp.CAPACITY:
            f = img[mask] / (cap * peak[tracer])
            eligible = 6 * f**2 * (1 - f) ** 2 * (448 / 2110)

            def cumulative(pixel):
                return np.bincount(idx, weights=pixel).cumsum() / slots

            donors = cumulative(f / 4)
            v = cumulative(eligible**2 / 2)
            old, new = (v * y * (1 - y) for y in (3 / 28, 5 / 28))
            direct = cumulative(eligible / 4)
            for n in sizes:
                direct_old, direct_new = (rate(direct * y, n) for y in (3 / 28, 5 / 28))
                a, b = rate(old, n), rate(new, n)
                for label, cut in (
                    ("single", rate(donors, n)),
                    ("middle_original", direct_old),
                    ("middle_repaired", direct_new),
                    ("middle_increment", direct_new - direct_old),
                    ("original", a),
                    ("repaired", b),
                    ("increment", b - a),
                ):
                    assert min(cut) >= 0 and min(np.diff(cut), default=0) >= -1e-13
                    key = f"{tracer}/{cap}/{n or 'infinite'}/{label}"
                    predictions.setdefault(key, {})[name] = integrated(stations, cut, r)
    panels = {key: predict(curves, inputs) for key, curves in predictions.items()}
    files = (
        source,
        HERE / "FLUCTUATION_CUTS.md",
        HERE / "source_profile.json",
        HERE / "record_profile_maps.json",
        HERE / "pair_galaxy_maps.json",
        profiles.galaxies.sparc.DATA / "table1.dat",
        profiles.galaxies.sparc.DATA / "table2.dat",
    )
    return dict(
        sha256={p.name: sha256(p.read_bytes()).hexdigest() for p in files},
        possible_addresses=slots,
        panels=panels,
        convexity=shape_catalog(inputs),
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--image-cache", type=Path, required=True)
    cli.add_argument("--write", action="store_true")
    opt = cli.parse_args()
    result = run(opt.image_cache)
    dst = HERE / "fluctuation_cuts.json"
    if opt.write:
        dst.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        rp.galaxy.check(json.loads(dst.read_bytes()), result)
    print(
        {
            label: {dom: v["failures"] for dom, v in tab.items()}
            for label, tab in result["convexity"].items()
        }
    )
