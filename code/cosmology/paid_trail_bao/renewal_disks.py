"""Ordinary-calibrated galaxy potential differences from the native renewal receiver."""

import argparse
from hashlib import sha256
import json
from math import comb
from pathlib import Path

import numpy as np

from . import clock_width as width
from . import closed_inventory as inventory
from . import passage_images
from . import record_profiles as rp

HERE = inventory.HERE


def integral_at(r, v, x):
    # Integral from the pixel's log radius through the final observed radius.
    tail = width.inward_integrals(r, v)
    j = np.minimum(np.searchsorted(r, x, side="right") - 1, len(r) - 2)
    assert min(j) >= 0
    u = np.log(x / r[j])
    slope = np.diff(v) / np.diff(np.log(r))
    return tail[j] - v[j] * u - slope[j] * u * u / 2


def acquire(cache):
    inputs, *_ = rp.galaxy.transfer.budget.load_inputs()
    peaks = json.loads((HERE / "record_profile_maps.json").read_bytes())[
        "peak_msun_pc2"
    ]
    result = {}
    for name, tracer, img, _, _ in rp.observed_cells(cache, inputs):
        tab = inputs[name]
        r = np.array(tab["radius_kpc"])
        need = max(0, int(np.ceil(r[-1] / rp.galaxy.WIDTH)) + 1 - len(img) // 2)
        image = np.pad(img, need)
        axis = (np.arange(len(image)) + 0.5 - len(image) // 2) * rp.galaxy.WIDTH
        radii = np.hypot(axis[:, None], axis[None, :])
        shell = np.searchsorted(r, radii)
        valid = (radii >= r[0]) & (radii <= r[-1])
        available = np.unique(shell[valid])
        inner = int(available[np.argmin(abs(r[available] - r[-1] / 2))])
        outer = int(available[-1])
        assert inner < outer
        where = [valid & (shell == j) for j in (inner, outer)]
        v, e, b = [
            np.array(tab[k]) for k in ("observed_kms", "error_kms", "baryons_kms")
        ]
        targets = {}
        for k, terms in dict(
            ordinary=b * b,
            observed=v * v - b * b,
            low=np.maximum(0, v - 2 * e) ** 2 - b * b,
            high=(v + 2 * e) ** 2 - b * b,
        ).items():
            values = [integral_at(r, terms, radii[mask]).mean() for mask in where]
            targets[k] = float(values[0] - values[1])
        shellstats = {}
        for cap in rp.CAPACITY:
            moments = []
            for mask in where:
                p = image[mask] / (cap * peaks[tracer])
                assert np.all((0 <= p) & (p <= 1))
                moments.append(
                    dict(
                        mean_filling=float(p.mean()),
                        inventory=[
                            float(np.mean(comb(4, k) * p**k * (1 - p) ** (4 - k)))
                            for k in range(5)
                        ],
                    )
                )
            shellstats[str(cap)] = moments
        if name not in result:
            result[name] = dict(
                annular_endpoints_kpc=[r[inner], r[outer]],
                targets=targets,
                address_counts=[int(mask.sum()) for mask in where],
                sources={},
            )
        else:
            np.testing.assert_allclose(
                list(targets.values()), list(result[name]["targets"].values())
            )
        result[name]["sources"][tracer] = shellstats
    return result


def compare(rows):
    def gain(names):
        a = np.array([rows[n]["bare"] for n in names])
        b = np.array([rows[n]["targets"]["ordinary"] for n in names])
        return max(0.0, float(a @ b / (a @ a)))

    names = sorted(rows)
    tests = {}
    intervals = []
    for name in names:
        row = rows[name]
        obs = row["targets"]
        other = gain([n for n in names if n != name])
        own = obs["ordinary"] / row["bare"] if row["bare"] > 0 else None
        possible = passage_images.intervals(
            np.array([row["extra"]]), np.array([obs["low"]]), np.array([obs["high"]])
        )[0]
        intervals.append((max(0.0, possible[0]), possible[1]))
        held = other * row["extra"]
        own_pred = None if own is None else own * row["extra"]
        tests[name] = dict(
            extra_clock_log=row["extra"],
            bare_clock_log=row["bare"],
            ordinary_gain_own=own,
            ordinary_gain_other=other,
            held_ordinary_kms2=other * row["bare"],
            held_extra_kms2=held,
            own_extra_kms2=own_pred,
            own_within=own_pred is not None and obs["low"] <= own_pred <= obs["high"],
            held_within=obs["low"] <= held <= obs["high"],
            positive_dark_sign=row["extra"] > 0,
        )
    best, k = passage_images.best_gain(intervals)
    return dict(
        galaxies=tests,
        held_within=sum(int(row["held_within"]) for row in tests.values()),
        own_within=sum(int(row["own_within"]) for row in tests.values()),
        positive_directions=sum(
            int(row["positive_dark_sign"]) for row in tests.values()
        ),
        best_any_dark_gain=int(best),
        favorable_dark_gain=float(k),
    )


def run(maps):
    from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture import (
        renewal_clock as lab,
    )

    _, classes = lab.lab.energy.pair.source.classes()
    counts = np.array([np.bincount(w[:, 1], minlength=3) for w in classes])
    distributions = counts / counts.sum(1)[:, None]
    clock = json.loads((lab.HERE / "renewal_clock.json").read_bytes())
    result = {}
    for t, reader in clock["readers"].items():
        period = 1 / np.mean(reader["exact"]["flags"], axis=1)
        typical = distributions @ period
        delta = (448 / len(classes[2])) / 14 * ((period[0] + period[2]) / 2 - period[1])
        panels = {}
        for name, tab in maps.items():
            for tracer, caps in tab["sources"].items():
                for cap, shells in caps.items():
                    p = np.array([s["mean_filling"] for s in shells])
                    stock = np.array([s["inventory"] for s in shells])
                    np.testing.assert_allclose(stock.sum(1), 1.0)
                    bare = period[0] + p * (period[2] - period[0])
                    added = np.log1p(delta * stock[:, 2] / (stock @ typical))
                    panels.setdefault(tracer + "/" + cap, {})[name] = dict(
                        targets=tab["targets"],
                        bare=float(np.log(bare[0] / bare[1])),
                        extra=float(added[0] - added[1]),
                    )
        result[t] = {key: compare(rows) for key, rows in panels.items()}
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                HERE / "RENEWAL_DISKS.md",
                HERE / "renewal_disk_maps.json",
                HERE / "source_profile.json",
                lab.HERE / "renewal_clock.json",
            )
        },
        native_middle_load_counts=counts.tolist(),
        readers=result,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--image-cache", type=Path)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    src = HERE / "renewal_disk_maps.json"
    if args.image_cache:
        obs = acquire(args.image_cache)
        if args.write:
            src.write_text(json.dumps(obs, indent=2, allow_nan=False) + "\n")
        else:
            rp.galaxy.check(json.loads(src.read_bytes()), obs)
    result = run(json.loads(src.read_bytes()))
    dst = HERE / "renewal_disks.json"
    if args.write:
        dst.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        rp.galaxy.check(json.loads(dst.read_bytes()), result)
    for t, v in result["readers"].items():
        print(
            t,
            {
                k: max(r[k] for r in v.values())
                for k in (
                    "own_within",
                    "held_within",
                    "positive_directions",
                    "best_any_dark_gain",
                )
            },
        )
