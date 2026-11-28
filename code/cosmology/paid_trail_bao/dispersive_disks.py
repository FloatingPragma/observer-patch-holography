"""Held galaxy rotation curves from distributed local dispersive clock responses."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import linprog

from . import thermal_mass as maps
from . import dispersive_work as work

HERE = work.HERE
LENGTH = 10.25


def disk_derivative(radii, bins, order):
    """R d(-profile)/dR per input annulus and per positive profile drop."""
    points, weights = leggauss(order)
    edges = np.r_[0.0, radii]
    shell = np.arange(bins + 1) * LENGTH
    assert 2 * radii[-1] <= shell[-1]
    result = np.zeros((len(radii), len(radii), bins))
    for idx, R in enumerate(radii):
        # Radial breaks at every place an angular shell starts or stops
        # intersecting the integration circle. R is itself a breakpoint.
        extra = np.r_[R + shell, R - shell, shell - R]
        breaks = np.unique(np.r_[edges, extra[(extra > 0) & (extra < edges[-1])]])
        half = np.diff(breaks) / 2
        rho = (breaks[:-1, None] + half[:, None] * (1 + points)).ravel()
        ann = np.searchsorted(radii, rho)
        radial_weight = (
            2
            * rho
            * np.outer(half, weights).ravel()
            / (edges[ann + 1] ** 2 - edges[ann] ** 2)
        )
        # Compute cuts from sin(theta/2) for accuracy close to R=rho.
        sine2 = (shell**2 - (R - rho[:, None]) ** 2) / (4 * R * rho[:, None])
        cuts = 2 * np.arcsin(np.sqrt(np.clip(sine2, 0, 1)))
        angle_weight = np.diff(cuts, axis=1) / 2
        angle = cuts[:, :-1, None] + angle_weight[:, :, None] * (1 + points)
        half_sine = np.sin(angle / 2) ** 2
        dist = np.sqrt(
            (R - rho[:, None, None]) ** 2 + 4 * R * rho[:, None, None] * half_sine
        )
        proj = (R - rho[:, None, None] + 2 * rho[:, None, None] * half_sine) / dist
        integrand = (proj @ weights) * angle_weight / np.pi
        np.add.at(result[idx], ann, (R / LENGTH) * radial_weight[:, None] * integrand)
    return result


def all_fields(inputs, kernels, records, order):
    bins = len(next(iter(kernels.values()))) - 1
    geometry = {
        n: disk_derivative(np.array(g["radius_kpc"]), bins, order)
        for n, g in inputs.items()
    }
    fields = {}
    for tag, capacities in records.items():
        for cap, menus in capacities.items():
            for menu, galaxies in menus.items():
                tensor = {
                    n: np.einsum("rjb,j->rb", geometry[n], np.diff(enclosed, prepend=0))
                    for n, enclosed in galaxies.items()
                }
                for k, profile in kernels.items():
                    direct = {n: mat @ (-np.diff(profile)) for n, mat in tensor.items()}
                    for sign in (1, -1):
                        fields[f"{tag}/{cap}/{menu}/{k}/{sign}"] = {
                            n: sign * f for n, f in direct.items()
                        }
    return fields


def compare(inputs, fields):
    names = sorted(inputs)
    observations = {n: np.array(inputs[n]["observed_kms"]) for n in names}
    bary = {n: np.array(inputs[n]["baryons_kms"]) ** 2 for n in names}
    errors = {n: np.array(inputs[n]["error_kms"]) for n in names}
    target = {n: observations[n] ** 2 - bary[n] for n in names}
    tests = {}
    chosen = {}
    for k, x in sorted(fields.items()):
        folds = {}
        total = full = undefined = 0
        cross = np.array([np.mean(x[n] * target[n]) for n in names])
        squares = np.array([np.mean(x[n] ** 2) for n in names])
        for i, held in enumerate(names):
            train = np.arange(len(names)) != i
            den = squares[train].sum()
            g = max(0, float(cross[train].sum() / den)) if den else 0.0
            # This score is exactly the training loss minus the common
            # baryons-only training loss. Held velocities never enter selection.
            loss = float(g * (g * den - 2 * cross[train].sum()) / train.sum())
            v2 = bary[held] + g * x[held]
            good = v2 >= 0
            valid = good & (
                abs(np.sqrt(np.maximum(v2, 0)) - observations[held]) <= 2 * errors[held]
            )
            hits = int(valid.sum())
            missing = int((~good).sum())
            folds[held] = [g, hits, missing]
            total += hits
            undefined += missing
            full += int(valid.all())
            if held not in chosen or (loss, k) < (
                chosen[held]["training_loss_change"],
                chosen[held]["setting"],
            ):
                chosen[held] = dict(
                    setting=k,
                    training_loss_change=loss,
                    gain=g,
                    within_two=hits,
                    undefined=missing,
                    model_squared_kms=v2.tolist(),
                    observed_kms=observations[held].tolist(),
                    error_kms=errors[held].tolist(),
                    radius_kpc=inputs[held]["radius_kpc"],
                )
        tests[k] = dict(
            within_two=total, undefined=undefined, full_curves=full, folds=folds
        )
    return tests, chosen


def mixture_certificates(inputs, first, final):
    """Test arbitrary nonnegative weights on the six kernels for each image law."""
    names = sorted(inputs)
    rings = [
        (name, i) for name in names for i in range(len(inputs[name]["radius_kpc"]))
    ]
    limits = [maps.rp.galaxy.transfer.budget.bounds(inputs[n]) for n in names]
    b = np.r_[
        np.concatenate([v[1] for v in limits]), -np.concatenate([v[0] for v in limits])
    ]
    groups = {}
    for key in sorted(final):
        tag, cap, menu, _, _, sign = key.split("/")
        groups.setdefault("/".join((tag, cap, menu, sign)), []).append(key)
    certs = {}
    for family, keys in groups.items():

        def constraints(fields):
            f = np.column_stack(
                [np.concatenate([fields[k][n] for n in names]) for k in keys]
            )
            return np.vstack((f, -f))

        matrices = [constraints(x) for x in (first, final)]
        union = np.hstack(matrices)
        scaled = union / np.abs(union).max(0)
        result = linprog(
            b / max(abs(b)),
            A_ub=-scaled.T,
            b_ub=np.full(scaled.shape[1], -1e-5),
            A_eq=np.ones((1, len(b))),
            b_eq=[1],
            bounds=(0, None),
            method="highs",
        )
        if not result.success:
            certs[family] = dict(certified=False)
            continue
        y = np.rint(result.x * 1e10).astype(np.int64)
        support = np.flatnonzero(y > 0)
        if not len(support):
            certs[family] = dict(certified=False)
            continue

        # These inequalities are checked in exact arithmetic on both finite
        # quadrature tables. The numerical optimizer only proposes the weights.
        def check(column):
            return sum(
                int(y[i]) * Fraction.from_float(float(column[i])) for i in support
            )

        margins = [check(column) for column in union.T]
        weighted_limit = check(b)
        certified = min(margins) >= 0 and weighted_limit < 0
        certs[family] = dict(
            certified=certified,
            profile_order=keys,
            weighted_limits_per_unit=float(weighted_limit / int(y.sum())),
            minimum_profile_margin_per_unit=float(min(margins) / int(y.sum())),
            nonnegative_weights=[
                [
                    *rings[i % len(rings)],
                    "upper" if i < len(rings) else "negated_lower",
                    int(y[i]),
                ]
                for i in support
            ],
        )
    return certs


def run():
    sourcefile = maps.LAB / "dispersive_clock.json"
    same = json.loads(sourcefile.read_bytes())
    mapfile = HERE / "thermal_mass_maps.json"
    mapped = json.loads(mapfile.read_bytes())
    inputs, *_ = maps.rp.galaxy.transfer.budget.load_inputs()
    profiles = {
        f"{species}/{tag.replace('/', '_')}": np.mean(entry["mixed_shell"], axis=1)
        for species, panels in same["results"].items()
        for tag, entry in panels["inactive"].items()
    }
    first = all_fields(inputs, profiles, mapped, 12)
    fields = all_fields(inputs, profiles, mapped, 24)
    initial, prior = compare(inputs, first)
    final, chosen = compare(inputs, fields)
    change = max(
        abs(initial[k]["within_two"] - r["within_two"]) for k, r in final.items()
    )
    rms_shift = max(
        abs(
            np.sqrt(np.maximum(0, r["model_squared_kms"]))
            - np.sqrt(np.maximum(0, prior[n]["model_squared_kms"]))
        ).max()
        for n, r in chosen.items()
    )
    best = max(r["within_two"] for r in final.values())
    baseline = sum(
        int(
            np.sum(
                abs(np.array(v["observed_kms"]) - v["baryons_kms"])
                <= 2 * np.array(v["error_kms"])
            )
        )
        for v in inputs.values()
    )
    certificates = mixture_certificates(inputs, first, fields)
    return dict(
        mixture_families=len(certificates),
        excluded_mixture_families=sum(
            int(v["certified"]) for v in certificates.values()
        ),
        mixture_certificates=certificates,
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (HERE / "DISPERSIVE_DISKS.md", mapfile, sourcefile)
        },
        length_kpc=LENGTH,
        galaxies=len(inputs),
        rings=sum(len(v["radius_kpc"]) for v in inputs.values()),
        baseline_within_two=baseline,
        settings=len(final),
        best_held=best,
        retrospectively_best=[k for k, r in final.items() if r["within_two"] == best],
        nested_within_two=sum(r["within_two"] for r in chosen.values()),
        nested_full_curves=sum(
            r["within_two"] == len(r["radius_kpc"]) for r in chosen.values()
        ),
        nested_undefined=sum(r["undefined"] for r in chosen.values()),
        quadrature=dict(
            orders=[12, 24],
            max_setting_hit_change=change,
            changed_selections=[
                n for n in chosen if chosen[n]["setting"] != prior[n]["setting"]
            ],
            max_nested_speed_change_kms=rms_shift,
        ),
        nested=chosen,
        fold_columns=["training_gain", "within_two", "undefined"],
        tests=final,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    output = run()
    path = HERE / "dispersive_disks.json"
    if args.write:
        path.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n")
    else:
        work.lab.flight.lab.contact.timing.lab.lab.cold.apparatus.check(
            json.loads(path.read_bytes()), output
        )
    print(
        output["nested_within_two"],
        output["nested_full_curves"],
        output["best_held"],
        output["quadrature"],
    )
