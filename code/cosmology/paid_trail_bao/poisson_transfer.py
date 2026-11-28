"""Compare galaxy orbits with the weak Poisson source and the strong localized preparation."""

import argparse
import json
from hashlib import sha256

import numpy as np

from . import free_rest
from . import kernel_budget as budget
from . import newton_transport as newton

HERE = budget.HERE


def scale_test(inputs, fields):
    folds = {}
    errors = []
    hits = 0
    for held, rows in inputs.items():
        train = sorted(set(inputs) - {held})

        def target(n):
            return (
                np.array(inputs[n]["observed_kms"]) ** 2
                - np.array(inputs[n]["baryons_kms"]) ** 2
            )

        coeff = max(
            0,
            sum(np.mean(fields[n] * target(n)) for n in train)
            / sum(np.mean(fields[n] ** 2) for n in train),
        )
        v = np.sqrt(np.array(rows["baryons_kms"]) ** 2 + coeff * fields[held])
        diff = abs(np.array(rows["observed_kms"]) - v)
        col = fields[held]
        assert np.all(col >= 0)
        lo, hi = budget.bounds(rows)
        pos = col > 0
        low = max(0, float(np.max(lo[pos] / col[pos])))
        high = float(np.min(hi[pos] / col[pos]))
        permitted = low <= high and np.all((lo[~pos] <= 0) & (hi[~pos] >= 0))
        inside = diff <= 2 * np.array(rows["error_kms"])
        hits += int(inside.sum())
        errors.extend(diff)
        folds[held] = dict(
            gain=coeff,
            forecast_kms=v.tolist(),
            within_two=int(inside.sum()),
            any_target_scale=bool(permitted),
            target_lower=low,
            target_upper=high,
        )
    return dict(
        galaxies=folds,
        rings=len(errors),
        within_two=hits,
        median_abs_kms=float(np.median(errors)),
        possible_full_curves=sum(x["any_target_scale"] for x in folds.values()),
    )


def run():
    inputs, charts, meta, masses = budget.load_inputs()
    lab = free_rest.lab
    origin = lab.HERE / "poisson_response.json"
    standard = lab.HERE / "mediated_release.json"
    rows = json.loads(origin.read_bytes())["endpoint"]["observations"]["endpoint"][
        "65536"
    ]
    tape = json.loads(standard.read_bytes())
    shape = free_rest.shell_field(rows, tape)
    disk = {n: np.array(row["baryons_kms"]) ** 2 for n, row in inputs.items()}
    sphere = {}
    for n, row in inputs.items():
        r = np.array(row["radius_kpc"])
        enclosed = newton.cdf_matrix(r, np.array([0.0]), charts[n])[:, 0]
        sphere[n] = budget.local.profiles.galaxies.G * 1e9 * masses[n] * enclosed / r
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                origin,
                standard,
                HERE / "POISSON_TRANSFER.md",
                HERE / "source_profile.json",
                budget.local.profiles.maps.HERE / "leroy_table7.dat",
                budget.local.profiles.maps.HERE / "leroy_sample.dat",
                budget.local.profiles.galaxies.sparc.DATA / "table1.dat",
            )
        },
        disk=scale_test(inputs, disk),
        sphere=scale_test(inputs, sphere),
        shell_field=shape.tolist(),
        strong=budget.local.evaluate(
            inputs, charts, meta, (np.arange(len(shape)), shape)
        ),
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    result = run()
    path = HERE / "poisson_transfer.json"
    if args.write:
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), result)
    print(
        {
            k: {a: b for a, b in result[k].items() if a != "galaxies"}
            for k in ("disk", "sphere", "strong")
        }
    )
