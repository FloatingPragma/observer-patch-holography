"""Eliminate an arbitrary probe-charge-odd contribution to each material reading."""

import argparse
from hashlib import sha256
import json

from . import ordinary_service_distribution as ordinary

mat = ordinary.mat
HERE = ordinary.HERE


def even(panel):
    out = {k: v for k, v in panel.items() if k.startswith("clock/")}
    for name in mat.METERS:
        lo, hi = (mat.iv.mpf(panel[f"{name}/{sign}"]) for sign in (-1, 1))
        out[name] = mat.lab.six.bounds((lo + hi) / 2)
    return out


def run():
    mat.iv.dps = 35
    distribution = json.loads(
        (HERE / "ordinary_service_distribution.json").read_bytes()
    )
    base = json.loads((HERE / "density_material.json").read_bytes())
    out = {}
    for window, recorded in distribution["windows"].items():
        costs = mat.iv.mpf(base["windows"][window]["ordinary_energy"])
        calibration = even(recorded["held_ordinary"])
        columns = {k: even(v) for k, v in recorded["ordinary_columns"].items()}
        targets = {k: even(v) for k, v in recorded["target"].items()}
        groups = {
            "impulse": ["shift"],
            "material": list(mat.METERS),
            "joint": [*mat.METERS, *recorded["groups"]["clocks"]],
        }
        prices = {}
        results = {}
        for instrument, target in targets.items():
            prices[instrument] = {}
            for name in mat.METERS:
                denominator = mat.iv.mpf(calibration[name])
                prices[instrument][name] = (
                    None
                    if denominator.a <= 0 <= denominator.b
                    else mat.lab.six.bounds(
                        costs * mat.iv.mpf(target[name]) / denominator
                    )
                )
            results[instrument] = {
                k: ordinary.separate(columns, target, v) for k, v in groups.items()
            }
        out[window] = {
            "uniform_calibration": calibration,
            "ordinary_energy": mat.lab.six.bounds(costs),
            "uniform_energy_equivalents": prices,
            "comparisons": results,
            "groups": groups,
            "ordinary_columns": columns,
            "target": targets,
        }
    return {
        "windows": out,
        "inputs": {
            n: sha256((HERE / n).read_bytes()).hexdigest()
            for n in (
                "CHARGE_EVEN_MATERIAL.md",
                "ordinary_service_distribution.json",
                "density_material.json",
            )
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    record = run()
    path = HERE / "charge_even_material.json"
    if args.write:
        path.write_text(json.dumps(record, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), record)
    for window, r in record["windows"].items():
        print(
            window,
            r["uniform_energy_equivalents"]["direct"],
            {k: v["status"] for k, v in r["comparisons"]["direct"].items()},
        )
