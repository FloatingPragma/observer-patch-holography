"""Test the uniform quantum shell using the native exchange amplitudes."""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

from model.executable_model.refinement import load_refinement
from .measure import HERE, Routes


def modes(routes, radius):
    home = routes.range_ball(0, radius + 1)
    source = sum(d == radius - 1 for d in home.values())
    row = []
    for degree in (5, 6):
        bins = Counter(
            sum(home.get(j) == radius - 1 for j in routes.neighbours(y))
            for y, d in home.items()
            if d == radius and len(routes.neighbours(y)) == degree
        )
        if not bins:
            continue
        count = sum(bins.values())
        linear = sum(k * v for k, v in bins.items())
        squared = sum(k * k * v for k, v in bins.items())
        total = squared / (4 * source)
        matched = linear**2 / (4 * source * count)
        row.append(
            dict(
                degree=degree,
                output_ports=count,
                inward_multiplicities=dict(sorted(bins.items())),
                transferred_norm_coefficient=total,
                uniform_mode_coefficient=matched,
                fraction_orthogonal=(count * squared - linear * linear)
                / (count * squared),
            )
        )
    return dict(launch_ports=source, blocks=row)


def run(directory):
    panel = json.loads((HERE / "results.json").read_text())
    records = []
    for p in panel["panels"]:
        level = p["level"]
        car = load_refinement(directory / f"FEDERATION_LEVEL{level}.npz", level=level)
        routes = Routes(car)
        records.append(
            dict(
                level=level,
                cells=p["cells"],
                radius=p["radius"],
                **modes(routes, p["radius"]),
            )
        )
        routes.neighbours.cache_clear()
    return dict(
        protocol_sha256=hashlib.sha256(
            (HERE / "SHELL_PROTOCOL.md").read_bytes()
        ).hexdigest(),
        records=records,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    data = run(args.directory)
    with args.output.open("x") as f:
        json.dump(data, f, indent=2, allow_nan=False)
        f.write("\n")
    for r in data["records"]:
        print(
            r["cells"], [(s["degree"], s["fraction_orthogonal"]) for s in r["blocks"]]
        )
