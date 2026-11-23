"""Cross-check a held contribution to the quantum stress with small-ball area."""

import argparse
import json
from pathlib import Path

import numpy as np

from . import measure
from ..observer_field import local_state


def density_by_axis(tape, mass2):
    _, arms, *_ = local_state(tape["radar"])
    contributions = []
    for length, modes in zip(arms, tape["modes"], strict=True):
        qq, pp = [
            sum(measure.second_moment(row[key]) for row in modes) for key in ("q", "p")
        ]
        contributions.append(
            (pp / arms.prod() ** 2 + (mass2 + 1 / length**2) * qq) / (2 * measure.VOL0)
        )
    return np.array(contributions)


def run(primary, geometry):
    fine = primary["surveys"][-1]
    base = density_by_axis(fine["records"]["ground"][1], primary["recovered_mass2"])
    vacuum_balls = geometry["records"]["ground"]["spheres"]
    result = {}
    for name in ("z", "y", "even_x"):
        tangent = (
            density_by_axis(fine["records"][name][1], primary["recovered_mass2"]) - base
        )
        ordinary = float(tangent[:2].sum())
        reserved = float(tangent[2])
        comparisons = []
        for sphere, ground in zip(
            geometry["records"][name]["spheres"], vacuum_balls, strict=True
        ):
            curvature = sphere["estimated_Guu"] - ground["estimated_Guu"]
            excess = curvature / fine["calibration"]["coupling"] - ordinary
            comparisons.append(
                dict(
                    reference=sphere["reference"],
                    delta_curvature=curvature,
                    inferred_reserved=excess,
                    discrepancy=excess - reserved,
                )
            )
        result[name] = dict(
            ordinary_change=ordinary, reserved_change=reserved, comparisons=comparisons
        )
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(
        *(
            json.loads((measure.HERE / name).read_text())
            for name in ("results.json", "balls.json")
        )
    )
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    for label, row in result.items():
        print(
            label,
            row["ordinary_change"],
            row["reserved_change"],
            row["comparisons"][-1],
        )
