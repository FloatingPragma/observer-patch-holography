"""Cluster and scaling checks with independent preparation bits on 20 ancestors."""

import argparse
import json
from pathlib import Path

import numpy as np

from model.executable_model.refinement import load_refinement
from . import measure
from .replay import verify_replay


def support(directory, level):
    board = load_refinement(directory / f"FEDERATION_LEVEL{level}.npz", level=level)
    reader = measure.ProtectedRanges(board)
    d = reader.survey(0)
    r = 6 * 2 ** ((level - 1) // 2)
    group = np.arange(board.cells) // (4**level)
    remote = int(np.argmax(d))
    near_choices = np.flatnonzero((group == 0) & (d >= 3 * r) & (d <= 4 * r))
    origins = {"far": remote}
    if len(near_choices):
        origins["near"] = int(near_choices[0])
    paired = []
    for name, origin in origins.items():
        partner = reader.survey(origin)
        for multiple in (1, 2) if name == "far" else (1,):
            lo, hi = r * multiple / 2, r * multiple
            a = (d > lo) & (d <= hi)
            b = (partner > lo) & (partner <= hi)
            if np.any(a & b):
                raise ValueError("Correlated cluster test requires disjoint annuli")
            admissions = np.column_stack(
                [
                    np.bincount(group[member], minlength=20)
                    for member in (a, b, ~(a | b))
                ]
            )
            paired.append(
                dict(
                    partner=name,
                    origin=origin,
                    distance=int(d[origin]),
                    scale=multiple,
                    admissions=admissions.tolist(),
                )
            )
    return dict(
        level=level,
        cells=board.cells,
        calibration_range=int(d.max()),
        radius=r,
        settings=paired,
    )


def tagged_slots(table):
    """Ordered different-cell draws in each of A×A, A×B, ..., outside×outside."""
    table = np.array(table, dtype=float)
    n = table.sum()
    total = measure.slot_membership(table.sum(axis=0), False)
    same = (table.T @ table - np.diag(table.sum(axis=0))) / (n * (n - 1))
    return same, total - same


def counts(laws, same, different, arm, rng):
    chances = np.stack([same, different])
    bins = rng.multinomial(measure.SHOTS, chances.ravel()).reshape(2, 3, 3)
    record = {}
    for tag, a, b in np.ndindex(2, 3, 3):
        law = laws[arm if tag == 0 else 1]
        if a == 2 and b == 2:
            readings = [int(bins[tag, a, b])]
        else:
            distribution = (
                law.ravel() if max(a, b) < 2 else law.sum(axis=0 if a == 2 else 1)
            )
            readings = rng.multinomial(bins[tag, a, b], distribution).tolist()
        record[f"{tag}{a}{b}"] = readings
    return record


def experiment(geometries):
    laws, score, ref = measure.source_reader()
    panels = []
    for j, geometry in enumerate(geometries):
        arms = []
        for k, setting in enumerate(geometry["settings"]):
            same, different = tagged_slots(setting["admissions"])
            reads = [
                counts(
                    laws,
                    same,
                    different,
                    arm,
                    np.random.default_rng([measure.SEED, 2, j, k, arm]),
                )
                for arm in range(3)
            ]
            radius = (
                geometry["radius"] * setting["scale"] / geometry["calibration_range"]
            )
            results = {}
            for label, slots in (("home", ((0, 0),)), ("cross", ((0, 1), (1, 0)))):
                p = sum(same[x, y] for x, y in slots)

                def read(index):
                    return (
                        sum(
                            np.array(reads[index][f"0{x}{y}"]) @ score.ravel()
                            for x, y in slots
                        )
                        / measure.SHOTS
                    )

                expected = float(p * ref["score_when_both_shared_change"])
                m = float(np.max(abs(score)))
                logtail = np.log(2 * 384 / 0.01)
                width = (
                    2 * np.sqrt(2 * p * m**2 * logtail / measure.SHOTS)
                    + 8 * m * logtail / (3 * measure.SHOTS)
                    if p
                    else 0.0
                )
                results[label] = dict(
                    shared=float(read(2) - read(1)),
                    reversed=float(read(0) - read(1)),
                    expected=expected,
                    normalized=expected / radius**4,
                    simultaneous_radius=float(width),
                )
            arms.append(dict(readings=reads, contrasts=results))
        panels.append(dict(geometry=geometry, arms=arms))
    return dict(calibration=ref, panels=panels)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    path = measure.HERE / "ancestor.json"
    if args.directory:
        geometry = [support(args.directory, lvl) for lvl in (1, 3, 5, 7)]
    else:
        geometry = [row["geometry"] for row in json.loads(path.read_text())["panels"]]
    result = experiment(geometry)
    if args.output:
        with args.output.open("x") as out:
            json.dump(result, out, separators=(",", ":"), allow_nan=False)
            out.write("\n")
    else:
        verify_replay(result, json.loads(path.read_text()), name="readings")
    for panel in result["panels"]:
        for setting, arm in zip(
            panel["geometry"]["settings"], panel["arms"], strict=True
        ):
            h = arm["contrasts"]["home"]
            x = arm["contrasts"]["cross"]
            print(
                panel["geometry"]["cells"],
                setting["partner"],
                setting["scale"],
                h["normalized"],
                x["normalized"],
                x["shared"],
                x["simultaneous_radius"],
            )


if __name__ == "__main__":
    main()
