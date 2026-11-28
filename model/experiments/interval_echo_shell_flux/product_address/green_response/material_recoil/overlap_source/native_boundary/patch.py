"""Count written material boundaries and weigh the complete prepared board."""

import argparse
from fractions import Fraction as F
import json
from pathlib import Path

import numpy as np

from model.executable_model import CELL_SEAMS, federation
from model.executable_model.refinement import load_refinement
from .source import HERE, encode, run_case


def board(level):
    if level < 3:
        car = federation(level)
        return car.cells, np.asarray(car.external, int)
    data = (
        Path(__file__).resolve().parents[8]
        / "code"
        / "cosmology"
        / "record_abundance"
        / "FEDERATION_LEVEL3.npz"
    )
    car = load_refinement(data, level=3)
    return car.cells, car.external


def cases(level, external, cells):
    face = np.arange(cells) // 4**level
    faces = face[external // 12]
    near = sorted(
        set(int(y) for a, b in faces for x, y in ((a, b), (b, a)) if x == 0 and y != 0)
    )
    far = min(set(range(1, 20)) - set(near))
    return {
        "one face": face == 0,
        "adjacent faces": np.isin(face, [0, near[0]]),
        "separate faces": np.isin(face, [0, far]),
        "whole": face >= 0,
        "half by owner index": np.arange(cells) < cells // 2,
        "single owner": np.arange(cells) == 0,
    }


def weigh(occupied, external, local):
    flat = np.repeat(occupied.astype(int), 12)
    a, b = external.T
    assert len(set(external.flatten().tolist())) == len(external.flatten())
    cut = flat[a] != flat[b]
    count = int(cut.sum())
    # Twice each one-port occupation probability after the sharp write.
    two = 2 * flat
    two[a] = two[b] = flat[a] + flat[b]
    degree = np.full(len(flat), 5)
    degree[external.ravel()] += 1
    term_degree = F(8 * int(((6 - degree) * two).sum()), 2)
    internal = (
        np.array(CELL_SEAMS)[None, :, :]
        + (12 * np.arange(len(occupied)))[:, None, None]
    )
    left, right = internal.reshape(-1, 2).T
    mismatch_twice = two[left] + two[right] - two[left] * two[right]
    internal_after = F(17 * int(mismatch_twice.sum()), 4)
    external_before = local["ordinary_incident"] * count
    before = 8 * int(((6 - degree) * flat).sum()) + external_before
    retained = term_degree + internal_after + local["sharp"]["retained_energy"] * count
    erased = term_degree + internal_after + local["sharp"]["independent_energy"] * count
    assert erased - retained == local["sharp"]["gain"] * count
    n = int(occupied.sum())
    assert before == 72 * n + F(17, 2) * count and count <= 3 * n
    return dict(
        occupied_owners=n,
        active_written_seams=count,
        initial_whole_energy=before,
        retained_whole_energy=retained,
        independent_whole_energy=erased,
        degree_energy=term_degree,
        other_seam_energy_after_write=internal_after,
        loss_of_pairing_gain=erased - retained,
        gain_relative_to_initial=(erased - retained) / before if before else F(0),
    )


def thresholds(external, level, cells):
    # A deterministic three-layer input, independent of the collision or recoil.
    heights = 1 + np.arange(cells) // 4**level % 4
    owner = external // 12
    left, right = heights[owner].T
    shared = sum(int(((left >= t) != (right >= t)).sum()) for t in range(1, 5))
    assert shared == abs(left - right).sum()
    independent = sum(
        int(((left >= t) != (right >= u)).sum())
        for t in range(1, 5)
        for u in range(1, 5)
    )
    assert independent == (left * (4 - right) + right * (4 - left)).sum()
    return dict(
        shared_expected_cut=F(shared, 4),
        independent_expected_cut=F(independent, 16),
        expected_occupied=F(int(heights.sum()), 4),
    )


def run():
    local = run_case(8)[1]
    panels = []
    for level in range(4):
        size, edges = board(level)
        shapes = {
            name: weigh(occupied, edges, local)
            for name, occupied in cases(level, edges, size).items()
        }
        assert shapes["one face"]["active_written_seams"] == 3 * 2**level
        assert shapes["separate faces"]["active_written_seams"] == 6 * 2**level
        assert shapes["adjacent faces"]["active_written_seams"] == 4 * 2**level
        panels.append(
            dict(
                level=level,
                owners=size,
                patches=shapes,
                layers=thresholds(edges, level, size),
            )
        )
    return panels


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--write", action="store_true")
    args = ap.parse_args()
    result = json.loads(json.dumps(run(), default=encode))
    path = HERE / "patch.json"
    if args.write:
        path.write_text(json.dumps(result, indent=2) + "\n")
    else:
        assert json.loads(path.read_text()) == result
