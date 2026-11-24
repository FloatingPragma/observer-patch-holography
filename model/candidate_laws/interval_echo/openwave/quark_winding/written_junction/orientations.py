"""Rotated wall backgrounds and their measurable coherences on the written-field cell."""

import argparse
from itertools import permutations
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh

from model.executable_model import carrier, rotations
from model.candidate_laws.interval_echo import candidate, currents
from ..belt_walls import belt
from ...row_13_walls.energy import is_wall
from . import loop_colour


def turn(patch, fields, matter, g):
    new = np.zeros_like(fields)
    loads = np.zeros_like(matter)
    for i, j in enumerate(g):
        loads |= ((matter >> i) & 1) << j
    for (a, b), mask in patch.bit.items():
        target = patch.bit[tuple(sorted((g[a], g[b])))]
        toward_a = (fields & mask) != 0
        new[toward_a if g[a] < g[b] else ~toward_a] |= target
    return new, loads


def block(patch, fields, matter, degrees):
    ordered = np.argsort(fields)
    sorted_fields = fields[ordered]
    row, col, values = [], [], []
    n = (matter[:, None] >> np.arange(12)) & 1
    bg = -n.copy()

    def moves(ids, flip, weight):
        other = fields[ids] ^ flip
        pos = np.searchsorted(sorted_fields, other)
        assert np.array_equal(sorted_fields[pos], other)
        row.append(ids)
        col.append(ordered[pos])
        values.append(np.full(len(ids), weight))

    for (a, b), mask in patch.bit.items():
        head = (fields & mask) != 0
        bg[:, a] += head
        bg[:, b] += ~head
        ids = np.flatnonzero((n[:, a] != n[:, b]) & (n[:, a] == head))
        moves(ids, mask, 8)
    for mask, orientation in patch.rings:
        bits = fields & mask
        moves(
            np.flatnonzero((bits == orientation) | (bits == (mask ^ orientation))),
            mask,
            1,
        )
    assert np.all(bg == bg[0])
    adj = sparse.csr_matrix(
        (np.concatenate(values), (np.concatenate(row), np.concatenate(col))),
        shape=(len(n), len(n)),
    )
    assert (adj != adj.T).nnz == 0
    # Twice the positive Laplacian and the degree term for actual records 1 and 2.
    h = sparse.csgraph.laplacian(adj).tocsr() + sparse.diags(
        6 * int(candidate.KAPPA) * (n @ (6 - degrees)), dtype=np.int64
    )
    return h, bg[0].tolist()


def run():
    board = carrier.federation(0)
    degrees = board.degree[:12]
    idx = np.arange(12)
    rotations3 = [np.array(g) for g in rotations()]
    g = min(
        (
            tuple(g)
            for g in rotations3
            if np.all(g != idx)
            and np.array_equal(g[g[g]], idx)
            and np.array_equal(degrees[g], degrees)
        )
    )
    patch = loop_colour.Patch()
    order, seed = loop_colour.wall_source(patch)
    panels = {}
    wall_words = {
        sum(1 << path[i] for i in range(10) if w & (1 << i))
        for path in belt.belts().values()
        for w in range(1024)
        if is_wall(w)
    }
    for name, sites in [("minus", (0, 4, 7)), ("plus", (0, 1, 4, 7))]:
        start = sum(1 << order[i] for i in sites)
        f, m = loop_colour.component(patch, start, seed)
        original = (f, m)
        copies = []
        for _ in range(3):
            h, b = block(patch, f, m, degrees)
            broken, _ = block(patch, f, m, board.degree[12:24])
            copies.append((h, b, broken))
            f, m = turn(patch, f, m, g)
        assert np.array_equal(f, original[0]) and np.array_equal(m, original[1])
        assert len({tuple(b) for _, b, _ in copies}) == 3
        base = copies[0][0]
        miss = [(mat != base).nnz for mat, _, _ in copies]
        assert not any(miss)
        mask = np.isin(original[1], list(wall_words))
        lines = {}
        for key, which in [("lowest", "SA"), ("highest", "LA")]:
            hs = base.astype(float) * 0.5
            value, vec = eigsh(
                hs,
                k=1,
                which=which,
                v0=np.random.default_rng(121313).normal(size=len(m)),
                tol=1e-12,
            )
            state = vec[:, 0]
            residual = float(np.linalg.norm(hs @ state - value[0] * state))
            assert residual < 1e-8
            lines[key] = dict(
                energy=float(value[0]),
                wall_probability=float(state[mask] @ state[mask]),
                residual=residual,
            )
        panels[name] = dict(
            states_per_orientation=len(m),
            backgrounds=[b for _, b, _ in copies],
            degree_symmetric_mismatches=miss,
            changed_owner_mismatches=[(a != copies[0][2]).nnz for _, _, a in copies],
            lines=lines,
        )
    # A permutation of the three background labels is its own superselection
    # outcome. Compare coherences before and after measuring these labels.
    D = np.zeros((3, 3, 3))
    for a, b, c in permutations(range(3)):
        D[a, b, c] = (-1) ** ((a > b) + (b > c) + (a > c)) / np.sqrt(6)
    rho = np.outer(D.ravel(), D.ravel())
    mixed = np.diag(np.diag(rho))
    cycle = np.roll(np.eye(3), 1, axis=0)
    rotation = np.kron(np.kron(cycle, cycle), cycle)
    return currents._portable(
        dict(
            rotation=list(map(int, g)),
            panels=panels,
            largest_diagonal_difference=float(abs(np.diag(rho - mixed)).max()),
            determinant_rotated_bright=(1 + np.trace(rotation @ rho)) / 2,
            mixture_rotated_bright=(1 + np.trace(rotation @ mixed)) / 2,
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    result = run()
    with args.output.open("x") as f:
        f.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
