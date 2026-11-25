"""Compile the sharp crest-copy and shared comparison instruments on actual seams."""

from collections import Counter
from functools import cache
from itertools import product
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import expm

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from relay import (
    CELL,
    apply,
    federation,
    observations,
    sharp_state,
    topology,
    transaction,
    transpose,
)


@cache
def instructions(g, h):
    """Both instruments restore borrowed loads; the stream receiver then copies the measured bit."""
    route = transaction(g, h)
    copy = np.empty((2, 2, 2), int)
    signed = {}
    law = observations.native_law
    echo = law.echo()
    agreement = (np.arange(law.DIM) % len(law.SRC)) == law.SRC.index((1, 1))
    for x, y, live in product(range(2), repeat=3):
        left, right = sharp_state(x, 1), sharp_state(y, 1)
        original_left, original_right = left.copy(), right.copy()
        # The stream is offered even after the first crest has stopped the heat.
        entry = route["sender"][0]
        ret = route["receiver"]
        apply(left, entry)
        apply(right, ret)
        if live:
            left[g], right[h] = right[h], left[g]
        read = int(right[h] == 2)
        if live:
            left[g], right[h] = right[h], left[g]
        apply(left, entry[::-1])
        apply(right, ret[::-1])
        assert left == original_left and right == original_right
        if read != y:
            apply(right, transpose(0, 1))
        assert right == sharp_state(read, 1)
        copy[live, x, y] = read

        # Move the complementary left rail to g and write ONE external seam.
        # Each local owner can unwrite, read its end, and restore the same write.
        left, right = original_left.copy(), original_right.copy()
        labels = list(range(12))
        apply(labels, entry)
        complement = transpose(labels.index(0), g)
        apply(left, entry)
        apply(right, ret)
        if live:
            left[g], right[h] = right[h], left[g]
        apply(left, complement)
        apply(labels, complement)
        ref = min(CELL[g])
        left_read = transpose(labels.index(11), ref)
        raw = (left[g], right[h])
        col = law.SRC.index(raw)
        written = echo[:, col]
        accepted = np.isclose(
            np.sum(abs(written[agreement]) ** 2), 1, atol=1e-13, rtol=0
        )
        readings = []
        if accepted:
            for end in (0, 1):
                restored = echo.conj().T @ written
                np.testing.assert_allclose(
                    restored, np.eye(law.DIM)[:, col], atol=2e-13, rtol=0
                )
                # Sharp single-end echo against an owner-local vacancy.
                val = law.SRC[int(np.argmax(abs(restored)))][end]
                if end == 0:
                    before = left.copy()
                    apply(left, left_read)
                    assert (left[g], left[ref]) == (val, 0)
                    apply(left, left_read[::-1])
                    assert left == before
                else:
                    assert (right[h], right[route["reference"]]) == (val, 0)
                single_col = law.SRC.index((val, 0))
                local_write = echo[:, single_col]
                detected = np.sum(abs(local_write[agreement]) ** 2)
                bright = bool(np.isclose(detected, 1, atol=1e-13, rtol=0))
                branch = local_write * (agreement if bright else ~agreement)
                np.testing.assert_allclose(
                    echo.conj().T @ branch,
                    np.eye(law.DIM)[:, single_col],
                    atol=2e-13,
                    rtol=0,
                )
                readings.append(int(bright if end == 1 else not bright))
                np.testing.assert_allclose(echo @ restored, written, atol=2e-13, rtol=0)
        apply(left, complement[::-1])
        if live:
            left[g], right[h] = right[h], left[g]
        apply(left, entry[::-1])
        apply(right, ret[::-1])
        assert left == original_left and right == original_right
        assert readings == ([x, x] if live or not x else [])
        signed[live, x] = readings
    return copy, signed


def clock_gates():
    law = observations.native_law
    exchange, contrast = law.local_term_parts()
    h = exchange + observations.KAPPA * contrast
    idx = [law.SRC.index(word) for word in ((2, 0), (0, 2))]
    half = expm(-0.5j * np.pi * h)[:, idx]
    np.testing.assert_allclose(np.delete(half, idx, axis=0), 0, atol=2e-13)
    u = half[idx]
    flip = expm(-1j * np.pi * h)[np.ix_(idx, idx)]
    parked = expm(-1j * np.pi * observations.KAPPA * contrast)[np.ix_(idx, idx)]
    np.testing.assert_allclose(abs(u) ** 2, 0.5, atol=2e-13)
    for gate, expected in ((u @ u, [1, 0]), (flip, [1, 0]), (parked, [0, 1])):
        np.testing.assert_allclose(abs(gate) ** 2, np.eye(2)[:, expected], atol=2e-13)
    return u, tuple(np.argmax(abs(flip), axis=0)), tuple(np.argmax(abs(parked), axis=0))


def screen(level):
    board, groups, dist = topology(level)
    seams = {}
    for a, b in board.external:
        u, v = a // 12, b // 12
        seams[u, v] = a % 12, b % 12
        seams[v, u] = b % 12, a % 12
    # Audit the angular enrollment against actual intercell ports at every nesting.
    for fine in range(1, level + 1):
        fine_edges = federation(fine).external
        cuts = Counter(
            tuple(sorted((a // 48, b // 48)))
            for a, b in fine_edges
            if a // 48 != b // 48
        )
        coarse = {
            tuple(sorted((a // 12, b // 12))) for a, b in federation(fine - 1).external
        }
        assert set(cuts) == coarse and set(cuts.values()) == {2}
        for parent in range(federation(fine - 1).cells):
            inside = {
                (a // 12 % 4, b // 12 % 4)
                for a, b in fine_edges
                if a // 48 == b // 48 == parent
            }
            assert inside == {(0, 3), (1, 3), (2, 3)}
    neighbors = [[] for _ in range(board.cells)]
    for u, v in seams:
        neighbors[u].append(v)
    radii = dist[0]
    parent = [0] + [
        min(u for u in neighbors[v] if radii[u] + 1 == radii[v])
        for v in range(1, board.cells)
    ]
    # Evaluate a complete local port truth table at every gateway used, then
    # batch equal instructions on each depth without ever forwarding a fresh write.
    levels = []
    for radius in range(int(max(radii)), 0, -1):
        children = np.flatnonzero(radii == radius)
        parents = np.array([parent[v] for v in children])
        tables = np.stack(
            [
                instructions(*seams[u, v])[0]
                for u, v in zip(parents, children, strict=True)
            ]
        )
        levels.append((parents, children, tables))
    return board, levels, radii, seams, neighbors, parent


def stream(levels, n, horizon, launch, cut=False):
    states = np.zeros((horizon + 1, n), dtype=np.int8)
    states[0, 0] = int(launch == 0)
    for t in range(1, horizon + 1):
        state = states[t - 1].copy()
        for parents, children, tables in levels:
            state[children] = tables[
                np.arange(len(children)), int(not cut), state[parents], state[children]
            ]
        state[0] = int(t - launch in (0, 2, 4, 6))
        states[t] = state
    return states
