"""Escape witnesses for every mobile background and every three-wall sharp word."""

import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path

import numpy as np
import sympy as sp

from . import backgrounds


def rotate(w, length):
    return ((w << 1) & ((1 << length) - 1)) | (w >> (length - 1))


def neighbours(w, e, length):
    forward = (w >> 1) | ((w & 1) << (length - 1))
    return (w ^ forward) & ~(w ^ e)


def back_mask(w, selected, length):
    return sum(
        1 << j
        for j in range(length)
        if (w ^ ((1 << j) | (1 << ((j + 1) % length)))) in selected
        and ((w >> j) ^ (w >> ((j + 1) % length))) & 1
    )


def unresolved_rank(columns, selected, length):
    parents = defaultdict(set)
    for i, (w, e) in enumerate(columns):
        allowed = neighbours(w, e, length)
        for j in range(length):
            dest = w ^ ((1 << j) | (1 << ((j + 1) % length)))
            if allowed & (1 << j) and dest not in selected:
                parents[(dest, e ^ (1 << j))].add(i)
    alive = set(range(len(columns)))
    while alive:
        rows = [row & alive for row in parents.values()]
        peel = set().union(*(r for r in rows if len(r) == 1))
        if not peel:
            break
        alive -= peel
    if not alive:
        return 0, 0
    # All native allowed hops have the same amplitude. Rank over the rationals
    # resolves any interference not already eliminated by single-parent rows.
    ids = sorted(alive)
    mat = sp.Matrix([[int(i in r) for i in ids] for r in parents.values() if r & alive])
    rank = mat.rank() if mat.rows else 0
    return len(ids), len(ids) - rank


def panel(length, n):
    chart = backgrounds.arcs(length, n)
    selected = set(map(int, np.flatnonzero(np.isfinite(chart))))
    heads = np.arange(1 << length, dtype=np.int64)
    prev = rotate(heads, length)
    remaining = set(selected)
    counts = Counter()
    unresolved = defaultdict(list)
    unique_hist = Counter()
    while remaining:
        w = min(remaining)
        orbit = [w]
        v = rotate(w, length)
        while v != w:
            orbit.append(v)
            v = rotate(v, length)
        remaining.difference_update(orbit)
        # g=n-e+previous. Integer values -1 and 2 force a static boundary.
        mobile = ((~w & heads & ~prev) | (w & ~heads & prev)) == 0
        f = neighbours(w, heads, length)
        unique = np.zeros(len(heads), dtype=np.int8)
        exits = np.zeros_like(unique)
        for j in range(length):
            z = 1 << j
            dst = w ^ (z | (1 << ((j + 1) % length)))
            if dst in selected:
                continue
            enabled = (f & z) != 0
            back = neighbours(dst, heads ^ z, length) & back_mask(dst, selected, length)
            unique += enabled & (back == z)
            exits += enabled
        m = len(orbit)
        counts["sharp_joint_words"] += m * len(heads)
        counts["pin_free_words"] += m * int(mobile.sum())
        counts["frozen_words"] += m * int((f == 0).sum())
        counts["pin_free_frozen"] += m * int(((f == 0) & mobile).sum())
        counts["pin_free_without_exit"] += m * int(((exits == 0) & mobile).sum())
        for k, v in enumerate(np.bincount(unique[mobile])):
            unique_hist[k] += m * int(v)
        for e in heads[mobile & (unique == 0)]:
            v = int(e)
            for material in orbit:
                g = material ^ v ^ rotate(v, length)
                unresolved[g].append((material, v))
                v = rotate(v, length)
    algebra = Counter()
    for columns in unresolved.values():
        left, null = unresolved_rank(columns, selected, length)
        algebra["unpeeled_columns"] += left
        algebra["boundary_nullity"] += null
        algebra["sectors_with_kernel"] += int(null > 0)
    return dict(
        length=length,
        loads=n,
        material_patterns=len(selected),
        **counts,
        single_parent_histogram={
            str(k): v for k, v in sorted(unique_hist.items()) if v
        },
        backgrounds_requiring_further_elimination=len(unresolved),
        unpeeled_columns=algebra["unpeeled_columns"],
        boundary_nullity=algebra["boundary_nullity"],
        sectors_with_kernel=algebra["sectors_with_kernel"],
    )


def run():
    return [
        panel(length, length // 3 + add) for length in (12, 15, 18) for add in (0, 1)
    ]


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", required=True, type=Path)
    a = p.parse_args()
    r = run()
    with a.output.open("x") as f:
        f.write(json.dumps(r, indent=2, allow_nan=False) + "\n")
