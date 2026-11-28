"""Test a necessary integer detour on all port pairs of the populated carriers."""

import argparse
from dataclasses import replace
from hashlib import sha256
import json
import math

import numpy as np
from scipy.sparse.csgraph import dijkstra

from . import collective_geometry as geom
from . import collective_source as source

HERE = source.tidal.HERE


def possible_detours(m, full, single, other):
    # A single repair has its own local price and changes the global pool
    # everywhere else. Consider the cheapest single repair of this route.
    return {
        k: list(
            range(math.floor(k * full) + 1, math.floor(single + (k - 1) * other) + 1)
        )
        for k in range(1, m + 1)
        if math.floor(k * full) < math.floor(single + (k - 1) * other)
    }


def address_pairs(s, ends):
    pairs = np.array(sorted(ends), int).reshape(-1, 2)
    starts = np.unique(np.r_[s.center, pairs.ravel()])
    index = np.full(len(s.radii), -1, int)
    index[starts] = np.arange(len(starts))
    return replace(s, pairs=pairs, starts=starts, index=index)


def run_case(board, layout):
    m = len(layout)
    n = board.n
    all_orig = ["original"] * m
    solo = ["source", *all_orig[1:]]
    original = 1 / source.chances(n, all_orig)[0]
    full = 1 / source.chances(n, ["source"] * m)[0]
    one_cost = [1 / p for p in source.chances(n, solo)]
    single = one_cost[0]
    other = one_cost[-1] if m > 1 else original
    assert full <= single < original
    if m > 1:
        assert single < other < original and full < other
    gaps = possible_detours(m, full, single, other)
    answer = dict(
        ports=n,
        unordered_ranges=n * (n - 1) // 2,
        source_paths=[p for p, _ in layout],
        signed_costs=list(map(str, (original, full, single, other))),
        possible_integer_detours=gaps,
        exclusively_collective_ranges=0,
        exclusively_collective_exterior_chords=0,
    )
    if not gaps:
        return answer
    s = geom.survey(board, layout, True)
    limit = 2 * int(s.radii.max())
    scale, prices, grid = geom.encode_lengths([full] * m, limit)
    price = int(prices[0])
    # Counts are unique in this finite range. Otherwise every tied
    # representation would have to be retained.
    assert full.denominator > m
    table = {
        b * scale + k * price: (b, k) for b in range(limit + 1) for k in range(m + 1)
    }
    assert len(table) == (limit + 1) * (m + 1) and set(table) == set(grid)
    net = (s.blank * scale).tolil()
    for a, b in s.contacts:
        net[a, b] = net[b, a] = price
    net = net.tocsr()
    potential = []
    crossings = np.zeros(m + 1, int)
    mod = np.zeros(scale, np.int16)
    assert len({k * price % scale for k in range(m + 1)}) == m + 1
    for k in range(m + 1):
        mod[k * price % scale] = k
    permitted = np.zeros((m + 1, limit + 1), bool)
    for k, vals in gaps.items():
        for val in vals:
            if val <= limit:
                permitted[k, val] = True
    for offset in range(0, n, 96):
        starts = np.arange(offset, min(n, offset + 96))
        d = dijkstra(net, indices=starts, directed=False).astype(np.int64)
        k = mod[d % scale]
        unique = np.arange(n)[None, :] > starts[:, None]
        crossings += np.bincount(k[unique], minlength=m + 1)
        cand = unique & np.isin(k, list(gaps))
        if not cand.any():
            continue
        i, j = np.nonzero(cand)
        bare = dijkstra(s.blank, indices=starts, directed=False, unweighted=True)[
            i, j
        ].astype(int)
        b = (d[i, j] - k[i, j].astype(np.int64) * price) // scale
        diff = bare - b
        retain = (
            (diff >= 0) & (diff <= limit) & permitted[k[i, j], np.clip(diff, 0, limit)]
        )
        for start, end in zip(starts[i[retain]], j[retain], strict=True):
            potential.append((int(start), int(end)))
    assert sum(crossings) == n * (n - 1) // 2
    answer.update(
        shortest_route_paid_counts=crossings.tolist(),
        actual_route_candidates=len(potential),
    )
    if not potential:
        return answer
    cand = address_pairs(s, potential)
    base = geom.measure(cand, [original] * m)
    repaired = geom.measure(cand, [full] * m)
    quiet = base.length > repaired.length
    for j in range(m):
        run = geom.measure(cand, [single if i == j else other for i in range(m)])
        quiet &= run.length == base.length
    exact = cand.pairs[quiet]
    answer["exclusively_collective_ranges"] = len(exact)
    if not len(exact):
        return answer
    answer["range_example"] = exact[0].tolist()
    # An unchanged-release chord may use that changed range in either
    # orientation of either launch star. Every candidate is tested.
    launch = set()
    for a, b in exact:
        for v, w in ((a, b), (b, a)):
            for u in board.neighbours[w]:
                if (
                    u != v
                    and min(s.clearance[v], s.clearance[u]) >= 2
                    and s.radii[u] == s.radii[v]
                ):
                    launch.add(tuple(sorted((int(v), int(u)))))
    cand = address_pairs(s, launch)
    base = geom.measure(cand, [original] * m)
    repaired = geom.measure(cand, [full] * m)
    effect = (
        (repaired.curvature < base.curvature)
        & (repaired.length == base.length)
        & base.equal_radial
        & repaired.equal_radial
    )
    for j in range(m):
        row = geom.measure(cand, [single if i == j else other for i in range(m)])
        effect &= (
            (row.curvature == base.curvature)
            & (row.length == base.length)
            & row.equal_radial
        )
    answer["exclusively_collective_exterior_chords"] = int(effect.sum())
    if effect.any():
        i = int(np.flatnonzero(effect)[0])
        answer["chord_witness"] = dict(
            ports=cand.pairs[i].tolist(),
            change=str(repaired.curvature[i] - base.curvature[i]),
        )
    return answer


def run():
    out = {}
    for level in (0, 1, 2):
        board = source.tidal.lab.energy.pair.stock.observations.federation(level)
        for mode in ("lexicographic", "dispersed"):
            paths = source.packing(board, mode, board.cells // 4)
            for m in sorted({1, 2, 4, 8, 12, 16, 24, len(paths)}):
                if m <= len(paths):
                    out[f"{level}/{mode}/{m}"] = run_case(board, paths[:m])
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (HERE / "COLLECTIVE_SWITCH.md", HERE / "collective_tide.json")
        },
        cases=out,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    answer = run()
    dst = HERE / "collective_switch.json"
    raw = json.dumps(answer, indent=2, allow_nan=False) + "\n"
    if args.write:
        dst.write_text(raw)
    else:
        assert raw == dst.read_text()
    print(
        [
            (
                k,
                v["exclusively_collective_ranges"],
                v["exclusively_collective_exterior_chords"],
            )
            for k, v in answer["cases"].items()
            if v["possible_integer_detours"]
        ]
    )
