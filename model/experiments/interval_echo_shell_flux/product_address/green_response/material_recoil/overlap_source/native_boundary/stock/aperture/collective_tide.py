"""Measure the collective exterior tide of disjoint native overlap sources."""

import argparse
from fractions import Fraction as Q
from hashlib import sha256
import heapq
import json

import numpy as np
from scipy.stats import beta

from . import collective_geometry as geom
from . import collective_source as source

tidal = source.tidal
HERE = tidal.HERE


def interval_acquisition(s, costs, pair, rng):
    """Independent end-to-end counters on the shortest mathematical route.
    Retain the paid contacts, the number of completed bits and the erasures.
    """
    N = 2**34
    chord, truth, paths = geom.measure(s, costs, pair)
    adj = [
        [(int(k), Q(1), -1) for k in s.blank[i].indices]
        for i in range(s.blank.shape[0])
    ]
    for j, ((a, b), c) in enumerate(zip(s.contacts, costs, strict=True)):
        adj[a].append((b, c, j))
        adj[b].append((a, c, j))

    def routes(start):
        d = {start: Q(0)}
        parent = {}
        todo = [(Q(0), start)]
        while todo:
            x, u = heapq.heappop(todo)
            if x != d[u]:
                continue
            for v, c, j in adj[u]:
                y = x + c
                if v not in d or y < d[v]:
                    d[v] = y
                    parent[v] = (u, j)
                    heapq.heappush(todo, (y, v))
        return d, parent

    trees = {u: routes(u) for u in pair}
    routes_out = []
    alpha = 0.01 / (4 * (len(paths) + 1) * max(1, len(costs)))
    lo = []
    hi = []
    for (u, v), exact in zip([tuple(pair), *paths], [chord, *truth], strict=True):
        dist, parent = trees[u]
        assert dist[v] == exact
        walk = v
        blank = 0
        paid = []
        lower = upper = 0.0
        while walk != u:
            walk, j = parent[walk]
            if j < 0:
                blank += 1
            else:
                p = float(1 / costs[j])
                fail = int(rng.negative_binomial(N, p))
                lower += 1 / beta.ppf(1 - alpha, N, fail) if fail else 1.0
                upper += 1 / beta.ppf(alpha, N, fail + 1)
                paid.append([j, fail])
        lo.append(lower + blank)
        hi.append(upper + blank)
        routes_out.append([int(u), int(v), blank, paid])
    # K=(sum cross² - degree)/L²-degree. The numerator may be signed.
    degree = len(paths)
    raw = [sum(x * x for x in v[1:]) - degree for v in (lo, hi)]
    ext = [z / d**2 - degree for z in raw for d in (lo[0], hi[0])]
    predicted = (sum(x * x for x in truth) - degree) / chord**2 - degree
    assert min(ext) - 1e-12 <= float(predicted) <= max(ext) + 1e-12
    return dict(
        exact=str(predicted),
        lower=min(ext),
        upper=max(ext),
        routes=routes_out,
        completed_bits=N,
    )


def counts(v):
    return dict(
        inward=int(np.sum(v < 0)),
        outward=int(np.sum(v > 0)),
        silent=int(np.sum(v == 0)),
    )


def run_panel(board, layout, external_only):
    m = len(layout)
    s = geom.survey(board, layout, external_only)
    original = ("original",) * m
    repaired = ("source",) * m

    def prices(arms, extra=0):
        return [1 / p for p in source.chances(board.n, arms, extra)]

    base = prices(original)
    all_src = prices(repaired)
    one = prices(original, m)
    ref, src, ordinary = (geom.measure(s, x) for x in (base, all_src, one))
    # A double adds zero singles. The partial heads have two or more singles
    # when present, preserving their entire reachable-label set.
    diff = src.curvature - ref.curvature
    ord_diff = ordinary.curvature - ref.curvature
    common = (ref.length == src.length) & (ref.length == ordinary.length)
    radial = ref.equal_radial & src.equal_radial & ordinary.equal_radial
    joint = np.full(len(diff), Q(0), object)
    locked_joint = joint.copy()
    all_quiet = np.ones(len(diff), bool)
    locked_quiet = all_quiet.copy()
    for idx in range(m):
        tags = list(original)
        tags[idx] = "source"
        live = geom.measure(s, prices(tags)) if m > 1 else src
        locked = (
            geom.measure(s, [all_src[j] if j == idx else base[j] for j in range(m)])
            if m > 1
            else src
        )
        for read, cum, quiet in (
            (live, joint, all_quiet),
            (locked, locked_joint, locked_quiet),
        ):
            delta = read.curvature - ref.curvature
            cum += delta
            quiet &= delta == 0
            common &= read.length == ref.length
            radial &= read.equal_radial
    connected = diff - joint
    locked_connected = diff - locked_joint
    U, V = s.pairs.T
    clear = np.minimum(s.clearance[U], s.clearance[V])
    control = common & radial

    def panel(mask):
        n = int(mask.sum())
        if not n:
            return dict(chords=0)

        def avg(v):
            return str(sum(v[mask], Q(0)) / n)

        return dict(
            chords=n,
            repair=counts(diff[mask]),
            ordinary_singles=counts(ord_diff[mask]),
            ordinary_doubles=counts(np.zeros(n, int)),
            departure_from_addition=counts(connected[mask]),
            locked_link_departure_from_addition=counts(locked_connected[mask]),
            exclusively_collective=int(np.sum(mask & all_quiet & (diff < 0))),
            exclusively_collective_locked=int(np.sum(mask & locked_quiet & (diff < 0))),
            average_change=avg(diff),
            average_ordinary=avg(ord_diff),
            average_single_repairs=avg(joint),
        )

    # Shell and clearance statistics use the same fully matched cohort.
    out = dict(
        paths=[list(p) for p, _ in layout],
        home_to_gate_routes=source.preparation_routes(board, layout),
        seams=[seams for _, seams in layout],
        transparent_center=s.center,
        nearest_source_ranges=geom.nearest_sources(s, base),
        mean_renewal_costs=[[str(c) for c in a] for a in (base, all_src, one)],
        all=panel(np.ones(len(diff), bool)),
        common_initial_chord=panel(common),
        matched_radii_and_chord=panel(control),
        shells={
            str(r): panel(control & (s.radii[U] == r)) for r in np.unique(s.radii[U])
        },
        clearance={str(r): panel(control & (clear == r)) for r in np.unique(clear)},
    )
    cand = np.flatnonzero(control & (diff < 0))
    if len(cand):
        # Prefer a purely collective witness, then maximum clearance,
        # then the strongest departure from individual addition.
        j = max(
            cand,
            key=lambda i: (
                bool(all_quiet[i]),
                bool(locked_quiet[i]),
                int(clear[i]),
                abs(connected[i]),
                -int(i),
            ),
        )
        out["witness"] = dict(
            pair=s.pairs[j].tolist(),
            clearance=int(clear[j]),
            radius=int(s.radii[U[j]]),
            single_sum=str(joint[j]),
            total=str(diff[j]),
            mixed=str(connected[j]),
            locked_single_sum=str(locked_joint[j]),
            ordinary_change=str(ord_diff[j]),
            initial_chord=str(ref.length[j]),
        )
    return out


def run():
    source.check_gate()
    banks = {}
    cases = {}
    for level in (0, 1, 2):
        board = tidal.lab.energy.pair.stock.observations.federation(level)
        lim = 4 if level == 0 else 8
        for mode in ("lexicographic", "dispersed"):
            chosen = source.packing(board, mode, lim)
            for m in (1, 2, 4, 8):
                if m > lim:
                    continue
                key = f"{level}/{mode}/{m}"
                if m > len(chosen):
                    cases[key] = dict(enrolled=False, available=len(chosen))
                    continue
                layout = chosen[:m]
                out = run_panel(board, layout, level == 2)
                cases[key] = dict(enrolled=True, **out)
                if mode == "lexicographic":
                    for tag in ("original", "source"):
                        banks[f"{level}/{m}/{tag}"] = source.pooled_draw(
                            board.n,
                            [tag] * m,
                            0,
                            0,
                            6000 + level * 100 + m * 3 + int(tag == "source"),
                        )
                # Acquire a held route bank for a pair of sources on every scale.
                if mode == "dispersed" and m == 2 and "witness" in out:
                    s = geom.survey(board, layout, level == 2)
                    pair = out["witness"]["pair"]
                    rng = np.random.default_rng(430626 + level)
                    arms = [
                        interval_acquisition(s, list(map(Q, prices)), pair, rng)
                        for prices in out["mean_renewal_costs"][:2]
                    ]
                    out = cases[key]["witness"]
                    out["route_banks"] = arms
                    out["joint99_change_interval"] = [
                        arms[1]["lower"] - arms[0]["upper"],
                        arms[1]["upper"] - arms[0]["lower"],
                    ]
    old = json.loads((HERE / "exterior_tide.json").read_bytes())
    small = cases["0/lexicographic/1"]
    assert small["mean_renewal_costs"] == [
        [old["chart_costs"][k]] for k in ("original", "source", "single")
    ]
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (HERE / "COLLECTIVE_TIDE.md", HERE / "exterior_tide.json")
        },
        local_moments={
            k: list(map(str, stats)) for k, stats in source.records()[-1].items()
        },
        gate_acquisitions=banks,
        carriers=cases,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    result = run()
    dst = HERE / "collective_tide.json"
    raw = json.dumps(result, indent=2, allow_nan=False) + "\n"
    if args.write:
        dst.write_text(raw)
    else:
        assert dst.read_text() == raw
    for name, row in result["carriers"].items():
        if row["enrolled"]:
            print(
                name,
                row["matched_radii_and_chord"]["chords"],
                row["matched_radii_and_chord"].get("exclusively_collective", 0),
            )
