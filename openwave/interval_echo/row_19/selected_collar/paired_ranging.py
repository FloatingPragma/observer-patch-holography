"""Measure the cost of transferring a pair record and survey its interior shells."""

import argparse
from hashlib import sha256
from itertools import combinations, combinations_with_replacement, product
import json
from math import acos, pi, sqrt

import networkx as nx
import numpy as np
from scipy.linalg import expm
from scipy.sparse.csgraph import shortest_path
from scipy.spatial import cKDTree

from . import midpoint_bulk as bulk
from . import word_routes as routes

HERE = bulk.HERE
N = 2**54
SEED = 28122861
SHELL_EDGES = [(i, j) for i, j in combinations(range(6), 2) if i // 2 != j // 2]
FACES = list(product((0, 1), (2, 3), (4, 5)))


def gram_signature(hop):
    # Three times each dot product in Q(sqrt(5)), in graph-hop order.
    rational = np.array([3, 0, 1, -1, 0, -3])[hop]
    radical = np.array([0, 1, 0, 0, -1, 0])[hop]
    return np.column_stack((rational, radical))


def midpoint_classes(points, hop):
    sig = gram_signature(hop)
    groups = {}
    for a, b in combinations_with_replacement(range(len(points)), 2):
        groups.setdefault(tuple(sig[a] + sig[b]), []).append((a, b))
    pools = list(groups.values())
    means = np.array([(points[g[0][0]] + points[g[0][1]]) / 2 for g in pools])
    for p, group in zip(means, pools, strict=True):
        for i, j in group:
            np.testing.assert_allclose((points[i] + points[j]) / 2, p, atol=1e-13)
    return pools, means


def queries(pools, edges):
    rows, first = [], []
    for i, j in edges:
        first.append(len(rows))
        rows.extend((a, b, c, d) for a, b in pools[i] for c, d in pools[j])
    return np.array(rows).T, np.array(first)


def deliver(dist, query):
    (a, b, c, d), first = query
    costs = np.minimum(
        np.maximum(dist[a, c], dist[b, d]), np.maximum(dist[a, d], dist[b, c])
    )
    return np.stack([costs[first], np.minimum.reduceat(costs, first)])


def angle_test(graphs, positions):
    checks = []
    for graph, points in zip(graphs, positions, strict=True):
        tree = cKDTree(points)
        hop = nx.floyd_warshall_numpy(graph, nodelist=range(len(points))).astype(int)
        rows = []
        for i in range(3):
            n, t, u = np.roll(np.eye(3), i, axis=0)
            for radius in (0.5, 0.75, sqrt(15) / 4):
                for name, (x, y) in (
                    ("radial", (0 * n, radius * n)),
                    ("transverse", (radius * n, radius * t)),
                ):
                    w = sqrt(1 - x @ x) * u
                    z = sqrt(1 - y @ y) * u
                    true = np.array([x + w, x - w, y + z, y - z])
                    _, ids = tree.query(true)
                    a, b, c, d = ids

                    def theta(v):
                        return acos(float(np.clip(v, -1, 1)))

                    ideal = theta(x @ y + sqrt(1 - x @ x) * sqrt(1 - y @ y))
                    dots = points @ points.T if len(points) <= 20 else None

                    def dist(p, q):
                        return theta(
                            dots[p, q] if dots is not None else points[p] @ points[q]
                        )

                    finite = min(
                        max(dist(a, c), dist(b, d)), max(dist(a, d), dist(b, c))
                    )
                    angular_error = np.max(
                        np.arccos(np.clip((points[ids] * true).sum(axis=1), -1, 1))
                    )
                    clock = min(max(hop[a, c], hop[b, d]), max(hop[a, d], hop[b, c]))
                    assert abs(finite - ideal) <= 2 * angular_error + 1e-12
                    rows.append(
                        {
                            "plane": i,
                            "kind": name,
                            "depth": radius,
                            "endpoints": ids.tolist(),
                            "hemisphere": ideal,
                            "finite_angular_range": finite,
                            "surface_quantization": angular_error,
                            "protected_time": 7 * int(clock),
                        }
                    )
        ratios = [r["protected_time"] / r["hemisphere"] for r in rows]
        checks.append(
            {
                "stations": len(points),
                "rows": rows,
                "single_speed_minimax_relative_error": (max(ratios) - min(ratios))
                / (max(ratios) + min(ratios)),
            }
        )
    return checks


def triangle_test(dist, pools):
    q = queries(pools, list(combinations(range(len(pools)), 2)))
    times = deliver(dist, q)
    final = {}
    for tag, vec in zip(("fixed", "best"), times, strict=True):
        n = len(pools)
        d = np.zeros((n, n), int)
        lo, hi = np.triu_indices(n, 1)
        d[lo, hi] = d[hi, lo] = vec
        worst, count, witness = 0, 0, None
        for j in range(n):
            gap = d - d[:, [j]] - d[[j], :]
            count += int(np.count_nonzero(gap[lo, hi] > 0))
            i, k = np.unravel_index(np.argmax(gap), gap.shape)
            if gap[i, k] > worst:
                worst = int(gap[i, k])
                witness = {
                    "classes": [int(i), j, int(k)],
                    "direct": int(d[i, k]),
                    "legs": [int(d[i, j]), int(d[j, k])],
                }

                def choice(u, v):
                    return min(
                        (max(dist[a, c], dist[b, h]), [a, b], [c, h])
                        for a, b in pools[u]
                        for c, h in [z for pair in pools[v] for z in (pair, pair[::-1])]
                    )

                _, p, m = choice(i, j)
                _, m2, t = choice(j, k)
                witness["leg_encodings"] = [p, m, m2, t]
                a, b = m
                c, h = m2
                witness["switch"] = int(
                    min(max(dist[a, c], dist[b, h]), max(dist[a, h], dist[b, c]))
                )
        final[tag] = {
            "violated_triangles": count,
            "maximum_deficit": worst,
            "witness": witness,
        }
    return final


def shells(means):
    radii = np.linalg.norm(means, axis=1)
    ranks = [
        z
        for z in sorted(set(radii.round(8)))
        if np.isclose(radii, z, atol=1e-8).sum() == 30
    ]
    assert len(ranks) == 3
    middle = np.flatnonzero(np.isclose(radii, ranks[1]))

    def direction(i):
        return means[i] / radii[i]

    axes = next(
        t
        for t in combinations(middle, 3)
        if all(abs(direction(a) @ direction(b)) < 1e-12 for a, b in combinations(t, 2))
    )
    bearings = np.array([sign * direction(a) for a in axes for sign in (1, -1)])
    plans = []
    for r in ranks:
        ids = cKDTree(means).query(r * bearings)[1]
        assert np.max(abs(means[ids] - r * bearings)) < 1e-8
        center = int(np.argmin(radii))
        plans.append(
            [(center, int(j)) for j in ids]
            + [(int(ids[i]), int(ids[j])) for i, j in SHELL_EDGES]
        )
    return ranks, bearings, plans


def chamber(edge):
    lengths = np.zeros((7, 7), dtype=edge.dtype)
    for (i, j), length in zip(
        [(6, i) for i in range(6)] + SHELL_EDGES, edge, strict=True
    ):
        lengths[i, j] = lengths[j, i] = length
    area, vol = 0.0, 0.0
    determinants = []
    for face in FACES:
        base = lengths[6, list(face)] ** 2
        d2 = lengths[np.ix_(face, face)] ** 2
        gram = (base[:, None] + base[None, :] - d2) / 2
        det = np.linalg.det(gram)
        determinants.append(det.real)
        if np.linalg.eigvalsh(gram.real)[0] <= 0:
            return None
        a, b, c = [lengths[i, j] for i, j in combinations(face, 2)]
        area += np.sqrt((a + b + c) * (-a + b + c) * (a - b + c) * (a + b - c)) / 4
        vol += np.sqrt(det) / 6
    return area, vol, min(determinants)


def fixed_volume(edge, change):
    ref = edge.mean(axis=0)
    h = 1e-12
    base = chamber(ref)
    if base is None:
        return {"valid": False}

    def derivative(vector):
        calc = chamber(ref.astype(complex) + 1j * h * vector)
        return np.imag(calc[:2]) / h

    source = derivative(change.mean(axis=0))
    radial = derivative(edge[1] - edge[0])
    fraction = -source[1] / radial[1]
    return {
        "valid": True,
        "area": base[0],
        "volume": base[1],
        "minimum_chamber_determinant": base[2],
        "source_dA_dV": source.tolist(),
        "radius_coin_dA_dV": radial.tolist(),
        "coin_slope_at_fixed_volume": fraction,
        "fixed_volume_area_slope": source[0] + fraction * radial[0],
    }


def bounded_chambers(edge, slope):
    iv = routes.iv
    side = {
        tuple(sorted(k)): (x, y)
        for k, x, y in zip(
            [(6, i) for i in range(6)] + SHELL_EDGES, edge, slope, strict=True
        )
    }
    darea = dvol = iv.mpf(0)
    for face in FACES:
        pairs = list(combinations(face, 2))
        radial = [side[i, 6] for i in face]
        transverse = [side[k] for k in pairs]
        length = [v for v, _ in transverse]
        deriv = [v for _, v in transverse]
        factors = [sum(length), *[sum(length) - 2 * h for h in length]]
        slopes = [sum(deriv), *[sum(deriv) - 2 * h for h in deriv]]
        if min(float(v.a) for v in factors) <= 0:
            return None
        area = iv.sqrt(factors[0] * factors[1] * factors[2] * factors[3]) / 4
        darea += area * sum(d / t for d, t in zip(slopes, factors, strict=True)) / 2
        x = [v**2 for v, _ in radial]
        dx = [2 * v * d for v, d in radial]
        y = []
        dy = []
        for (i, j), (v, d) in zip(combinations(range(3), 2), transverse, strict=True):
            y.append((x[i] + x[j] - v**2) / 2)
            dy.append((dx[i] + dx[j] - 2 * v * d) / 2)
        # Opposite off-diagonal entries for vertices 0,1,2.
        y = y[::-1]
        dy = dy[::-1]
        det = (
            x[0] * x[1] * x[2]
            + 2 * y[0] * y[1] * y[2]
            - sum(x[i] * y[i] ** 2 for i in range(3))
        )
        if float(det.a) <= 0 or float((x[0] * x[1] - y[2] ** 2).a) <= 0:
            return None
        ddet = 0
        for i, j, k in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
            ddet += dx[i] * (x[j] * x[k] - y[i] ** 2) + 2 * dy[i] * (
                y[j] * y[k] - x[i] * y[i]
            )
        dvol += ddet / (12 * iv.sqrt(det))
    return darea, dvol


def certify_chambers(tapes, hist, ref, sign):
    iv = routes.iv
    rows = []
    for arm in (ref, 3, 4):
        pair = []
        for shell in (0, 2):
            edge = []
            for i in range(18 * shell, 18 * (shell + 1)):
                sums = sum(int(n) * k for k, n in enumerate(tapes[arm, i]))
                support = np.flatnonzero(hist[arm, i])
                delta = int(support[-1] - support[0]) * iv.sqrt(
                    iv.ln(2000000) / (2 * N)
                )
                edge.append(iv.mpf(sums) / N + iv.mpf([-1, 1]) * delta)
            pair.append(edge)
        rows.append(pair)
    near, far = rows[0]
    edge = [(a + b) / 2 for a, b in zip(near, far, strict=True)]
    radial = [b - a for a, b in zip(near, far, strict=True)]
    slope = [
        sign * (a + b - c - d) / 144
        for a, b, c, d in zip(*rows[1], *rows[2], strict=True)
    ]
    nat = bounded_chambers(edge, slope)
    comp = bounded_chambers(edge, radial)
    if nat is None or comp is None or 0 in comp[1]:
        return None
    coin = -nat[1] / comp[1]
    return {
        "coin_slope": routes.bounds(coin),
        "area_slope": routes.bounds(nat[0] + coin * comp[0]),
    }


def local_swap():
    law = routes.volume.candidate.native_law
    ex, co = law.local_term_parts()
    h = ex + 8 * co
    blank = np.eye(law.DIM)[:, : law.N_SRC]
    target = blank[:, [law.SRC.index((b, a)) for a, b in law.SRC]]
    np.testing.assert_allclose(expm(-1j * pi * h) @ blank, target, atol=3e-13, rtol=0)
    return {
        "checked_blank_load_pairs": [list(k) for k in law.SRC],
        "ordinary_home_mode": json.loads(
            (HERE / "settled_midpoint_force.json").read_bytes()
        )["ordinary_mode"]["home_amplitudes_before_normalization"],
    }


def standards(board):
    src = routes.volume.radar_source_and_hand(board)[-1]
    laws, _ = routes.shared_variance.source_words(src)
    weights = np.rint(2592 * laws).astype(np.int64)
    delta = weights[2] - weights[1]
    turned = weights.reshape(3, 9, 3, 3).transpose(0, 1, 3, 2).reshape(3, 9, 9)
    np.testing.assert_array_equal(turned[2] - turned[1], -delta)
    symmetric = turned + turned.transpose(0, 2, 1)
    assert not np.any(symmetric % 2)

    weights = np.concatenate(
        [
            weights,
            [72 * np.maximum(delta, 0), 72 * np.maximum(-delta, 0)],
            symmetric // 2,
        ]
    )
    assert np.all(weights.sum(axis=(1, 2)) == 2592)
    return weights, turned


def run():
    routes.iv.dps = 35
    board = bulk.federation(0)
    graphs, positions = bulk.geometry._atlas(
        [board, bulk.federation(1), bulk.federation(2)]
    )
    points = positions[0]
    hops = nx.floyd_warshall_numpy(graphs[0], nodelist=range(board.cells)).astype(int)
    pools, means = midpoint_classes(points, hops)
    radii, bearings, plans = shells(means)
    query = queries(pools, [p for plan in plans for p in plan])
    code, _ = routes.codebook()
    cal, _ = routes.calibrate(code)
    make = routes.graph(board, cal)

    def arrival(adj):
        return shortest_path(adj, indices=np.arange(0, board.n, 12), unweighted=True)[
            :, ::12
        ].astype(int)

    blank = arrival(make(-1, -1, 0, 0)[0])
    tri = {
        name: triangle_test(d, pools)
        for name, d in (("absence", blank), ("protected", 7 * hops))
    }
    weights, _ = standards(board)
    # A histogram counts every distinct source placement and both word tags.
    tags = [*combinations(range(board.cells), 2)]
    bins = np.arange(40)
    hist = np.zeros((len(weights), 2, len(query[1]), len(bins)), np.int64)
    for i, j in tags:
        vals = np.stack(
            [
                deliver(arrival(make(i, j, k + 1, w + 1)[0]), query)
                for k, w in np.ndindex(9, 9)
            ]
        )
        assert vals.max() < len(bins)
        for mode in range(2):
            offset = len(bins) * np.arange(len(query[1]))
            linear = (vals[:, mode, :] + offset).ravel()
            for arm, weights_arm in enumerate(weights):
                counts = np.bincount(
                    linear,
                    weights=np.repeat(weights_arm.ravel(), len(offset)),
                    minlength=len(offset) * len(bins),
                ).reshape(len(offset), -1)
                hist[arm, mode] += counts.astype(np.int64)
    total = len(tags) * 2592
    assert np.all(hist.sum(axis=-1) == total)
    rng = np.random.default_rng(SEED)
    max_time = np.max(np.nonzero(hist)[-1])
    hist = hist[..., : max_time + 1]
    bins = bins[: max_time + 1]
    summary = {}
    for mode, name in enumerate(("fixed", "best")):
        numerator = hist[:, mode] @ bins
        edges = numerator.reshape(len(weights), 3, 18) / total
        tapes = np.array(
            [[rng.multinomial(N, p / total) for p in arm] for arm in hist[:, mode]]
        )
        held = (tapes @ bins).reshape(len(weights), 3, 18) / N
        slope = edges[2] - edges[1]
        result = {}
        for tag, reference, sign in (("native", 0, 1), ("reversed_second", 5, -1)):
            ref = edges[reference]
            result[tag] = {
                "count_enclosure": certify_chambers(
                    tapes, hist[:, mode], reference, sign
                ),
                "per_shell": [
                    list(chamber(v)) if chamber(v) is not None else None for v in ref
                ],
                "compensated": fixed_volume(ref[[0, 2]], sign * slope[[0, 2]]),
                "from_counts": fixed_volume(
                    held[reference, [0, 2]],
                    sign * (held[3, [0, 2]] - held[4, [0, 2]]) / 72,
                ),
            }
        summary[name] = {
            "area_volumes": result,
            "time_sums": numerator.tolist(),
            "counts": tapes.tolist(),
            "maximum_absolute_source_time_change": float(abs(slope).max()),
        }
    return {
        "seed": SEED,
        "shots_per_edge_arm": N,
        "per_arm_exact_denominator": total,
        "maximum_elapsed_appointments": int(max_time),
        "classes": pools,
        "radii": radii,
        "bearings": bearings.tolist(),
        "plans": plans,
        "ordinary": local_swap(),
        "angular_surveys": angle_test(graphs, positions),
        "triangle_tests": tri,
        "shells": summary,
        "inputs": {
            f: sha256((HERE / f).read_bytes()).hexdigest()
            for f in (
                "PAIRED_RANGING.md",
                "three_regions.json",
                "settled_midpoint_force.json",
                "word_routes.json",
            )
        },
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = json.loads(json.dumps(run(), allow_nan=False))
    path = HERE / "paired_ranging.json"
    if args.write:
        path.write_text(json.dumps(result, indent=1) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), result)
    print(result["triangle_tests"])
    for name, summary in result["shells"].items():
        print(name, {k: z["compensated"] for k, z in summary["area_volumes"].items()})
