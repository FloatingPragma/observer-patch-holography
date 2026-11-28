"""Thermal source, material echo ranges and entropy at exactly fixed marker volume."""

import argparse
from hashlib import sha256
import json
from itertools import product

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import dijkstra, shortest_path

from . import thermal_pair as lab
from . import exchange_geometry
from . import thermal_ordinary as standard

HERE = lab.HERE
RADII = (2, 4, 6, 8, 9, 10, 11, 12)
N = 2**55
BATCH = 2**47
LOG_RISK = np.log(2 * 10**6 / 0.01)


def chart(board, edges):
    adjacency = board.adjacency().tolil().astype(float)
    for a, b in edges:
        adjacency[a, b] = adjacency[b, a] = 0
    live = adjacency.tocsr()
    hops = shortest_path(live, unweighted=True, directed=False).astype(int)
    span = hops[:, np.ravel(edges)].max(axis=1)
    center = int(np.flatnonzero(span == min(span))[0])
    return live, hops, center


def thermal(board, edges, reduce_sectors=lab.motion.binary_responses):
    h, z, i, j, _ = lab.motion.hamiltonian(board, edges)
    fixed = int(
        np.flatnonzero(
            lab.motion.cold.distances(board, sorted(set(np.ravel(edges)))) == 1
        )[0]
    )
    keep = (i != fixed) & (j != fixed)
    remap = np.full(len(i), -1, int)
    remap[keep] = np.arange(keep.sum())
    sites = remap[z]
    i, j = i[keep], j[keep]
    h = h[keep][:, keep]
    adj = board.adjacency().toarray().astype(int)
    near = adj[fixed, i] + adj[fixed, j]
    placed = h - sparse.diags((2 * lab.motion.K + 1) * near, format="csr", dtype=float)
    distance = shortest_path(h != 0, directed=False, indices=sites[0], unweighted=True)[
        sites[1]
    ]
    dist = {}
    error = {}
    for name, mat, degree in (("blank", h, 14), ("ordinary", placed, 18)):
        w, e = lab.kernel.baseline(mat, degree=degree)
        inc, cert = lab.kernel.single_updates(mat, sites)
        joint_bound = lab.kernel.mixed_bounds(mat, sites, distance)[
            "upper_each_diagonal"
        ]
        sectors = {}
        relative = 0.0
        for a, b in product(lab.STRENGTHS, repeat=2):
            x = inc[0].get(a, 0) + inc[1].get(b, 0)
            weights = w + x
            eps = e + 2 * cert["tail"] + x * 2 * cert["rounding_relative"] + joint_bound
            rel = max(eps / weights)
            relative = max(relative, float(2 * rel / (1 - rel)))
            sectors[f"{-a},{-b}"] = weights / weights.sum()
        dist[name] = reduce_sectors(sectors)
        error[name] = relative + 1e-12
    old = json.loads((HERE / "pair_contact.json").read_bytes())
    table = np.asarray(old["accepted_weight"])
    p = table.ravel() / table.sum()
    q = np.outer(table.sum(0), table.sum(1)).T.ravel() / table.sum() ** 2
    return (
        i,
        j,
        fixed,
        dict(
            native=p[:, None] * dist["blank"],
            ordinary=p[:, None] * dist["ordinary"],
            product=q[:, None] * dist["blank"],
        ),
        error,
    )


def ranges(live, center, i, j):
    # Charge the target vertex when entered. The fixed cost of clearing the
    # origin is removed on every path, including the bare range calibration.
    out = np.empty((len(i), live.shape[0]), np.int16)
    mat = live.copy()
    cols = mat.indices
    for k, (u, v) in enumerate(zip(i, j, strict=True)):
        mat.data[:] = 1 + (cols == u).astype(np.int8) + (cols == v).astype(np.int8)
        out[k] = dijkstra(mat, indices=center, directed=True)
    bare = shortest_path(live, indices=center, directed=True, unweighted=True)
    assert np.all((bare <= out) & (out <= bare + 2))
    return out


def certify_routing(board, live, edges):
    swap = exchange_geometry.contact_swap()
    for background in product((0, 1), repeat=3):
        for protect in (False, True):
            # Original occupants have their own empty side rail. A separately
            # tagged double is inserted after parking the sending occupant.
            state = list(background) + [0] * 3
            cleared = [k for k, a in enumerate(background) if a and (k == 0 or protect)]

            def pulse(a, b):
                word = 3 * state[a] + state[b]
                result = int(np.argmax(swap[:, word]))
                state[a], state[b] = divmod(result, 3)

            for k in cleared:
                pulse(k, k + 3)
            state[0] = 2
            for a, b in ((0, 1), (1, 2), (2, 1), (1, 0)):
                pulse(a, b)
            assert state[0] == 2
            state[0] = 0
            for k in reversed(cleared):
                pulse(k, k + 3)
            assert state == [*background, 0, 0, 0]
    # In every admitted written sector, all material port words stay fixed.
    lab.contact.blocks()
    assert len(board.seams) - live.nnz // 2 == len(edges)


def certify_paths(live, center, i, j, dist):
    swap = exchange_geometry.contact_swap().argmax(0)
    n = live.shape[0]
    tested = 0
    for k in (
        0,
        len(i) // 2,
        len(i) - 1,
        int(np.flatnonzero((i == center) | (j == center))[0]),
    ):
        for protect in (False, True):
            graph = live.copy()
            occupied = {int(i[k]), int(j[k])}
            # Explicit integer transition cost, including the count for a parked single.
            graph.data[:] = 1 + np.isin(graph.indices, list(occupied)).astype(
                int
            ) * int(protect)
            travel, pred = dijkstra(graph, indices=center, return_predecessors=True)
            for dest in (int(i[k]), int(j[k]), int(np.argmax(travel))):
                path = [dest]
                while path[-1] != center:
                    path.append(int(pred[path[-1]]))
                path.reverse()
                route = list(zip(path[:-1], path[1:], strict=True))
                initial = np.zeros(2 * n, int)
                initial[list(occupied)] = 1
                state = initial.copy()
                steps = 0

                def pulse(a, b):
                    nonlocal steps
                    state[a], state[b] = divmod(int(swap[3 * state[a] + state[b]]), 3)
                    steps += 1

                parked = [a for a in path if a in occupied and (a == center or protect)]
                for a in parked:
                    pulse(a, a + n)
                state[center] = 2
                for a, b in route:
                    pulse(a, b)
                assert state[dest] == 2
                for a, b in reversed(route):
                    pulse(b, a)
                state[center] = 0
                for a in reversed(parked):
                    pulse(a, a + n)
                np.testing.assert_array_equal(state, initial)
                assert steps - 2 * (center in occupied) == 2 * travel[dest]
                if protect:
                    assert travel[dest] == dist[k, dest]
                tested += 1
    return tested


def matched_positions(weights, rng, rel):
    p = weights["native"].sum(1)
    q = weights["product"].sum(1)
    sharp = weights["native"] / p[:, None]
    low = sharp.min(0)
    success = low.sum()
    accept = low / sharp
    assert np.all((0 < accept) & (accept <= 1))
    conditional = low / success
    output = {}
    for name, source in (("native", p), ("product", q)):
        before = binomial_hist(rng, weights[name])
        kept = rng.binomial(before, accept)
        after = (weights[name] * accept) / success
        np.testing.assert_allclose(after, source[:, None] * conditional, rtol=3e-14)
        output[name] = dict(
            accepted=kept.sum(1).tolist(), rejected=(before - kept).sum(1).tolist()
        )
    native_positions = weights["native"].sum(0)
    diff = weights["product"].sum(0) - native_positions
    tv = float(abs(diff).sum() / 2)
    # Positive vs negative part of the signed difference attains the bound.
    assert abs(diff[diff > 0].sum() - tv) < 1e-13
    charge = float(-(q - p) @ np.log(p))
    return dict(
        success=float(success),
        relative_kernel_error=rel,
        counts=output,
        common_head_entropy_tangent=charge,
        unselected_position_tv=band(tv, rel + 1e-12).tolist(),
        near_match_tv_bound=2 * rel / (1 - rel),
        readout="one acceptance flag per attempt; head order 00, 01, 10, 11",
    )


def binomial_hist(rng, p):
    flat = p.ravel()
    pnorm = flat / flat.sum()
    total = np.zeros(flat.shape, dtype=np.int64)
    for _ in range(N // (8 * BATCH)):
        total += rng.multinomial(BATCH, pnorm, size=8).sum(0)
    assert total.sum() == N
    return total.reshape(p.shape)


def geometry(order, volume, i, j, edges, u, v, internal):
    n = order.shape[1]
    inside = np.zeros(order.shape, dtype=bool)
    np.put_along_axis(inside, order[:, :volume], True, axis=1)
    assert np.all(inside.sum(1) == volume)
    seen = np.stack([inside[:, a] & inside[:, b] for a, b in edges])
    surface = []
    for select in (internal, ~internal):
        surface.append(
            np.count_nonzero(inside[:, u[select]] ^ inside[:, v[select]], axis=1)
        )
    _, ball = np.unique(np.packbits(inside, axis=1), axis=0, return_inverse=True)
    # Store only the visible identities, sorted; no index of an unseen load leaks.
    row = np.arange(len(i))
    a = np.where(inside[row, i], i, n)
    b = np.where(inside[row, j], j, n)
    position = np.minimum(a, b) * (n + 1) + np.maximum(a, b)
    codes = np.array(
        [
            (
                (
                    ball * 9
                    + 3 * np.where(seen[0], h // 2, 2)
                    + np.where(seen[1], h % 2, 2)
                )
                * (n + 1) ** 2
                + position
            )
            for h in range(4)
        ]
    )
    _, record = np.unique(codes, return_inverse=True)
    return np.array(surface), ball, record.reshape(codes.shape), seen.all(axis=0)


def difference(a, b):
    return np.asarray(a) - np.asarray(b)[::-1]


def band(v, e):
    return np.array([v - e, v + e])


def read_meters(weights, ball, record, surface, counts, area_counts, rel_error, delta):
    empty, ordinary, held = (weights[k] for k in ("native", "ordinary", "product"))
    shape = record.shape

    def coarse(q, ix):
        return np.bincount(ix.ravel(), q.ravel())

    ball_idx = np.broadcast_to(ball, shape)
    old = coarse(empty, record)
    old_ball = coarse(empty, ball_idx)
    support = empty > 0
    # The derivative of entropy on this face requires common support.
    assert np.array_equal(ordinary > 0, support) and np.array_equal(held > 0, support)
    assert np.all(old[record[support]] > 0) and old_ball.min() > 0

    def scores(ref):
        vals = coarse(ref, record)
        balls = coarse(ref, ball_idx)
        assert np.all(vals[record[support]] > 0) and min(balls) > 0
        # Impossible head/flag records have zero weight on every arm.
        full = -np.log(np.where(vals > 0, vals, 1)[record])
        conditional = full + np.log(balls[ball_idx])
        return conditional, full

    score, full = scores(empty)
    ref_count = counts["reference"] / N
    observed_score, _ = scores(ref_count)
    # Bernstein for a Bernoulli cell: given model p>=pmin, simultaneously
    # bound relative count error at every cell. Grouping preserves this bound.
    pmin = float(empty[support].min())
    rel = np.sqrt(2 * LOG_RISK / (N * pmin)) + 2 * LOG_RISK / (3 * N * pmin)
    # Grouped positive sums may accumulate one rounding error per input cell.
    # Include that loss twice for a conditional log probability.
    rounding = 4 * empty.size * np.finfo(float).eps
    ref_err = 2 * (-np.log1p(-rel) - np.log1p(-rel_error["blank"]) + rounding)

    def entropy(q):
        a = coarse(q, record)
        b = coarse(q, ball_idx)
        pos = a > 0
        return -a[pos] @ np.log(a[pos]) + b @ np.log(b)

    tests = {}
    for tag, q, e0 in (
        ("ordinary", ordinary, rel_error["ordinary"]),
        ("product", held, rel_error["blank"]),
    ):
        diff = q - empty
        tv = abs(diff).sum() / 2

        def assess(meter, obs_meter=None, reference_error=0.0):
            if obs_meter is None:
                obs_meter = meter
            span = np.ptp(meter)
            span_obs = np.ptp(obs_meter)
            # Each exact equilibrium differs from the numerical probability law
            # by relative e on every atom. Meter centers may therefore be subtracted.
            numerical = (e0 + rel_error["blank"] + 2 * rounding) * span / 2
            model = float(np.sum(diff * meter))
            bank = counts if meter.ndim == 2 else area_counts
            obs = float(np.sum((bank[tag] - bank["native"]) / N * obs_meter))
            shot = span_obs * np.sqrt(2 * LOG_RISK / N)
            answer = dict(
                model=model,
                from_counts=band(obs, shot + numerical + reference_error).tolist(),
            )
            if tag == "product":
                # Average separate sharp-head references and use the measured signed source.
                sharp = empty / empty.sum(1, keepdims=True)
                means = np.sum(sharp * meter, axis=1)
                mixed = means @ [1.0, -1.0, -1.0, 1.0]
                # For the source score the finite heat approximation also changes the
                # reference logarithm; both its numerator and denominator are bounded.
                score_error = (
                    2 * (-np.log1p(-rel_error["blank"]) + rounding)
                    if meter.ndim == 2
                    else 0.0
                )
                eps = (
                    2 * span * (rel_error["blank"] + rounding) + 4 * score_error + 1e-12
                )
                answer["from_charge"] = lab.probe.mix_bands(
                    delta, band(mixed, eps)
                ).tolist()
            return answer

        source = assess(
            score, observed_score, (2 * tv + e0 + rel_error["blank"]) * ref_err
        )
        measured = [assess(a) for a in (*surface, surface.sum(0))]
        tests[tag] = dict(
            entropy=source,
            full_entropy_tangent=float(np.sum(diff * full)),
            full_step_conditional_entropy=float(entropy(q) - entropy(empty)),
            area=measured,
            tv_full_records=float(tv),
        )
    return tests


def compare(tests):
    cal, held = tests["ordinary"], tests["product"]
    a, b = (x["area"][-1]["from_counts"] for x in (cal, held))
    c, d = (x["entropy"]["from_counts"] for x in (cal, held))
    # eta * delta A + delta S = 0 at exactly matched counting volume.
    coeff = standard.quotient([-c[1], -c[0]], a)
    predicted = None if coeff is None else standard.quotient([-d[1], -d[0]], coeff)
    contrast = difference(lab.probe.mix_bands(a, d), lab.probe.mix_bands(b, c))
    return dict(
        ordinary_entropy_per_area=coeff,
        predicted_held_area=predicted,
        determinant=contrast.tolist(),
        common_gain_excluded=bool(min(contrast) > 0 or max(contrast) < 0),
    )


def run():
    board = lab.motion.cold.federation(0)
    edges = lab.contact.stock.carrier()[1][:2]
    live, hops, center = chart(board, edges)
    certify_routing(board, live, edges)
    i, j, fixed, weights, relative = thermal(board, edges)
    rng = np.random.default_rng(696277)
    counts = {
        k: binomial_hist(rng, weights[tag])
        for k, tag in (
            ("reference", "native"),
            ("native", "native"),
            ("ordinary", "ordinary"),
            ("product", "product"),
        )
    }
    area_counts = {k: binomial_hist(rng, v.sum(axis=0)) for k, v in weights.items()}
    delta = lab.probe.log_factor(
        json.loads((HERE / "pair_contact.json").read_bytes())["source"]
    )
    dist = ranges(live, center, i, j)
    route_checks = certify_paths(live, center, i, j, dist)
    matching = matched_positions(weights, rng, relative["blank"])
    direct = {}
    direct_order = np.broadcast_to(np.argsort(hops[center], kind="stable"), dist.shape)
    u, v = sparse.triu(live, 1).nonzero()
    internal = u // 12 == v // 12
    empty = hops[center]
    result = {}
    for reverse in (False, True):
        asc = np.argsort(dist[:, ::-1] if reverse else dist, axis=1, kind="stable")
        order = board.n - 1 - asc if reverse else asc
        rows = {}
        for r in (*RADII, int(empty.max())):
            volume = int(np.count_nonzero(empty <= r))
            surface, ball, record, seen = geometry(
                order, volume, i, j, edges, u, v, internal
            )
            protected = read_meters(
                weights, ball, record, surface, counts, area_counts, relative, delta
            )
            if not reverse:
                ds, db, dr, dseen = geometry(
                    direct_order, volume, i, j, edges, u, v, internal
                )
                direct[str(r)] = dict(
                    contains_heads=bool(dseen.all()),
                    tests=read_meters(
                        weights, db, dr, ds, counts, area_counts, relative, delta
                    ),
                )
            kth = dist[np.arange(len(i)), order[:, volume - 1]]
            transparent = (empty[u] <= r) ^ (empty[v] <= r)
            rows[str(r)] = dict(
                volume=volume,
                radius_range=list(map(int, [kth.min(), kth.max()])),
                fully_contains_source_everywhere=bool(seen.all()),
                source_coverage={
                    k: float(np.sum(v * seen)) for k, v in weights.items()
                },
                direct_area=[
                    int(sum(transparent & internal)),
                    int(sum(transparent & ~internal)),
                ],
                distinct_ball_masks=len(np.unique(ball)),
                tests=protected,
                comparison=compare(protected),
            )
        result["descending" if reverse else "ascending"] = rows
    import itertools

    def determinants(rows):
        from mpmath import iv

        strongest = None
        for triple in itertools.combinations(
            [
                (r, v[k])
                for r, v in (
                    (r, row["tests"])
                    for r, row in rows.items()
                    if row["fully_contains_source_everywhere"]
                    and sum(row["direct_area"]) > 0
                )
                for k in ("ordinary", "product")
            ],
            3,
        ):
            if not any(t is rows[r]["tests"]["product"] for r, t in triple):
                continue
            matrix = [
                [iv.mpf(t["area"][j]["from_counts"]) for j in (0, 1)]
                + [iv.mpf(t["entropy"]["from_counts"])]
                for _, t in triple
            ]
            a, b, c = matrix
            val = (
                a[0] * (b[1] * c[2] - b[2] * c[1])
                + a[1] * (b[2] * c[0] - b[0] * c[2])
                + a[2] * (b[0] * c[1] - b[1] * c[0])
            )
            lo, hi = lab.interval_output(val)
            margin = max(lo, -hi)
            if strongest is None or margin > strongest[0]:
                strongest = (
                    margin,
                    dict(
                        rows=[
                            [
                                r,
                                "ordinary"
                                if t is rows[r]["tests"]["ordinary"]
                                else "product",
                            ]
                            for r, t in triple
                        ],
                        signed_interval=[lo, hi],
                    ),
                )
        return strongest[1]

    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                HERE / "THERMAL_RADAR.md",
                HERE / "thermal_pair.json",
                HERE / "thermal_ordinary.json",
            )
        },
        center=center,
        ordinary_clamp=fixed,
        source_edges=[list(e) for e in edges],
        live_endpoint_distances=empty[np.ravel(edges)].tolist(),
        per_arm=N,
        numerical_relative=relative,
        count_sha256={
            k: sha256(v.astype("<i8").tobytes()).hexdigest() for k, v in counts.items()
        },
        geometry_count_sha256={
            k: sha256(v.astype("<i8").tobytes()).hexdigest()
            for k, v in area_counts.items()
        },
        charge_cell_displacement=delta.tolist(),
        matching=matching,
        replayed_routes=route_checks,
        direct=direct,
        tests=result,
        two_area_classes={k: determinants(v) for k, v in result.items()},
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    out = run()
    path = HERE / "thermal_radar.json"
    if args.write:
        path.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
    else:
        lab.contact.timing.lab.lab.cold.apparatus.check(
            json.loads(path.read_bytes()), out
        )
    for key, rows in out["tests"].items():
        print(
            key,
            [
                (
                    r,
                    v["fully_contains_source_everywhere"],
                    v["comparison"]["common_gain_excluded"],
                )
                for r, v in rows.items()
            ],
        )
    print(out["two_area_classes"])
