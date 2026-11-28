"""Fixed-volume radar with the native home loads present throughout the survey."""

from itertools import product

import numpy as np
from scipy.sparse.csgraph import dijkstra

from . import closed_balls
from . import exchange_geometry
from . import sector_boundary_clock


def routes(graph, words, homes, center, max_volume):
    materials, groups = np.unique(words[:, : len(homes)], axis=0, return_inverse=True)
    swap = exchange_geometry.contact_swap().argmax(0)
    table = graph.toarray() != 0
    n = len(table)
    trips = 0
    distances = []
    for loads in materials:
        state = np.zeros(n, int)
        state[homes] = loads
        occupied = set(np.flatnonzero(state))
        assert center not in occupied
        live = graph.copy()
        live.data[:] = 1 + (state[live.indices] > 0)
        dist, pred = dijkstra(live, indices=center, return_predecessors=True)
        distances.append(dist.astype(int))
        # The positive weights give a lower bound for any route preserving
        # each intermediate resident. Supply actual unused ordinary neighbors
        # as storage for every shortest route needed by this survey.
        horizon = np.sort(dist)[max_volume - 1]
        for target in np.flatnonzero(dist <= horizon):
            if target == center:
                continue
            path = [int(target)]
            while path[-1] != center:
                path.append(int(pred[path[-1]]))
            path = path[::-1]
            on = sorted(set(path) & occupied)
            options = [
                list(set(np.flatnonzero(table[i])) - set(path) - occupied) for i in on
            ]
            rails = next(
                rails for rails in product(*options) if len(set(rails)) == len(rails)
            )
            for load in (1, 2):
                current = state.copy()

                def pulse(i, j):
                    assert table[i, j]
                    current[i], current[j] = divmod(
                        int(swap[3 * current[i] + current[j]]), 3
                    )

                for i, j in zip(on, rails, strict=True):
                    pulse(i, j)
                current[center] = load
                for i, j in zip(path[:-1], path[1:], strict=True):
                    pulse(i, j)
                assert current[target] == load
                for i, j in zip(path[-2::-1], path[:0:-1], strict=True):
                    pulse(j, i)
                assert current[center] == load
                current[center] = 0
                for i, j in zip(on, rails, strict=True):
                    pulse(i, j)
                assert (
                    np.array_equal(current, state)
                    and len(path) - 1 + len(rails) == dist[target]
                )
                trips += 1
    return np.array(distances), groups, trips


def evaluate(
    words,
    homes,
    seams,
    board,
    center,
    sources,
    ref,
    clock,
    warming,
    observed,
    stat,
    acquired,
):
    graph, hops, _ = closed_balls.radar.chart(board, seams)
    u, v = np.nonzero(np.triu(graph.toarray(), 1))
    volume0 = int(sum(hops[center] < 2))
    max_volume = int(sum(hops[center] < 4))
    dist, g, trips = routes(graph, words, homes, center, max_volume)
    material = np.unique(words[:, :4], axis=0)
    durations = sector_boundary_clock.rates(graph, material, homes, u, v)
    before, after, ordinary_route = sector_boundary_clock.move_double(
        graph, homes, words, acquired
    )
    empty = words[:, 1] == 0
    six_meter = None
    rows = {}
    markers = np.arange(board.n)
    selection = [
        volume0 - 1,
        volume0,
        volume0 + 1,
        int(sum(hops[center] < 3)),
        max_volume,
    ]
    for desc in (False, True):
        order = markers[::-1] if desc else markers
        rank = np.argsort(
            order[np.argsort(dist[:, order], axis=1, kind="stable")], axis=1
        )
        for volume in selection:
            inside = rank < volume
            cut = inside[:, u] ^ inside[:, v]
            area = cut.sum(1)[g]
            if not desc and volume == volume0:
                six_meter = cut.sum(1).tolist()
            # Six times the rates in empty-link units are exact integers.
            capacity = (cut * (6 // durations)).sum(1)[g]
            assert np.all(6 % durations == 0)
            mask = np.column_stack(
                (inside[:, homes], *[inside[:, x] & inside[:, y] for x, y in seams])
            )[g]
            raw = np.column_stack(
                (np.packbits(inside[g], axis=1), np.where(mask, words, -1))
            )
            _, ids, count = np.unique(
                raw, axis=0, return_inverse=True, return_counts=True
            )
            probs = np.bincount(ids, weights=ref)
            assert min(probs) > 0
            meter = -np.log(probs[ids])
            heat_area = float(warming @ area)
            heat_tangent = float(warming @ meter)
            eta = None if abs(heat_area) < 1e-12 else -heat_tangent / heat_area
            move_entropy = float(meter[after] - meter[before])
            move_area = int(area[after] - area[before])
            move_cap = int(capacity[after] - capacity[before]) / 6
            move_eta = -move_entropy / move_area if move_area else None
            rate_eta = -move_entropy / move_cap if move_cap else None
            per_bin = np.array(
                [np.unique(area[ids == i])[0] for i in range(len(count))]
            )
            assert np.array_equal(per_bin[ids], area)
            band = np.column_stack([np.bincount(ids, weights=z) for z in clock.T])
            values = -np.log(band[ids, ::-1])
            assert np.all((values[:, 0] <= meter) & (meter <= values[:, 1]))
            by_rate = {}
            for div, (d, err) in sources.items():
                da, ds = float(d @ (area - area[0])), float(d @ (meter - meter[0]))
                dcap = float(d @ (capacity - capacity[0])) / 6
                result = dict(
                    area_change=da,
                    rate_weighted_area_change=dcap,
                    moved_double_residual=None
                    if move_eta is None
                    else ds + move_eta * da,
                    moved_double_rate_residual=None
                    if rate_eta is None
                    else ds + rate_eta * dcap,
                    retained_sharp_tangent_nats=ds,
                    ordinary_stationarity_residual=None
                    if eta is None
                    else ds + eta * da,
                )
                if div in observed:
                    # These specific held checks are downstream of the native
                    # acquisition. Use independent bounds simultaneous over the
                    # five volumes, two orderings, two rates and both meters.
                    width = stat[div] * np.ptp(area)
                    ad = float(observed[div] @ (area - area[0]))
                    q = observed[div] @ values.mean(1)
                    eps = (
                        stat[div] * (values[:, 1].max() - values[:, 0].min())
                        + abs(observed[div]) @ (np.ptp(values, axis=1) / 2)
                        + 1e-10
                    )
                    result.update(
                        area_from_counts=[ad - width, ad + width],
                        entropy_from_counts=[q - eps, q + eps],
                    )
                by_rate[str(div)] = result
            binary = None
            if len(count) == 2 and per_bin[0] != per_bin[1]:
                # At a common reference any probability perturbation has one
                # coordinate, so a single ratio is valid for all preparations.
                odds = np.log(band[1] / band[0, ::-1])
                slopes = np.sort(odds / (per_bin[1] - per_bin[0]))
                assert len(set(ids[empty])) == len(set(ids[~empty])) == 1
                empty_gain = int(area[empty][0] - area[~empty][0])
                binary = dict(
                    extra_boundary_links_per_empty_home=empty_gain,
                    area_by_record=per_bin.tolist(),
                    probabilities=probs.tolist(),
                    entropy_per_area_bounds=slopes.tolist(),
                )
            rows[f"{'descending' if desc else 'ascending'}/{volume}"] = dict(
                distinct_record_outputs=len(count),
                distinct_geometries=len(np.unique(inside, axis=0)),
                area_range=[int(area.min()), int(area.max())],
                rate_weighted_area_range=[
                    int(capacity.min()) / 6,
                    int(capacity.max()) / 6,
                ],
                ordinary_load_inside_before_after=[
                    int((words[i, :4] * mask[i, :4]).sum()) for i in (before, after)
                ],
                moved_double=dict(
                    entropy_change=move_entropy,
                    area_change=move_area,
                    rate_weighted_area_change=move_cap,
                    eta=move_eta,
                    rate_eta=rate_eta,
                ),
                ordinary_area_change=heat_area,
                ordinary_entropy_tangent=heat_tangent,
                ordinary_eta=eta,
                binary_identity=binary,
                rates=by_rate,
            )
    return dict(
        center=center,
        six_marker_surface_by_inventory=six_meter,
        volumes=selection,
        native_material_round_trips=trips,
        ordinary_transfer=ordinary_route,
        ordinary_words=[before, after],
        links_measured_per_inventory=len(u),
        distinct_home_inventories=len(dist),
        surveys=rows,
    )
