"""Insert and remove physical rim plugs, measuring successive aperture admission probabilities."""

from collections import Counter
from itertools import combinations, permutations
from math import comb, log, sqrt
import numpy as np
from model.candidate_laws.interval_echo.clock import Registers
from model.executable_model import federation
from .. import clock_ball, measure as region


def run(rng, shots, stocks):
    car = federation(0)
    ports, contacts, border, _ = region.apparatus(car, [0])
    reach = clock_ball.shortest_chains(ports, contacts, clock_ball.CENTER)
    inner = [p for p in ports if len(reach[p]) < 2]
    rim = [p for p in ports if len(reach[p]) == 2]
    words = list(combinations(sorted(inner + rim), 2))
    gateway = min(border)
    blank = border[gateway]
    feeds = [p for p in list(range(12, car.n)) if p != blank][: len(rim)]
    reg = Registers(car, seed=31059)
    owner = set(range(car.cells))
    moves = {p: reg.route(p, gateway, {0}) for p in inner + rim if p != gateway}
    routes = {(f, p): reg.route(f, p, owner) for f in feeds for p in rim}

    def swap(route):
        for a, b in route:
            reg.swap(a, b, owner)

    def occupied(p):
        path = moves.get(p, ())
        swap(path)
        bit = reg.flag(gateway, blank, owner)
        swap(reversed(path))
        return bit

    records = []
    for w in words:
        reg.loads.fill(0)
        reg.loads[list(w) + feeds] = 2
        ref = reg.loads.copy()
        k = 2 - sum(occupied(p) for p in inner)
        stop = len(rim) - k
        outcomes = Counter()
        for order in permutations(rim):
            parked = []
            for feed, p in zip(feeds[:stop], order):
                if occupied(p):
                    break
                route = routes[feed, p]
                swap(route)
                parked.append(route)
            outcomes[len(parked)] += 1
            for route in reversed(parked):
                swap(reversed(route))
            assert np.array_equal(reg.loads, ref)
        n = sum(outcomes.values())
        assert outcomes[stop] * comb(len(rim), k) == n
        records.append(
            dict(
                residents=w,
                k=k,
                stop=stop,
                permutations_by_inserted_prefix=[outcomes[s] for s in range(stop + 1)],
            )
        )
    # Independent per-class calibration and a fresh incoming classification for each stock.
    representatives = {row["k"]: row for row in records}
    epsilon = sqrt(log(2 * 2000 / 0.01) / (2 * shots))
    scores, cal = {}, {}
    for k, row in representatives.items():
        table = np.array(row["permutations_by_inserted_prefix"])
        trials = rng.multinomial(shots, table / table.sum()).tolist()
        measured = trials[-1] / shots
        scores[k] = [-log(min(1, measured + epsilon)), -log(measured - epsilon)]
        cal[k] = dict(capped_prefix_counts=trials, work_interval_nats=scores[k])
    arms = {}
    for name, stock in stocks.items():
        distribution = np.array(
            [sum(sum(p in rim for p in w) == k for w in stock) for k in range(3)]
        ) / len(stock)
        counts = rng.multinomial(shots, distribution).tolist()
        bounds = [
            sum(
                max(0, n / shots - epsilon) * scores[k][0] for k, n in enumerate(counts)
            ),
            sum(
                min(1, n / shots + epsilon) * scores[k][1] for k, n in enumerate(counts)
            ),
        ]
        arms[name] = dict(
            class_counts=counts,
            mean_work_nats=bounds,
            negative_log_pooled_success=-log(
                sum(distribution[k] / comb(5, k) for k in range(3))
            ),
        )

    def change(a, b):
        return [
            arms[b]["mean_work_nats"][0] - arms[a]["mean_work_nats"][1],
            arms[b]["mean_work_nats"][1] - arms[a]["mean_work_nats"][0],
        ]

    return dict(
        rim=rim,
        inner=inner,
        feeders=feeds,
        reader_gateway=gateway,
        blank=blank,
        permutations_per_word=native_orders(records),
        native_gates=dict(reg.gates),
        closed_word_histograms=records,
        standards=cal,
        arms=arms,
        measured_ordinary_nats=change("zero", "one"),
        measured_held_nats=change("split", "spread"),
    )


def native_orders(rows):
    counts = {sum(row["permutations_by_inserted_prefix"]) for row in rows}
    assert counts == {120}
    return counts.pop()
