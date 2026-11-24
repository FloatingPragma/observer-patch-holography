"""Three prepared walls on native linked belts and free breakup on one comparison ring."""

import argparse
from collections import Counter, defaultdict
from itertools import permutations, product, combinations
import json
from math import factorial
from pathlib import Path

import networkx as nx
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.contrast_memory_spherical_v1.sectors import BinarySector
from model.candidate_laws.interval_echo import candidate, currents
from model.executable_model import carrier
from . import belt


def path():
    board = carrier.federation(0)
    graph = nx.Graph((u // 12, v // 12) for u, v in board.external)
    owners = min(
        (a, b, c)
        for b in graph
        for a, c in combinations(sorted(graph[b]), 2)
        if not graph.has_edge(a, c)
    )
    bridges = [e for e in board.external if all(v // 12 in owners for v in e)]
    assert len(bridges) == 2
    orders = []
    for owner in owners:
        ends = {u % 12 for e in bridges for u in e if u // 12 == owner}
        orders.append(
            next(order for order in belt.belts().values() if ends <= set(order))
        )
    address = [12 * owner + p for owner, order in zip(owners, orders) for p in order]
    contacts = [(address.index(a), address.index(b)) for a, b in bridges]
    edges = [
        (10 * j + a, 10 * j + b) for j in range(3) for a, b in belt.ring(10).seams
    ] + contacts
    assert {tuple(sorted((address[a], address[b]))) for a, b in edges} == {
        (a, b) for a, b in board.seams if a in address and b in address
    }
    return dict(
        owners=owners, addresses=address, bridges=bridges, local_bridges=contacts
    ), edges


def mask(sites):
    return sum(1 << int(i) for i in sites)


def hop_list(edges):
    return [(1 << a, 1 << b, (1 << a) | (1 << b)) for a, b in edges]


def neighbours(word, moves):
    for a, b, flip in moves:
        if bool(word & a) != bool(word & b):
            yield word ^ flip


def diagonal(word, moves, kappa):
    return sum(
        0.5 * (bool(word & a) != bool(word & b))
        - 2 * kappa * (bool(word & a) and bool(word & b))
        for a, b, _ in moves
    )


def closed_ball(source, moves, kappa, allowed):
    """Full-H measurement bounds via the residual outside three exact hop layers."""
    depth = 3
    layers = [[source]]
    known = {source}
    norms = []
    for _ in range(depth + 1):
        next_counts = Counter()
        max_out = 0
        for word in layers[-1]:
            outgoing = [v for v in neighbours(word, moves) if v not in known]
            max_out = max(max_out, len(outgoing))
            next_counts.update(outgoing)
        norms.append(np.sqrt(max_out * max(next_counts.values())) / 2)
        layers.append(sorted(next_counts))
        known.update(next_counts)
    states = sum(layers[: depth + 1], [])
    index = {s: i for i, s in enumerate(states)}
    rows, cols = [], []
    for i, s in enumerate(states):
        for v in neighbours(s, moves):
            if v in index:
                rows.append(index[v])
                cols.append(i)
    h = sparse.csr_matrix(
        (np.full(len(rows), -0.5), (rows, cols)), shape=(len(states), len(states))
    )
    h += sparse.diags([diagonal(s, moves, kappa) for s in states])
    start = np.zeros(len(states))
    start[0] = 1
    inside = np.array([s in allowed for s in states])
    reads = []
    for t in (1 / 32, 1 / 16, 1 / 8):
        psi = expm_multiply(-1j * t * h, start)
        # Unitary Duhamel at each new layer gives the product and ordered simplex volume.
        error = np.prod(norms) * t ** (depth + 1) / factorial(depth + 1)
        p = np.linalg.norm(psi[~inside]) ** 2
        reads.append(
            dict(
                time=t,
                ball_outside_three_walls=p,
                amplitude_error_upper=error,
                full_outside_bounds=[
                    max(0.0, np.sqrt(p) - error) ** 2,
                    min(1.0, (np.sqrt(p) + error) ** 2),
                ],
            )
        )
    return dict(
        layer_sizes=list(map(len, layers)), cross_block_norm_upper=norms, reads=reads
    )


def native():
    geometry, edges = path()
    moves = hop_list(edges)
    bridges = hop_list(geometry["local_bridges"])
    local = {}
    for n in (3, 4):
        sec = BinarySector(belt.ring(10), n)
        local[n] = {
            mask(sec.states[i]): belt.wall(sec.states[i], 10)[2]
            for i in belt.selected(sec)
        }
    result = {}
    for labels in ((3, 3, 4), (3, 4, 4)):
        orders = sorted(set(permutations(labels)))
        all_walls = {}
        for order in orders:
            for parts in product(*(local[n] for n in order)):
                word = sum(s << (10 * i) for i, s in enumerate(parts))
                centres = [local[n][s] for n, s in zip(order, parts)]
                distance = sum(
                    abs(float(belt.signed_distance(centres[p // 10], p % 10, 10)))
                    for edge in geometry["local_bridges"]
                    for p in edge
                )
                all_walls[word] = (order, distance)
        exit_counts = Counter(
            v
            for word in all_walls
            for v in neighbours(word, moves)
            if v not in all_walls
        )
        unique_exits = [
            sum(exit_counts[v] == 1 for v in neighbours(word, moves))
            for word in all_walls
        ]
        panels = []
        for order in orders:
            records = defaultdict(list)
            for word, (types, dist) in all_walls.items():
                if order != types:
                    continue
                # Disjoint bridge pulses have independent stay/swap amplitudes on a sharp word.
                outcomes = {word: 1.0}
                for bridge in bridges:
                    nxt = defaultdict(float)
                    for w, weight in outcomes.items():
                        v = next(neighbours(w, [bridge]), w)
                        nxt[w] += weight / 2
                        nxt[v] += weight / 2
                    outcomes = nxt
                stay = sum(weight for v, weight in outcomes.items() if v in all_walls)
                at_pi = word
                for bridge in bridges:
                    at_pi = next(neighbours(at_pi, [bridge]), at_pi)
                matches = sum(bool(word & a) != bool(word & b) for a, b, _ in bridges)
                records[dist].append(
                    (
                        word,
                        diagonal(word, bridges, candidate.KAPPA),
                        (candidate.KAPPA + 0.5) * matches,
                        stay,
                        int(at_pi in all_walls),
                    )
                )
            sets = []
            for dist, rows in sorted(records.items()):
                values = np.array([r[1:] for r in rows])
                sets.append(
                    dict(
                        distance=dist,
                        count=len(rows),
                        connected_energy_values=sorted(set(values[:, 0])),
                        bare_contact_values=sorted(set(values[:, 1])),
                        mean_connected_energy=values[:, 0].mean(),
                        pulse_half_survival=values[:, 2].mean(),
                        pulse_pi_survival=values[:, 3].mean(),
                    )
                )
            flight = {}
            if order == orders[0]:
                for key, dist in [("near", min(records)), ("far", max(records))]:
                    source = min(r[0] for r in records[dist])
                    flight[key] = dict(
                        source_bits=source,
                        **closed_ball(source, moves, candidate.KAPPA, all_walls),
                    )
            panels.append(
                dict(order=order, placements_by_distance=sets, free_flight=flight)
            )
        result[str(sum(labels))] = dict(
            placements=len(all_walls),
            whole_excess=sum(labels) - 10,
            reached_complement=len(exit_counts),
            min_unique_exits=min(unique_exits),
            squared_leakage_norm_lower=min(unique_exits) / 4,
            orderings=panels,
        )
    return dict(geometry=geometry, sectors=result)


def gap_arc(sites, length, ones, fours):
    gaps = [(a, (b - a) % length) for a, b in zip(sites, (*sites[1:], sites[0]))]
    if sorted(g for _, g in gaps) != sorted(
        [1] * ones + [4] * fours + [3] * (len(sites) - 3)
    ):
        return None
    centres = sorted((a + g / 2) % length for a, g in gaps if g != 3)
    return length - max(
        (b - a) % length for a, b in zip(centres, (*centres[1:], centres[0]))
    )


def breakup(length, ones):
    n = length // 3 + ones - 1
    sites = list(combinations(range(length), n))
    words = list(map(mask, sites))
    index = {w: i for i, w in enumerate(words)}
    moves = hop_list(belt.ring(length).seams)
    row, col = [], []
    for i, w in enumerate(words):
        for v in neighbours(w, moves):
            row.append(index[v])
            col.append(i)
    hopping = sparse.csr_matrix(
        (np.full(len(row), -0.5), (row, col)), shape=(len(words), len(words))
    )
    arcs = np.array([gap_arc(s, length, ones, 3 - ones) for s in sites], dtype=float)
    valid = np.isfinite(arcs)
    low, high = np.min(arcs[valid]), np.max(arcs[valid])
    launch = [int(np.flatnonzero(arcs == v)[0]) for v in (low, high)]
    state = np.zeros((len(words), 2), complex)
    state[launch, [0, 1]] = 1
    panels = []
    for kappa in (candidate.KAPPA, 0.0):
        h = hopping + sparse.diags([diagonal(w, moves, kappa) for w in words])
        psi = state.copy()
        last = 0
        for t in (0, 1, 4, 16):
            if t != last:
                psi = expm_multiply(-1j * (t - last) * h, psi)
            last = t
            prob = abs(psi) ** 2
            np.testing.assert_allclose(prob.sum(axis=0), 1, atol=2e-11)
            panels.append(
                dict(
                    time=t,
                    kappa=kappa,
                    three_wall=prob[valid].sum(axis=0),
                    arc_weights={
                        str(a): prob[arcs == a].sum(axis=0)
                        for a in np.unique(arcs[valid])
                    },
                    breakup=prob[~valid].sum(axis=0),
                )
            )
    return dict(
        length=length,
        loads=n,
        dimension=len(words),
        three_wall_words=int(valid.sum()),
        closest_arc=low,
        farthest_arc=high,
        sources=[sites[i] for i in launch],
        reads=panels,
    )


def run():
    return currents._portable(
        dict(
            native=native(),
            single_string=[
                breakup(length, ones) for length in (12, 15) for ones in (1, 2)
            ],
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    data = run()
    encoded = json.dumps(data, indent=2, allow_nan=False) + "\n"
    with args.output.open("x") as output:
        output.write(encoded)
