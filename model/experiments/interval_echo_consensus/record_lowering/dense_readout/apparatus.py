"""Compile addressed publication/receiver routes on a registered Interval Echo board."""

from collections import deque
from functools import cache
from itertools import product

import networkx as nx

from model.executable_model.carrier import CELL_SEAMS


SITES = tuple(product(range(4), repeat=3))


def squared_distance(i, j):
    return sum((a - b) ** 2 for a, b in zip(SITES[i], SITES[j], strict=True))


@cache
def cell_layout(gateways):
    """Two independent words, a reader and its reference, all off the gateways."""
    graph = nx.Graph(CELL_SEAMS)
    edges = [pair for pair in CELL_SEAMS if set(pair).isdisjoint(gateways)]
    for public, old, meter in product(edges, repeat=3):
        occupied = set(public + old + meter)
        if len(occupied) != 6:
            continue
        transit = set(graph) - occupied
        if nx.is_connected(graph.subgraph(transit)) and all(
            transit.intersection(graph[p]) for p in [public[0], old[0], meter[0]]
        ):
            return (*public, *old, *meter)
    raise ValueError("No connected protected layout")


def mount(board):
    sites = []
    for owner in range(len(SITES)):
        offset = 12 * owner
        gates = tuple(i for i in range(12) if board.degree[offset + i] == 6)
        sites.append(tuple(offset + a for a in cell_layout(gates)))
    return sites


def paths(board, sites, stale=False):
    """Shortest protected paths with sorted-neighbour breadth-first ties."""
    reserved = {p for site in sites for p in site}
    result = {}
    for receiver, site in enumerate(sites):
        head = site[4]
        toward = {head: None}
        queue = deque([head])
        while queue:
            p = queue.popleft()
            for q in sorted(board.neighbours[p]):
                if q not in toward and q not in reserved:
                    toward[q] = p
                    queue.append(q)
        distance = {head: 0}
        for p, parent in toward.items():
            if parent is not None:
                distance[p] = distance[parent] + 1
        for source, record in enumerate(sites):
            rail = record[2 if stale else 0]
            q = min(
                (v for v in board.neighbours[rail] if v in toward),
                key=lambda v: (distance[v], v),
            )
            path = [rail, q]
            while q != head:
                q = toward[q]
                path.append(q)
            result[receiver, source] = tuple(path)
    return result


def compile_layer(sites, routes):
    """Reserve all hardware of a read through inverse transport; disjoint reads overlap."""
    free = {}
    plan = []
    for receiver in range(len(sites)):
        for source in range(len(sites)):
            if squared_distance(receiver, source) > 4:
                continue
            path = routes[receiver, source]
            ref = sites[receiver][5]
            ports = {*path, ref}
            start = max(free.get(p, 0) for p in ports)
            end = start + 2 * (len(path) - 1) + 3
            for port in ports:
                free[port] = end
            plan.append(
                dict(
                    receiver=receiver,
                    source=source,
                    start=start,
                    end=end,
                    ref=ref,
                    path=path,
                )
            )
    return plan


def clock_trial(routes):
    """First detector slots on independent cohorts, relative to simultaneous launch."""
    # Transport d edges, write the local pointer, then project.
    arrivals = [
        [len(routes[r, s]) + 1 for s in range(len(SITES))] for r in range(len(SITES))
    ]

    def evaluate(receivers):
        cases = [
            (arrivals[r][s], squared_distance(r, s) <= 4, r, s)
            for r in receivers
            for s in range(len(SITES))
        ]
        best = min(
            (sum((t <= k) != inside for t, inside, _, _ in cases), k)
            for k in {0, *(c[0] for c in cases)}
        )
        late = max((t, r, s) for t, inside, r, s in cases if inside)
        early = min((t, r, s) for t, inside, r, s in cases if not inside)
        return dict(
            best_disagreements=best[0],
            earliest_best_deadline=best[1],
            latest_inside=late,
            earliest_outside=early,
        )

    return dict(
        ordered_pair_slots=arrivals,
        pooled=evaluate(range(len(SITES))),
        per_receiver=[evaluate([r]) for r in range(len(SITES))],
    )


def owner_obstruction(board):
    owners = nx.Graph((a // 12, b // 12) for a, b in board.external)
    assert all(degree == 3 for _, degree in owners.degree())
    planar, embedding = nx.check_planarity(owners)
    assert planar
    embedding.check_structure()
    # Nine vertical edges connect three disjoint rows to three disjoint columns.
    rows = [[(i, j, 0) for j in range(3)] for i in range(3)]
    cols = [[(i, j, 1) for i in range(3)] for j in range(3)]
    vertical = [[rows[i][j], cols[j][i]] for i, j in product(range(3), repeat=2)]
    minor = nx.Graph()
    for a, b in vertical:
        assert sum(abs(x - y) for x, y in zip(a, b, strict=True)) == 1
        minor.add_edge(("r", a[0]), ("c", b[1]))
    assert len(set(sum(rows + cols, []))) == 18 and not nx.check_planarity(minor)[0]
    return dict(
        planar_owner_graph=True,
        row_branches=rows,
        column_branches=cols,
        connecting_edges=vertical,
    )
