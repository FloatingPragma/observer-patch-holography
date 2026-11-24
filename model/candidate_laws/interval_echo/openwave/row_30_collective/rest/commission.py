"""Prepare a grown closed ring by addressed operations on a full federation."""

from hashlib import sha256

import numpy as np

from model.executable_model import federation
from ..experiment import grow
from ..formation import edge, path_to


def prepare_ring(level, seed):
    board = federation(level)
    rings = [route for closed, route in grow(board, seed) if closed]
    ring = max(rings, key=len)
    n = len(ring)
    seams = {edge(a, b) for a, b in board.seams}
    dimers = {edge(a, b) for a, b in ring}
    links = {edge(ring[i - 1][1], ring[i][0]) for i in range(n)}
    sites = {p for pair in ring for p in pair}
    assert len(sites) == 2 * n and len(links) == n and len(dimers) == n
    assert dimers.isdisjoint(links) and (dimers | links) <= seams
    assert all(board.degree[a] == board.degree[b] for a, b in ring)
    # Park one sharp vacancy per dimer outside the ring during the connector
    # writes. Every connector carries its own two-load input from the stock.
    available = sorted(set(range(board.n)) - sites)
    assert len(available) >= n
    stock = np.ones(board.n, dtype=int)
    stock[available[:n]] = 0
    for a, b in links:
        stock[a], stock[b] = 2, 0
    total = int(stock.sum())
    blank_contacts = seams - links
    neighbors = [[] for _ in stock]
    for a, b in sorted(blank_contacts):
        neighbors[a].append(b)
        neighbors[b].append(a)
    circuit = sha256()
    totals = {}

    def step(name, a, b):
        circuit.update(f"{name} {a} {b}\n".encode("ascii"))
        totals[name] = totals.get(name, 0) + 1

    def swap(a, b):
        assert edge(a, b) in blank_contacts
        stock[a], stock[b] = stock[b], stock[a]
        step("swap", a, b)

    for a, b in sorted(links):
        assert (stock[a], stock[b]) == (2, 0)
        route = path_to(b, neighbors, lambda i: board.degree[i] == 5)
        outward = list(zip(route, route[1:]))
        step("half-contact", a, b)
        for x, y in outward:
            step("swap", x, y)
        step("private-phase", route[-1], route[-1])
        for x, y in reversed(outward):
            step("swap", x, y)
        step("echo-conditional", a, b)
        stock[a], stock[b] = 1, 1
    wanted = {a for a, b in ring}
    for port in sorted(wanted):
        if stock[port] == 0:
            continue
        route = path_to(port, neighbors, lambda i: i not in wanted and stock[i] == 0)
        edges = list(zip(route, route[1:]))
        for a, b in edges + list(reversed(edges[:-1])):
            swap(a, b)
    mates = {p: q for pair in ring for p, q in (pair, pair[::-1])}
    incident = {e for e in seams if sites.intersection(e)}
    for a, b in sorted(incident - dimers - links):
        for p in (a, b):
            if stock[p] == 0:
                swap(p, mates[p])
        assert (stock[a], stock[b]) == (1, 1)
        step("echo-idle", a, b)
    for a, b in ring:
        if stock[a] == 1:
            swap(a, b)
    assert set(np.flatnonzero(stock == 0)) == wanted and stock.sum() == total
    assert all(stock[p] == 1 for p in available)
    return ring, dict(
        level=level,
        seed=seed,
        board_cells=board.cells,
        ring=[list(pair) for pair in ring],
        ring_cells=len({i // 12 for i in sites}),
        rings=len(rings),
        load_total=total,
        services=totals,
        circuit_sha256=circuit.hexdigest(),
    )
