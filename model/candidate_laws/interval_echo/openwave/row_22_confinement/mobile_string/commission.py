"""Addressed, number-preserving preparation on the three actual graph deposits."""

from collections import Counter
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from model.executable_model import federation
from ...row_30_collective.formation import edge, native, path_to
from .material import route


def material_stock(level, n):
    board = federation(level)
    sequence = route(level)[: 2 * n]
    dimers = list(zip(sequence[::2], sequence[1::2]))
    links = {edge(a, b) for a, b in zip(sequence[1:-1:2], sequence[2::2])}
    seams = {edge(a, b) for a, b in board.seams}
    assert len(set(sequence)) == 2 * n
    assert len(links) == n - 1 and links <= seams
    mates = {a: b for pair in dimers for a, b in (pair, pair[::-1])}
    blanks = {edge(a, b) for a, b in dimers}
    assert len(blanks) == n and blanks <= seams and blanks.isdisjoint(links)
    others = set(range(board.n)) - set(sequence)
    paths = [[] for _ in range(board.n)]
    for a, b in sorted(seams - links):
        paths[a].append(b)
        paths[b].append(a)
    # Initially a sorted two-load stock per written connector, and
    # one sharp vacancy outside the strand per future one-load dimer.
    stock = np.ones(board.n, dtype=int)
    stock[sorted(others)[:n]] = 0
    for a, b in links:
        stock[a], stock[b] = 2, 0
    conserved = int(stock.sum())
    digest = sha256()
    ops = Counter()

    def record(name, a, b):
        ops[name] += 1
        digest.update(f"{name} {a} {b}\n".encode())

    def swap(a, b):
        assert b in paths[a]
        stock[a], stock[b] = stock[b], stock[a]
        record("swap", a, b)

    for a, b in sorted(links):
        assert (stock[a], stock[b]) == (2, 0)
        outward = path_to(b, paths, lambda i: board.degree[i] == 5)
        chain = list(zip(outward, outward[1:]))
        record("half-contact", a, b)
        for x, y in chain:
            record("swap", x, y)
        record("degree-phase", outward[-1], outward[-1])
        for x, y in reversed(chain):
            record("swap", x, y)
        record("echo-conditional", a, b)
        stock[a], stock[b] = 1, 1
    # One hole on the degree-five end of every dimer gives lower everywhere.
    holes = {p for p in sequence if board.degree[p] == 5}
    assert len(holes) == n
    for dest in sorted(holes):
        if stock[dest]:
            walk = path_to(dest, paths, lambda i: i not in holes and stock[i] == 0)
            edges = list(zip(walk, walk[1:]))
            for a, b in edges + edges[-2::-1]:
                swap(a, b)
    assert set(np.flatnonzero(stock == 0)) == holes
    for a, b in sorted(seams - links - blanks):
        if a not in mates and b not in mates:
            continue
        for site in (a, b):
            if stock[site] == 0:
                swap(site, mates[site])
        assert (stock[a], stock[b]) == (1, 1)
        record("echo-idle", a, b)
    for dest in sorted(holes):
        if stock[dest]:
            swap(dest, mates[dest])
    assert stock.sum() == conserved and set(np.flatnonzero(stock == 0)) == holes
    return dict(
        level=level,
        dimers=dimers,
        board_ports=board.n,
        total_load=conserved,
        operations=dict(ops),
        ordered_address_sha256=digest.hexdigest(),
    )


def report():
    return {
        "max_local_compiler_residual": native(),
        "stocks": [material_stock(level, n) for level in range(3) for n in (16, 14)],
    }


if __name__ == "__main__":
    path = Path(__file__).with_name("commission.json")
    path.write_text(json.dumps(report(), indent=1) + "\n")
