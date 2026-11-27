"""Close the unequal-degree strand using only deposited native seams."""

from collections import Counter
from hashlib import sha256

import numpy as np

from model.candidate_laws.interval_echo.openwave.row_30_collective.formation import (
    edge,
    path_to,
)
from model.executable_model import federation


def first_circuit(board, size):
    bridge = {
        (a // 12, b // 12): a for pair in board.external for a, b in (pair, pair[::-1])
    }
    adjacent = {i: sorted(b for a, b in bridge if a == i) for i in range(board.cells)}

    def material(owners):
        strands = []
        for i, owner in enumerate(owners):
            first = bridge[owner, owners[i - 1]]
            last = bridge[owner, owners[(i + 1) % size]]
            inners = sorted(
                (a, b)
                for a in board.neighbours[first]
                for b in board.neighbours[last]
                if board.degree[a] == board.degree[b] == 5 and b in board.neighbours[a]
            )
            if not inners:
                return None
            a, b = inners[0]
            strands.extend([(first, a), (b, last)])
        return strands

    def extend(path):
        if len(path) == size:
            if path[0] in adjacent[path[-1]] and path[1] < path[-1]:
                return material(path)
            return None
        for owner in adjacent[path[-1]]:
            if owner > path[0] and owner not in path:
                found = extend([*path, owner])
                if found:
                    return found
        return None

    for start in adjacent:
        found = extend([start])
        if found:
            return found
    return None


def commission(board, ring):
    """Write connector, idle and blank collars using the conserved sharp stock."""
    n = len(ring)
    seams = {edge(a, b) for a, b in board.seams}
    blanks = {edge(a, b) for a, b in ring}
    links = {edge(ring[i - 1][1], ring[i][0]) for i in range(n)}
    sites = {p for dimer in ring for p in dimer}
    assert len(sites) == 2 * n and len(blanks) == len(links) == n
    assert not links & blanks and (blanks | links) <= seams
    vacancies = {p for p in sites if board.degree[p] == 5}
    assert len(vacancies) == n and all(
        sum(p in vacancies for p in pair) == 1 for pair in ring
    )
    assert all((a in vacancies) == (b in vacancies) for a, b in links)
    outside = sorted(set(range(board.n)) - sites)
    paths = [[] for _ in range(board.n)]
    for a, b in sorted(seams - links):
        paths[a].append(b)
        paths[b].append(a)
    stock = np.ones(board.n, dtype=int)
    stock[outside[:n]] = 0
    for a, b in sorted(links):
        stock[a], stock[b] = 2, 0
    digest = sha256()
    steps = Counter()

    def record(name, a, b):
        digest.update(f"{name} {a} {b}\n".encode())
        steps[name] += 1

    def swap(a, b):
        assert b in paths[a]
        stock[a], stock[b] = stock[b], stock[a]
        record("swap", a, b)

    for a, b in sorted(links):
        assert (stock[a], stock[b]) == (2, 0)
        route = path_to(b, paths, lambda p: board.degree[p] == 5)
        chain = list(zip(route, route[1:]))
        record("half-contact", a, b)
        for x, y in chain:
            swap(x, y)
        record("degree-phase", route[-1], route[-1])
        for x, y in reversed(chain):
            swap(x, y)
        record("write-conditional", a, b)
        stock[a], stock[b] = 1, 1
    for dest in sorted(vacancies):
        if stock[dest]:
            route = path_to(dest, paths, lambda p: p not in vacancies and stock[p] == 0)
            chain = list(zip(route, route[1:]))
            for a, b in chain + chain[-2::-1]:
                swap(a, b)
    mates = {p: q for x, y in ring for p, q in ((x, y), (y, x))}
    for a, b in sorted(seams - links - blanks):
        if a not in sites and b not in sites:
            continue
        for p in (a, b):
            if stock[p] == 0:
                swap(p, mates[p])
        assert (stock[a], stock[b]) == (1, 1)
        record("write-idle", a, b)
    for p in sorted(vacancies):
        if stock[p]:
            swap(p, mates[p])
    assert set(np.flatnonzero(stock == 0)) == vacancies and np.all(stock[outside] == 1)
    return dict(
        ring=ring,
        connectors=sorted(links),
        owners=[a // 12 for a, b in ring],
        private=sorted(vacancies),
        exterior_unit=next(i for i in outside if board.degree[i] == 5),
        degree_signs=[1 if a in vacancies else -1 for a, b in ring],
        operations=dict(steps),
        address_digest=digest.hexdigest(),
        total_load=int(stock.sum()),
    )


def charts():
    cases = []
    for level in range(3):
        board = federation(level)
        for size in (5, 6):
            ring = first_circuit(board, size)
            cases.append(
                dict(
                    cells=board.cells,
                    requested_owners=size,
                    stock=commission(board, ring) if ring else None,
                )
            )
    return cases
