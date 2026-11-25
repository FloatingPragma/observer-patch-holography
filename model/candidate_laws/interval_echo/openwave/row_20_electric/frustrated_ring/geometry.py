"""Select closed physical routes and sharp source comparators by native adjacency."""

from itertools import permutations, product

import networkx as nx

from model.executable_model.carrier import federation


def attach(board, pairs, marks):
    used = {p for pair in pairs for p in pair}
    choices = [sorted(set(board.neighbours[pairs[i][0]]) - used) for i in marks]
    for pins in product(*choices):
        if len(set(pins)) != len(marks):
            continue
        busy = used.union(pins)
        refs = [sorted(set(board.neighbours[p]) - busy) for p in pins]
        for standards in product(*refs):
            if len(set(standards)) == len(marks):
                return list(pins), list(standards)
    return None


def owner_rings(board):
    owners = nx.Graph()
    gates = {}
    for a, b in board.external:
        i, j = a // 12, b // 12
        gates[i, j], gates[j, i] = a, b
        owners.add_edge(i, j)
    # A canonical five-owner circuit through owner zero.
    loops = set()
    for a in owners[0]:
        for path in nx.all_simple_paths(owners, a, 0, cutoff=4):
            if len(path) == 5:
                interior = tuple(path[:-1])
                loops.add((0,) + min(interior, interior[::-1]))
    return sorted(loops), gates


def layouts():
    board = federation(0)
    graph = nx.Graph(board.seams)
    small = None
    for cell in range(board.cells):
        ports = range(12 * cell, 12 * cell + 12)
        inside = graph.subgraph(ports)
        for a in ports:
            for b in sorted(inside[a]):
                for path in sorted(nx.all_simple_paths(inside, b, a, cutoff=5)):
                    if len(path) != 6:
                        continue
                    route = [a, *path[:-1]]
                    pairs = list(zip(route[::2], route[1::2], strict=True))
                    if any(board.degree[x] != board.degree[y] for x, y in pairs):
                        continue
                    marks = (0, 1, 2)
                    endpoints = attach(board, pairs, marks)
                    if endpoints:
                        small = describe(graph, pairs, marks, endpoints)
                        break
                if small:
                    break
            if small:
                break
        if small:
            break

    loops, gates = owner_rings(board)
    for owners in loops:
        options = []
        for i, owner in enumerate(owners):
            first = gates[owner, owners[i - 1]]
            last = gates[owner, owners[(i + 1) % len(owners)]]
            private = {
                p for p in range(12 * owner, 12 * owner + 12) if board.degree[p] == 5
            }
            threads = nx.all_simple_paths(
                graph.subgraph(private | {first, last}), first, last, cutoff=5
            )
            options.append(sorted(x for x in threads if len(x) in (4, 6)))
        # k=(2,1,1,1,1), in one of five rotations.
        for doubled in range(5):
            sets = [
                [x for x in opt if len(x) == (6 if i == doubled else 4)]
                for i, opt in enumerate(options)
            ]
            for pieces in product(*sets):
                ports = [p for piece in pieces for p in piece]
                for offset in range(1, len(ports), 2):
                    ordered = ports[offset:] + ports[:offset]
                    pairs = list(zip(ordered[::2], ordered[1::2], strict=True))
                    marks = (0, 3, 6, 2)
                    endpoints = attach(board, pairs, marks)
                    if endpoints and small:
                        return small, describe(graph, pairs, marks, endpoints)
    raise ValueError("No admissible pentagonal-owner route with eleven native dimers")


def describe(graph, pairs, marks, endpoints):
    pins, standards = endpoints
    return dict(
        pairs=pairs,
        marks=marks,
        pins=pins,
        standards=standards,
        routes={
            f"{a}{b}": min(nx.all_shortest_paths(graph, pins[a], pins[b]))
            for a, b in permutations(range(len(pins)), 2)
            if a < b
        },
    )
