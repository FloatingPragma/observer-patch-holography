"""Local splitters, routed stocks, two-ended writes and binary load meters."""

from functools import cache
from itertools import combinations, product
from math import pi

import networkx as nx
import numpy as np

from model.executable_model import federation
from ....seam_word_meter.apparatus import pulse
from ....seam_word_meter.chart import exchange_word
from ...repair_outputs.apparatus import INPUTS, law, writers


@cache
def local_calibration():
    writer = law.echo()
    balanced = np.arange(law.DIM) % 9 == INPUTS.index((1, 1))
    meter = []
    for val in range(3):
        col = np.eye(law.DIM)[:, INPUTS.index((val, 1))]
        yes = writer.T.conj() @ (balanced * (writer @ col))
        bit = int(round(np.linalg.norm(yes) ** 2))
        np.testing.assert_allclose(yes, bit * col, atol=1e-12, rtol=0)
        meter.append(bit)
    if meter != [0, 1, 0]:
        raise ValueError("One-load standard cannot read the middle load")
    flip = pulse(pi)
    for i, pair in enumerate(INPUTS):
        target = np.eye(9)[:, INPUTS.index(pair[::-1])]
        np.testing.assert_allclose(flip[:, i], target, atol=4e-12, rtol=0)
    splitters = {}
    for duration, label in ((pi / 2, "half"), (pi / 3, "third")):
        col = pulse(duration)[:, INPUTS.index((1, 0))]
        chances = [
            sum(abs(a) ** 2 for a, (_, dest) in zip(col, INPUTS) if dest == n)
            for n in range(3)
        ]
        numerator = round(4 * chances[1])
        np.testing.assert_allclose(
            chances, [1 - numerator / 4, numerator / 4, 0], atol=1e-12, rtol=0
        )
        splitters[label] = numerator
    tables = {}
    for mode, columns in writers().items():
        table = np.empty((2, 2, 2), dtype=int)
        for a, b in product(range(2), repeat=2):
            col = columns[:, INPUTS.index((2 * a, 2 * b))].reshape(-1, 3, 3)
            chance = (abs(col) ** 2).sum(axis=0)
            ports = np.unravel_index(np.argmax(chance), chance.shape)
            np.testing.assert_allclose(chance[ports], 1, atol=1e-12, rtol=0)
            table[a, b] = [meter[v] for v in ports]
        tables[mode] = table
    return splitters, tables


def cycle(graph):
    paths = []
    inner = graph.subgraph(set(graph) - {0})
    for left, right in combinations(sorted(graph[0]), 2):
        for path in nx.all_shortest_paths(inner, left, right):
            paths.extend([(0, *path), (0, *reversed(path))])
    return min(paths, key=lambda p: (len(p), p))


def route(graph, a, b, avoid=()):
    sub = graph.subgraph(set(graph) - set(avoid))
    path = min(nx.all_shortest_paths(sub, a, b))
    steps = exchange_word(path)
    labels = {p: p for p in graph}
    for x, y in steps:
        labels[x], labels[y] = labels[y], labels[x]
    if any(labels[k] != (b if k == a else a if k == b else k) for k in graph):
        raise ValueError("Route disturbed a spectator")
    return path


def owner_layout(board, owner, ends):
    graph = nx.Graph()
    for p in range(12 * owner, 12 * (owner + 1)):
        graph.add_edges_from((p, q) for q in board.neighbours[p] if q // 12 == owner)
    edges = sorted(
        tuple(sorted(e)) for e in graph.edges if not set(e).intersection(ends)
    )
    coins = next(
        (a, b) for a, b in combinations(edges, 2) if not set(a).intersection(b)
    )
    free = sorted(set(graph) - set(ends) - set(sum(coins, ())))
    spares = free[:2]
    standard = free[2]
    swaps = [route(graph, p, q) for p, q in zip(ends, spares, strict=True)]
    instruments = []
    for p in [coins[0][1], coins[1][1], *ends]:
        q = min(graph[p])
        instruments.append(
            dict(input=p, comparison=q, route=route(graph, standard, q, avoid=(p,)))
        )
    return dict(
        owner=owner,
        ports=ends,
        coins=coins,
        spares=spares,
        standard=standard,
        swaps=swaps,
        readers=instruments,
    )


@cache
def board_layout(level):
    board = federation(level)
    cross = {}
    for p, q in board.external:
        a, b = p // 12, q // 12
        if (a, b) in cross:
            raise ValueError("Nonunique owner contact")
        cross[a, b] = (p, q)
        cross[b, a] = (q, p)
    owners = cycle(nx.Graph(cross.keys()))
    links = [cross[a, b] for a, b in zip(owners, owners[1:] + owners[:1], strict=True)]
    configs = [
        owner_layout(board, a, [links[i - 1][1], links[i][0]])
        for i, a in enumerate(owners)
    ]
    return dict(level=level, owners=list(owners), seams=links, stations=configs)


def load_table(station, shared):
    """Push the sharp data and spare rails through each conditional native swap."""
    rows = []
    for first, second in product((0, 1), repeat=2):
        state = {
            k: 0 for k in range(12 * station["owner"], 12 * (station["owner"] + 1))
        }
        state[station["standard"]] = 1
        for pair, bit in zip(station["coins"], [first, second], strict=True):
            state[pair[0]], state[pair[1]] = 1 - bit, bit
        for p in station["spares"]:
            state[p] = 2
        for bit, path in zip(
            (first, first if shared else second), station["swaps"], strict=True
        ):
            if bit:
                for a, b in exchange_word(path):
                    state[a], state[b] = state[b], state[a]
        rows.append([state[p] // 2 for p in station["ports"]])
    return np.array(rows, dtype=np.int8)
