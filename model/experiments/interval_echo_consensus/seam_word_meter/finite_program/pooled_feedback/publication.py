"""Two restored material reads of both finite terminal scores."""

from functools import cache
from math import pi

import networkx as nx
import numpy as np

from ... import apparatus as app, chart

WIDTH = 17


@cache
def layout():
    board = app.federation(0)
    pairs = [(12*owner + a, 12*owner + a + 1) for owner in range(7, 13) for a in range(0, 12, 2)][:2*WIDTH]
    if len(set(v for pair in pairs for v in pair)) != 4*WIDTH or any(b not in board.neighbours[a] for a, b in pairs):
        raise ValueError('missing disjoint score storage')
    graph = nx.Graph(board.seams)
    stations = app.layout()['stations']
    graph.remove_nodes_from([v for s in stations for v in s['standards']])
    paths = []
    for station in stations:
        target = station['head']
        all_paths = nx.single_source_shortest_path(graph, target)
        station_paths = []
        for pair in pairs:
            for port in pair:
                path = all_paths[port][::-1]
                labels = np.arange(board.n)
                for a, b in chart.exchange_word(path):
                    labels[a], labels[b] = labels[b], labels[a]
                want = np.arange(board.n)
                want[port], want[target] = want[target], want[port]
                np.testing.assert_equal(labels, want)
                station_paths.append(path)
        paths.append(station_paths)
    return dict(width=WIDTH, pairs=pairs, owners=[1, 3], paths=paths)


@cache
def bit_read(bit):
    state = np.eye(9)[:, app.LAW.SRC.index((0, 2))]
    if bit:
        state = app.pulse(pi) @ state
    word = []
    for left in (True, False):
        probabilities = [np.linalg.norm((np.kron(p, np.eye(3)) if left else np.kron(np.eye(3), p)) @ state)**2
                         for p in app.projectors()]
        sharp = int(np.argmax(probabilities))
        np.testing.assert_allclose(probabilities, np.eye(3)[sharp], atol=app.TOL)
        word.append(('01', '00', '10')[sharp])
    return word


def write_scores(sign, records):
    codes = ([sum(row['terminal_counts']) + row['terminal_counts'][2] - row['terminal_counts'][0]
              for row in records[-2:]] if sign else [0, 0])
    if len(codes) != 2 or any(not 0 <= v < 1 << WIDTH for v in codes):
        raise ValueError('unrepresentable terminal publication')
    flags = [[bit_read((value >> k) & 1) for k in range(WIDTH)] for value in codes]
    return [flags for _ in layout()['owners']]


def read_scores(flags, n, radius):
    def decode(word):
        bits = {('01', '10'): 0, ('10', '01'): 1}
        return sum(bits[tuple(pair)] << k for k, pair in enumerate(word))
    values = [(decode(word) - n) / n for word in flags]
    return [max(-1., max(values) - radius), min(1., min(values) + radius)]
