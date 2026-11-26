"""Independent short-time polynomial flight of two planted sharp material sources."""

from collections import defaultdict
import math

import networkx as nx

from model.executable_model import federation


def apply(board, state):
    """Act on the complete inventory sector, generating only supported words."""
    result = defaultdict(complex)
    for (plus, minus), amplitude in state.items():
        charges = {**dict.fromkeys(plus, 1), **dict.fromkeys(minus, -1)}
        edges = {tuple(sorted((a, b))) for a in charges for b in board.neighbours[a]}
        potential = 8 * sum(
            (charges.get(a, 0) - charges.get(b, 0)) ** 2 for a, b in edges
        )
        potential += 8 * sum(
            (6 - int(board.degree[a])) * ((1 + q) ** 2 - 1) for a, q in charges.items()
        )
        diagonal = potential - 220  # centre of the proven [128,312] energy window
        for a, b in edges:
            if charges.get(a, 0) == charges.get(b, 0):
                continue
            diagonal += 0.5

            def swapped(pair):
                return tuple(sorted(b if i == a else a if i == b else i for i in pair))

            result[swapped(plus), swapped(minus)] -= amplitude / 2
        result[plus, minus] += amplitude * diagonal
    return result


def release(plus, minus, denominator=4096):
    board = federation(0)
    owners = nx.Graph((a // 12, b // 12) for a, b in board.external)
    distance = nx.single_source_shortest_path_length(owners, 0)
    t = 1 / denominator
    word = (tuple(plus), tuple(minus))
    term = {word: 1}
    approximation = defaultdict(complex, term)
    for order in range(1, 5):
        term = {w: -1j * t * z / order for w, z in apply(board, term).items()}
        for w, z in term.items():
            approximation[w] += z
    moment = 0.0
    outer_probability = 0.0
    for (p, m), value in approximation.items():
        if all(i < 12 for i in (*p, *m)):
            continue
        squared = abs(value) ** 2
        outer_probability += squared
        moment += squared * (
            sum(distance[i // 12] ** 2 for i in p)
            - sum(distance[i // 12] ** 2 for i in m)
        )
    norm_error = (92 * t) ** 5 / math.factorial(5)
    absolute_error = (
        2
        * max(distance.values()) ** 2
        * (2 * math.sqrt(outer_probability) * norm_error + norm_error**2)
    )
    return dict(
        plus=plus,
        minus=minus,
        denominator=denominator,
        degree_four_words=len(approximation),
        signed_radius_per_squared_time=moment / t**2,
        truncation_radius_per_squared_time=absolute_error / t**2,
    )


def standards():
    # Material ports 1 and 3 are gateways. Ports 0 and 2 are ordinary.
    return [release([1, 3], [0, 2]), release([0, 2], [1, 3])]
