from collections import defaultdict
from fractions import Fraction as Q
from itertools import product

import networkx as nx
import numpy as np

from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture import (
    collective_geometry as geometry,
    collective_source as source,
    collective_switch as switch,
)


def test_pooled_read_uses_native_columns_and_the_joint_material_count():
    words, weight, _ = source.records()
    laws = []
    for name in ("original", "source"):
        law = defaultdict(Q)
        for word, p in zip(words, weight[name], strict=True):
            law[
                int(word[-1]), int(np.sum(word[:4] == 1)), int(np.sum(word[:4] == 2))
            ] += p
        laws.append(law)
    gate = source.tidal.fastest_relay().reshape(2, 11, 3, 3)
    total = Q(0)
    for (x, px), (y, py) in product(laws[0].items(), laws[1].items()):
        _, a, b = x
        m, c, d = y
        n1, n2 = a + c + 1, b + d
        pop = (240 - n1 - n2, n1, n2)
        heads = source.tidal.lab.stationary_heads(pop)[m]
        for h in heads:
            for load, count in enumerate(pop):
                total += (
                    px
                    * py
                    * count
                    * int(gate[0, h, 0, load])
                    / Q(32 * 240 * len(heads))
                )
    assert total == source.chances(240, ["original", "source"], 1)[1]
    a, b, c = source.ordinary_coefficients()
    for n, m in ((240, 1), (960, 8), (3840, 2)):
        substitute = a * n + b * m + c
        assert source.chances(n, ["source"] * m) == source.chances(
            n, ["original"] * m, substitute
        )


def test_weighted_release_matches_independent_rational_graph_on_full_launch_stars():
    board = source.tidal.lab.energy.pair.stock.observations.federation(0)
    s = geometry.survey(board, source.packing(board, "lexicographic", 2), False)
    cost = [1 / x for x in source.chances(board.n, ["source", "original"])]
    result = geometry.measure(s, cost)
    g = nx.Graph()
    for a, b in zip(*s.blank.nonzero(), strict=True):
        g.add_edge(int(a), int(b), cost=Q(1))
    for (a, b), price in zip(s.contacts, cost, strict=True):
        g.add_edge(a, b, cost=price)
    for j in range(0, len(s.pairs), max(1, len(s.pairs) // 12)):
        u, v = s.pairs[j]
        U = nx.single_source_dijkstra_path_length(g, u, weight="cost")
        V = nx.single_source_dijkstra_path_length(g, v, weight="cost")
        D = U[v]
        cross = sum(V[i] ** 2 for i in board.neighbours[u]) + sum(
            U[i] ** 2 for i in board.neighbours[v]
        )
        deg = len(board.neighbours[u]) + len(board.neighbours[v])
        assert (
            result.length[j] == D
            and result.curvature[j] == (cross - deg * (1 + D * D)) / D**2
        )
    assert geometry.nearest_sources(s, cost) == [10, 10]


def test_collective_switch_detects_a_pure_threshold_and_retains_ties():
    # Cyclic route with two tolls against five blank hops. A single
    # change exactly ties. Equality must be allowed on that side only.
    assert switch.possible_detours(2, Q(9, 4), Q(9, 4), Q(11, 4)) == {2: [5]}
    assert switch.possible_detours(2, Q(5, 2), Q(5, 2), Q(11, 4)) == {}
    assert switch.possible_detours(
        4, Q(13, 4) - Q(1, 100), Q(13, 4), Q(13, 4) + Q(1, 100)
    )[4] == [13]
