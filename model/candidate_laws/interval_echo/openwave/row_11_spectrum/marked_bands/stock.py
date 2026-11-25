"""Prepare every static side collar with its own supplied material pair."""

from functools import cache
from itertools import pairwise
from math import pi

import networkx as nx
import numpy as np
from scipy.linalg import expm

from model.executable_model import federation
from model.experiments.material_inertia.dark_collar import label
from . import material, paired
from ..conjugate_flight.measure import closed_walk, edge, seam_data


@cache
def write_column(reference_endpoint):
    """Blank (0,2) in ascending physical port order and a positive private dwell."""
    strength = 8
    law = material.LAW
    blank, private = seam_data(strength)
    ordered = [law.SRC.index(w) for w in ((0, 2), (2, 0))]
    half = expm(-1j * pi * blank[np.ix_(ordered, ordered)] / 2)[:, 0]
    # The second order carries the larger private load only at the first endpoint.
    dwell = (1 if reference_endpoint == 0 else 3) * pi / (8 * strength)
    local_cost = private[[3 * law.SRC[j][reference_endpoint] for j in ordered]]
    balanced = np.exp(-1j * strength * dwell * local_cost) * half
    stored = (law.echo()[:, ordered] @ balanced).reshape(11, 9)
    expected = label(2, 1)
    filled = law.SRC.index((1, 1))
    phase = np.vdot(expected, stored[:, filled])
    target = np.zeros_like(stored)
    target[:, filled] = phase * expected
    np.testing.assert_allclose(stored, target, atol=2e-13)
    np.testing.assert_allclose(abs(phase), 1, atol=2e-13)
    code = np.kron(expected[:, None], np.eye(9))
    exchange, contrast = law.local_term_parts()
    diagonal = np.diag([0.0 if word == (1, 1) else 32.0 for word in law.SRC])
    np.testing.assert_allclose(
        (exchange + strength * contrast) @ code, code @ diagonal, atol=2e-13
    )
    return dwell


def setting(which):
    if which == "ring":
        ports, marks = material.geometry()
        private = ports[1::3]
        return (
            federation(2),
            [ports],
            marks,
            private,
            material.enclosure(-1)["remote_store"],
            True,
        )
    if which == "paired":
        rails = paired.rail.RAILS
        return (
            federation(0),
            list(rails),
            paired.side_contacts(),
            paired.rail.preparation()["silent_private"],
            paired.store()["dimer"],
            False,
        )
    raise ValueError(f"Unknown bench {which}")


def verify(which, sign):
    """Replay native transpositions, deposits, private writes and every exterior seal."""
    if sign not in (-1, 1):
        raise ValueError("The rail carries a single vacancy or excess.")
    board, rails, marks, private, (store, rest), cycle = setting(which)
    graph = nx.Graph(board.seams)
    active = {edge(a, b) for rail in rails for a, b in pairwise(rail)}
    if cycle:
        active.add(edge(rails[0][-1], rails[0][0]))
    physical = set().union(*rails)
    assert len(set(sum(([p, q] for p, q in marks), []))) == 2 * len(marks)
    assert all(p in physical and q not in physical for p, q in marks)
    distances = nx.multi_source_dijkstra_path_length(graph, physical)
    reserved = {store, rest}
    stocks = []
    for a, b in board.seams:
        if distances[a] > 1 and distances[b] > 1 and not reserved.intersection((a, b)):
            stocks.append((a, b))
            reserved.update((a, b))
            if len(stocks) == len(marks):
                break
    assert len(stocks) == len(marks)
    word = np.ones(board.n, dtype=int)
    for a, b in stocks:
        word[a], word[b] = 0, 2
    start = rails[0][0]
    word[start], word[store] = 0, 2
    swaps = 0

    def swap(a, b):
        nonlocal swaps
        assert graph.has_edge(a, b)
        word[a], word[b] = word[b], word[a]
        swaps += 1

    def transpose(a, b):
        if a != b:
            for u, v in closed_walk(nx.shortest_path(graph, a, b)):
                swap(u, v)

    sequence = []
    for (stock0, stock2), (p, q) in zip(stocks, marks, strict=True):
        a, b = edge(p, q)
        assert word[a] == word[b] == 1
        transpose(stock0, a)
        transpose(stock2, b)
        assert (word[a], word[b]) == (0, 2)
        # All private references are blank until the side stock has been spent.
        reference = next(r for r in (p, q) if board.degree[r] == 5)
        ref_end = (a, b).index(reference)
        dwell = write_column(ref_end)
        word[a] = word[b] = 1
        graph.remove_edge(a, b)
        sequence.append(
            dict(
                stock=[stock0, stock2],
                seam=[a, b],
                reference=reference,
                dwell_pi_over_64=round(dwell * 64 / pi),
            )
        )
    assert np.flatnonzero(word != 1).tolist() == sorted([start, store])
    private_swaps = swaps
    for p in private:
        assert board.degree[p] == 5
        transpose(start, p)
        assert word[p] == 0
        # seam_data checks that echo on (0,0) preserves both endpoints and
        # leaves the entire nine-word material action zero on the written line.
        transpose(start, p)
    if sign == 1:
        transpose(start, store)
    assert (word[start], word[store]) == (1 + sign, 1 - sign)
    domain = physical | {store, rest}
    cuts = sorted(
        e
        for e in graph.edges
        if domain.intersection(e)
        and edge(*e) not in active
        and edge(*e) != edge(store, rest)
    )
    for a, b in cuts:
        relay = None
        if start in (a, b):
            relay = (start, rails[0][1])
        if store in (a, b):
            assert relay is None
            relay = (store, rest)
        if relay is not None:
            swap(*relay)
        assert (word[a], word[b]) == (1, 1)
        graph.remove_edge(a, b)
        if relay is not None:
            swap(*relay)
    # Dark contacts vanish. Written side contacts are diagonal for every
    # endpoint load and therefore add no edge to the material transition graph.
    assert all(set(nx.node_connected_component(graph, r[0])) == set(r) for r in rails)
    assert set(nx.node_connected_component(graph, store)) == {store, rest}
    assert all(not any(nx.has_path(graph, q, r[0]) for r in rails) for _, q in marks)
    assert np.flatnonzero(word != 1).tolist() == sorted([start, store])
    assert (word[start], word[store]) == (1 + sign, 1 - sign)
    return dict(
        bench=which,
        sign=sign,
        stock_pairs=len(stocks),
        initial_vacancies=len(stocks) + 1,
        initial_excesses=len(stocks) + 1,
        stored_even_side_marks=sequence,
        private_zero_writes=len(private),
        silent_exterior_contacts=len(cuts),
        separate_dimer=[store, rest],
        unspent_vacancies=1,
        unspent_excesses=1,
        swaps_during_side_deposits=private_swaps,
        total_full_swaps=swaps,
    )


if __name__ == "__main__":
    import json

    print(
        json.dumps(
            [verify(b, s) for b in ("ring", "paired") for s in (-1, 1)], indent=2
        )
    )
