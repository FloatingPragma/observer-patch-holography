"""Transfer the native vacancy columns over an unfiltered geometric cycle basis."""

import argparse
from collections import Counter
import hashlib
import json
from itertools import pairwise

import networkx as nx
import numpy as np
from scipy.linalg import eigh_tridiagonal

from model.candidate_laws.interval_echo.openwave.quark_winding.belt_walls import belt
from model.executable_model import federation
from .measure import HERE, SHOTS, source


PROTOCOL = "2c3a7c250a72d38f4c258c8b27d0bb3be458099787b778b05bc72a7d6e4bb587"


def cycles(graph):
    tree = nx.bfs_tree(graph, 0, sort_neighbors=sorted)
    paths = nx.single_source_shortest_path(tree, 0)
    count = Counter()
    representatives = {}
    for a, b in sorted(graph.edges()):
        u, v = paths[a], paths[b]
        common = next(
            (i for i, (x, y) in enumerate(zip(u, v)) if x != y), min(len(u), len(v))
        )
        cycle = u[: common - 1 : -1] + v[common - 1 :]
        if len(cycle) < 4 or len(cycle) % 2:
            continue
        assert len(set(cycle)) == len(cycle)
        index = cycle.index(min(cycle))
        ordered = cycle[index:] + cycle[:index]
        key = min(ordered, [ordered[0], *ordered[:0:-1]])
        n = len(key)
        count[n] += 1
        if n not in representatives or key < representatives[n]:
            representatives[n] = key
    return count, representatives


def qualify(board, graph, route):
    n = len(route)
    assert all(graph.has_edge(v, route[(i + 1) % n]) for i, v in enumerate(route))
    distances = nx.multi_source_dijkstra_path_length(graph, route)
    donor = max(graph, key=lambda v: (distances[v], -v))
    allowed = graph.copy()
    allowed.remove_node(donor)
    marks = [(route[i], route[(i + 1) % n]) for i in range(1, n, 2)]
    allowed.remove_edges_from(marks)
    refs = sorted(
        p for p in graph if p not in route and p != donor and board.degree[p] == 5
    )
    ok = distances[donor] > 2 and bool(refs) and nx.is_biconnected(allowed)
    return dict(
        size=n,
        compensator=donor,
        minimum_distance=distances[donor],
        off_route_phase=refs[0] if refs else None,
        qualified=ok,
    )


def microscopic(route, kappa=8, charge=-1):
    if charge not in (-1, 1):
        raise ValueError("Choose a single vacancy or excess.")
    n = len(route)
    heads = 0
    words = [(route[0], heads)]
    for x in range(n + 1):
        i = x % n
        if i % 2:
            heads ^= 1 << (i // 2)
        words.append((route[(x + 1) % n], heads))
    index = {word: i for i, word in enumerate(words)}
    assert len(words) == len(index)
    h = np.zeros((len(words), len(words)), dtype=int)
    tables = source.measure.blocks(kappa)
    loads = source.measure.LAW.SRC
    for i, (a, b) in enumerate(pairwise([*route, route[0]])):
        marked = i % 2
        term = tables[(1, 1) if marked else "blank"]
        for col, (hole, record) in enumerate(words):
            coin = (((record >> (i // 2)) & 1) ^ int(charge > 0)) if marked else 0
            start = 9 * coin + loads.index(
                (1 + charge * int(hole == a), 1 + charge * int(hole == b))
            )
            filled = 9 * coin + loads.index((1, 1))
            h[col, col] -= term[filled, filled]
            for end in np.flatnonzero(term[:, start]):
                new = loads[end % 9]
                position = (
                    a if new[0] == 1 + charge else b if new[1] == 1 + charge else hole
                )
                history = (
                    record ^ ((int(end // 9) ^ coin) << (i // 2)) if marked else record
                )
                h[index[position, history], col] += term[end, start]
    return h, np.array([route.index(p) for p, _ in words])


def physical(psi, positions, n):
    p = np.array(
        [np.bincount(positions, weights=abs(row) ** 2, minlength=n) for row in psi]
    )
    np.testing.assert_allclose(p.sum(axis=1), 1, atol=1e-12)
    return p


def compare(route, rng):
    n = len(route)
    h, positions = microscopic(route)
    word, _ = source.lifted("01" * (n // 2))
    diag, hop = source.coefficients(word, "even")
    model = np.diag(diag) + np.diag(hop, 1) + np.diag(hop, -1)
    mismatch = int(np.max(abs(h - model)))
    assert mismatch == 0
    native, a = np.linalg.eigh(h / 8)
    predicted, b = eigh_tridiagonal(np.array(diag) / 8, np.array(hop) / 8)
    t = [1, 4, 12]
    first = 2 * (n // 4)
    observed = physical(
        (np.exp(-1j * np.outer(t, native)) * a[first]) @ a.T, positions, n
    )
    forecast = physical(
        (np.exp(-1j * np.outer(t, predicted)) * b[first]) @ b.T, positions, n
    )
    shots = np.array([rng.multinomial(SHOTS, p / p.sum()) for p in observed])
    err = float(np.max(abs(shots / SHOTS - forecast)))
    return dict(
        route=route,
        reachable_words=len(h),
        maximum_column_error=mismatch,
        maximum_port_error=err,
        passed=err < 0.012,
        counts=shots.tolist(),
    )


def run():
    assert hashlib.sha256((HERE / "GEOMETRY.md").read_bytes()).hexdigest() == PROTOCOL
    rng = np.random.default_rng(304144)
    census, flights, belts = [], [], []
    for level in range(3):
        board = federation(level)
        graph = nx.Graph(board.seams)
        count, representatives = cycles(graph)
        results = {
            n: qualify(board, graph, route)
            for n, route in sorted(representatives.items())
        }
        lengths = [n for n, rec in results.items() if rec["qualified"]]
        for n in sorted({min(lengths), max(lengths)}) if lengths else []:
            flights.append(dict(level=level, **compare(representatives[n], rng)))
        census.append(
            dict(
                level=level,
                total_even=sum(count.values()),
                counts_by_length=dict(sorted(count.items())),
                stock_by_length=results,
            )
        )
    for poles, route in belt.belts().items():
        checks = qualify(board, graph, list(route))
        flight = compare(list(route), rng) if checks["qualified"] else {}
        belts.append(dict(poles=poles, **checks, **flight))
    return (
        json.dumps(
            dict(protocol_sha256=PROTOCOL, census=census, flights=flights, belts=belts),
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    out = run()
    path = HERE / "geometry.json"
    if parser.parse_args().write:
        path.write_text(out)
    else:
        assert path.read_text() == out, "Geometry replay mismatch"
    print("Verified winding cycle survey")
