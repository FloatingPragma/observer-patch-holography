"""Registered bridge preparations, reversible ternary readers and raw rail keys."""

from collections import deque
from functools import cache
from math import pi

import networkx as nx
import numpy as np

from model.candidate_laws.interval_echo import candidate
from model.executable_model.carrier import federation
from ..seam_word_meter.apparatus import projectors, pulse
from ..seam_word_meter.chart import exchange_word
from ..seam_word_meter.finite_program.machine import READ_CODES


@cache
def rail_swaps():
    """Certify port transfer probabilities from the full native blank columns."""
    ex, contrast = candidate.native_law.local_term_parts()
    src = candidate.native_law.SRC
    pair = [src.index((1, 0)), src.index((0, 1))]
    h = ex + candidate.KAPPA * contrast
    kept = h[:, pair].copy()
    kept[pair, :] = 0
    np.testing.assert_allclose(kept, 0, rtol=0, atol=3e-13)
    np.testing.assert_allclose(
        h[np.ix_(pair, pair)],
        (candidate.KAPPA + 0.5) * np.eye(2) - np.array([[0, 0.5], [0.5, 0]]),
        rtol=0,
        atol=3e-13,
    )
    times = {"parked": 0, "third": pi / 3, "half": pi / 2, "full": pi}
    numerators = {"parked": 0, "third": 1, "half": 2, "full": 4}
    swap = np.zeros((9, 9))
    for j, (a, b) in enumerate(src):
        swap[src.index((b, a)), j] = 1
    np.testing.assert_allclose(pulse(pi), swap, atol=4e-12, rtol=0)
    for mode, time in times.items():
        unitary = pulse(time)
        for j in (0, 1):
            prob = abs(unitary[:, pair[j]]) ** 2
            expected = np.zeros(9)
            expected[pair[1 - j]] = numerators[mode] / 4
            expected[pair[j]] = 1 - numerators[mode] / 4
            np.testing.assert_allclose(prob, expected, rtol=0, atol=3e-12)
    # Deleting off-diagonal exchange keeps the same contact, diagonal and alphabet.
    no_exchange = np.diag(np.diag(h[np.ix_(pair, pair)]))
    assert no_exchange[1, 0] == 0
    return numerators


def layout(level, size=None):
    board = federation(level)
    n = (6, 9, 12)[level] if size is None else size
    owners = board.n // 12
    auditors = (owners - 2, owners - 1)
    graph = nx.Graph()
    graph.add_nodes_from(range(owners))
    cross = [tuple(sorted((a, b))) for a, b in board.seams if a // 12 != b // 12]
    graph.add_edges_from((a // 12, b // 12) for a, b in cross)
    chosen = []
    queue = deque([0])
    while len(chosen) < n:
        u = queue.popleft()
        if u in chosen or u in auditors:
            continue
        chosen.append(u)
        queue.extend(sorted(graph[u]))
    inside = sorted(e for e in cross if all(v // 12 in chosen for v in e))
    cut = sorted(e for e in cross if sum(v // 12 in chosen for v in e) == 1)
    boundary = tuple(sorted(cut[0], key=lambda v: v // 12 not in chosen))
    outside = boundary[1]
    assert outside // 12 not in auditors
    assert len({v for e in [*inside, boundary] for v in e}) == 2 * (len(inside) + 1)
    owner_edges = [(chosen.index(a // 12), chosen.index(b // 12)) for a, b in inside]
    region = nx.Graph()
    region.add_nodes_from(range(n))
    region.add_edges_from(owner_edges)
    bridges = {frozenset(e) for e in nx.bridges(region)}
    skip = next(
        ([i] for i, e in enumerate(owner_edges) if frozenset(e) in bridges), None
    )
    if skip is None:
        skip = [i for i, e in enumerate(owner_edges) if n - 1 in e]
    return dict(
        level=level,
        owners=chosen,
        auditors=list(auditors),
        inside=inside,
        boundary=boundary,
        skip=skip,
    )


def meter_response(setting="native"):
    effects = projectors()
    cols = {}
    for v in range(3):
        hits = [abs(m[v, v]) ** 2 for m in effects]
        np.testing.assert_allclose(hits, np.eye(3)[v], rtol=0, atol=1e-12)
        normal = READ_CODES[int(np.argmax(hits))]
        cols[v] = normal if setting == "native" else normal & 2
    return cols


def calibrate(config):
    """Both stations route each addressed load to distinct heads and restore all spectators."""
    board = federation(config["level"])
    graph = nx.Graph(board.seams)
    ports = [12 * c + j for c in config["owners"] for j in range(12)] + [
        config["boundary"][1]
    ]
    heads = [12 * c for c in config["auditors"]]
    references = [[h + 1, h + 5] for h in heads]
    reserved = set(heads + [p for rr in references for p in rr])
    stations = []
    for h, refs in zip(heads, references, strict=True):
        assert all(graph.has_edge(h, r) for r in refs)
        sub = graph.subgraph(set(graph) - (reserved - {h}))
        paths = nx.single_source_shortest_path(sub, h)
        reads = []
        for p in ports:
            path = list(reversed(paths[p]))
            schedule = exchange_word(path)
            local = list(range(board.n))
            for a, b in schedule:
                local[a], local[b] = local[b], local[a]
            assert local[h] == p and local[p] == h and all(local[r] == r for r in refs)
            for a, b in reversed(schedule):
                local[a], local[b] = local[b], local[a]
            assert local == list(range(board.n))
            reads.append(dict(port=p, path=path))
        stations.append(
            dict(
                head=h,
                references=refs,
                paths=reads,
                sharp=[meter_response()[v] for v in range(3)],
                shuttered=[meter_response("shutter")[v] for v in range(3)],
            )
        )
    return stations


def raw_word(config, transfers, deposit):
    loads = {12 * c + j: 0 for c in config["owners"] for j in range(12)}
    for j, (a, b) in enumerate(config["inside"]):
        loads[(a, b)[(transfers >> j) & 1]] = 1
    a, b = config["boundary"]
    loads[a], loads[b] = int(deposit), 1 - int(deposit)
    return [loads[p] for p in loads]


def read_words(loads, stations):
    return [[station["sharp"][v] for v in loads] for station in stations]


def decode(raw, n):
    if len(raw) != 2 or raw[0] != raw[1] or len(raw[0]) != 12 * n + 1:
        raise ValueError("Incomplete double read")
    cells = [READ_CODES.index(k) for k in raw[0]]
    w = sum((sum(cells[12 * i : 12 * i + 12]) % 2) << i for i in range(n))
    return w, cells[-1]
