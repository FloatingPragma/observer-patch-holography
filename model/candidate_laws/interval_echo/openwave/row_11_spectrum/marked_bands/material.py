"""Material preparation and bounds for a charged insertion into a filled federation."""

import json
from functools import cache
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from model.experiments.interval_echo_angular.charged_orientation.live_trimer.measure import (
    blocks,
)
from ..conjugate_flight.measure import seam_data

HERE = Path(__file__).resolve().parent
LAW = candidate.native_law


def formation_certificate():
    """Compare each seam with the SAME seam memory in unit-filled material."""
    filled = LAW.SRC.index((1, 1))
    spectrum = {}
    for sector, h in blocks().items():
        heads = len(h) // 9
        host = np.kron(h[filled::9, filled::9], np.eye(9, dtype=int))
        diff = h - host
        defects = np.tile([int(a != 1) + int(b != 1) for a, b in LAW.SRC], heads)
        assert np.array_equal(h @ host, host @ h)
        assert not np.any(diff[:, defects == 0])
        assert not np.any(h * (defects[None, :] != defects[:, None]))
        # All roots and the absence of any other roots are checked over integers.
        roots = (
            {
                "blank": (0, 8, 9, 32, 33),
                (0, 1): (0,),
                (1, -1): (-9, 0),
                (1, 1): (-8, 0),
                (2, -1): (0, 33),
                (2, 1): (0, 32),
            }
        )[sector]
        ident = np.eye(len(h), dtype=object)
        rem = ident
        for root in roots:
            rem = rem @ (diff.astype(object) - 8 * root * ident)
        assert not rem.any()
        spectrum[str(sector)] = sorted(
            set(np.rint(np.linalg.eigvalsh(diff) / 8).astype(int).tolist())
        )
    return dict(
        native_columns=99,
        matched_host_commutes=True,
        relative_seam_levels=spectrum,
        physical_degree_limit=6,
        absolute_energy_per_defect_ceiling=198,
    )


@cache
def geometry():
    route_file = HERE.parents[1] / "electron_field/lepton_family/bands/route.json"
    ports = json.loads(route_file.read_text())["ports"]
    board = federation(2)
    ring = set(ports)
    assert len(ring) == len(ports) == 144
    assert list(board.degree[ports]) == [6, 5, 6] * 48
    assert all(b in board.neighbours[a] for a, b in zip(ports, ports[1:] + ports[:1]))
    marks = [
        (p, min(q for q in board.neighbours[p] if q not in ring and q // 12 == p // 12))
        for p in ports[1::3]
    ]
    assert len({q for _, q in marks}) == 48
    return ports, marks


def enclosure(sign):
    ports, marks = geometry()
    board = federation(2)
    graph = nx.Graph(board.seams)
    graph.remove_edges_from(marks)
    distance = nx.multi_source_dijkstra_path_length(graph, ports)
    a, b = next((a, b) for a, b in board.seams if min(distance[a], distance[b]) > 3)
    edges = {tuple(sorted(pair)) for pair in zip(ports, ports[1:] + ports[:1])}
    marked = {tuple(sorted(e)) for e in marks}
    pair = (a, b)
    domain = set(ports) | {a, b}
    cuts = sorted(
        e
        for e in board.seams
        if domain.intersection(e) and e not in edges and e not in marked and e != pair
    )
    word = np.ones(board.n, dtype=int)
    word[ports[0]], word[a] = 0, 2

    def swap(i, j):
        assert graph.has_edge(i, j)
        word[i], word[j] = word[j], word[i]

    # One stock vacancy visits every middle port to write a private (0,0) memory.
    for i in range(len(ports)):
        assert word[ports[i]] == 0
        swap(ports[i], ports[(i + 1) % len(ports)])
    if sign == 1:
        path = nx.shortest_path(graph, ports[0], a)
        bonds = list(zip(path, path[1:]))
        for i, j in bonds + bonds[-2::-1]:
            swap(i, j)
    for i, j in cuts:
        relay = None
        for src, dst in [(ports[0], ports[1]), (a, b)]:
            if src in (i, j):
                relay = (src, dst)
                swap(src, dst)
        assert word[i] == word[j] == 1
        graph.remove_edge(i, j)
        if relay:
            swap(*relay)
    assert word[ports[0]] == 1 + sign and word[a] == 1 - sign
    assert set(nx.node_connected_component(graph, a)) == {a, b}
    assert set(nx.node_connected_component(graph, ports[0])) == set(ports)
    return dict(remote_store=[a, b], silent_seams=len(cuts), private_zero_writes=48)


def programmed_mark(kappa, orient_at=8):
    """A native (0,2)+(2,0) echo prepares the +2/-2, even-coin side memory."""
    blank, _ = seam_data(kappa)
    pair = [LAW.SRC.index(w) for w in ((0, 2), (2, 0))]
    half = expm(-1j * np.pi / 2 * blank[np.ix_(pair, pair)])[:, 0]
    # Private contrast with fixed zero reference on the first endpoint.
    half[1] *= np.exp(-4j * kappa * np.pi / (8 * orient_at))
    out = (LAW.echo()[:, pair] @ half).reshape(11, 9)
    filled = LAW.SRC.index((1, 1))
    head = out[:, filled]
    assert np.linalg.norm(head) > 1 - 1e-13
    code = np.kron(head[:, None], np.eye(9))
    exchange, contrast = LAW.local_term_parts()
    term = exchange + kappa * contrast
    acted = term @ code
    energies = (code.conj().T @ acted).diagonal().real
    np.testing.assert_allclose(acted, code * energies, atol=2e-13)
    np.testing.assert_allclose(
        energies, [0 if w == (1, 1) else 4 * kappa for w in LAW.SRC], atol=2e-13
    )
    return energies[LAW.SRC.index((0, 1))] - energies[filled]


def route_matrix(kappa, sign):
    ports, marks = geometry()
    n = len(ports)
    h = np.zeros((n, n))
    pair = [LAW.SRC.index(w) for w in ((1 + sign, 1), (1, 1 + sign))]
    term = seam_data(kappa)[0][np.ix_(pair, pair)]
    for i in range(n):
        ix = [i, (i + 1) % n]
        h[np.ix_(ix, ix)] += term
    shift = programmed_mark(kappa)
    for port, _ in marks:
        h[ports.index(port), ports.index(port)] += shift
    return h


def synthesize(target, actual_k=8, source_k=8):
    """Emit a one-charge state by scheduled adjacent exchanges and marked-port phases."""
    state = np.zeros(len(target), complex)
    state[0] = 1
    n = len(state)

    def phase(j, angle):
        # Address the nearest marked middle port by a swap and its inverse.
        mid = (j // 3) * 3 + 1
        if mid != j:
            state[[j, mid]] = state[[mid, j]]
        hold = (-angle % (2 * np.pi)) / (4 * source_k)
        state[mid] *= np.exp(-4j * actual_k * hold)
        if mid != j:
            state[[j, mid]] = state[[mid, j]]

    tail = np.sqrt(np.cumsum(abs(target[::-1]) ** 2))[::-1]
    for j in range(n - 1):
        time = 2 * np.arctan2(tail[j + 1], abs(target[j]))
        c, s = np.cos(time / 2), 1j * np.sin(time / 2)
        state[j : j + 2] = (
            np.exp(-1j * (actual_k + 0.5) * time)
            * np.array([[c, s], [s, c]])
            @ state[j : j + 2]
        )
        phase(j, (source_k + 0.5) * time)
        phase(j + 1, (source_k + 0.5) * time - np.pi / 2)
    for j, z in enumerate(target):
        phase(j, np.angle(z))
    return state
