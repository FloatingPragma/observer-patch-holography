"""Blank native contacts, restored sign readers and a sharp one-unit message."""

from functools import cache
from itertools import combinations
from math import pi

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from model.experiments.interval_echo_gauge.exterior_register import material
from model.experiments.interval_echo_gauge.matter_class.compiler import compile_unitary

report = material.shared_report
MODES = (0, 5, 11)


@cache
def layout():
    graph, base = report.layout()
    board = federation(0)
    stations = [base["alice"], base["bob"]]
    rails = [[12 * s + p for p in MODES] for s in stations]
    standards = [base["readers"][name] for name in ("alice", "bob")]
    guards = {
        p: v for pairs in standards for (_, p), v in zip(pairs, (0, 2), strict=True)
    }
    route_graph = graph.subgraph(set(graph) - guards.keys())
    sources = []
    for s in range(board.n // 12):
        pair = [12 * s, 12 * s + 5]
        if s not in stations and set(pair).isdisjoint(guards) and graph.has_edge(*pair):
            sources.append(pair)
        if len(sources) == 3:
            break
    paths = [
        [
            report.path(route_graph, p, rails[station][j])
            for station, p in enumerate(pair)
        ]
        for j, pair in enumerate(sources)
    ]
    perm = np.arange(board.n)
    for pair in paths:
        for path in pair:
            report.transpose(perm, path)
    for j, pair in enumerate(sources):
        assert [perm[rail[j]] for rail in rails] == pair
    for rail in rails:
        assert all(board.degree[p] == 5 for p in rail)
        assert all(graph.has_edge(*edge) for edge in combinations(rail, 2))
    # The three reported bits reuse the same fixed-charge pair. Both transport
    # paths and the local encode swap avoid every parked data rail and standard.
    bus = route_graph.subgraph(set(route_graph) - set(sum(rails, [])))
    ends = []
    for s in stations:
        home = bus.subgraph(range(12 * s, 12 * (s + 1)))
        ends.append(min(tuple(sorted(edge)) for edge in home.edges))
    message = [report.path(bus, a, b) for a, b in zip(*ends, strict=True)]
    return dict(
        stations=stations,
        rails=rails,
        guards=guards,
        standards=standards,
        sources=sources,
        source_routes=paths,
        bus_ends=ends,
        bus_routes=message,
    )


def raw_read(stock, station, addresses):
    graph, _ = report.layout()
    geo = layout()
    owner = geo["stations"][station]
    local = graph.subgraph(range(12 * owner, 12 * (owner + 1)))
    routes = [
        [(report.path(local, p, g), r) for g, r in geo["standards"][station]]
        for p in addresses
    ]
    return material.read(stock, routes, calibrated_echoes())


@cache
def calibrated_echoes():
    return report.certify()


@cache
def word_read(station, m, mask):
    geo = layout()
    stock = np.zeros(federation(0).n, dtype=int)
    for p, v in geo["guards"].items():
        stock[p] = v
    for i, p in enumerate(geo["rails"][station][:m]):
        stock[p] = (mask >> i) & 1
    return raw_read(stock, station, geo["rails"][station][:m])


@cache
def announce(bit):
    geo = layout()
    near, far = geo["bus_ends"]
    stock = np.zeros(federation(0).n, dtype=int)
    for p, v in geo["guards"].items():
        stock[p] = v
    stock[near[0]] = 1
    before = stock.copy()
    if bit:
        report.transpose(stock, near)
    for path in geo["bus_routes"]:
        report.transpose(stock, path)
    flags = raw_read(stock, 1, far)
    for path in reversed(geo["bus_routes"]):
        report.transpose(stock, path)
    if bit:
        report.transpose(stock, near)
    np.testing.assert_array_equal(stock, before)
    return flags


@cache
def generator(kappa):
    law = candidate.native_law
    ex, co = law.local_term_parts()
    h = ex + kappa * co
    idx = [law.SRC.index((a, b)) for b in range(2) for a in range(2)]
    cols = np.eye(law.DIM)[:, idx]
    out = cols.T @ h @ cols
    np.testing.assert_allclose(h @ cols, cols @ out, atol=2e-13)
    np.testing.assert_allclose(
        expm(-1j * pi * out), np.eye(4)[[0, 2, 1, 3]], atol=2e-13
    )
    # A private zero service has the contrast-squared part on its one active rail.
    zero_ref = [law.SRC.index((x, 0)) for x in range(3)]
    np.testing.assert_allclose(
        co[:, zero_ref], np.eye(law.DIM)[:, zero_ref] * [0, 1, 4], atol=2e-13
    )
    return out


@cache
def contact(m, i, j, time, kappa):
    gate = expm(-1j * time * generator(kappa))
    dim = 1 << m
    matrix = np.zeros((dim, dim), complex)
    for x in range(dim):
        b = ((x >> i) & 1) + 2 * ((x >> j) & 1)
        outside = x & ~((1 << i) | (1 << j))
        for k in range(4):
            y = outside | ((k & 1) << i) | ((k >> 1) << j)
            matrix[y, x] = gate[k, b]
    return matrix


def phase(m, i, angle, kappa):
    # Times are fixed using the native eight-unit standard, including at κ=0.
    time = ((-angle) % (2 * pi)) / candidate.KAPPA
    return np.diag([np.exp(-1j * kappa * time * ((x >> i) & 1)) for x in range(1 << m)])


def programmed(target, kappa):
    m = len(target)
    result = np.eye(1 << m, dtype=complex)
    for kind, i, j, theta in compile_unitary(target):
        if kind == "phase":
            result = phase(m, i, theta, kappa) @ result
            continue
        # Refocusing leaves a pure adjacent exchange rotation in every inventory.
        # Subtracting a contact uses its exact 2π period at both tested laws.
        z = phase(m, i, pi, kappa)
        result = phase(m, i, pi / 2, kappa) @ result
        result = contact(m, i, j, theta / 2 % (2 * pi), kappa) @ result
        result = z @ contact(m, i, j, -theta / 2 % (2 * pi), kappa) @ z @ result
        result = phase(m, i, -pi / 2, kappa) @ result
    return result
