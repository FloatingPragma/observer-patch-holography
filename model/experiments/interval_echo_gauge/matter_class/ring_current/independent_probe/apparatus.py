"""Matched blank exchanges between two material stocks and an independent probe."""

from functools import cache
from math import pi

import networkx as nx
import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from ..apparatus import geometry, sharp_calibration

RAILS = 13
REF = 12


@cache
def layout():
    board = federation(0)
    graph = nx.Graph(board.seams)
    sources = [list(range(12)) + [24], list(range(156, 168)) + [25]]
    probe = list(range(12, 24)) + [26]
    standard = 239
    graph.remove_node(standard)
    paths = []
    for stock in sources:
        for home, destination in zip(stock, probe, strict=True):
            path = nx.shortest_path(graph, home, destination)
            edges = list(zip(path[:-2], path[1:-1]))
            tags = list(range(board.n))
            for a, b in edges:
                tags[a], tags[b] = tags[b], tags[a]
            if [tags[x] for x in path[-2:]] != [home, destination]:
                raise ArithmeticError("Contact has the wrong physical endpoints")
            for a, b in reversed(edges):
                tags[a], tags[b] = tags[b], tags[a]
            if tags != list(range(board.n)):
                raise ArithmeticError("Swap route displaced a spectator")
            paths.append(path)
    # All phases, analyzers and final occupation queries are routed by the
    # same transpositions. The double standard is parked at a fixed endpoint.
    if not nx.is_connected(graph) or not any(board.degree[x] == 5 for x in graph):
        raise ArithmeticError("Gate and degree reference are inaccessible")
    return dict(sources=sources, probe=probe, double_standard=standard, contacts=paths)


@cache
def generator(rate):
    law = candidate.native_law
    exchange, contrast = law.local_term_parts()
    raw = rate * exchange + candidate.KAPPA * contrast
    blank = raw[: len(law.SRC), : len(law.SRC)]
    np.testing.assert_allclose(raw[len(blank) :, : len(blank)], 0, atol=1e-14)
    swap = [law.SRC.index((v, u)) for u, v in law.SRC]
    np.testing.assert_allclose(
        expm(-1j * pi / rate * blank), np.eye(9)[swap], atol=5e-13
    )
    ids = [law.SRC.index(w) for w in ((0, 0), (1, 0), (0, 1), (1, 1))]
    edges = blank[np.ix_(ids, ids)]
    np.testing.assert_allclose(blank[:, ids], np.eye(9)[:, ids] @ edges, atol=1e-13)
    return edges


@cache
def contact(rate, time):
    return expm(-1j * time * generator(rate))


@cache
def source(q, sign, rate):
    ring = geometry()[0][q]
    target = np.exp(2j * pi * sign * np.arange(5) / 5) / np.sqrt(5)
    u = np.eye(RAILS, dtype=complex)
    for step, (a, b) in enumerate(zip(ring[:-1], ring[1:], strict=True)):
        t = 2 * np.arccos(1 / np.sqrt(5 - step)) / rate
        u[[a, b]] = contact(rate, t)[1:3, 1:3] @ u[[a, b]]
    # A single occupied private reference charges kappa per exchange radian.
    for a, z in zip(ring, target, strict=True):
        dwell = ((np.angle(u[a, ring[0]]) - np.angle(z)) % (2 * pi)) / candidate.KAPPA
        u[a] *= np.exp(-1j * candidate.KAPPA * dwell)
    s = u[:, ring[0]]
    np.testing.assert_allclose(s[ring], target, rtol=0, atol=2e-14)
    np.testing.assert_allclose(u.conj().T @ s, np.eye(RAILS)[:, ring[0]], atol=3e-14)
    return s


@cache
def half(rate):
    g = contact(rate, pi / (2 * rate))[1:3, 1:3]
    return g / (g[0, 0] / abs(g[0, 0]))


@cache
def initial(rate):
    """Columns j are independently stocked probe and reference superpositions."""
    out = np.zeros((RAILS, REF), complex)
    splitter = half(rate)
    out[:REF] = np.eye(REF) * splitter[0, 0]
    out[REF] = -1j * splitter[1, 0]
    return out


@cache
def terminal(rate, quadrature):
    phase = quadrature * pi / 2
    return half(rate) @ np.diag([1, np.exp(-1j * (phase + pi / 2))])


@cache
def encounter(q, sign, rate, time):
    s = source(q, sign, rate)
    g = contact(rate, time)
    a, b = g[1, 1] * g[2, 2], g[1, 2] * g[2, 1]
    shared = g[0, 0] * g[3, 3]
    census = np.diag(abs(s) ** 2)
    ray = np.outer(s, s.conj())
    keep = a * np.eye(RAILS) + b * ray + (shared - a - b) * census
    bunch = abs(g[1, 1] * g[1, 2]) ** 2 * (np.eye(RAILS) + ray - 2 * census)
    return keep, bunch


def offer(inputs, q, sign, rate, time):
    """Failure probabilities and unnormalized retained probe on every input column."""
    keep, bunch = encounter(q, sign, rate, time)
    retained = keep @ inputs
    weight = (abs(inputs) ** 2).sum(axis=0)
    one_end = (inputs.conj() * (bunch @ inputs)).sum(axis=0).real
    wrong = weight - 2 * one_end - (abs(retained) ** 2).sum(axis=0)
    return retained, np.stack([one_end, one_end, wrong], axis=-1)


def banks(order, marks, signs, rate, time, quadrature, pairs):
    wave = initial(rate)
    failed = [np.zeros((REF, 3)), np.zeros((REF, 3))]
    for stage, station in enumerate(order):
        wave, failed[stage] = offer(wave, marks[station], signs[station], rate, time)
    i, j = pairs.T
    two = terminal(rate, quadrature) @ np.array([wave[i, j], wave[REF, j]])
    counts = abs(two.T) ** 2
    rest = (abs(wave[:, j]) ** 2).sum(axis=0) - counts.sum(axis=1)
    pmf = np.column_stack([failed[0][j], failed[1][j], counts, rest])
    if pmf.min() < -2e-13:
        raise ArithmeticError("Encounter returned a negative event weight")
    np.testing.assert_allclose(pmf.sum(axis=1), 1, atol=1e-12, rtol=0)
    return np.maximum(0, pmf)


def standards(rate):
    sharp_calibration()
    # Two independently selected reference phases at both interferometer exits.
    splitter = half(rate)
    outcomes = []
    for quadrature in (0, 1):
        for sign in (-1, 1):
            incoming = splitter[:, 0] * [1, -1j * sign * (1j**quadrature)]
            outcomes.append(abs(terminal(rate, quadrature) @ incoming) ** 2)
    return np.array(outcomes)
