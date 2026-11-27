"""Single-ring sources, free native release, and signed material bond analyzers."""

from functools import cache
from math import pi
from types import SimpleNamespace

import networkx as nx
import numpy as np
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo import candidate
from model.executable_model import cell, federation, rotations
from model.experiments.interval_echo_space.shell_source import (
    local_source as native_contact_source,
    prepare,
)
from model.experiments.interval_echo_null_balance.directional import (
    native_analyzer as native_half_exchange,
)
from ..compiler import execute as addressed


@cache
def geometry():
    graph = nx.Graph(cell().seams)
    distance = dict(nx.all_pairs_shortest_path_length(graph))
    anti = tuple(next(j for j, d in distance[i].items() if d == 3) for i in range(12))
    neighbors = sorted(graph[0])
    ring = neighbors[:1]
    while len(ring) < 5:
        ring.append(min((set(graph[ring[-1]]) & set(neighbors)) - set(ring)))
    cycles = []
    for i in range(12):
        mount = next(g for g in rotations() if g[0] == i)
        moved = [mount[j] for j in ring]
        start = moved.index(min(moved))
        cycles.append(moved[start:] + moved[:start])
    return cycles, anti, [[distance[i][j] for j in range(12)] for i in range(12)]


local_source = cache(native_contact_source)
native_analyzer = cache(native_half_exchange)


def board(level):
    return cell() if level == -1 else federation(level)


@cache
def source(level, owner, port, winding):
    car = board(level)
    ring = geometry()[0][port]
    target = np.zeros(car.n, complex)
    target[12 * owner + np.array(ring)] = np.exp(
        2j * pi * winding * np.arange(5) / 5
    ) / np.sqrt(5)
    return prepare_load(car, target)


def prepare_load(car, target, root=None):
    wave, _ = prepare(SimpleNamespace(car=car, active=car.seams), target, root=root)
    return wave


@cache
def free(level, rate):
    car = board(level)
    _, edge = local_source(rate, candidate.KAPPA)
    # The degree reference adds the missing 8 per port to the six-contact total.
    native = candidate.generator(car)
    kin = car.laplacian() * (-edge[0, 1])
    if rate == 1:
        np.testing.assert_allclose(
            (native - kin).diagonal(), 6 * candidate.KAPPA, atol=3e-13
        )
        off = native - kin
        off.setdiag(0)
        if np.linalg.norm(off.data) >= 1e-12:
            raise ArithmeticError("Native bond off-diagonal changed")
    return kin.tocsr()


def evolve(level, owner, port, winding, time, rate=1):
    # For half rate we match prep and read pulses to their calibrated rates:
    # double every exchange dwell, rephase the resulting extra diagonal.
    start = source(level, owner, port, winding)
    return expm_multiply(-1j * time * free(level, rate), start) if time else start


@cache
def detector(level, owner, edge, sorted_parity):
    a, b = edge
    anti = geometry()[1]
    unitary = np.eye(12)
    if sorted_parity:
        # The routed native rotation maps |a>-|a'> into rail a.
        car = board(level)
        mount = next(
            g
            for g in rotations()
            if all(car.degree[12 * owner + g[k]] == 5 for k in (0, 5))
        )
        chart = np.argsort(mount)
        word = [("real", chart[p], chart[anti[p]], -pi / 2) for p in (a, b)]
        unitary = addressed(word)[np.ix_(chart, chart)]
        for p in (a, b):
            want = (np.eye(12)[p] - np.eye(12)[anti[p]]) / np.sqrt(2)
            np.testing.assert_allclose(unitary[p], want, atol=4e-13)
    return native_analyzer() @ unitary[[a, b]]


def endpoint_probabilities(level, owner, local, edge, sorted_parity):
    ab = abs(detector(level, owner, edge, sorted_parity) @ local) ** 2
    return np.r_[ab, max(0.0, 1.0 - sum(ab))]


@cache
def sharp_calibration():
    # A fresh double standard is coupled at readout. Restored native
    # 11 flags distinguish occupancy zero (bright) from one (dark).
    law = candidate.native_law
    bright = np.tile([w == (1, 1) for w in law.SRC], len(law.LABELS) + 1)
    unitary = law.echo()
    for value in (0, 1):
        for w in ((value, 2), (2, value)):
            sharp = np.eye(law.DIM)[:, law.SRC.index(w)]
            recovered = unitary.T @ ((bright == (value == 0)) * (unitary @ sharp))
            np.testing.assert_allclose(recovered, sharp, rtol=0, atol=2e-13)
    return {"vacant": 1, "occupied": 0}
