"""Local ternary occupation and signed coincidence counters on both receiver seams."""

from functools import cache
from itertools import product

import networkx as nx
import numpy as np
from scipy.linalg import expm

from model.executable_model import federation
from model.experiments.interval_echo_gauge.arrival_ladder.measure import LAW, read_loads
from .geometry import layout
from ....apparatus import fixture
from .....live_trimer.measure import blocks

SENSORS = {"opposite": (207, 206), "second_owner": (106, 96)}


@cache
def verify():
    ring, _, leaves, _, _ = layout()
    raw = fixture()
    board = federation(0)
    graph = nx.Graph(board.seams)
    echo = LAW.echo()
    unit_index = LAW.SRC.index((1, 1))
    unit_mark = echo[:, unit_index].reshape(11, 9)[:, unit_index]
    eye = np.eye(9)
    swap = expm(-1j * np.pi * blocks()["blank"] / 8)
    np.testing.assert_allclose(
        swap, eye[[LAW.SRC.index(x[::-1]) for x in LAW.SRC]], atol=3e-13
    )
    assert read_loads() == [0, 1, 2]

    # A single echo on the actual blank receiving bond writes the signed
    # difference. Equal-sign exceptional pairs cannot occur on this source.
    for x, y in product(range(3), repeat=2):
        i = LAW.SRC.index((x, y))
        written = echo[:, i].reshape(11, 9)
        labels = [False, *(d == x - y for _, d in LAW.LABELS)]
        assert abs(np.linalg.norm(written[labels]) ** 2 - 1) < 2e-13
        np.testing.assert_allclose(echo @ written.ravel(), np.eye(99)[i], atol=3e-14)
    # A unit shutter and exterior unit load allow a blank/nonblank trial. A
    # nonunit endpoint leaves its mark and both loads exactly unchanged.
    for n in range(3):
        for w in [(n, 1), (1, n)]:
            i = LAW.SRC.index(w)
            shutter = np.kron(unit_mark, eye[i])
            prepared = echo @ shutter
            np.testing.assert_allclose(
                np.vdot(prepared[:9], prepared[:9]).real, int(n == 1), atol=2e-14
            )
            np.testing.assert_allclose(
                prepared, np.eye(99)[i] if n == 1 else shutter, atol=2e-14
            )

    excluded = set(ring) | set(leaves) | set(raw["ports"]) | set(raw["donor"])
    checks = []
    for owner, (a, b) in SENSORS.items():
        assert graph.has_edge(a, b) and a // 12 == b // 12
        stations = []
        for u in (a, b):
            v = min(
                x
                for x in graph[u]
                if x not in excluded and x // 12 == u // 12 and board.degree[x] == 5
            )
            stations.append(v)
        # Opening a hatch requires a unit load at the chosen receiver site.
        # If both trials fail, d=+2/-2 on the blank receiver bond resolves both
        # loads directly. Otherwise route the two recorded loads, one at a time,
        # through the open hatch to its ordinary external private reference.
        for x, y in product(range(3), repeat=2):
            if x == y != 1:
                continue
            if x != 1 and y != 1:
                assert (x, y) == ((2, 0) if x - y == 2 else (0, 2))
                continue
            v = int(x != 1)
            for j in (0, 1):
                original = [x, y, 1]
                value = original.copy()
                if v != j:
                    value[v], value[j] = value[j], value[v]
                value[v], value[2] = value[2], value[v]
                assert read_loads()[value[2]] == original[j]
                value[v], value[2] = value[2], value[v]
                if v != j:
                    value[v], value[j] = value[j], value[v]
                assert value == original
        checks.append(dict(name=owner, receiver=[a, b], private_stations=stations))
    return checks
