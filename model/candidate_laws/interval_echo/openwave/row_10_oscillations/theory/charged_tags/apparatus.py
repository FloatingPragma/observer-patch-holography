"""Certify the two blank hatches, conservative guard swaps and positive-time inverse."""

from functools import cache
import networkx as nx
import numpy as np
from scipy.linalg import expm
from model.candidate_laws.interval_echo.openwave.row_09_neutral.arrival_current.apparatus import (
    certify,
    install,
    station,
    transp,
    SENSORS,
)
from model.experiments.interval_echo_angular.charged_orientation.three_arm_spin.cyclic_charge.neutral_ladder.measure import (
    gate,
    ladder,
    WORDS,
)
from model.executable_model import federation


@cache
def bench(sensor):
    certify()
    assert ladder.read_loads() == [0, 1, 2]
    setup = install(sensor, 2, 1)
    (u, v), (x, a, g, d), flag, read, slots, graph = station(sensor)
    board = federation(0)
    graph = graph.copy()
    graph.remove_edges_from([(x, a), (g, a)])
    for edge in [(u, x), (v, a)]:
        assert tuple(sorted(edge)) in {tuple(sorted(e)) for e in board.seams}
        graph.add_edge(*edge)
    load = np.ones(board.n, int)
    for i, z in setup["mounted_nonunit"]:
        load[i] = z
    assert all(load[i] == 1 for i in (u, v, x, a, d))
    # These two unused mounted stocks prepare the input, not a substituted pair.
    zero, double = slots[11], slots[10]
    assert (load[zero], load[double]) == (0, 2)

    def transfer(path, word):
        before = word.copy()
        for left, right in transp(path):
            word[left], word[right] = word[right], word[left]
        s, t = path[0], path[-1]
        before[s], before[t] = before[t], before[s]
        np.testing.assert_array_equal(word, before)

    stocks = {}
    routes = {}
    for name, st in (("02", (zero, double)), ("20", (double, zero))):
        ready = load.copy()
        routes[name] = [nx.shortest_path(graph, i, j) for i, j in zip(st, (u, v))]
        for path in routes[name]:
            transfer(path, ready)
        stocks[name] = ready
        assert tuple(ready[[u, v]]) == tuple(map(int, name))
        assert tuple(ready[[x, a, d, g]]) == (1, 1, 1, 2)
        assert ready[flag[1]] == 0
    # The spare flag rail remains sharp at zero in both input preparations.
    # Internal transits are exact SWAPs, hence hold even on a coherent ladder.
    path = nx.shortest_path(graph, g, flag[1])
    for name, ready in stocks.items():
        for rail in ((u, x), (v, a)):
            transfer(rail, ready)
        transferred = ready.copy()
        transfer(path, ready)
        assert ready[g] == 0 and ready[flag[1]] == 2
        transfer(path, ready)
        np.testing.assert_array_equal(ready, transferred)
    return dict(
        sensor=sensor,
        guard_route=path,
        input_routes=routes,
        hatches=setup["hatches"],
        workstation=setup["workstation"],
        zero_stock=flag[1],
    )


@cache
def inverse(angle, guard, k):
    """Reverse the actual subpulses, using their separate periods at both contrasts."""
    h = ladder.free_seam(k)
    energies = ladder.stored_energies(k)
    echo = ladder.LAW.echo()
    np.testing.assert_allclose(expm(-2j * np.pi * h), np.eye(9), atol=2e-13)
    np.testing.assert_allclose(
        np.exp(-2j * np.pi * energies), np.ones((3, 3)), atol=2e-13
    )
    np.testing.assert_allclose(echo @ echo, np.eye(99), atol=2e-13)
    swap = expm(-1j * np.pi * h)
    np.testing.assert_allclose(
        swap,
        np.eye(9)[:, [ladder.LAW.SRC.index(w[::-1]) for w in ladder.LAW.SRC]],
        atol=2e-13,
    )
    dwell = angle / (8 * ladder.K)
    adjust = (-angle / 2 % (2 * np.pi)) / (4 * ladder.K)

    # An absent pulse is skipped; all other complements have positive duration.
    def reverse(t):
        return 2 * np.pi - t if t else 0

    raw = np.array(ladder.LAW.SRC)
    phases = np.zeros((3, 3), complex)
    for x in range(3):
        for a in range(3):
            col = ladder.LAW.SRC.index((x, guard))
            z = echo[:, col] * np.exp(-1j * reverse(adjust) * energies[x, a])
            wait = np.exp(
                -1j * reverse(dwell) * (energies[raw[:, 0], a] + energies[raw[:, 1], a])
            )
            out = echo @ (z * np.tile(wait, 11))
            np.testing.assert_allclose(np.delete(out, col), 0, atol=3e-13)
            phases[x, a] = out[col]
    begin = expm(-1j * reverse(3 * np.pi / 2) * h)
    end = expm(-1j * reverse(np.pi / 2) * h)
    undo = (
        np.kron(np.eye(3), end)
        @ np.diag(phases.ravel().repeat(3))
        @ np.kron(np.eye(3), begin)
    )
    np.testing.assert_allclose(
        undo @ gate(angle, guard, k), np.eye(len(WORDS)), atol=4e-13
    )
    return undo


def stations():
    return [bench(sensor) for sensor in SENSORS]
