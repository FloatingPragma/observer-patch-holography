"""Two material terminal analyzers for all three arm characters and their conjugate basis."""

from functools import cache
import networkx as nx
import numpy as np

from model.executable_model import federation
from ..apparatus import fixture
from ..material_interferometer.apparatus import Terminal, inverse
from ..measure import ONE

FOURIER = np.exp(2j * np.pi * np.outer(np.arange(3), np.arange(3)) / 3) / np.sqrt(3)
BASES = {"arm": np.eye(3), "character": FOURIER}


def trivial(t):
    split = [
        (t.hub, a, d)
        for (a, c, g), d in zip(
            t.arms, (2 * np.arcsin(1 / np.sqrt(3)), np.pi / 2, np.pi), strict=True
        )
    ]
    laps = [
        edge
        for a, c, g in t.arms
        for edge in ((a, g, np.pi / 8), (g, c, np.pi), (c, a, np.pi))
    ]
    return split + laps + inverse(split)


@cache
def material(side):
    lab = Terminal(side)
    assert lab.zero_flags()
    maps = {}
    for name, basis in BASES.items():
        emits = []
        for j in range(3):
            a, c, g = lab.arms[j]
            if name == "character":
                gates = trivial(lab) if j == 0 else lab.recipe(1 if j == 1 else -1)
            else:
                # Write only chosen head, then use a blank spoke for the same
                # 1/3 amplitude-squared aperture as the three-way split.
                gates = [
                    (lab.hub, a, np.pi),
                    (a, g, np.pi / 8),
                    (g, c, np.pi),
                    (c, a, np.pi),
                    (a, lab.hub, np.pi),
                    (lab.hub, a, 2 * np.arccos(1 / np.sqrt(3))),
                ]
            prepared = lab.run(gates)
            out = prepared[0].copy()
            undo = prepared.copy()
            for a, b, dt in inverse(gates):
                lab.gate(undo, a, b, dt)
            np.testing.assert_allclose(
                undo, np.outer(np.eye(10)[0], np.eye(8)[0]), atol=8e-13
            )
            target = np.zeros(8, complex)
            target[ONE] = basis[:, j] / np.sqrt(3)
            phase = np.vdot(target, out) * 3
            np.testing.assert_allclose(abs(phase), 1, atol=1e-12)
            np.testing.assert_allclose(out, target * phase, atol=8e-13)
            emits.append(out)
        maps[name] = np.column_stack(emits)
    return maps


def loan():
    # After the independent B / A charge check the remaining guide is filled.
    # Open one incident unit-pair shutter beside the emitting hub, swap in a
    # held probe from a blank exterior station, and use it for the donor read.
    f = fixture()
    board = federation(0)
    graph = nx.Graph(board.seams)
    guide = set(f["ports"])
    excluded = guide | set(f["donor"])
    outside = set(graph) - excluded
    transit = nx.Graph(f["blank"])
    receipts = []
    for key in ("hub", "tip"):
        hub = f[key]
        candidates = []
        remote = f["tip"] if key == "hub" else f["hub"]
        walk = transit.copy()
        walk.remove_node(remote)
        for p in guide - {remote}:
            for q in set(graph[p]) & outside:
                for r in sorted(set(graph[q]) & outside - {q}):
                    inside = nx.shortest_path(walk, p, hub)
                    candidates.append([r, q, *inside])
        trail = min(candidates, key=lambda p: (len(p), p))
        # Hold the countercharge on a disjoint port. The shortest stock
        # route visits only blanks after the unit-pair shutter opens.
        assert all(graph.has_edge(x, y) for x, y in zip(trail, trail[1:]))
        balance = min(outside - set(trail))
        # Test the second probe both with the original charge at the other
        # hub and with it parked coherently on these local storage ports.
        storage = [c for _, c, _ in Terminal(key).arms]
        pairs = list(zip(trail, trail[1:]))
        for spectator in [remote, *storage]:
            records = np.ones(board.n, int)
            records[[spectator, trail[0]]] = 0
            records[[f["donor"][0], balance]] = 2
            assert records[trail[1]] == records[trail[2]] == 1
            expected = records.copy()
            expected[hub], expected[trail[0]] = expected[trail[0]], expected[hub]
            for i, j in [*pairs, *reversed(pairs[:-1])]:
                records[i], records[j] = records[j], records[i]
            np.testing.assert_array_equal(records, expected)
            assert records[hub] == records[spectator] == 0
        receipts.append(
            dict(
                terminal=hub,
                closed_edge=trail[1:3],
                borrowed_probe_route=trail,
                held_balance=balance,
            )
        )
    return receipts


@cache
def extract(side):
    """Place an unknown one-hot head word on three real vacancy positions."""
    lab = Terminal(side)
    split = [
        (lab.hub, a, t)
        for (a, c, g), t in zip(
            lab.arms, (2 * np.arcsin(1 / np.sqrt(3)), np.pi / 2, np.pi), strict=True
        )
    ]
    route = split.copy()
    for a, c, g in lab.arms:
        route.extend(((a, c, np.pi), (c, g, np.pi), (a, g, np.pi / 8), (a, c, np.pi)))
    c_ports = [c for _, c, _ in lab.arms]
    g_ports = [g for _, _, g in lab.arms]
    storage = [lab.ports.index(c) for c in c_ports]
    rejected = [lab.ports.index(g) for g in g_ports]
    emission = []
    for m in range(8):
        final = lab.run(route, m)
        residual = final.copy()
        residual[storage + rejected] = 0
        np.testing.assert_allclose(residual, 0, atol=5e-13)
        # Reading all g occupations as filled selects the storage c subspace.
        # A second probe's sharp zero flags act only on h,a,g and never on c.
        emission.append(final[storage, 0])
    matrix = np.column_stack(emission)
    phases = np.diag(matrix[:, ONE]) * np.sqrt(3)
    np.testing.assert_allclose(abs(phases), 1, atol=5e-13)
    standard = np.zeros((3, 8), complex)
    standard[:, ONE] = np.diag(phases) / np.sqrt(3)
    np.testing.assert_allclose(matrix, standard, atol=5e-13)
    assert lab.zero_flags()
    return matrix, c_ports, g_ports
