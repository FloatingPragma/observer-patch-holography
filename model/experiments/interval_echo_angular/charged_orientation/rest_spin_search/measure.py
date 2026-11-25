"""Stationary charged pentagon doublets on a connected, two-owner Interval Echo bench."""

import json
from functools import cache
from itertools import pairwise
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.linalg import expm
from scipy.stats import beta

from model.executable_model import federation
from model.executable_model.symmetry import rotations
from model.candidate_laws.interval_echo.openwave.row_11_spectrum.conjugate_flight import (
    measure as flight,
)

HERE = Path(__file__).resolve().parent
K = int(flight.candidate.KAPPA)
SHOTS = 2**20
EDGES = [
    *((i, (i + 1) % 5) for i in range(5)),
    *((i, 5) for i in range(5)),
    (5, 6),
    *((7 + i, 7 + (i + 1) % 5) for i in range(5)),
    *((6, 7 + i) for i in range(5)),
]
ORDER = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)


CUT = {(4, 0), (1, 5), (11, 7), (6, 8)}


def fixture(released=False):
    board = federation(0)
    hubs = board.external[0]
    rims = []
    actions = []
    for hub in hubs:
        owner = hub // 12
        ring = sorted(p for p in board.neighbours[hub] if p // 12 == owner)
        p = ring[0]
        action = next(
            r
            for r in rotations()
            if r[hub % 12] == hub % 12
            and 12 * owner + r[p % 12] in set(ring).intersection(board.neighbours[p])
        )
        rim = []
        for _ in range(5):
            rim.append(p)
            p = 12 * owner + action[p % 12]
        assert p == rim[0] and set(rim) == set(ring)
        rims.append(rim)
        actions.append(list(action))
    ports = rims[0] + list(hubs) + rims[1]
    active = {
        tuple(sorted((ports[i], ports[j])))
        for i, j in EDGES
        if not released or (i, j) not in CUT
    }
    assert active <= set(board.seams)
    assert all(p in board.neighbours[q] for p, q in pairwise(ports))
    assert all(tuple(sorted((a, b))) in active for a, b in pairwise(ports))
    region = set(ports)
    shutters = {e for e in board.seams if region.intersection(e)} - active
    tip = ports[1]
    ref = min(
        p
        for p in board.neighbours[tip]
        if p not in region and p // 12 == tip // 12 and board.degree[p] == 5
    )
    lane = tuple(sorted((tip, ref)))
    far_tip = ports[8]
    far_ref = min(
        p
        for p in board.neighbours[far_tip]
        if p not in region and p // 12 == far_tip // 12 and board.degree[p] == 5
    )
    assert lane in shutters and tuple(sorted((far_tip, far_ref))) in shutters
    store = next(
        e
        for e in board.seams
        if all(
            p not in region | {ref, far_ref}
            and not (region | {ref, far_ref}).intersection(board.neighbours[p])
            for p in e
        )
    )
    store_cuts = {e for e in board.seams if set(store).intersection(e) and e != store}
    w = np.ones(board.n, int)
    w[store[0]], w[ports[0]] = 2, 0
    for e in store_cuts:
        relay = store[0] in e
        if relay:
            w[store[0]], w[store[1]] = w[store[1]], w[store[0]]
        assert all(w[p] == 1 for p in e)
        if relay:
            w[store[0]], w[store[1]] = w[store[1]], w[store[0]]
    graph = nx.Graph(board.seams)
    graph.remove_edges_from(store_cuts)
    for p in ports:
        if board.degree[p] == 5:
            for a, b in flight.closed_walk(nx.shortest_path(graph, ports[0], p)):
                w[a], w[b] = w[b], w[a]
            assert w[p] == 0
            for a, b in flight.closed_walk(nx.shortest_path(graph, ports[0], p)):
                w[a], w[b] = w[b], w[a]
    for e in shutters - {lane}:
        relay = ports[0] in e
        if relay:
            w[ports[0]], w[hubs[0]] = w[hubs[0]], w[ports[0]]
        assert all(w[p] == 1 for p in e)
        if relay:
            w[ports[0]], w[hubs[0]] = w[hubs[0]], w[ports[0]]
    graph.remove_edges_from(shutters)
    assert set(nx.node_connected_component(graph, hubs[0])) == region
    assert set(nx.node_connected_component(graph, store[0])) == set(store)
    assert flight.meter() == [0, 1, 2]
    return dict(
        ports=ports,
        owner_actions=actions,
        open_bridge=list(hubs),
        private_erased=[p for p in ports if board.degree[p] == 5],
        reference_lane=lane,
        other_reference_lane=[far_tip, far_ref],
        shutter_count=len(shutters),
        stock=store,
        stock_shutters=len(store_cuts),
    )


@cache
def native():
    blank, _ = flight.seam_data(K)
    ids = [flight.LAW.SRC.index(v) for v in ((0, 1), (1, 0))]
    two = blank[np.ix_(ids, ids)]
    np.testing.assert_array_equal(np.rint(2 * two).astype(int), [[17, -1], [-1, 17]])
    h = np.zeros((12, 12))
    for a, b in EDGES:
        h[np.ix_([a, b], [a, b])] += two
    return h, two


def contact(a, b, t):
    u = np.eye(12, dtype=complex)
    u[np.ix_([a, b], [a, b])] = expm(-1j * t * native()[1])
    return u


@cache
def transpose(a, b):
    path = list(range(a, b + 1)) if a <= b else list(range(a, b - 1, -1))
    u = np.eye(12, dtype=complex)
    for x, y in flight.closed_walk(path):
        u = contact(x, y, np.pi) @ u
    result = np.eye(12)
    result[[b, a]] = result[[a, b]]
    np.testing.assert_allclose(u, result, atol=1e-12)
    return u


def phase(site, angle, station=1):
    # Open reference lane meets rim 1. Reroute any one coordinate onto it.
    relay = transpose(site, station)
    blank, private = flight.seam_data(K)
    swap = expm(-1j * np.pi * blank)
    dwell = (angle % (2 * np.pi)) / K
    diag = np.diag(
        [np.exp(-1j * K * private[3 * b] * dwell) for _, b in flight.LAW.SRC]
    )
    acted = swap.conj().T @ diag @ swap
    indices = [flight.LAW.SRC.index((i, 1)) for i in range(2)]
    unit = acted[indices[1], indices[1]]
    np.testing.assert_allclose(
        acted[:, indices] / unit,
        np.eye(9)[:, indices] @ np.diag([np.exp(1j * angle), 1]),
        atol=1e-12,
    )
    d = np.eye(12, dtype=complex)
    d[station, station] = acted[indices[0], indices[0]] / unit
    return relay.conj().T @ d @ relay


def emit(target):
    u = np.eye(12, dtype=complex)
    for j in range(11):
        tail = np.linalg.norm(target[j + 1 :])
        if tail > 1e-13 or abs(target[j]) > 1e-13:
            u = contact(j, j + 1, 2 * np.arctan2(tail, abs(target[j]))) @ u
        if abs(target[j]) > 1e-13:
            u = phase(j, np.angle(target[j]) - np.angle(u[j, 0])) @ u
    if abs(target[-1]) > 1e-13:
        u = phase(11, np.angle(target[-1]) - np.angle(u[-1, 0])) @ u
    np.testing.assert_allclose(u[:, 0], target, atol=1e-12)
    return u


def shutter(state, station=1, hub=5):
    """Use an empty hub as a buffer while writing or unwriting a filled side seam."""
    assert abs(state[hub]) < 1e-10
    if station == 1:
        relay = transpose(1, 5)
    else:
        relay = transpose(8, 6)
    parked = relay @ state
    assert abs(parked[station]) < 1e-10
    law = flight.LAW
    col = np.eye(99)[:, law.SRC.index((1, 1))]
    head = law.echo() @ col
    np.testing.assert_allclose(law.echo() @ head, col, atol=2e-13)
    for term in law.local_term_parts():
        np.testing.assert_allclose(term @ head, 0, atol=1e-13)
    np.testing.assert_allclose(relay @ parked, state, atol=1e-12)


def mode(m, side=0):
    x = np.zeros(12, complex)
    x[7 * side : 7 * side + 5] = np.exp(2j * np.pi * m * np.arange(5) / 5) / np.sqrt(5)
    return x


def bond_current(state, rng):
    # A calibrated half swap maps imaginary coherence into a population difference.
    hits = []
    for j in range(5):
        i = (j + 1) % 5
        y = contact(j, i, np.pi / 2) @ state
        p = abs(y) ** 2
        p /= sum(p)
        n = rng.multinomial(SHOTS, p)
        hits.append([int(n[j]), int(n[i]), int(SHOTS - n[j] - n[i])])
        np.testing.assert_allclose(
            (p[i] - p[j]) / 2, np.imag(state[j].conjugate() * state[i]), atol=1e-12
        )
    return hits


def half_mirror():
    # Rotate both disjoint pairs by the same angle while spectator rim 0 stays fixed.
    u = np.eye(12, dtype=complex)
    for a, b in [(1, 4), (2, 3)]:
        if (a, b) in EDGES or (b, a) in EDGES:
            move = np.eye(12)
            other = b
        else:
            other = (a + 1) % 5
            move = transpose(other, b)
        u = move.conj().T @ contact(a, other, np.pi / 2) @ move @ u
    return u


def rotation():
    turn = np.eye(12, dtype=complex)
    for j in reversed(range(4)):
        turn = transpose(j, j + 1) @ turn
    for j in range(5):
        np.testing.assert_allclose(turn[:, j], np.eye(12)[:, (j + 1) % 5], atol=2e-12)
    return turn


def cost_meter(state, rng):
    # The centrality cost is diagonal in signed position. Other arms are
    # independent copies read on each native bond in its swap-parity basis.
    position = rng.multinomial(SHOTS, abs(state) ** 2 / np.vdot(state, state).real)
    degree = np.bincount(np.ravel(EDGES), minlength=12)
    odds = []
    for a, b in EDGES:
        # The source is localized in the left wheel. Empty remote contacts
        # give identically zero; the phase lane is required only on the left.
        if abs(state[a]) + abs(state[b]) < 1e-13:
            odds.append(0)
            continue
        read = contact(a, b, np.pi / 2) @ phase(b, np.pi / 2) @ state
        expected = abs(state[a] - state[b]) ** 2 / 2
        np.testing.assert_allclose(abs(read[a]) ** 2, expected, atol=1e-12)
        odds.append(int(rng.binomial(SHOTS, np.clip(abs(read[a]) ** 2, 0, 1))))
    # Two full lines, with a single joint Bonferroni family. The 24-unit
    # potential at each rim differs from a gateway by exactly 24.
    tail = 0.01 / (2 * 2 * (len(EDGES) + 1))

    def intervals(k):
        low = 0 if k == 0 else beta.ppf(tail, k, SHOTS - k + 1)
        high = 1 if k == SHOTS else beta.isf(tail, k + 1, SHOTS - k)
        return np.array([low, high])

    band = 24 + 24 * intervals(int(position[5] + position[6]))
    for hits in odds:
        band += intervals(hits)
    return dict(
        port_counts=position.tolist(),
        odd_bond_counts=odds,
        simultaneous99=band.tolist(),
        measured=float((K * (degree @ position) + sum(odds)) / SHOTS),
    )


def supplied_reference():
    # A supplied two-component spin under the SAME 72-degree material angle.
    spin = expm(-1j * (2 * np.pi / 5) * np.diag([0.5, -0.5]))
    reference = np.array([1.0, 0.0])
    marks = {}
    for steps in (0, 1, 2, 5, 10):
        amplitude = (
            1 + reference @ np.linalg.matrix_power(spin, steps) @ reference
        ) / 2
        marks[str(steps)] = round(float(abs(amplitude) ** 2), 12)
    return marks


def run():
    from .release import pair_meter

    bench = fixture()
    h = native()[0]
    rng = np.random.default_rng(192516)
    turn = rotation()
    mixer = half_mirror()
    records = []
    for m in (1, 2):
        left, right = mode(m), mode(m, 1)
        opposite = mode(-m)
        e = 25.5 - np.cos(2 * np.pi * m / 5)
        u = emit(left)
        for v in (left, right, opposite):
            np.testing.assert_allclose(h @ v, e * v, atol=1e-12)
        shutter(u[:, 0])
        after = expm(-27j * h) @ u[:, 0]
        shutter(after)
        np.testing.assert_allclose(
            abs(u.conj().T @ after) ** 2, np.eye(12)[0], atol=1e-12
        )
        rotated = mixer @ after
        receiver_axes = pair_meter(left, opposite, tuple(range(5)))
        p = abs((receiver_axes @ rotated)[:2]) ** 2
        back = abs((u.conj().T @ mixer.conj().T @ rotated)[0]) ** 2
        np.testing.assert_allclose(p, [0.5, 0.5], atol=1e-12)
        np.testing.assert_allclose(
            abs(left.conj() @ mixer.conj().T @ rotated) ** 2, 1, atol=2e-12
        )
        ramsey = (left + right) / np.sqrt(2)
        f = emit(ramsey)
        shutter(f[:, 0])
        scan = []
        for steps in (0, 1, 2, 5, 10):
            y = np.linalg.matrix_power(turn, steps) @ expm(-7j * h) @ f[:, 0]
            shutter(y)
            actual = abs((f.conj().T @ y)[0]) ** 2
            expect = (1 + np.cos(2 * np.pi * m * steps / 5)) / 2
            np.testing.assert_allclose(actual, expect, atol=2e-12)
            scan.append(
                dict(
                    turns=steps,
                    return_hits=int(rng.binomial(SHOTS, np.clip(actual, 0, 1))),
                    probability=round(float(actual), 12),
                )
            )
        # Break symmetry with one ordinary (vacancy -8) private reference,
        # available at physical rim 0 of the first wheel.
        assert bench["ports"][0] in bench["private_erased"]
        asym = h.copy()
        asym[0, 0] -= 8
        pert = asym @ left
        variance = np.linalg.norm(pert - (left.conj() @ pert) * left) ** 2
        np.testing.assert_allclose(variance, 256 / 25 * 1.0, atol=1e-10)
        records.append(
            dict(
                m=m,
                formation=e,
                variance=0,
                sharp_charge=-1,
                cost_meter=cost_meter(after, rng),
                clockwise=bond_current(after, rng),
                counterclockwise=bond_current(opposite, rng),
                interferometer=scan,
                rest_return_hits=SHOTS,
                mirror_plus_minus=rng.multinomial(SHOTS, p / p.sum()).tolist(),
                mirror_inverse_hits=int(rng.binomial(SHOTS, np.clip(back, 0, 1))),
                broken_chiral_variance=round(float(variance), 12),
            )
        )
    # The azimuthal unit phase kick converts the +1 winding into the bright
    # symmetric wave; the identical kick sends -1 to -2.
    kick = np.eye(12, dtype=complex)
    for j in range(5):
        kick = phase(j, -2 * np.pi * j / 5) @ kick
    release = []
    for t in (0, 16, 128, 512, 2048, 8192):
        probs = []
        angular = []
        for m in (1, -1):
            src = kick @ mode(m)
            shutter(src)
            x = expm(-1j * t * h) @ src
            probs.append(float(np.sum(abs(x[6:]) ** 2)))
            angular.append(
                float(sum(abs(mode(q, 1).conj() @ x) ** 2 for q in (1, -1, 2, -2)))
            )
        assert max(angular) < 1e-19
        release.append(
            dict(
                time=t,
                regional_probabilities=probs,
                incoming_chiral_charge_arrivals=[
                    int(rng.binomial(SHOTS, np.clip(p, 0, 1))) for p in probs
                ],
            )
        )
    return dict(
        fixture=bench,
        starts=SHOTS,
        stationary=records,
        open_bridge_kick=release,
        supplied_spinor_reference=supplied_reference(),
    )


if __name__ == "__main__":
    (HERE / "results.json").write_text(json.dumps(run(), indent=2) + "\n")
