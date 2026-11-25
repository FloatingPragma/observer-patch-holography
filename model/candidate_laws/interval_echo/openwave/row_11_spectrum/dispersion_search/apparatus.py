"""Native adjacent controls and material interferometers on a terminated guide."""

from functools import cache

import numpy as np
from scipy.linalg import eigh_tridiagonal

from model.experiments.interval_echo_angular.charged_orientation.winding_memory import (
    mixed_masses,
)

MOTIFS = (4, 6, 8, 12, 16)
WIDTHS = (96, 144)


def motif(period):
    return "01" * ((period - 2) // 2) + "00"


class Guide:
    def __init__(self, n, period):
        self.n = n
        self.period = period
        self.word = (motif(period) * (n // period))[:-1]
        d, e = mixed_masses.coefficients(self.word, "even")
        self.diag, self.off = np.array(d) / 8, np.array(e) / 8
        self.levels, self.columns = eigh_tridiagonal(self.diag, self.off)
        # Fix the otherwise arbitrary phase of each real eigenvector at a
        # detector address reflected from the first interior well.
        phase_port = n - 1 - (period + period // 2)
        for column in self.columns.T:
            value = column[phase_port]
            if abs(value) < 1e-12:
                value = column[np.argmax(abs(column))]
            column *= np.sign(value)
        self.w = np.abs(self.off)
        self.potential = self.diag - np.r_[self.w, 0] - np.r_[0, self.w]
        assert set(self.potential) == {0, 8, 16}
        # A full addressed pulse gives +SWAP on a blank and -SWAP on a mark.
        self.swaps = np.array([1 if mark == "0" else -1 for mark in self.word])
        self.sources = [self.word.index(c) for c in "01"]

    @cache
    def transported_sign(self, ref, a, b):
        """Signed physical permutation taking the chosen reference bond to (a,b)."""
        assert a < b
        j = self.sources[ref]
        sites = [j, j + 1]
        sign = 1
        # Move the outermost targets first to keep the rails ordered.
        order = (1, 0) if a >= j + 1 else (0, 1)
        for k in order:
            dest = (a, b)[k]
            start = sites[k]
            low, high = sorted((start, dest))
            assert sites[1 - k] not in range(low, high + 1)
            sign *= int(np.prod(self.swaps[low:high]))
            sites[k] = dest
        assert sites == [a, b]
        return sign

    def pair_phase(self, a, b, angle):
        """Two commuting transported contacts give exp(i angle) on exactly this pair."""
        a, b = sorted((a, b))
        sb = self.transported_sign(0, a, b)
        sm = self.transported_sign(1, a, b)
        rate = 8.5 + 0.5 * sb / sm
        tb = (-angle / rate) % (2 * np.pi)
        tm = (-tb * sb / (8 * sm)) % (np.pi / 4)
        c = -0.5 * sb * tb - 4 * sm * tm
        assert abs(np.sin(c)) < 2e-13
        return np.exp(-1j * (8.5 * tb - 4 * tm)) * np.cos(c)

    def phase(self, state, a, angle):
        b, c = ((a + 1) % self.n, (a + 2) % self.n)
        for i, j, phi in ((a, b, angle / 2), (a, c, angle / 2), (b, c, -angle / 2)):
            state[[i, j]] *= self.pair_phase(i, j, phi)

    def prepare(self, target):
        """Eliminate the tail left to right by physical adjacent contacts and phases."""
        state = np.zeros(self.n, complex)
        state[0] = 1
        tail = np.sqrt(np.cumsum(abs(target[::-1]) ** 2))[::-1]
        for j in range(self.n - 1):
            theta = np.arctan2(tail[j + 1], abs(target[j]))
            t = theta / self.w[j]
            center = 8.5 if self.word[j] == "0" else -4.0
            z = state[j] * np.exp(-1j * t * center)
            state[j : j + 2] = [z * np.cos(theta), 1j * z * np.sin(theta)]
            self.phase(state, j, np.angle(target[j]) - np.angle(state[j]))
        self.phase(state, self.n - 1, np.angle(target[-1]) - np.angle(state[-1]))
        np.testing.assert_allclose(state, target, atol=3e-12)
        return state

    def read_pair(self, state, a, b):
        """Charge at b after a physical half splitter counts the real antisymmetric rail."""
        s = self.transported_sign(0, a, b)
        # A transported blank half swap has relative matrix I+i*s*X.
        work = state.astype(complex, copy=True)
        self.phase(work, b, -s * np.pi / 2)
        work = work[[a, b]]
        probe = (work + 1j * s * work[::-1]) / np.sqrt(2)
        return abs(probe) ** 2


def geometry():
    ports, _ = mixed_masses.material.geometry()
    board = mixed_masses.material.federation(2)
    mixed_masses.contacts()
    assert all(
        ports[j + 1] in board.neighbours[ports[j]] for j in range(len(ports) - 1)
    )
    return ports, board


def certify_chart(ports, board, guide, sign):
    """Native lifted columns, host preparation and sharp charge/head reader on this path."""
    import networkx as nx
    from scipy.linalg import expm
    from model.experiments.interval_echo_angular.charged_orientation.live_trimer.measure import (
        LAW,
        blocks,
    )

    path = list(map(int, ports))
    active = [tuple(sorted((a, b))) for a, b in zip(path, path[1:])]
    marks = [e for e, kind in zip(active, guide.word, strict=True) if kind == "1"]
    graph = nx.Graph(board.seams)
    dist = nx.multi_source_dijkstra_path_length(graph, set(path))
    donor = max(graph, key=lambda s: (dist[s], -s))
    assert dist[donor] > 2
    bypass = graph.copy()
    bypass.remove_node(donor)
    bypass.remove_edges_from(marks)
    assert nx.is_biconnected(bypass)
    outside = {
        e
        for x in path
        for y in board.neighbours[x]
        if (e := tuple(sorted((x, y)))) in graph.edges
    } - set(active)
    quiet = graph.copy()
    quiet.remove_edges_from(outside)
    assert set(nx.node_connected_component(quiet, path[0])) == set(path)
    assert guide.word[0] == "0"
    # The defect can park at port 1 by an addressed blank swap while writing any
    # exterior contact of port 0. All other exterior contacts have unit pairs.
    assert all(path[1] not in cut for cut in outside if path[0] in cut)

    coin0 = int(sign > 0)
    raw = blocks()
    h = np.zeros((guide.n, guide.n))
    for x in range(guide.n):
        for a in (x - 1, x):
            if a < 0 or a >= guide.n - 1:
                continue
            live = guide.word[a] == "1"
            local = raw[1, 1] if live else raw["blank"]
            coin = (coin0 ^ int(a < x)) if live else 0
            loads = tuple(1 + sign * int(i == x) for i in (a, a + 1))
            start = 9 * coin + LAW.SRC.index(loads)
            host = 9 * coin + LAW.SRC.index((1, 1))
            h[x, x] -= local[host, host] / 8
            for end in np.flatnonzero(local[:, start]):
                after = LAW.SRC[end % 9]
                row = a + after.index(1 + sign)
                if live:
                    assert end // 9 == (coin0 ^ int(a < row))
                h[row, x] += local[end, start] / 8
    expected = np.diag(guide.diag) + np.diag(guide.off, 1) + np.diag(guide.off, -1)
    np.testing.assert_array_equal(h, expected)
    # A fresh sharp-load copy, followed by an addressed full marked pulse,
    # reads the direction on arrival by testing whether the defect moves back.
    act = expm(-1j * np.pi / 8 * raw[1, 1] / 8)
    rr = LAW.SRC.index((1, 1 + sign))
    left = LAW.SRC.index((1 + sign, 1))
    reader = []
    for c in (0, 1):
        reader.append(float(sum(abs(act[left::9, 9 * c + rr]) ** 2)))
    np.testing.assert_allclose(reader, [coin0, 1 - coin0], atol=1e-13)
    return dict(
        marks=len(marks),
        exterior_shutters=len(outside),
        donor=int(donor),
        private_zero_writes=sum(int(board.degree[p] == 5) for p in path),
        column_count=guide.n,
        maximum_sharp_charge_dimension=1,
        reverse_transfer_by_arrived_head=[round(a, 12) for a in reader],
    )
