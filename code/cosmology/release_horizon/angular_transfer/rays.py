"""Compile the signed port steps and their three-component payload on the registered federation."""

from itertools import combinations
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from radial_transport.carrier_routes import Transport, coords
from relay import CELL, align, apply
from model.executable_model.carrier import (
    cell,
    icosahedron_from_coordinates,
    isomorphism,
)
from model.executable_model.symmetry import rotations

PORTS = (2, 3, 4)
XYZ = np.array([coords(i) for i in range(125)])


def frame():
    board = cell()
    vertices, edges = icosahedron_from_coordinates()
    vertices = vertices[isomorphism(board.seams, edges)]
    phi = (1 + np.sqrt(5)) / 2
    assert np.all(
        np.isclose(abs(vertices), 0)
        | np.isclose(abs(vertices), 1)
        | np.isclose(abs(vertices), phi)
    )
    # Evaluation Z[phi] -> F5 is well defined since 3²=3+1 mod 5.
    delta = np.rint(
        np.where(np.isclose(abs(vertices), phi), 3 * np.sign(vertices), vertices)
    ).astype(int)
    faces = [
        tuple(sorted(t))
        for t in combinations(range(12), 3)
        if all(b in CELL[a] for a, b in combinations(t, 2))
    ]
    rays = [(i,) for i in range(12)] + faces
    # Check that the spatial quotient carries the same sixty frame actions
    # as the angular ports. The inverse below is verified in F5.
    basis_ids = next(
        list(ii)
        for ii in combinations(range(12), 3)
        if int(round(np.linalg.det(delta[list(ii)]))) % 5
    )
    basis = delta[basis_ids]
    det = int(round(np.linalg.det(basis)))
    inverse = (np.rint(np.linalg.inv(basis) * det).astype(int) * pow(det, -1, 5)) % 5
    np.testing.assert_array_equal(inverse @ basis % 5, np.eye(3, dtype=int))
    for rot in rotations():
        image = delta[list(rot)]
        action = inverse @ image[basis_ids] % 5
        np.testing.assert_array_equal(delta @ action % 5, image % 5)
        assert {tuple(sorted(rot[p] for p in word)) for word in rays} == set(rays)
    directions = np.stack([vertices[list(word)].sum(axis=0) for word in rays])
    directions /= np.linalg.norm(directions, axis=1)[:, None]
    return delta, rays, directions


def command_table(net):
    # Read a requested port of a one-load command register through a real
    # external seam, against a vacant local reference; restore the command.
    a, b = next(
        (a, b) for a, b in net.board.external if a // 12 >= 125 and b // 12 >= 125
    )
    g, h = a % 12, b % 12
    lookup = []
    for p in range(12):
        ops = align([p], [g])
        row = []
        for token in range(12):
            left = [2 * int(k == token) for k in range(12)]
            right = [0] * 12
            apply(left, ops)
            left[g], right[h] = right[h], left[g]
            assert right[min(CELL[h])] == 0
            row.append(net.pointer[right[h]])
            left[g], right[h] = right[h], left[g]
            apply(left, ops[::-1])
            assert left[token] == 2 and sum(left) == 2 and sum(right) == 0
        lookup.append(row)
    np.testing.assert_array_equal(lookup, np.eye(12, dtype=int))
    return np.asarray(lookup), [a, b]


def cycles(image):
    available = set(range(len(image)))
    while available:
        start = min(available)
        orbit = []
        x = start
        while x in available:
            orbit.append(x)
            available.remove(x)
            x = int(image[x])
        assert x == start
        yield orbit


class Rays:
    def __init__(self):
        net = Transport()
        self.net = net
        allowed, self.gateway = command_table(net)
        delta, self.words, self.directions = frame()
        self.delta = delta
        opposite = np.argmin(self.directions[:12] @ self.directions[:12].T, axis=1)
        self.positive = np.flatnonzero(np.arange(12) < opposite).tolist()
        self.steps = {}
        self.cost = {}
        pairs = {}
        for i in self.positive:
            dest = (XYZ + delta[i]) % 5 @ np.array([25, 5, 1])
            pairs[i] = [(orbit[0], v) for orbit in cycles(dest) for v in orbit[1:]]
        edges = sorted({e for batch in pairs.values() for e in batch})
        for port in PORTS:
            net.compile(edges, port)
        for i in self.positive:
            labels = list(range(net.board.n))
            length = 0
            for a, b in pairs[i]:
                for p in PORTS:
                    ops = net.ops[12 * a + p, 12 * b + p]
                    apply(labels, ops)
                    length += len(ops)
            dest = np.argsort(labels)  # source physical port -> observed arrival
            target = (XYZ + delta[i]) % 5 @ np.array([25, 5, 1])
            expected = np.arange(net.board.n)
            for p in PORTS:
                expected[12 * np.arange(125) + p] = 12 * target + p
            np.testing.assert_array_equal(dest, expected)
            self.steps[i] = dest
            self.steps[int(opposite[i])] = np.asarray(labels)
            self.cost[i] = length
        self.commands = [int(np.flatnonzero(allowed[:, p])[0]) for p in range(12)]
        self.opposite = opposite.tolist()

    def destination(self, owner, word):
        rail = 12 * owner + np.array(PORTS)
        for port in word:
            rail = self.steps[self.commands[port]][rail]
        np.testing.assert_array_equal(rail % 12, PORTS)
        assert len(set(rail // 12)) == 1
        return int(rail[0] // 12)

    def witnesses(self):
        # Terminal address calibration on separate sharp preparations.
        # The second reader borrows a queried rail through the owner's seam.
        links = {}
        for a, b in self.net.board.external:
            for u, v in ((a, b), (b, a)):
                links.setdefault(u // 12, (u, v))
        observers = []
        for owner in range(125):
            u, v = links[owner]
            g, h = u % 12, v % 12
            observers.append(v // 12)
            for port in PORTS:
                self.net.read_ops[port] = align([port, 11], [port, min(CELL[port])])
                move = align([port], [g])
                for token in (*PORTS, None):
                    local = [2 * int(k == token) for k in range(12)]
                    self.net.word[12 * owner : 12 * owner + 12] = local
                    first = self.net.read(owner, port)
                    apply(local, move)
                    remote = [0] * 12
                    local[g], remote[h] = remote[h], local[g]
                    assert remote[min(CELL[h])] == 0
                    second = self.net.pointer[remote[h]]
                    local[g], remote[h] = remote[h], local[g]
                    apply(local, move[::-1])
                    assert local == self.net.word[
                        12 * owner : 12 * owner + 12
                    ] and second == first == int(port == token)
            self.net.word[12 * owner : 12 * owner + 12] = [0] * 12
        return observers

    def audit(self):
        rows = []
        for i, j in combinations(self.positive, 2):
            a, b = self.steps[i], self.steps[j]
            rows.append([i, j, int(np.count_nonzero(a[b] != b[a]))])
        assert all(row[-1] == 0 for row in rows)
        left, right = self.positive[:2]
        fault = np.arange(self.net.board.n)
        for p in PORTS:
            fault[p], fault[12 + p] = 12 + p, p
        bad = fault[self.steps[left][fault]]
        other = self.steps[right]
        affected = np.flatnonzero(bad[other] != other[bad])
        assert len(affected) > 0
        arrivals = []
        for word in self.words:
            inverse = [self.opposite[p] for p in reversed(word)]
            emitter = self.destination(0, inverse)
            assert self.destination(emitter, word) == 0
            observer = self.destination(0, word)
            arrivals.append([emitter, observer])
        return dict(
            deltas=self.delta.tolist(),
            positive=self.positive,
            opposites=self.opposite,
            gateway=self.gateway,
            address_witnesses=self.witnesses(),
            primitive_swaps=[self.cost[i] for i in self.positive],
            pair_disagreements=rows,
            broken_pair=[left, right],
            broken_owners=sorted(set((affected // 12).tolist())),
            aperture_words=self.words,
            arrivals=arrivals,
            distinct_destinations=len(set(b for _, b in arrivals)),
        )
