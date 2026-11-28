"""A written collar with independent material routes on two actual outside faces."""

from collections import defaultdict, deque
from itertools import pairwise

import networkx as nx
import numpy as np
from scipy import sparse


def choose_collar(path):
    fan = path.cells[-1]
    enrolled = set(fan)
    original = {frozenset((fan[a], fan[b])) for a, b in path.field.edges}
    faces = sorted(
        {
            tuple(sorted(face))
            for face in path.board.faces()
            if enrolled.intersection(face)
        }
    )
    choices = []
    for first in faces:
        shared = frozenset(enrolled.intersection(first))
        outside = set(first) - enrolled
        if shared not in original or len(outside) != 1:
            continue
        for second in faces:
            outer = set(second) - enrolled
            if (
                len(outer) == 2
                and len(set(first) & set(second)) == 2
                and outside <= outer
            ):
                choices.append((first, second, tuple(sorted(outer))))
    first, second, probe = min(choices)
    assert not set(probe) & enrolled
    board_graph = nx.Graph(path.board.seams)
    ranges = [nx.shortest_path_length(board_graph, fan[0], p) for p in probe]
    return [first, second], probe, ranges, len(choices)


class Collar:
    def __init__(self, path, probe_sign):
        self.source_sign = path.sign
        self.probe_sign = probe_sign
        self.fan = path.cells[-1]
        self.hub = self.fan[0]
        faces, self.probe, self.ranges, self.eligible = choose_collar(path)
        self.faces = faces
        self.original = [(self.fan[a], self.fan[b]) for a, b in path.field.edges]
        self.inside_mask = (1 << len(self.original)) - 1
        internal_sets = {frozenset(pair) for pair in self.original}
        outside_seams = set()
        for triangle in faces:
            outside_seams.update(
                tuple(sorted(edge))
                for edge in pairwise((*triangle, triangle[0]))
                if frozenset(edge) not in internal_sets
            )
        self.added = sorted(outside_seams)
        assert len(self.added) == 4 and self.probe in self.added
        self.edges = self.original + self.added
        self.edge_index = {frozenset(edge): j for j, edge in enumerate(self.edges)}
        self.source_moves = defaultdict(list)
        for bit, (a, b) in enumerate(self.original):
            self.source_moves[a].append((b, 1 << bit, True))
            self.source_moves[b].append((a, 1 << bit, False))
        self.probe_mask = 1 << self.edge_index[frozenset(self.probe)]
        triangles = [tuple(self.fan[p] for p in tri) for tri in path.field.rings]
        self.rings = []
        for k, triangle in enumerate([*triangles, *faces]):
            kind = (
                "fan_rings"
                if k < len(triangles)
                else ("connector" if k == len(triangles) else "outer_ring")
            )
            mask = clockwise = 0
            for p, q in pairwise((*triangle, triangle[0])):
                bit = self.edge_index[frozenset((p, q))]
                mask |= 1 << bit
                clockwise |= int(self.edges[bit][0] == p) << bit
            self.rings.append((kind, mask, clockwise))
        for kind, mask, _ in self.rings:
            if kind == "fan_rings":
                assert not mask & ~self.inside_mask
            if kind == "outer_ring":
                assert not mask & self.inside_mask
            if kind == "connector":
                assert mask & self.inside_mask and mask & ~self.inside_mask
        assert not self.probe_mask & self.inside_mask
        self.start_words = path.words_at[0]
        starts = [
            (self.hub, self.probe[0], w | (background << len(self.original)))
            for background in range(16)
            for w in self.start_words
        ]
        self.basis = list(starts)
        self.ids = {s: i for i, s in enumerate(starts)}
        self.terms = self.build(path.board.degree)

    def moves(self, s):
        source, probe, word = s
        for target, bit, first in self.source_moves[source]:
            if (bool(word & bit) == first) == (self.source_sign == 1):
                yield "source_hop", 4, (target, probe, word ^ bit)
        if ((probe == self.probe[0]) == bool(word & self.probe_mask)) == (
            self.probe_sign == 1
        ):
            target = self.probe[1] if probe == self.probe[0] else self.probe[0]
            yield "probe_hop", 4, (source, target, word ^ self.probe_mask)
        for kind, mask, clockwise in self.rings:
            if word & mask in (clockwise, mask ^ clockwise):
                yield kind, 0.5, (source, probe, word ^ mask)

    def build(self, degrees):
        terms = defaultdict(list)
        queue = deque(self.basis)
        while queue:
            s = queue.popleft()
            i = self.ids[s]
            for kind, rate, target in self.moves(s):
                if target not in self.ids:
                    self.ids[target] = len(self.basis)
                    self.basis.append(target)
                    queue.append(target)
                terms[kind].append((i, self.ids[target], rate))
        dimension = len(self.basis)
        named = {}
        for key, values in terms.items():
            i, j, rate = zip(*values, strict=True)
            off = sparse.coo_matrix(
                (-np.asarray(rate), (i, j)), shape=(dimension,) * 2
            ).tocsr()
            assert (off != off.T).nnz == 0
            named[key] = off - sparse.diags(
                np.asarray(off.sum(axis=1)).ravel(), dtype=float
            )
        assert set(named) == {
            "source_hop",
            "probe_hop",
            "fan_rings",
            "connector",
            "outer_ring",
        }
        named["source_degree"] = sparse.diags(
            [8 * (6 - int(degrees[x])) for x, _, _ in self.basis],
            format="csr",
            dtype=float,
        )
        named["probe_degree"] = sparse.diags(
            [8 * (6 - int(degrees[x])) for _, x, _ in self.basis],
            format="csr",
            dtype=float,
        )
        return named

    def generators(self):
        h = self.terms
        neutral = (
            h["probe_hop"]
            + h["fan_rings"]
            + h["connector"]
            + h["outer_ring"]
            + h["probe_degree"]
        )
        mobile = neutral + h["source_hop"] + h["source_degree"]
        at_hub = np.array([x == self.hub for x, _, _ in self.basis])
        np.testing.assert_array_equal(h["source_degree"].diagonal()[at_hub], 0)
        held = (h["source_degree"] + h["source_hop"]).diagonal()[at_hub]
        assert np.all(held == held[0])
        assert neutral[at_hub][:, ~at_hub].nnz == 0
        # Without the connector the interior and exterior supports are disjoint.
        # The exact basis may retain several total head sectors, but no legal
        # operator on one factor depends on the other factor's word or load.
        return {
            "mobile": mobile,
            "neutral": neutral,
            "mobile_cut": mobile - h["connector"],
            "neutral_cut": neutral - h["connector"],
        }, float(held[0])

    def verify_probe(self, board):
        """The separate material meter includes its native degree detuning."""
        degrees = [int(board.degree[port]) for port in self.probe]
        assert all(d in (5, 6) for d in degrees)
        for i, (_, position, word) in enumerate(self.basis):
            expected = (position == self.probe[0]) == bool(word & self.probe_mask)
            row = self.terms["probe_hop"].getrow(i)
            if expected == (self.probe_sign == 1):
                target = self.probe[1] if position == self.probe[0] else self.probe[0]
                partner = self.ids[self.basis[i][0], target, word ^ self.probe_mask]
                assert row.nnz == 2 and row[0, i] == 4 and row[0, partner] == -4
            else:
                assert row.nnz == 0
        return degrees
