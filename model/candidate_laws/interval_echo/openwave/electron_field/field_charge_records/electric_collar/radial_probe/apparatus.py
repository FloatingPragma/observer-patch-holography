"""Construct a reciprocal radial probe on the four-fan receiver."""

from collections import defaultdict
from itertools import pairwise, product

import networkx as nx
import numpy as np
from scipy import sparse

from ..apparatus import Collar


def seams(triangle):
    return {tuple(sorted(edge)) for edge in pairwise((*triangle, triangle[0]))}


def geometry(path):
    fan = path.cells[-1]
    owner = fan[0] // 12
    vertices = set(range(owner * 12, (owner + 1) * 12))
    original = {tuple(sorted((fan[a], fan[b]))) for a, b in path.field.edges}
    outside = vertices - set(fan)
    anchor = tuple(sorted(fan[i] for i in path.field.rings[0] if i != 0))
    faces = sorted(
        {
            tuple(sorted(f))
            for f in path.board.faces()
            if set(f) <= vertices and set(f) & outside
        }
    )
    graph = nx.Graph(path.board.seams)
    distance = nx.single_source_shortest_path_length(graph, fan[0])
    chains = [(f,) for f in faces if anchor in seams(f)]
    while chains:
        choices = []
        for chain in chains:
            for a, b in seams(chain[-1]):
                if (
                    {a, b} <= outside
                    and distance[a] != distance[b]
                    and path.board.degree[a] == path.board.degree[b]
                ):
                    probe = tuple(sorted((a, b), key=distance.get))
                    choices.append((chain, probe))
        if choices:
            chain, probe = min(choices)
            return chain, probe, [distance[p] for p in probe], len(choices)
        chains = [
            (*chain, face)
            for chain in chains
            for face in faces
            if face not in chain and (seams(face) & seams(chain[-1])) - original
        ]
    raise ValueError("No equal-degree radial seam at the chosen face")


class RadialCollar(Collar):
    def __init__(self, path, sign, pair):
        self.source_sign = path.sign
        self.probe_sign = sign
        self.fan = path.cells[-1]
        self.hub = self.fan[0]
        self.faces, self.radial, self.ranges, self.eligible = geometry(path)
        self.probe = tuple(sorted(self.radial))
        self.original = [(self.fan[a], self.fan[b]) for a, b in path.field.edges]
        self.inside_bits = len(self.original)
        self.inside_mask = (1 << self.inside_bits) - 1
        originals = {tuple(sorted(edge)) for edge in self.original}
        self.added = sorted(
            set().union(*(seams(face) for face in self.faces)) - originals
        )
        self.edges = self.original + self.added
        self.edge_index = {frozenset(edge): j for j, edge in enumerate(self.edges)}
        self.source_moves = defaultdict(list)
        for bit, (a, b) in enumerate(self.original):
            self.source_moves[a].append((b, 1 << bit, True))
            self.source_moves[b].append((a, 1 << bit, False))
        self.probe_mask = 1 << self.edge_index[frozenset(self.probe)]
        triangles = [tuple(self.fan[p] for p in tri) for tri in path.field.rings]
        self.rings = []
        for index, triangle in enumerate([*triangles, *self.faces]):
            mask = clockwise = 0
            for p, q in pairwise((*triangle, triangle[0])):
                bit = self.edge_index[frozenset((p, q))]
                mask |= 1 << bit
                clockwise |= int(self.edges[bit][0] == p) << bit
            if index < len(triangles):
                kind = "fan_rings"
                assert not mask & ~self.inside_mask
            elif index == len(triangles):
                kind = "connector"
                assert mask & self.inside_mask and mask & ~self.inside_mask
            else:
                kind = "outer_ring"
                assert not mask & self.inside_mask
            self.rings.append((kind, mask, clockwise))
        self.backgrounds = 1 << len(self.added)
        self.start_words = pair
        starts = [
            (self.hub, location, w | (background << self.inside_bits))
            for location, background, w in product(
                self.radial, range(self.backgrounds), pair
            )
        ]
        self.starts = len(starts)
        self.basis = list(starts)
        self.ids = {s: i for i, s in enumerate(starts)}
        self.terms = self.build(path.board.degree)
        self.probe_degrees = self.verify_probe(path.board)
        assert self.probe_degrees[0] == self.probe_degrees[1]
        assert self.ranges == [2, 3]

    def static_difference(self, hamiltonians, scalar):
        """Clamp both material locations and retain all native diagonals."""
        for port in self.radial:
            positions = np.array(
                [s == self.hub and x == port for s, x, _ in self.basis]
            )
            a = hamiltonians["mobile"][positions][:, positions]
            b = hamiltonians["neutral"][positions][:, positions]
            assert (a - b - scalar * sparse.eye(a.shape[0])).nnz == 0
        return 0
