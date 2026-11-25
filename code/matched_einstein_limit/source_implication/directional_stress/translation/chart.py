"""Find short, oriented translations of actual decorated cell neighborhoods."""

from collections import Counter
from functools import lru_cache

import networkx as nx

from model.executable_model.carrier import cell, federation
from model.executable_model.symmetry import rotations


@lru_cache(maxsize=1)
def relabelings():
    # Include orientation-reversing relabelings of the microscopic port graph.
    # The shift itself preserves the oriented cyclic order of macro-cell edges.
    graph = nx.Graph(cell().seams)
    antipode = {
        i: max(
            nx.single_source_shortest_path_length(graph, i),
            key=nx.single_source_shortest_path_length(graph, i).get,
        )
        for i in range(12)
    }
    perms = list(rotations())
    perms += [tuple(antipode[v] for v in p) for p in perms]
    assert len(set(perms)) == 120
    return sorted(perms)


@lru_cache(maxsize=None)
def match(src, dst):
    return next((p for p in relabelings() if tuple(p[i] for i in src) == dst), None)


class Atlas:
    def __init__(self, level):
        self.board = federation(level)
        self.pin = {}
        for a, b in self.board.external:
            self.pin[a // 12, b // 12] = a % 12
            self.pin[b // 12, a // 12] = b % 12
        self.graph = nx.Graph(list(self.pin))
        good, embedding = nx.check_planarity(self.graph)
        assert good
        self.order = {i: list(embedding.neighbors_cw_order(i)) for i in self.graph}
        assert all(len(v) == 3 for v in self.order.values())
        self.dist = {
            i: dict(d) for i, d in nx.all_pairs_shortest_path_length(self.graph)
        }

    def extension(self, src, dst, known, nodes, labels, inside):
        """Map all three incident edges at this pair, using a registered first edge."""
        v, w = known
        a, b = self.order[src], self.order[dst]
        offset = b.index(w) - a.index(v)
        images = [b[(i + offset) % 3] for i in range(3)]
        perm = match(
            tuple(self.pin[src, j] for j in a), tuple(self.pin[dst, j] for j in images)
        )
        if perm is None:
            return False
        for u, t in zip(a, images):
            if (
                self.dist[u][t] != 2
                or (u in nodes and nodes[u] != t)
                or (t in inside and inside[t] != u)
            ):
                return False
        labels[src] = perm
        for u, t in zip(a, images):
            nodes[u] = t
            inside[t] = u
        return True

    def balls(self, src, dst, k):
        """A fixed direction of a two-macro-hop shift. Return its witnessed ball radii."""
        if self.dist[src][dst] != 2:
            return []
        nodes, inside, labels = {src: dst}, {dst: src}, {}
        witness = []
        queue = [(src, self.order[src][0], self.order[dst][k])]
        for radius in range(8):
            next_ring = []
            for u, ref, img in queue:
                if not self.extension(u, nodes[u], (ref, img), nodes, labels, inside):
                    return witness
                next_ring += [
                    (v, u, nodes[u])
                    for v in self.order[u]
                    if self.dist[src][v] == radius + 1
                ]
            witness.append(
                dict(
                    radius=radius,
                    center=src,
                    target=dst,
                    orientation=k,
                    cells=sorted(labels),
                    images=[nodes[i] for i in sorted(labels)],
                    lifts=[labels[i] for i in sorted(labels)],
                )
            )
            unique = {}
            for row in next_ring:
                unique[row[0]] = row
            queue = list(unique.values())
            if not queue:
                return witness
        return witness

    def patch(self, src, dst, k):
        """Connected extension, retaining matched vertices and cutting failed continuations."""
        nodes, inside, labels = {src: dst}, {dst: src}, {}
        queue = [(src, self.order[src][0], self.order[dst][k])]
        tried = set()
        for u, ref, img in queue:
            if u in tried:
                continue
            tried.add(u)
            if not self.extension(u, nodes[u], (ref, img), nodes, labels, inside):
                continue
            queue.extend(
                (v, u, nodes[u]) for v in sorted(self.order[u]) if v not in tried
            )
        return dict(
            center=src,
            target=dst,
            orientation=k,
            cells=sorted(labels),
            images=[nodes[i] for i in sorted(labels)],
            lifts=[labels[i] for i in sorted(labels)],
        )

    def census(self):
        counts = Counter()
        best = {}
        for src in sorted(self.graph):
            for dst, d in self.dist[src].items():
                if d != 2:
                    continue
                for k in range(3):
                    for rec in self.balls(src, dst, k):
                        radius = rec["radius"]
                        counts[radius] += 1
                        best.setdefault(radius, rec)
        return counts, best
