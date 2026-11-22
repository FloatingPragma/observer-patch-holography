"""Wedge reconstruction. The only inputs are a simple screen and paired public hits."""

from dataclasses import dataclass
import heapq
from itertools import combinations


def need(ok, message):
    if not ok:
        raise ValueError(message)


class Screen:
    def __init__(self, n, edges):
        need(type(n) is int and n >= 0, "invalid port count")
        self.n = n
        self.edges = []
        incident = [[] for _ in range(n)]
        for i, edge in enumerate(edges):
            need(isinstance(edge, (list, tuple)) and len(edge) == 2, "invalid edge")
            a, b = edge
            need(type(a) is int and type(b) is int and 0 <= a < b < n, "invalid endpoints")
            self.edges.append((a, b))
            incident[a].append(i)
            incident[b].append(i)
        need(len(self.edges) == len(set(self.edges)), "parallel address")
        self.wedges = {}
        for center, items in enumerate(incident):
            for i, j in combinations(items, 2):
                leaves = sorted((set(self.edges[i]) | set(self.edges[j])) - {center})
                self.wedges[i, j] = (center, *leaves)


@dataclass
class Recovered:
    edges: list
    higher: list
    order: list


def reconstruct(screen, rows):
    need(isinstance(rows, list), "invalid transcript")
    remaining = set(screen.wedges)
    parents = [set() for _ in screen.edges]
    children = [[] for _ in screen.edges]
    for row in rows:
        need(isinstance(row, list) and len(row) == 4 and all(type(x) is int for x in row),
             "invalid transcript row")
        i, j, a, b = row
        need((i, j) in screen.wedges, "nonincident or noncanonical wedge")
        need((i, j) in remaining, "duplicate wedge")
        remaining.remove((i, j))
        need((a, b) in ((0, 1), (1, 0)), "not exactly one addressed hit")
        winner, loser = (i, j) if a else (j, i)
        parents[loser].add(winner)
        children[winner].append(loser)
    need(not remaining, "missing wedge")

    # Among simultaneously available addresses, take boundary seams first.
    def tie(i):
        u, v = screen.edges[i]
        return (u // 12 == v // 12, i)

    counts = list(map(len, parents))
    pending = [tie(i) for i, n in enumerate(counts) if n == 0]
    heapq.heapify(pending)
    order = []
    while pending:
        _, i = heapq.heappop(pending)
        order.append(i)
        for j in children[i]:
            counts[j] -= 1
            if counts[j] == 0:
                heapq.heappush(pending, tie(j))
    need(len(order) == len(screen.edges), "cyclic priority observations")
    return Recovered(screen.edges, [sorted(s) for s in parents], order)


def separator(screen, left, right):
    for i, j in screen.wedges:
        if (i in left.higher[j]) != (i in right.higher[j]):
            center, a, b = screen.wedges[i, j]
            return dict(pair=[i, j], double=center, vacancies=[a, b],
                        default_load=1, initial_collars="blank", read=i)
    return None
