"""Two disjoint contact calendars with fully expanded physical addresses."""

from functools import cache
from itertools import combinations

import networkx as nx
import numpy as np

from ..seam_word_meter import apparatus as seam
from ..seam_word_meter.chart import exchange_word
from ..seam_word_meter.finite_program.machine import digit_instruments
from .compiler import compile_table


@cache
def calibrated():
    read, written = digit_instruments()
    echo = seam.LAW.echo()
    bright = np.array([ab == (1, 1) for ab in seam.LAW.SRC] * (seam.LAW.DIM // 9))
    for value in range(3):
        j = seam.LAW.SRC.index((value, 0))
        projected = echo[:, j] * (bright == bool(read[value, 0]))
        np.testing.assert_allclose(
            echo.conj().T @ projected, np.eye(seam.LAW.DIM)[:, j], atol=seam.TOL, rtol=0
        )
    for (a, b), (out, _) in written.items():
        j = seam.LAW.SRC.index((a, b))
        vector = echo[:, j]
        selected = np.tile([ab == out for ab in seam.LAW.SRC], (seam.LAW.DIM // 9))
        np.testing.assert_allclose(
            echo.conj().T @ (vector * selected),
            np.eye(seam.LAW.DIM)[:, j],
            atol=seam.TOL,
            rtol=0,
        )
        if read[out[0], 0] != int(a == b == 2):
            raise ValueError("The native joint detector does not implement conjunction")
    return {v: read[v, 0] for v in range(3)}, {k: v[0] for k, v in written.items()}


def domains(board):
    owners = nx.Graph((a // 12, b // 12) for a, b in board.external)
    triples = [
        g
        for g in combinations(range(board.cells), 3)
        if sum(owners.has_edge(a, b) for a, b in combinations(g, 2)) >= 2
    ]
    first = triples[0]
    return first, next(g for g in triples if set(g).isdisjoint(first))


class Lane:
    def __init__(self, board, owners):
        self.owners = owners
        ports = {12 * c + j for c in owners for j in range(12)}
        self.graph = nx.Graph(
            (a, b) for a, b in board.seams if a in ports and b in ports
        )
        ends = [min(owners), max(owners)]
        self.heads, self.refs = [], []
        for owner in ends:
            # Reserve local service ports away from the three gateway ports.
            local = range(12 * owner, 12 * owner + 12)
            free = {p for p in local if board.degree[p] == 5}
            h, neighbours = next(
                (p, sorted(set(self.graph[p]) & free))
                for p in sorted(free)
                if len(set(self.graph[p]) & free) >= 3
            )
            self.heads.append(h)
            self.refs.append(neighbours[0])
            if owner == ends[0]:
                self.writer = (h, neighbours[1])
        self.service = set(self.writer) | set(self.heads) | set(self.refs)
        data = self.graph.subgraph(ports - self.service)
        if not nx.is_connected(data):
            raise ValueError("Blocked transit lane")
        # Unit edge weights are uniquely tie-broken by their binary ranks.
        matching = nx.Graph()
        for j, (a, b) in enumerate(sorted(data.edges)):
            if a // 12 == b // 12:
                matching.add_edge(a, b, weight=1 << j)
        self.pairs = tuple(
            sorted(
                tuple(sorted(edge))
                for edge in nx.max_weight_matching(matching, maxcardinality=True)
            )
        )
        if len(self.pairs) < 9:
            raise ValueError("Insufficient isolated data dimers")
        self.word = []
        self.reads = []

    @cache
    def route(self, source, target):
        kept = set(self.graph) - (self.service - {source, target})
        path = nx.shortest_path(self.graph.subgraph(kept), source, target)
        return tuple(exchange_word(path))

    def moves(self, edges):
        self.word.extend(("swap", a, b) for a, b in edges)

    def probe(self, pair):
        for name in ("probe_write", "project", "probe_undo"):
            self.word.append((name, *pair))

    def gate(self, literals, dest):
        addresses = [r for r, _ in literals] + [dest]
        if (
            len(literals) > 2
            or len(set(addresses)) != len(addresses)
            or any(
                type(r) is not int or r not in range(len(self.pairs)) for r in addresses
            )
            or any(b not in (0, 1) for _, b in literals)
        ):
            raise ValueError(
                "A primitive needs distinct enrolled control and target addresses"
            )
        if not literals:
            self.word.append(("flip", *self.pairs[dest]))
            return
        ports = [self.pairs[r][1 - value] for r, value in literals]
        restore = []
        for p, target in zip(ports, self.writer, strict=False):
            route = self.route(p, target)
            self.moves(route)
            restore.extend(route)
        if len(ports) == 2:
            self.word.append(("write", *self.writer))
        self.probe((self.heads[0], self.refs[0]))
        if len(ports) == 2:
            self.word.append(("unwrite", *self.writer))
        self.moves(reversed(restore))
        self.word.append(("shutter", *self.pairs[dest]))

    def table(self, sources, target, truth, scratch=(8,)):
        for args, dest in compile_table(sources, (target,), truth, scratch):
            self.gate(args, dest)

    def survey(self, label, registers):
        for at, (head, ref) in enumerate(zip(self.heads, self.refs, strict=True)):
            for reg in registers:
                where = []
                for rail in self.pairs[reg]:
                    steps = self.route(rail, head)
                    self.moves(steps)
                    self.probe((head, ref))
                    where.append(len(self.word) - 2)
                    self.moves(reversed(steps))
                self.reads.append((label, at, reg, *where))

    def program(self, stale=False):
        copy = (0, 1)
        self.survey("start", range(9))
        self.table((0, 1, 2), 3, tuple(int(i.bit_count() >= 2) for i in range(8)))
        self.table((3,), 4, copy)
        self.survey("cached", (3, 4, 8))
        self.word.append(("intervene", *self.pairs[3]))
        self.survey("edited", (3, 4, 8))
        self.table((3,), 5, copy)
        self.gate((), 6)
        self.survey("sealed", (3, 4, 5, 6))
        self.word.append(("wait", self.heads[0], self.heads[0]))
        self.table((4 if stale else 5, 2), 7, (0, 1, 1, 0))
        self.survey("end", range(9))
        return self

    def address_map(self):
        return dict(
            owners=self.owners,
            heads=self.heads,
            zeros=self.refs,
            writer=self.writer,
            pairs=self.pairs,
            word=self.word,
            reads=self.reads,
        )


class LocalState:
    def __init__(self, parked):
        self.read, self.written = calibrated()
        self.parked = parked
        self.saved = self.pending = None
        self.lamp = None

    def step(self, loads, op, edit):
        name, a, b = op
        if name == "write":
            if self.saved is not None:
                raise ArithmeticError("Occupation overwrite")
            self.saved = (loads[a], loads[b])
            loads[a], loads[b] = self.written[self.saved]
        elif name == "unwrite":
            if (loads[a], loads[b]) != self.written[self.saved]:
                raise ArithmeticError("Missing inverse write")
            loads[a], loads[b] = self.saved
            self.saved = None
        elif name == "probe_write":
            if self.pending is not None or loads[b] != 0:
                raise ArithmeticError("Probe or zero standard invalid")
            self.pending = (a, b, loads[a])
        elif name == "project":
            x, y, v = self.pending
            if (a, b) != (x, y):
                raise ArithmeticError("Pointer selected the wrong contact")
            self.lamp = 0 if self.parked and self.saved is not None else self.read[v]
            return self.lamp
        elif name == "probe_undo":
            if self.pending != (a, b, loads[a]):
                raise ArithmeticError("Invalid projected unwrite")
            self.pending = None
        elif name in ("swap", "flip", "shutter", "intervene"):
            enable = (
                self.lamp if name == "shutter" else edit if name == "intervene" else 1
            )
            if enable:
                loads[a], loads[b] = loads[b], loads[a]
            if name in ("shutter", "intervene"):
                return enable
        elif name != "wait":
            raise ValueError(name)
        return None
