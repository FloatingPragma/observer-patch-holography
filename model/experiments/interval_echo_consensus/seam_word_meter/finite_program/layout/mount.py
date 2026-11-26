"""Pack binary shelves and route native endpoint exchanges on labelled cells."""

from collections import deque
from collections.abc import Sequence
from functools import cache
from itertools import combinations
from operator import index

import numpy as np
from model.executable_model import carrier
from ...chart import exchange_word
from ..machine import Bench, digit_instruments

TEMPLATE = tuple(
    tuple(b if a == i else a for a, b in carrier.CELL_SEAMS if i in (a, b))
    for i in range(12)
)
PORTS = frozenset(range(12))


@cache
def matching(vertices):
    if not vertices:
        return ()
    a = min(vertices)
    for b in TEMPLATE[a]:
        if b in vertices:
            rest = matching(vertices - {a, b})
            if rest is not None:
                return ((a, b),) + rest
    return None


@cache
def station(writer, gateways):
    """Choose a head and two adjacent references avoiding all three gateways."""
    available = PORTS - frozenset(gateways)
    for head in sorted(available):
        for zero, two in combinations(sorted(set(TEMPLATE[head]) & available), 2):
            data = matching(PORTS - {head, zero, two, writer})
            if data is not None:
                return head, (zero, two), data
    raise ValueError("No disjoint station and eight-port matching")


def local_path(start, end, blocked):
    parent = {start: None}
    todo = deque([start])
    while todo:
        cur = todo.popleft()
        if cur == end:
            result = []
            while cur is not None:
                result.append(cur)
                cur = parent[cur]
            return tuple(reversed(result))
        for nxt in TEMPLATE[cur]:
            if nxt not in parent and nxt not in blocked:
                parent[nxt] = cur
                todo.append(nxt)
    raise ValueError("Blocked local station route")


class Shelves(Sequence):
    def __init__(self, mount, count, width):
        self.mount, self.count, self.width = mount, count, width

    def __len__(self):
        return self.count

    def __getitem__(self, i):
        if isinstance(i, slice):
            return tuple(self[j] for j in range(*i.indices(self.count)))
        i = index(i)
        if i < 0:
            i += self.count
        if not 0 <= i < self.count:
            raise IndexError(i)
        start = 2 + i * self.width
        return tuple(self.mount.dimer(j) for j in range(start, start + self.width))


class Mount:
    """Stores the three gateways per owner and an owner spanning tree."""

    def __init__(self, board):
        self.cells = int(board.cells)
        if self.cells < 2:
            raise ValueError("The two stations require distinct owners")
        if isinstance(board, carrier.Carrier):
            inner = [
                (a // 12, a % 12, b % 12) for a, b in board.seams if a // 12 == b // 12
            ]
            template = set(carrier.CELL_SEAMS)
            if (
                board.n != 12 * self.cells
                or len(inner) != 30 * self.cells
                or len(set(inner)) != 30 * self.cells
                or any((a, b) not in template for _, a, b in inner)
            ):
                raise ValueError("Incomplete labelled internal seams")
        elif tuple(getattr(board, "cell_seams", ())) != carrier.CELL_SEAMS:
            raise ValueError(
                "The supplied board needs the labelled twelve-port template"
            )
        edges = np.asarray(board.external)
        c = self.cells
        if (
            edges.shape != (3 * c // 2, 2)
            or not np.issubdtype(edges.dtype, np.integer)
            or np.any(edges < 0)
            or np.any(edges >= 12 * c)
            or np.any(edges[:, 0] // 12 == edges[:, 1] // 12)
            or np.unique(edges).size != edges.size
        ):
            raise ValueError("Invalid labelled external three-port matching")
        if not np.all(np.bincount((edges // 12).ravel(), minlength=c) == 3):
            raise ValueError("Each owner needs three different gateways")
        both = np.concatenate((edges, edges[:, ::-1]))
        self.incident = both[np.argsort(both[:, 0], kind="stable")].reshape(c, 3, 2)
        self.parent = np.full(c, -1, dtype=np.int64)
        self.depth = np.zeros(c, dtype=np.int64)
        self.up = np.full((c, 2), -1, dtype=np.int64)
        self.parent[0] = 0
        todo = deque([0])
        while todo:
            a = todo.popleft()
            for start, end in self.incident[a]:
                b = int(end // 12)
                if self.parent[b] == -1:
                    self.parent[b] = a
                    self.depth[b] = self.depth[a] + 1
                    self.up[b] = (end, start)
                    todo.append(b)
        if np.any(self.parent < 0):
            raise ValueError("Disconnected owner graph")
        outside = (edges[:, 0] != self.up[edges[:, 0] // 12, 0]) & (
            edges[:, 1] != self.up[edges[:, 1] // 12, 0]
        )
        if not np.any(outside):
            raise ValueError("No non-tree writer available")
        candidates = np.sort(edges[outside], axis=1)
        self.writer = tuple(map(int, candidates[np.argmin(candidates[:, 0])]))
        self.heads, self.references, self.special = [], [], {}
        for port in self.writer:
            owner, writer = divmod(port, 12)
            gates = tuple(int(p % 12) for p in self.incident[owner, :, 0])
            head, refs, data = station(writer, gates)
            self.heads.append(12 * owner + head)
            self.references.append(tuple(12 * owner + j for j in refs))
            self.special[owner] = data
        self.heads, self.references = tuple(self.heads), tuple(self.references)
        self.protected = frozenset(
            (*self.writer, *self.heads, *self.references[0], *self.references[1])
        )
        lengths = np.full(c, 6, dtype=np.int64)
        for owner in self.special:
            lengths[owner] = 4
        self.stops = np.cumsum(lengths)
        self.capacity = int(self.stops[-1]) - 2
        self.flags = tuple(self.dimer(j) for j in (0, 1))

    def number_of_nodes(self):
        return 12 * self.cells

    def has_edge(self, a, b):
        i, x = divmod(a, 12)
        j, y = divmod(b, 12)
        if i == j:
            return y in TEMPLATE[x]
        return any(int(u) == a and int(v) == b for u, v in self.incident[i])

    def dimer(self, j):
        """Index into all nonservice pairs, with pairs zero and one reserved as flags."""
        if not 0 <= j < self.stops[-1]:
            raise IndexError(j)
        owner = int(np.searchsorted(self.stops, j, side="right"))
        start = int(self.stops[owner - 1]) if owner else 0
        a, b = self.special.get(owner, matching(PORTS))[j - start]
        return 12 * owner + a, 12 * owner + b

    def path(self, source, end):
        if not (
            0 <= source < self.number_of_nodes() and 0 <= end < self.number_of_nodes()
        ):
            raise ValueError("Port outside this board")
        a, b = source // 12, end // 12
        left, right = [], []
        while a != b:
            if self.depth[a] >= self.depth[b]:
                left.append(tuple(map(int, self.up[a])))
                a = int(self.parent[a])
            else:
                right.append(tuple(map(int, self.up[b][::-1])))
                b = int(self.parent[b])
        blocked = self.protected - {source, end}
        result = [source]

        def within(target):
            current = result[-1]
            owner = current // 12
            if target // 12 != owner:
                raise ValueError("Invalid intra-owner route")
            forbidden = {i % 12 for i in blocked if i // 12 == owner}
            result.extend(
                12 * owner + j
                for j in local_path(current % 12, target % 12, forbidden)[1:]
            )

        for p, q in (*left, *reversed(right)):
            within(p)
            result.append(q)
        within(end)
        return tuple(result)

    def bench(self, *, registers, axes, width):
        if (
            any(type(n) is not int or n < 0 for n in (registers, axes, width))
            or width == 0
        ):
            raise ValueError("Integer dimensions and positive material width required")
        if registers * axes * width > self.capacity:
            raise ValueError("Insufficient data dimers")
        return RoutedBench(self, registers * axes, width)


class RoutedBench(Bench):
    """Use the existing native write, both reads and all stores with compact routing."""

    def __init__(self, mount, count, width):
        self.graph = mount
        self.words = Shelves(mount, count, width)
        self.width = width
        self.writer, self.heads, self.references = (
            mount.writer,
            mount.heads,
            mount.references,
        )
        self.flags = mount.flags
        self.read_table, self.write_table = digit_instruments()

    def route(self, source, end):
        return exchange_word(self.graph.path(source, end))
