"""Sharp local records and addressed transport on a finite Interval Echo board."""

from collections import deque
from functools import cache

import numpy as np

from model.candidate_laws.interval_echo.clock import DIGITS, EXTRA_DIGITS, Registers
from model.executable_model import federation
from ..measure import rails, tree

PORTS = (
    "north",
    "south",
    *(f"upper_{i}" for i in range(5)),
    *(f"lower_{i}" for i in range(5)),
)
DATA, INBOX, OUTBOX = DIGITS[:1], DIGITS[1:], EXTRA_DIGITS

# In each five-cycle the successive transpositions with its first slot
# advance the initial item by one place around the cycle.
ROTATE = [(2, i) for i in range(3, 7)] + [(7, i) for i in range(8, 12)]
REFLECT = [(3, 6), (4, 5), (7, 11), (8, 10)]
ANTIPODE = [(0, 1)] + [(2 + i, 7 + (i + 2) % 5) for i in range(5)]


def primary_read(word):
    """The sole primary instrument requires one valid digit at every named station."""
    if len(word) != 12 or any(v not in (0, 1) for v in word):
        return None
    count = sum(word)
    if count == 1:
        return PORTS[word.index(1)] + "/write"
    if count == 11:
        return PORTS[word.index(0)] + "/check"
    return None


class Checkpoint:
    def __init__(self, level):
        self.board = federation(level)
        self.meter = Registers(self.board, 0)
        self.edges = tree(self.board)
        self.owners = tuple(child for _, child in self.edges[:12])
        self.cells = frozenset(range(self.board.cells))
        self.addresses = {
            (cell, register): rails(cell, register)
            for cell in self.cells
            for register in (DATA, INBOX, OUTBOX)
        }
        self._route = cache(self._compile_route)
        self._data_rails = np.array(
            [self.addresses[c, DATA][0][:2] for c in self.owners]
        )
        self._idle = np.ones(self.board.n, dtype=bool)
        self._idle[self._data_rails.ravel()] = False

    def prepare(self, label):
        self.meter.loads.fill(0)
        self.meter.gates.clear()
        port, orientation = label.split("/")
        for cell in self.cells:
            self.meter.initialize_digits(rails(cell, DATA + INBOX + OUTBOX))
        for name, cell in zip(PORTS, self.owners, strict=True):
            self.put(cell, DATA, int((name == port) != (orientation == "check")))

    def put(self, cell, register, value):
        self.meter.put(self.addresses[cell, register], value, {cell})

    def get(self, cell, register):
        return self.meter.get(self.addresses[cell, register], {cell})

    def reads(self):
        # An impaired continuation can move a reference rail. Retain the actual
        # local echo in that control, as at the calibrated checkpoint.
        outcomes = []
        for cell in self.owners:
            first, _, reference = self.addresses[cell, DATA][0]
            outcomes.append(int(self.meter.flag(first, reference, {cell})))
        return tuple(outcomes)

    def idle_snapshot(self):
        return self.meter.loads[self._idle].copy()

    def audit_closure(self, baseline):
        # Internal check of the subspace used by the concatenation proof. The
        # separately recorded detector words above never consume this check.
        rails = self.meter.loads[self._data_rails]
        if not np.array_equal(baseline, self.meter.loads[self._idle]) or not np.all(
            rails.sum(axis=1) == 2
        ):
            raise ValueError(
                "the addressed continuation left the prepared sharp workspace"
            )

    def _compile_route(self, left, right, region):
        return [
            (a, b, {a // 12, b // 12}) for a, b in self.meter.route(left, right, region)
        ]

    def exchange(self, left, right, region, parked=False):
        for a, b, owners in self._route(left, right, region):
            if not parked or len(owners) == 1:
                self.meter.swap(a, b, owners)

    def apply(self, setting, parked=False):
        r, s, a, f = (int(field[1:]) for field in setting.split("_"))
        for left, right in REFLECT * s + ROTATE * r + ANTIPODE * a:
            x, y = self.owners[left], self.owners[right]
            for rail in range(2):
                self.exchange(
                    self.addresses[x, DATA][0][rail],
                    self.addresses[y, DATA][0][rail],
                    self.cells,
                    parked,
                )
        if f:
            for cell in self.owners:
                first, second, _ = self.addresses[cell, DATA][0]
                self.meter.swap(first, second, {cell})

    def rooted(self, root):
        neighbors = {cell: [] for cell in self.cells}
        for left, right in self.edges:
            neighbors[left].append(right)
            neighbors[right].append(left)
        seen, queue, edges = {root}, deque([root]), []
        while queue:
            left = queue.popleft()
            for right in sorted(neighbors[left]):
                if right in seen:
                    continue
                seen.add(right)
                queue.append(right)
                edges.append((left, right))
        return edges

    def relay(self, omit=False):
        """A fresh local send word carries the chosen station's measured digit at each round."""
        received = {cell: [] for cell in self.owners}
        for turn, sender in enumerate(self.owners):
            for cell in self.cells:
                self.put(cell, INBOX, 0)
                self.put(cell, OUTBOX, 0)
            self.put(sender, INBOX, self.get(sender, DATA) + 1)
            for index, (parent, child) in enumerate(self.rooted(sender)):
                payload = self.get(parent, INBOX)
                if payload == 0 or (omit and turn == index == 0):
                    continue
                self.put(parent, OUTBOX, payload)
                for p, q in zip(
                    self.addresses[parent, OUTBOX],
                    self.addresses[child, INBOX],
                    strict=True,
                ):
                    for i in range(2):
                        self.exchange(p[i], q[i], frozenset((parent, child)))
            for cell in received:
                value = self.get(cell, INBOX)
                received[cell].append(None if value == 0 else value - 1)
        return list(received.values())
