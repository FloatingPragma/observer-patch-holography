"""Parallel truth-table columns of the sharp physical rail experiment."""

from functools import cache

from model.executable_model import CELL_SEAMS

from ...boolean_nodes import material

PORTS = 12
WORDS = 1 << PORTS
MARKS = min(tuple(sorted(edge)) for edge in CELL_SEAMS)
AUDITORS = (18, 19)


def transpositions(mount):
    if sorted(mount) != list(range(PORTS)):
        raise ValueError("Invalid mount permutation")
    seen = set()
    for root in range(PORTS):
        if root in seen:
            continue
        seen.add(root)
        other = mount[root]
        while other != root:
            yield root, other
            seen.add(other)
            other = mount[other]


class Board:
    def __init__(self):
        self.device = material.Machine(0, 20)
        self.slots = [pair for row in self.device.rows for pair in row[:2]][:36]
        if len({r for pair in self.slots for r in pair}) != 72:
            raise ValueError("Overlapping physical pairs")

    def load(self, words):
        # Bit k is one sharp ternary preparation, with loads 0 and 2 encoded 0 and 1.
        full = (1 << len(words)) - 1
        planes = [
            sum(((word >> i) & 1) << k for k, word in enumerate(words))
            for i in range(PORTS)
        ]
        planes.extend(full if i == mark else 0 for mark in MARKS for i in range(PORTS))
        loads = [0] * self.device.car.n
        for row in self.device.rows:
            loads[row[2][0]] = full
        for value, (a, b) in zip(planes, self.slots, strict=True):
            loads[a], loads[b] = value, full ^ value
        return loads

    @cache
    def program(self, mount, parked):
        pulses = []
        for left, right in transpositions(mount):
            for offset in (0, PORTS, 2 * PORTS):
                for start, end in zip(
                    self.slots[offset + left], self.slots[offset + right], strict=True
                ):
                    pulses.extend(self.device.move(start, end))
        return tuple((a, b) for a, b in pulses if not parked or a // 12 == b // 12)

    @cache
    def query(self, source, station):
        reference = self.device.rows[station][2][0]
        a, b = self.device.rows[station][3]
        route = self.device.move(source, a) + self.device.move(
            source if reference == a else reference, b
        )
        return route, a, b

    @staticmethod
    def swaps(loads, steps):
        for a, b in steps:
            loads[a], loads[b] = loads[b], loads[a]

    def read(self, loads, station, mode, full, pairs=None):
        # Compile the measured single-seam table once. OR disjoint truth-table columns.
        response = material.balanced_pointer(mode)
        raw = []
        for dim in self.slots if pairs is None else pairs:
            pair = []
            for rail in dim:
                route, a, b = self.query(rail, station)
                self.swaps(loads, route)
                x, y = loads[a], loads[b]
                output = 0
                for u, v in response:
                    if response[u, v]:
                        output |= (x if u else x ^ full) & (y if v else y ^ full)
                self.swaps(loads, reversed(route))
                pair.append(output)
            raw.append(pair)
        return raw

    def acquire(self, initial, mount, width, parked=False, mode="native"):
        loads = initial.copy()
        program = self.program(mount, parked)
        self.swaps(loads, program)
        before = loads.copy()
        raw = [self.read(loads, at, mode, (1 << width) - 1) for at in AUDITORS]
        if loads != before or raw[0] != raw[1]:
            raise ValueError("Reader disagreement or destructive echo")
        return raw, len(program)
