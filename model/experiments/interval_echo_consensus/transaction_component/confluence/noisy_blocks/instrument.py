"""Packed sharp branches on the count-addressed material triangles and a reset dimer."""

from functools import cache
from itertools import product
from math import pi

import numpy as np

from ....boolean_nodes.material import balanced_pointer
from ....seam_word_meter import apparatus as seam
from ..apparatus import COMMIT, FAR, MASKS, hardware, layout

NOMINAL = tuple((side, m) for m in range(4) for side in range(2))
BOUNDARY = (0, 1)


@cache
def pulse_calibration():
    start = seam.LAW.SRC.index((2, 0))
    moved = seam.LAW.SRC.index((0, 2))
    state = seam.pulse(pi / 3)[:, start]
    for j in range(9):
        np.testing.assert_allclose(
            abs(state[j]) ** 2, {start: 0.75, moved: 0.25}.get(j, 0), atol=seam.TOL
        )
    for j in (start, moved):
        projected = np.zeros(9, complex)
        projected[j] = state[j]
        restored = seam.pulse(pi) @ projected if j == moved else projected
        np.testing.assert_allclose(np.delete(restored, start), 0, atol=seam.TOL)
    for (a, b), click in balanced_pointer("native").items():
        if click != int(a != b):
            raise ArithmeticError("Calibrated extremal XOR pointer changed")
    return {
        "clock": list(hardware().rows[14][0]),
        "latch": [list(hardware().rows[14][1]), list(hardware().rows[13][0])],
        "source": [2, 0],
        "stay_cross_weights": [3, 1],
        "denominator": 4,
        "mode_weights": [12, 3, 1],
        "mode_names": ["correct", "hold", "wrong"],
    }


def packed(bits):
    """Bit j denotes retained sharp leaf j, in lexicographic order."""
    return int.from_bytes(
        np.packbits(np.asarray(bits, dtype=np.uint8), bitorder="little").tobytes(),
        "little",
    )


def classes(bits, full):
    result = [0] * 4
    for choices in product((0, 1), repeat=3):
        included = full
        for bit, value in zip(bits, choices, strict=True):
            included &= bit if value else full ^ bit
        result[sum(choices)] |= included
    return result


class Bench:
    def __init__(self, leaves, initial, hidden=1):
        self.device = hardware()
        self.full = (1 << len(leaves)) - 1
        self.loads = [0] * self.device.car.n
        for row in self.device.rows:
            self.loads[row[2][0]] = self.full
        for side, (b, word) in enumerate(zip(BOUNDARY, initial, strict=True)):
            p = layout(side)
            self.set_pair(p["boundary"], b)
            rank = sum(x != b for x in word)
            for k, bank in enumerate(p["banks"]):
                for pair, bit in zip(bank, word, strict=True):
                    self.set_pair(pair, bit if k == rank else hidden)
                self.loads[p["selector"][k]] = self.full if k == rank else 0
        self.clock = tuple(pulse_calibration()["clock"])
        self.set_pair(self.clock, 1)
        self.latch = tuple(tuple(pair) for pair in pulse_calibration()["latch"])
        for pair in self.latch:
            self.set_pair(pair, 0)
        self.columns = np.array(leaves, dtype=np.uint8)
        self.reads = []

    def set_pair(self, pair, value):
        a, b = pair
        self.loads[a] = self.full if value else 0
        self.loads[b] = self.full ^ self.loads[a]

    def swap(self, route, enabled):
        for a, b in route:
            change = (self.loads[a] ^ self.loads[b]) & enabled
            self.loads[a] ^= change
            self.loads[b] ^= change

    @cache
    def route(self, ports, station):
        first, second = ports
        a, b = self.device.rows[station][3]
        return self.device.move(first, a) + self.device.move(
            first if second == a else second, b
        )

    def compare(self, ports, station, active):
        route = self.route(tuple(ports), station)
        a, b = self.device.rows[station][3]
        self.swap(route, active)
        click = active & (self.loads[a] ^ self.loads[b])
        self.swap(reversed(route), active)
        self.reads.append([station, *ports, format(active, "x"), format(click, "x")])
        return click

    def rail(self, r, station, active):
        return self.compare((r, self.device.rows[station][2][0]), station, active)

    def read_bit(self, pair, station, active):
        lo, hi = (self.rail(r, station, active) for r in pair)
        if lo ^ hi != active:
            raise ArithmeticError("Unbalanced dimer encountered")
        return hi

    def state(self, side, station):
        p, full = layout(side), self.full
        high = [full ^ self.rail(port, station, full) for port in p["selector"]]
        b = self.read_bit(p["boundary"], station, full)
        word = [0] * 3
        for k, active in enumerate(high):
            if active:
                for j, pair in enumerate(p["banks"][k]):
                    word[j] |= self.read_bit(pair, station, active)
        if classes([x ^ b for x in word], full) != high:
            raise ArithmeticError("Selector failed independent mismatch census")
        return b, word, high

    def quarter(self, crossed, enabled):
        # Conditional sharp Kraus branches of the calibrated unrevealed clock.
        # The dwell and projection endpoints are derived from the native matrix.
        self.swap([self.clock], crossed)
        stopped = self.read_bit(self.clock, COMMIT, enabled)
        received = enabled ^ stopped
        if received != crossed:
            raise ArithmeticError(
                "Clocked material reading differs from its projected outcome"
            )
        self.swap([self.clock], received)
        if self.read_bit(self.clock, COMMIT, enabled) != enabled:
            raise ArithmeticError("Clock did not return to its initial material state")
        return received

    def draw(self, step):
        labels = self.columns[:, step]
        first = self.quarter(packed(labels != 0), self.full)
        second = self.quarter(packed(labels == 2), first)
        for pair, requested in zip(self.latch, (first, second), strict=True):
            present = self.read_bit(pair, COMMIT, self.full)
            self.swap([pair], present ^ requested)

    def latched(self):
        first, second = (self.read_bit(pair, COMMIT, self.full) for pair in self.latch)
        if second & ~first:
            raise ArithmeticError("Invalid material fault word")
        return self.full ^ first, first ^ second, second

    def offer(self, side, m, mode):
        p, full = layout(side), self.full
        b, old, ranks = self.state(side, 6 * side)
        if self.state(side, COMMIT) != (b, old, ranks):
            raise ArithmeticError("Stale snapshot")
        target = old.copy()
        for i in MASKS[m]:
            for station in (6 * side, COMMIT):
                flag = 0
                for k, active in enumerate(ranks):
                    if active:
                        flag |= self.compare(
                            (p["boundary"][0], p["banks"][k][i][0]), station, active
                        )
                if flag != (b ^ old[i]):
                    raise ArithmeticError("Native material guard disagrees")
            good, hold, bad = mode
            target[i] = (good & b) | (hold & old[i]) | (bad & (full ^ b))
        new = classes([v ^ b for v in target], full)
        for k, active in enumerate(new):
            if active:
                for bit, pair in zip(target, p["banks"][k], strict=True):
                    current = self.read_bit(pair, COMMIT, active)
                    self.swap([pair], active & (current ^ bit))
        for i, high in enumerate(ranks):
            for j, low in enumerate(new):
                if i != j and high & low:
                    self.swap(
                        self.device.move(p["selector"][i], p["selector"][j]), high & low
                    )
        if self.state(side, FAR) != (b, target, new):
            raise ArithmeticError("Publication did not reach the independent receiver")

    def snapshot(self):
        result = []
        for side in (0, 1):
            b, word, rank = self.state(side, FAR)
            result.append(
                [
                    format(b, "x"),
                    *[format(x, "x") for x in word],
                    *[format(x, "x") for x in rank],
                ]
            )
        return result


def acquire(order, initial, held=False, reverse_faults=False, selected=None):
    leaves = (
        list(product(range(3), repeat=(1 if held else len(order))))
        if selected is None
        else selected
    )
    bench = Bench(leaves, initial)
    steps = [bench.snapshot()]
    snapshot_ranges = [[0, len(bench.reads)]]
    rng_records = []
    for k, (side, m) in enumerate(order):
        start = len(bench.reads)
        if k == 0 or not held:
            bench.draw(0 if held else k)
        mode = bench.latched()
        if reverse_faults:
            mode = mode[0], mode[2], mode[1]
        rng_records.append([start, len(bench.reads), *[format(x, "x") for x in mode]])
        bench.offer(side, m, mode)
        start = len(bench.reads)
        steps.append(bench.snapshot())
        snapshot_ranges.append([start, len(bench.reads)])
    return dict(
        order=order,
        initial=initial,
        held=held,
        reverse_faults=reverse_faults,
        leaves=len(leaves),
        state=steps,
        snapshot_ranges=snapshot_ranges,
        selector_records=rng_records,
        reads=bench.reads,
    )
