"""Material endpoint reads and swaps compiled from registered Interval Echo seam columns."""

from collections import deque
from pathlib import Path
import json
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from background_clock.material import clock_gates, observations, screen
from relay import CELL, align, apply, gates

CONFIG = json.loads((HERE / "protocol.json").read_text())
SIZE = CONFIG["max_coordinate"] + 1
AXES = (SIZE**2, SIZE, 1)
POINTS = range(SIZE**3)


def coords(site):
    x, yz = divmod(site, SIZE**2)
    y, z = divmod(yz, SIZE)
    return x, y, z


def grid_edges():
    return [
        (s, s + step)
        for s in POINTS
        for k, step in enumerate(AXES)
        if coords(s)[k] < SIZE - 1
    ]


class Transport:
    """Compile logical routes, including injection, on a single shared array of ports."""

    def __init__(self):
        gates()
        self.flip = clock_gates()[1]
        self.board, self.levels, self.radii, *_ = screen(CONFIG["carrier_level"])
        self.word = [0] * self.board.n
        self.ops = {}
        self.pointer = {}
        law = observations.native_law
        echo = law.echo()
        mask = np.arange(law.DIM) % len(law.SRC) == law.SRC.index((1, 1))
        for value in (0, 2):
            col = law.SRC.index((value, 0))
            branch = echo[:, col] * (mask if value else ~mask)
            np.testing.assert_allclose(
                echo.conj().T @ branch, np.eye(law.DIM)[:, col], rtol=0, atol=2e-13
            )
            self.pointer[value] = int(bool(value))
        self.read_ops = {p: align([p, 11], [p, min(CELL[p])]) for p in (2, 4, 6, 8, 9)}
        for off, src in zip(
            CONFIG["box_offsets"], CONFIG["injection_owners"], strict=True
        ):
            local_edges = [(src, off)] + [(off + a, off + b) for a, b in grid_edges()]
            for port in (CONFIG["packet_port"], CONFIG["range_port"]):
                self.compile(local_edges, port)
        self.slot = max(map(len, self.ops.values()))
        transit = {p // 12 for ops in self.ops.values() for edge in ops for p in edge}
        self.clock_owner = min(set(range(self.board.cells)) - transit)
        self.word[12 * self.clock_owner + 6] = 2
        self.service_clock = 0

    def compile(self, edges, port):
        by_start = {}
        for a, b in edges:
            by_start.setdefault(12 * a + port, []).append(12 * b + port)
        for start, targets in by_start.items():
            pred = {start: start}
            queue = deque([start])
            missing = set(targets)
            while missing:
                node = queue.popleft()
                for nxt in self.board.neighbours[node]:
                    if nxt not in pred:
                        pred[nxt] = node
                        missing.discard(nxt)
                        queue.append(nxt)
            for end in targets:
                path = []
                node = end
                while node != start:
                    path.append((pred[node], node))
                    node = pred[node]
                path.reverse()
                ops = path + path[:-1][::-1]
                # Distinct labels certify the full permutation, including borrowed payloads.
                labels = {p: p for edge in path for p in edge}
                apply(labels, ops)
                assert labels[start] == end and labels[end] == start
                assert all(k == v for k, v in labels.items() if k not in (start, end))
                self.ops[start, end] = self.ops[end, start] = ops

    def read(self, owner, port):
        base = 12 * owner
        local = self.word[base : base + 12]
        apply(local, self.read_ops[port])
        assert local[min(CELL[port])] == 0
        # The echo / projection / inverse column above restores this complete input.
        return self.pointer[local[port]]

    def hop(self, a, b, port, cut=False, metered=False):
        ops = self.ops[12 * a + port, 12 * b + port]
        if not metered:
            if not cut:
                apply(self.word, ops)
            return
        # Independent sharp clock on a station outside every transit path. The
        # parked tail has the same clock exposures as active transfer slots.
        clock = 12 * self.clock_owner + 6
        last = self.read(self.clock_owner, 6)
        for k in range(self.slot):
            if k < len(ops) and not cut:
                x, y = ops[k]
                self.word[x], self.word[y] = self.word[y], self.word[x]
            self.word[clock], self.word[clock + 1] = (
                self.word[clock + 1],
                self.word[clock],
            )
            reading = self.read(self.clock_owner, 6)
            self.service_clock += int(last != reading)
            last = reading
