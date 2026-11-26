"""Sharp two-bank comparisons extended to four directly measured count addresses."""

from functools import cache
from itertools import combinations

from ...prepared_quorum.transport import read_rail
from ..board import hardware
from ..contract import require

COMMIT = 13
FAR = 15
MASKS = ((0, 1), (0, 2), (1, 2), (0, 1, 2))


@cache
def layout(side):
    machine = hardware()
    cells = machine.car
    selector_cell = 16 + side
    local = set(range(12 * selector_cell, 12 * (selector_cell + 1)))
    for one, two in combinations(sorted(local), 2):
        if two not in cells.neighbours[one]:
            continue
        opposite = sorted(
            local & set(cells.neighbours[one]) & set(cells.neighbours[two])
        )
        if len(opposite) == 2 and opposite[1] not in cells.neighbours[opposite[0]]:
            selector = (opposite[0], one, two, opposite[1])
            break
    else:
        raise ValueError("Missing adjacent native triangles")
    pairs = [machine.rows[6 * side + i][slot] for i in range(6) for slot in range(2)]
    return dict(
        banks=tuple(tuple(pairs[3 * k : 3 * k + 3]) for k in range(4)),
        boundary=machine.rows[12][side],
        selector=selector,
    )


def mismatch(b, word):
    return sum(bit != b for bit in word)


class Experiment:
    def __init__(self, boundaries, words, hidden=0):
        self.m = hardware()
        self.loads = self.m.fresh([0] * 16)
        self.raw = []
        self.observations = []
        self.attempts = []
        for side, (b, word) in enumerate(zip(boundaries, words, strict=True)):
            ports = layout(side)
            self.initialize(ports["boundary"], b)
            rank = mismatch(b, word)
            for k, bank in enumerate(ports["banks"]):
                for j, pair in enumerate(bank):
                    self.initialize(pair, word[j] if k == rank else hidden)
                self.loads[ports["selector"][k]] = 2 * int(k == rank)

    def initialize(self, pair, bit):
        left, right = pair
        self.loads[left], self.loads[right] = 2 * bit, 2 * (1 - bit)

    def query(self, port, station):
        _, _, reference, pointer = self.m.rows[station]
        click = read_rail(self.m, self.loads, port, reference[0], pointer)
        self.raw.append([station, port, int(click)])
        return click

    def read_bit(self, pair, station):
        pair_read = tuple(self.query(port, station) for port in pair)
        require(pair_read in ((1, 0), (0, 1)), "Unbalanced dimer")
        return pair_read[1]

    def state(self, side, station, purpose):
        start = len(self.raw)
        p = layout(side)
        selector = [1 - self.query(port, station) for port in p["selector"]]
        require(sum(selector) == 1, "Selector is not one-hot")
        k = selector.index(1)
        b = self.read_bit(p["boundary"], station)
        word = tuple(self.read_bit(pair, station) for pair in p["banks"][k])
        require(mismatch(b, word) == k, "Count disagrees with independent word")
        self.observations.append(
            [purpose, side, station, start, len(self.raw), b, list(word), k]
        )
        return b, word, k

    def guard(self, side, k, mask, station):
        ports = layout(side)
        head = self.m.rows[station][3]
        return [
            int(
                self.m.probe(
                    self.loads,
                    (ports["boundary"][0], ports["banks"][k][i][0]),
                    head,
                    "native",
                )
            )
            for i in mask
        ]

    def write(self, pair, x):
        if self.read_bit(pair, COMMIT) != x:
            u, v = pair
            self.loads[u], self.loads[v] = self.loads[v], self.loads[u]

    def offer(self, side, mask, park=False):
        b, old, k = self.state(side, 6 * side, "prepare")
        before = self.guard(side, k, mask, 6 * side)
        public = self.state(side, COMMIT, "guard")
        after = self.guard(side, k, mask, COMMIT)
        require(
            public == (b, old, k)
            and before == after == [int(old[i] != b) for i in mask],
            "Independent endpoint guard mismatch",
        )
        target = tuple(b if j in mask else x for j, x in enumerate(old))
        lower = k - sum(after)
        local = [side, list(mask), before, after, k, lower]
        if not any(after):
            self.attempts.append([*local, "stutter", []])
            return
        ports = layout(side)
        for pair, x in zip(ports["banks"][lower], target, strict=True):
            self.write(pair, x)
        require(
            self.state(side, FAR, "stage")[1] == old, "Staging leaked into public word"
        )
        require(self.state(side, COMMIT, "revalidate") == public, "Stale proposal")
        require(
            tuple(self.read_bit(pair, COMMIT) for pair in ports["banks"][lower])
            == target,
            "Wrong staged word",
        )
        route = self.m.move(ports["selector"][k], ports["selector"][lower])
        if k - lower <= 2:
            require(len(route) == 1, "Expected direct selector contact")
        orig = self.loads.copy()
        if not park:
            for u, v in route:
                self.loads[u], self.loads[v] = self.loads[v], self.loads[u]
            u, v = ports["selector"][k], ports["selector"][lower]
            orig[u], orig[v] = orig[v], orig[u]
        require(orig == self.loads, "Selector transit changed another load")
        read = self.state(side, FAR, "published" if not park else "parked")
        require(read == (public if park else (b, target, lower)), "Wrong publication")
        self.attempts.append([*local, "parked" if park else "published", list(route)])

    def checkpoint(self, name):
        words = [self.state(i, FAR, name) for i in range(2)]
        click = int(
            self.m.probe(
                self.loads,
                tuple(layout(i)["boundary"][0] for i in range(2)),
                self.m.rows[FAR][3],
                "native",
            )
        )
        require(
            click == int(words[0][0] != words[1][0]),
            "Independent protected parity mismatch",
        )
        return [[b, list(word), k] for b, word, k in words], click

    def exercise(self, offers):
        start = self.checkpoint("initial")
        for side, key, park in offers:
            self.offer(side, MASKS[key], park)
        end = self.checkpoint("terminal")
        return dict(
            initial=start,
            end=end,
            offers=self.attempts,
            observations=self.observations,
            raw=self.raw,
        )


def calibration():
    board = Experiment((0, 1), ((0, 0, 0), (1, 1, 1)))
    records = []
    for side in range(2):
        ports = layout(side)
        for pair in [ports["boundary"], *(p for bank in ports["banks"] for p in bank)]:
            for bit in (0, 1):
                board.initialize(pair, bit)
                snapshot = board.loads[:]
                require(
                    all(
                        board.read_bit(pair, station) == bit
                        for station in (6 * side, COMMIT, FAR)
                    ),
                    "Bit calibration failed",
                )
                require(snapshot == board.loads, "Calibration read disturbed stock")
        for high in range(4):
            for low in range(high):
                for j, port in enumerate(ports["selector"]):
                    board.loads[port] = 2 * int(j == high)
                old = board.loads[:]
                route = board.m.move(ports["selector"][high], ports["selector"][low])
                for a, b in route:
                    board.loads[a], board.loads[b] = board.loads[b], board.loads[a]
                a, b = ports["selector"][high], ports["selector"][low]
                old[a], old[b] = old[b], old[a]
                out = [1 - board.query(p, FAR) for p in ports["selector"]]
                require(
                    board.loads == old and out == [int(i == low) for i in range(4)],
                    "Route calibration failed",
                )
                records.append([side, high, low, list(route), out])
    return dict(raw=board.raw, selectors=records)
