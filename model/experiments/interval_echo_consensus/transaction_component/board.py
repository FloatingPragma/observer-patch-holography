"""Persistent sharp dimers with selectors, four registered observation stations."""

from functools import cache

from ..boolean_nodes.material import Machine
from ..prepared_quorum.transport import read_rail
from .contract import require


@cache
def hardware():
    return Machine(0, 16)


def allocation(side):
    m = hardware()
    pairs = [tuple(m.rows[5 * side + i][j]) for i in range(5) for j in range(2)]
    return dict(
        b=pairs[0],
        g=pairs[1],
        bank0=pairs[2:5],
        bank1=pairs[5:8],
        p=pairs[8],
        s=pairs[9],
    )


class Board:
    def __init__(self, states, hidden=0):
        self.m = hardware()
        self.loads = self.m.fresh([0] * 16)
        self.records = []
        self.raw = []
        for side, u in enumerate(states):
            a = allocation(side)
            for name, bit in zip(
                ("b", "g", "p", "s"), (u[0], u[1], u[5], u[6]), strict=True
            ):
                self.initialize(a[name], bit)
            for j, bit in enumerate(u[2:5]):
                self.initialize(a["bank0"][j], bit)
                self.initialize(a["bank1"][j], hidden)

    def initialize(self, pair, value):
        a, b = pair
        self.loads[a], self.loads[b] = 2 * value, 2 * (1 - value)

    def read(self, pair, station):
        _, _, std, head = self.m.rows[station]
        clicks = tuple(read_rail(self.m, self.loads, p, std[0], head) for p in pair)
        self.raw.append([station, pair[0], 2 * clicks[0] + clicks[1]])
        require(clicks in ((1, 0), (0, 1)), "unbalanced read")
        return clicks[1]

    def write(self, pair, value, station):
        if self.read(pair, station) != value:
            a, b = pair
            self.loads[a], self.loads[b] = self.loads[b], self.loads[a]

    def snapshot(self, side, station, purpose):
        a = allocation(side)
        # Isolated observation appointment; the selector actually chooses the
        # bank. The public reader cannot inspect inactive values by this call.
        p = self.read(a["p"], station)
        bank = a[f"bank{p}"]
        u = (
            self.read(a["b"], station),
            self.read(a["g"], station),
            *(self.read(pair, station) for pair in bank),
            p,
            self.read(a["s"], station),
        )
        self.records.append([purpose, side, station, "".join(map(str, u))])
        return u

    def local(self, side):
        return 5 * side

    def survey(self):
        return [self.snapshot(i, self.local(i), "prepare") for i in range(2)]

    def audit(self, side):
        return self.snapshot(side, 15, "far")
