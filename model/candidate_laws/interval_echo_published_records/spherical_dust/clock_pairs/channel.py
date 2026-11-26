"""Two local sign records from a pre-shared, finite lag-table source."""

from dataclasses import dataclass
from decimal import Decimal, localcontext
from fractions import Fraction
from functools import cached_property
from math import floor, isqrt, log, log1p


VACUUM_CORRELATION = Fraction(1, 2)
VACUUM_SLOPE = log(3) / 2


@dataclass(frozen=True)
class LagGrid:
    bits: int
    max_lag: int
    subphases: int = 16

    def __post_init__(self):
        if self.bits < 1 or not 0 <= self.max_lag < self.bins or self.subphases < 1:
            raise ValueError("The finite binary lag labels must lie strictly inside (-1,1)")

    @property
    def bins(self):
        return 1 << self.bits

    def fine_tick(self, clock, separation):
        if separation <= 0:
            raise ValueError("The material pair needs positive comoving separation")
        return floor(Fraction(clock) * (self.bins*self.subphases) / Fraction(separation))


def entropy_remainder(x):
    """(1+x) log(1+x) - x, including weak source weights."""
    if abs(x) < 1e-4:
        return x*x*(0.5 + x*(-1/6 + x*(1/12 - x/20)))
    return (1+x)*log1p(x) - x


@dataclass(frozen=True)
class ClockPair:
    grid: LagGrid
    weight: float
    a_max: float = 1.0

    def __post_init__(self):
        if self.a_max <= 0:
            raise ValueError("A positive clock-scale envelope is required")
        if (2*self.grid.max_lag+1)*self.class_width > self.tickets//4:
            raise ValueError("Refine the material pair weight to fit the finite lag classes")

    @cached_property
    def ticket_bits(self):
        w = abs(Fraction(self.weight))*Fraction(self.a_max)
        # Per-class rounding relative to the leading entropy response is O(b²).
        return max(52, w.denominator.bit_length() - w.numerator.bit_length() + 2*self.grid.bits + 16)

    @cached_property
    def tickets(self):
        return 1 << self.ticket_bits

    @cached_property
    def coupling_ticks(self):
        w = abs(Fraction(self.weight))*Fraction(self.a_max)
        with localcontext() as ctx:
            ctx.prec = 30 + self.ticket_bits
            x = Decimal(w.numerator*self.tickets) / w.denominator
            return int(2*x / Decimal(3).ln())

    def class_ticks(self, offset):
        if abs(offset) > self.grid.max_lag:
            return 0
        b2 = self.grid.bins**2
        return self.coupling_ticks*b2 // (b2-offset**2)

    @cached_property
    def class_width(self):
        return self.class_ticks(self.grid.max_lag)

    def gate(self, scale):
        a = Fraction(scale)/Fraction(self.a_max)
        if not 0 <= a <= 1:
            raise ValueError("Local clock scale exceeds the public envelope")
        return isqrt(self.tickets**2*a.numerator//a.denominator)

    def shift(self, clock_a, clock_b, separation, scale_a=1, scale_b=1):
        ta, tb = (self.grid.fine_tick(t, separation) for t in (clock_a, clock_b))
        n = self.grid.subphases
        lo, rem = divmod(tb-ta, n)
        matched = (n-rem)*self.class_ticks(lo) + rem*self.class_ticks(lo+1)
        sign = -1 if self.weight >= 0 else 1
        return Fraction(sign*matched*self.gate(scale_a)*self.gate(scale_b), n*self.tickets**3)

    def correlation(self, clock_a, clock_b, separation, scale_a=1, scale_b=1):
        return float(VACUUM_CORRELATION + self.shift(clock_a, clock_b, separation, scale_a, scale_b))

    def entropy_change(self, clock_a, clock_b, separation, scale_a=1, scale_b=1):
        d = float(self.shift(clock_a, clock_b, separation, scale_a, scale_b))
        return -VACUUM_SLOPE*d - self.divergence(d)

    def entropy_from_tick_differences(self, tick_difference, scale_a=1, scale_b=1):
        """Batch of the same dyadic channel, given differences of registered fine ticks."""
        import numpy as np

        ticks = np.asarray(tick_difference)
        if ticks.dtype.kind not in ("i", "u"):
            raise ValueError("Register each endpoint's fine tick before taking differences")
        n = self.grid.subphases
        lo, rem = np.divmod(ticks, n)

        def mapped(values, operation):
            values = np.asarray(values)
            unique, inverse = np.unique(values, return_inverse=True)
            table = np.array([operation(item.item() if isinstance(item, np.generic) else item)/self.tickets
                              for item in unique])
            return table[inverse].reshape(values.shape)

        matched = (1-rem/n)*mapped(lo, self.class_ticks) + (rem/n)*mapped(lo+1, self.class_ticks)
        sign = -1 if self.weight >= 0 else 1
        d = sign*matched*mapped(scale_a, self.gate)*mapped(scale_b, self.gate)

        def remainder(x):
            return np.where(np.abs(x) < 1e-4, x*x*(0.5+x*(-1/6+x*(1/12-x/20))),
                            (1+x)*np.log1p(x)-x)

        return -VACUUM_SLOPE*d - (1.5*remainder(d/1.5)+0.5*remainder(-d/0.5))/2

    @staticmethod
    def divergence(d):
        c = float(VACUUM_CORRELATION)
        return ((1+c)*entropy_remainder(d/(1+c)) + (1-c)*entropy_remainder(-d/(1-c))) / 2

    def issue(self, rng, separation, *, shared_lag=True):
        """Draw a copy of the clock material before either endpoint is queried."""
        ticket = rng.getrandbits(self.ticket_bits)
        class_ticket = ticket-3*self.tickets//4
        offset = 0
        if ticket < self.tickets//2:
            kind = "base"
        elif class_ticket < 0:
            kind = "white"
        else:
            kind = "independent"
            if self.class_width:
                block, position = divmod(class_ticket, self.class_width)
                offset = block-self.grid.max_lag
                if position < self.class_ticks(offset) and shared_lag:
                    kind = "lag"
        return PairTape(self, separation, kind, offset, rng.randrange(self.grid.subphases), rng)


class PairTape:
    """Deferred independent tape entries; rereading an addressed entry returns its bit."""

    def __init__(self, source, separation, kind, offset, phase, rng):
        self.source, self.separation = source, separation
        self.kind, self.offset, self.phase, self.rng = kind, offset, phase, rng
        self.signs, self.gates = {}, {}

    def sign(self, key):
        if key not in self.signs:
            self.signs[key] = 2*self.rng.getrandbits(1)-1
        return self.signs[key]

    def read(self, side, clock, scale=1):
        if side not in ("a", "b"):
            raise ValueError("Unknown endpoint")
        tick = self.source.grid.fine_tick(clock, self.separation)
        gate = self.source.gate(scale)
        if self.kind == "base":
            return self.sign("base")
        if self.kind == "white":
            return self.sign((side,tick,gate))
        index = (tick+self.phase)//self.source.grid.subphases
        orientation = 1
        if self.kind == "lag":
            key = ("shared",index-self.offset if side == "b" else index)
            if side == "b" and self.source.weight >= 0:
                orientation = -1
        else:
            key = (side,index)
        coin = side,tick,gate
        if coin not in self.gates:
            self.gates[coin] = self.rng.getrandbits(self.source.ticket_bits+1)
        positive = self.gates[coin] < self.source.tickets+gate
        return self.sign(key)*orientation*(1 if positive else -1)
