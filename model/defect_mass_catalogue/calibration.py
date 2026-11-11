"""One positive GeV scale for a specified source clock, with exact intervals."""
from dataclasses import dataclass
from fractions import Fraction as F


def interval(value, *, positive=False):
    if len(value) != 2:
        raise ValueError('an interval needs two endpoints')
    lo, hi = map(F, value)
    if lo > hi or (positive and lo <= 0):
        raise ValueError('ordered positive interval required')
    return lo, hi


@dataclass(frozen=True)
class Anchor:
    source_clock: str
    source_energy: tuple
    target_GeV: tuple
    interpretation: str
    comparison_convention: str = "source-energy calibration"

    def scale_interval(self):
        x0, x1 = interval(self.source_energy, positive=True)
        y0, y1 = interval(self.target_GeV, positive=True)
        return y0/x1, y1/x0


def common_scale(anchors):
    anchors = tuple(anchors)
    if not anchors:
        raise ValueError('at least one explicit empirical anchor is required')
    if not anchors[0].source_clock or len({a.source_clock for a in anchors}) != 1:
        raise ValueError('a common declared source and clock are required')
    if len({a.comparison_convention for a in anchors}) != 1:
        raise ValueError("incompatible target mass conventions require an explicit matching map")
    ranges = [a.scale_interval() for a in anchors]
    lo, hi = max(x[0] for x in ranges), min(x[1] for x in ranges)
    return None if lo > hi else (lo, hi)


def convert_energy(source_clock, expected_clock, energy, scale):
    if not source_clock or source_clock != expected_clock:
        raise ValueError('cannot transfer an energy scale between unmatched clocks')
    e0, e1 = interval(energy)
    s0, s1 = interval(scale, positive=True)
    products = [e*s for e in [e0,e1] for s in [s0,s1]]
    return min(products), max(products)


def ratio(a, b):
    a0,a1=interval(a,positive=True); b0,b1=interval(b,positive=True)
    return a0/b1,a1/b0


def strings(pair):
    return None if pair is None else [str(x) for x in pair]
