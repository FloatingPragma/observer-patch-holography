"""Completed finite censuses of three preparations, with independently decodable meters."""
from fractions import Fraction as F
from math import lcm

from ..body import archive_bodies
from ..survey import mass_tally
from .channel import VACUUM_CORRELATION
from .regions import Band


def census():
    """Reduced word counts for a complete sweep of the finite preparation tickets."""
    arms = {}
    for name, body in archive_bodies().items():
        tickets = [mass_tally(F(m).limit_denominator(33), 128) for m in (body.ordinary, body.collar)]
        stock = sum(F((hi-lo)*64, 7*n) for n, lo, hi in tickets)
        pair = Band(1).source(stock)
        # At simultaneous clocks the two gates open with product mean gain².
        # Enumerating the fixed class ticket and two independent signed gate
        # tickets gives this signed excess number of agreements for lag zero.
        n, g = pair.tickets, pair.gate(1)
        corr = F(1, 2)-F(pair.class_ticks(0)*g*g, n**3)
        if corr != VACUUM_CORRELATION + pair.shift(0, 0, F(1, 3)):
            raise ArithmeticError("class and gate census disagrees with the local-sign source")
        word_rates = ((1+corr)/4, (1-corr)/4, (1-corr)/4, (1+corr)/4)
        exterior = body.exterior_records()
        gap = exterior["received"][1]-exterior["received"][0]
        arms[name] = dict(mass=stock, receiver_gap=gap, word_rates=word_rates, correlation=corr)
    common_period = lcm(*(p.denominator for arm in arms.values() for p in arm["word_rates"]))
    for arm in arms.values():
        arm["word_counts"] = tuple(int(common_period*p) for p in arm["word_rates"])
        del arm["word_rates"]
    return common_period, arms


def verify():
    period, rows = census()
    if period & (period-1) or len({r["receiver_gap"] for r in rows.values()}) != 3:
        raise ArithmeticError("the balanced period or exterior clock no longer separates the stocks")
    for key in ("mass", "correlation", "word_counts"):
        if len({r[key] for r in rows.values()}) != 3:
            raise ArithmeticError(f"the completed {key} no longer identifies its stock")
    correlations = [rows[name]["correlation"] for name in ("held", "ordinary", "blank")]
    if not 0 < correlations[0] < correlations[1] < correlations[2] < 1:
        raise ArithmeticError("joint Shannon entropy ceased to distinguish the three stock censuses")
    for row in rows.values():
        if sum(row["word_counts"]) != period:
            raise ArithmeticError("incomplete census")
    return period, rows
