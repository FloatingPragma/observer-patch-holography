"""Exact ceiling on inclusive home detection before flight on row 30 material."""

import argparse
from fractions import Fraction

import numpy as np

from tools.strict_json import load
from . import flight
from .__main__ import ARCHIVE


def ceiling(write=Fraction(1, 4)):
    """At zero coupling, count every first vacancy on the two actual pulse seams."""
    if not isinstance(write, Fraction) or not 0 <= write <= 1:
        raise ValueError("Expected an exact write probability")
    rows = []
    for x in (False, True):
        for z in (False, True):
            if z:
                # No contact can move X: its two loads share the X seam.
                first = [Fraction(0)]*3
            elif x:
                # The remaining b-c contact couples X to Y. A 2π/3 exchange
                # leaves X with probability cos²(π/3)=1/4 after each read.
                first = [Fraction(3, 4**i) for i in (1, 2, 3)]
            else:
                first = [Fraction(4, 9**i) for i in (1, 2, 3)]
            prior = (write if x else 1-write) * (write if z else 1-write)
            rows.append((x, z, prior, first, sum(first)**8))
    return rows, 1 - sum(prior*failed for _, _, prior, _, failed in rows)


def verify():
    record = load(ARCHIVE)
    rows, limit = ceiling()
    # Recompute all of the pulse branches on the physical seam reductions,
    # including the two ways of obtaining a sharp accepted input.
    for x, z, _, first, failed in rows:
        _, prep = flight.acquired(0, x, z)
        np.testing.assert_allclose(prep["first_vacancy_per_attempt"], list(map(float, first)),
                                   atol=2e-13, rtol=0)
        np.testing.assert_allclose(prep["exhaustion"] + prep["failed_with_reset_unavailable"],
                                   float(failed), atol=2e-13, rtol=0)
    board = flight.federation(0)
    _, home = flight.open_pair.geometry()
    panel = record["ensemble"]
    if (type(record["units"]) is not int or record["units"] != flight.UNITS
            or record["pair_outcomes"] != [*map(list, home), "quiet"]
            or (panel["strength"], panel["exposure"]) != (0, 7)):
        raise ValueError("Wrong panel units, appointment or home reader")
    fair = flight.transfer.donor.native()
    failures = {(x, z): failed for x, z, _, _, failed in rows}
    for seed, bank in zip(flight.SEEDS, panel["banks"], strict=True):
        if type(bank["seed"]) is not int or bank["seed"] != seed:
            raise ValueError("Missing or substituted material seed")
        written, packed = flight.material(board, seed, fair)
        if packed != bank["coins_hex"] or type(bank["writes"]) is not int or len(written) != bank["writes"]:
            raise ValueError("Material history differs from the source")
        failed = failures[flight.transfer.X in written, flight.transfer.Z in written]
        count = bank["source_exhaustion"]
        reset = bank["reset_unavailable"]
        bins = bank["acquired"]
        if (type(count) is not int or count != round(flight.UNITS*failed)
                or type(reset) is not int or reset != 0 or len(bins) != 31
                or any(type(n) is not int or not 0 <= n <= flight.UNITS for n in bins)
                or abs(sum(bins)-flight.UNITS) > len(bins)/2 or bins[-1] < count-1):
            raise ArithmeticError("Source losses disagree with the coin words or the inclusive read")
    if limit >= Fraction(94, 100):
        raise ArithmeticError("The inherited preparation does not exclude the benchmark")
    return limit


def panel_ceiling(counts):
    """Empirical ceiling on a complete, equally weighted, fixed finite donor panel."""
    if len(counts) != 4 or any(type(n) is not int or n < 0 for n in counts) or not sum(counts):
        raise ValueError("Expected four nonnegative counts with nonzero total")
    rows, _ = ceiling()
    failures = {(x, z): f for x, z, _, _, f in rows}
    cases = [(False, False), (True, False), (False, True), (True, True)]
    lost = sum(n*failures[cut] for n, cut in zip(counts, cases, strict=True))
    return 1 - lost/sum(counts)


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--panel", nargs=4, type=int, metavar=("BOTH_OPEN", "ONLY_Z", "ONLY_X", "NEITHER"),
                     help="compute the ceiling for these material counts at the fixed address")
    args = cli.parse_args()
    result = verify() if args.panel is None else panel_ceiling(args.panel)
    print(f"At zero coupling, every free exposure: P(home) <= {result} = {float(result):.12f}")
    if args.panel is not None:
        print("This panel excludes 94%." if result < Fraction(94, 100)
              else "Acquisition alone does not decide persistence on this panel.")
