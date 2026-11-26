"""Check whether one common delay permits exact Euclidean three-space ranges."""
from decimal import Decimal
from fractions import Fraction
from itertools import combinations
import json
from pathlib import Path

from sympy import Matrix, Poly, Rational, Symbol, gcd


def delay_polynomial(ranges):
    """A Euclidean 3D fit requires a common zero of these anchored minors."""
    delay = Symbol("delay")
    def square(i, j):
        if i == j:
            return 0
        return (Rational(str(ranges[tuple(sorted((i, j)))])) - delay)**2

    def minor(stations):
        gram = Matrix([[(square(0, i) + square(0, j) - square(i, j))/2
                        for j in stations] for i in stations])
        return Poly(gram.det(method="domain-ge"), delay)

    return gcd(minor((1, 2, 3, 4)), minor((1, 2, 3, 5)))


def archived_ranges():
    path = Path(__file__).with_name("station_chart.json")
    data = json.loads(path.read_text(), parse_float=Decimal)
    held = list(combinations((4, 5, 6), 2))
    train = [pair for pair in combinations(range(7), 2) if pair not in held]
    for board in data["boards"]:
        population = board["population"]
        values = population["training_ranges"] + population["held_echo"]
        yield board["cells"], dict(zip(train + held, values, strict=True))


def archived_degrees():
    return {cells: delay_polynomial(ranges).degree() for cells, ranges in archived_ranges()}


# Post-hoc coordinates in thousandths of a hand unit, in archived station order.
# The delay 533/100 and these integer triples suffice to replay the fit.
MILLI_COORDINATES = {
    20: [(0, 0, 0), (9384, 0, 0), (-4111, 9867, 0), (-2388, -2544, -10337),
         (11855, 9165, -1), (5535, 15232, -1908), (3365, 2069, -15903)],
    80: [(0, 0, 0), (9304, 0, 0), (-5725, 5801, 0), (-2119, -7197, 6699),
         (14115, 4181, 5953), (-4750, 14074, 4360), (-5177, -2768, 14112)],
    320: [(0, 0, 0), (9279, 0, 0), (-4203, 7864, 0), (-2637, -4542, -8159),
          (16152, 4551, -2141), (-889, 16695, -3857), (-4375, 1204, -15084)],
}
COMMON_DELAY = Fraction(533, 100)


def ranges_within_tolerance(ranges, milli_coordinates, delay, tolerance):
    """Test all measured ranges by rational bounds on squared Euclidean length."""
    if not 0 <= tolerance < 1 or delay < 0:
        raise ValueError("a relative tolerance in [0, 1) and a nonnegative delay are required")
    for (i, j), reading in ranges.items():
        reading = Fraction(reading)
        if reading <= 0:
            raise ValueError("positive measured ranges are required")
        lower = max(0, (1 - tolerance)*reading - delay)
        upper = (1 + tolerance)*reading - delay
        squared = sum((x-y)**2 for x, y in zip(milli_coordinates[i], milli_coordinates[j], strict=True))
        if upper < 0 or not (1000*lower)**2 <= squared <= (1000*upper)**2:
            return False
    return True


if __name__ == "__main__":
    print("Exact-fit polynomial degrees:", archived_degrees())
    for cells, ranges in archived_ranges():
        print(cells, "cells, shared-delay ranges within 5%:",
              ranges_within_tolerance(ranges, MILLI_COORDINATES[cells], COMMON_DELAY, Fraction(1, 20)))
