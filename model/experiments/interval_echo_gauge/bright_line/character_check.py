"""Exclude every determinant character from an antipodal pulse and pair return counts."""
from fractions import Fraction
from itertools import product
import json
from math import exp, pi
from pathlib import Path

from scipy.stats import beta

from model.executable_model import cell
from .certificate import pair_bright


def check():
    record = json.loads(Path(__file__).with_name("records.json").read_text())
    board = cell()
    port = 1
    near = set(board.neighbours[port]) | {port}
    within_two = near.union(*(set(board.neighbours[i]) for i in near))
    antipode, = set(range(12)) - within_two
    star = [i for i, edge in enumerate(board.seams) if port in edge or antipode in edge]

    # The antipodal pulse has first one-load trace 2*(1/4). The determinant
    # of the product is therefore i, also without a product approximation.
    turns = Fraction(1, 2) * (Fraction(1, 4) + Fraction(1, 4))
    assert turns == Fraction(1, 4)
    low, high = pair_bright()
    assert Fraction(1, 2) < low < high < 1  # pair phase strictly between pi/2 and pi

    n = record["shots"]
    steps = record["symmetric_slices"]
    length = 12 * (pi / 2) * .81
    bias = exp(length / steps) * length**3 / (3 * steps**2)
    tail = .01 / 6  # three pair preparations, two tails
    intervals = []
    for arm in ("native", "sham", "phase_anchor"):
        hist = record["pair"][arm]["launches"][str(port)]["counts"]
        if len(hist) != 31 or sum(hist) != n or min(hist) < 0:
            raise ValueError("incomplete seam and quiet histogram")
        hits = sum(hist[i] for i in star)
        lower = beta.ppf(tail, hits, n - hits + 1) if hits else 0.
        upper = beta.isf(tail, hits + 1, n - hits) if hits < n else 1.
        delta = bias if arm == "phase_anchor" else 0.
        intervals.append((max(0., lower - delta), min(1., upper + delta)))
    if intervals[1][0] <= intervals[2][1]:
        raise ValueError("unresolved phase calibration")
    candidates = [1 + 2*(read - sham)/(sham - anchor)
                  for read, sham, anchor in product(*intervals)]
    lo, hi = min(candidates), max(candidates)
    if not -1 < lo <= hi < 0:
        raise ValueError("the measured phase fails to exclude all powers of i")
    print(f"Pair cosine in [{lo:.7f}, {hi:.7f}] at simultaneous 99% coverage")
    print("All integer powers of the measured service's exact determinant have real part -1, 0, or 1.")


if __name__ == "__main__":
    check()
