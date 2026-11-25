"""Fair source orientations and a connected thermal interaction on equal-energy pins."""

import argparse
from fractions import Fraction
from itertools import product
import json
from math import inf, nextafter

import numpy as np

from . import matched_energy, measure


def outward(lower, upper):
    return [nextafter(float(lower), -inf), nextafter(float(upper), inf)]


def orientations(board, pins):
    rows = []
    for q, reverse_a, reverse_b in product((1, -1), (False, True), (False, True)):
        first = pins[:2][::-1] if reverse_a else pins[:2]
        graph, h, va, vb = measure.operators(board, (*first, *pins[2:]), 8, q, reverse_b)
        if np.any(va * vb):
            raise ValueError("the two potential neighborhoods overlap")
        if any(np.min(h.diagonal() + i*va + j*vb - graph.sum(axis=1)) < 0
               for i, j in measure.SETTINGS):
            raise ValueError("the zero-mobile reference does not bound the heat traces")
        rows.append(dict(q=q, reverse_a=reverse_a, reverse_b=reverse_b,
                         heat=measure.thermal_trace(h, va, vb),
                         path_weight=measure.first_mixed_moment(graph, va, vb)["signed_path_weight"]))
    return rows


def summarize(panels):
    # Use the exact rational values of the outward floating-point bounds.
    summaries = []
    for p in panels:
        arms = p["group"]
        lower, upper = (sum(Fraction(r["heat"][side]) for r in arms) / 8 for side in (0, 1))
        weight = sum(Fraction(r["path_weight"]) for r in arms) / 8
        # The same bound works before or after mixing the conditional traces.
        error = Fraction(16 * 236**2, 10**72)
        summaries.append(dict(distance=p["distance"], orientations=arms,
                              mean_path_weight=str(weight),
                              mean_heat_trace=outward(lower, upper),
                              potential_per_activity=outward(-8*upper-error, -8*lower+error)))
    comparisons = []
    for previous, following in zip(summaries, summaries[1:]):
        r = previous["distance"]
        if following["distance"] != r + 1:
            raise ValueError("calibration compares consecutive station distances")
        a, b = map(Fraction, previous["potential_per_activity"])
        c, d = map(Fraction, following["potential_per_activity"])
        comparisons.append(dict(from_distance=r, to_distance=r+1, gain=outward(
            r*(r+1)*(c-b), r*(r+1)*(d-a))))
    return dict(kappa=8, beta="1/8", activity="1e-72", source_energy="233/2",
                panels=summaries, inverse_radius_calibrations=comparisons)


def run():
    board = measure.federation(0)
    panels = []
    for distance, _, original in measure.sources(board):
        if distance < 2:
            continue
        a = matched_energy.regular_pair(board, original[0])
        b = matched_energy.regular_pair(board, original[2])
        panels.append(dict(distance=distance, group=orientations(board, (*a, *b))))
    return summarize(panels)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    text = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (measure.HERE / "unoriented.json").write_text(text)
    else:
        print(text, end="")
