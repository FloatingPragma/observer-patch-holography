"""Distinguishability and the area range required by the local word experiment."""

import argparse
from fractions import Fraction
import json
from pathlib import Path

from . import word_entropy as words

HERE = Path(__file__).resolve().parent


def run():
    groups = words.probabilities()
    # dP/ds = D. The positive entries of D constitute the optimal sharp event.
    tv = sum(Fraction(count * max(d, 0), 2592) for _, d, count in groups)
    base = sum(Fraction(count * (p - d), 2592) for p, d, count in groups if d > 0)
    gain = words.entropy_and_slope(groups, -1)[1]
    if not 0 < base < 1 or not 0 < tv or gain.a <= 0:
        raise ValueError("the source has no positive regional contrast")
    pair_limit = Fraction(729, 256)
    panels = []
    for row in json.loads((HERE / "results.json").read_text())["midpoint_surveys"]:
        n, total = row["annulus_stations"], row["cells"]
        pair = Fraction(n * (n - 1), total * (total - 1))
        panels.append(dict(
            level=row["level"],
            eligible_probability=str(pair),
            optimal_event_probability_at_reference=str(pair * base),
            optimal_event_slope=str(pair * tv),
            entropy_slope_nats=words.endpoints(words.interval(pair) * gain),
        ))
    return dict(
        conditional_positive_event_reference_probability=str(base),
        conditional_total_variation_slope=str(tv),
        calibrated_area_range_times_copies_required=words.endpoints(gain / words.interval(tv)),
        one_copy_optimal_event_slope_per_fourth_power=str(tv * pair_limit),
        finite_panels=panels,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    out = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if parser.parse_args().write:
        (HERE / "area_read.json").write_text(out)
    else:
        print(out, end="")
