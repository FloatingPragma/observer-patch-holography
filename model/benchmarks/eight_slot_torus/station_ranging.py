"""Read the six-station ranging and pulse-control records of the eight-slot torus benchmark."""

import argparse
import json
from fractions import Fraction
from math import isqrt

import numpy as np

from .ranging_geometry import classify, distance_bands


def _source_geometry(pulsed, bands):
    """Resolve the zero-error Euclidean question for the symmetric source means.

    Signed spatial axis permutations fix the pulse and permute the station
    labels and slots together. The equal-slot launch, latch and local time
    queries have the same action after this permutation. This leaves just
    the opposite and adjacent classes of ordered pairs.
    """
    from .six_station import exact_menus

    def source_mean(j):
        radar, _ = exact_menus(0, j, pulsed)
        both = radar[1:,1:]
        total = sum(both.flat)
        if not total:
            return None
        ticks = sum((u+v+2)*p for (u,v),p in np.ndenumerate(both))
        return ticks / (2*total)

    opposite, adjacent = source_mean(1), source_mean(2)
    if opposite is None or adjacent is None:
        return {"status": "undefined conditional range"}
    # The centered Gram matrix acts on three antisymmetric vectors with
    # eigenvalue b²/2 and on two pair-symmetric contrasts with eigenvalue
    # a²-b²/2. The constant vector is in its kernel.
    excess = adjacent**2 - opposite**2/2
    dimension = 5 if excess > 0 else 3 if excess == 0 else None
    missed = []
    for (i,j), band in bands.items():
        d = opposite if i//2 == j//2 else adjacent
        if band is None or not band[0] <= d <= band[1]:
            missed.append([i,j])
    # Opposite pairs along two axes give the zero-sum contrast (1,1,-1,-1).
    # At distortion t this requires (1-t)b <= sqrt(2)(1+t)a. For b >=
    # sqrt(2)a the regular octahedron realizes equality on every pair.
    optimum = None
    if excess <= 0:
        grid = 10**12
        root = isqrt(2*grid**2)
        def bound(n):
            return max(0, (opposite - Fraction(n,grid)*adjacent)/(opposite + Fraction(n,grid)*adjacent))
        optimum = [str(bound(root+1)), str(bound(root))]
    return {"opposite_range": str(opposite), "adjacent_range": str(adjacent),
            "opposite_ticks": float(opposite), "adjacent_ticks": float(adjacent),
            "pair_symmetric_gram_eigenvalue": str(excess), "minimum_euclidean_dimension": dimension,
            "indefinite": excess < 0, "optimal_euclidean_distortion_enclosure": optimum,
            "source_ranges_outside_or_missing_from_count_bands": missed}


def outcomes(histogram):
    """Partition a single directed cohort by its independently recorded clicks."""
    return {
        "both_missed": int(histogram[0, 0]),
        "receiver_missed_sender_clicked": int(histogram[0, 1:].sum()),
        "receiver_clicked_sender_missed": int(histogram[1:, 0].sum()),
        "both_clicked": int(histogram[1:, 1:].sum()),
    }


def tag_outcomes(histogram, bit):
    """Compare the two time-bin letters with the launch bit after the read."""
    if bit not in (0, 1) or histogram.shape != (13, 13):
        raise ValueError("two six-tick slots and a binary launch required")
    correct = slice(1+6*bit, 7+6*bit)
    wrong = slice(7-6*bit, 13-6*bit)
    return dict(outcomes(histogram),
                wrong_receiver_letter=int(histogram[wrong, :].sum()),
                wrong_sender_letter=int(histogram[:, wrong].sum()),
                correct_both_letters=int(histogram[correct, correct].sum()),
                returned_receiver_letter=int(histogram[1:7,1:7].sum() + histogram[7:,7:].sum()),
                changed_receiver_letter=int(histogram[1:7,7:].sum() + histogram[7:,1:7].sum()),
                sender_correct_without_receiver=int(histogram[0, correct].sum()))


def report(arrays, detail=False):
    ranges = arrays["range_time_counts"]
    tags = arrays["tag_time_counts"]
    settings = []
    for pulse, bands in enumerate(distance_bands(ranges)):
        geometry = classify(bands, ((0,1), (2,3), (4,5)))
        pairs = [(i,j) for i in range(6) for j in range(6) if i != j]
        result = {
            "mark_pulse": ["sham", "pi"][pulse],
            "ranges": {f"{i}-{j}": [str(x) for x in band] if band is not None else None
                       for (i,j), band in bands.items()},
            "joint_confidence_failure_at_most": geometry["confidence_failure_probability"],
            "distortion_bracket": geometry["distortion_bracket"],
            "geometry": [
                {key: row[key] for key in ("relative_tolerance", "six_station_verdict", "excluded_subsets",
                                           "any_dimension_obstruction", "octahedron")}
                for row in geometry["thresholds"]],
            "range_totals": outcomes(ranges[pulse].sum(axis=(0,1))),
            "exact_source_comparison": _source_geometry(pulse, bands),
            "tag_totals": {},
        }
        for control, name in enumerate(("active", "parked", "detector_disabled")):
            result["tag_totals"][name] = [
                tag_outcomes(tags[pulse,:,:,control,bit].sum(axis=(0,1)), bit)
                for bit in (0,1)]
        if detail:
            result["geometry_certificate"] = geometry
            result["directed"] = [
                {"source": i, "receiver": j, "range": outcomes(ranges[pulse,i,j]),
                 "tags": {name: [tag_outcomes(tags[pulse,i,j,c,b], b) for b in (0,1)]
                          for c,name in enumerate(("active","parked","detector_disabled"))}}
                for i,j in pairs]
        settings.append(result)
    return settings


def main():
    from .six_station import verify_archive

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--details", action="store_true", help="include all directed reads and Gram certificates")
    parser.add_argument("--seed-replay", action="store_true", help="repeat the seeded acquisition as well")
    args = parser.parse_args()
    arrays = verify_archive(seed_replay=args.seed_replay)
    print(json.dumps(report(arrays, args.details), indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
