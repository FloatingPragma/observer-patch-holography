"""Exhaust the balanced address partitions for the three-receiver colour read."""

import argparse
from collections import Counter
from fractions import Fraction
from itertools import combinations
import json
import math
from pathlib import Path

from model.candidate_laws.interval_echo import currents
from model.executable_model import cell


def partitions():
    ports = set(range(12))
    for others in combinations(range(1, 12), 3):
        first = (0, *others)
        rest = sorted(ports.difference(first))
        for middle in combinations(rest[1:], 3):
            second = (rest[0], *middle)
            yield first, second, tuple(sorted(set(rest).difference(second)))


def lower_in_256ths(addresses, antipodes):
    # Bound |P_ab| in eighths by 2 for antipodes and 1 for other distinct ports.
    # The off-diagonal effect coefficients are 3 P_ab/8 or -3 P_ab P_cd/4.
    radii = [0] * len(addresses)
    for (i, x), (j, y) in combinations(enumerate(addresses), 2):
        off = [2 if antipodes[a] == b else 1 for a, b in zip(x, y, strict=True) if a != b]
        if len(off) == 1:
            weight = 12 * off[0]
        elif len(off) == 2:
            weight = 3 * off[0] * off[1]
        else:
            continue
        radii[i] += weight
        radii[j] += weight
    return 108 - max(radii)


def run():
    car = cell()
    faces = car.faces()
    antipodes, _ = currents._geometry(car)
    counts, bounds = Counter(), Counter()
    worst = None
    for groups in partitions():
        slot = {p: j for j, group in enumerate(groups) for p in group}
        addresses = [sorted(f, key=slot.__getitem__) for f in faces if len({slot[p] for p in f}) == 3]
        counts[len(addresses)] += 1
        if addresses:
            lower = lower_in_256ths(addresses, antipodes)
            bounds[lower] += 1
            if worst is None or lower < worst['bound_in_256ths']:
                worst = dict(groups=groups, accepted_addresses=addresses, bound_in_256ths=lower)
    maximum = max(counts)
    # det(P restricted to three mutually neighbouring ports) / 3! .
    weight = (1 + 1 / math.sqrt(5)) / 960
    return dict(partitions=sum(counts.values()), face_count_histogram=dict(sorted(counts.items())),
                lower_bound_histogram_in_256ths=dict(sorted(bounds.items())), worst=worst,
                common_frame_singlet_probability_per_face=weight,
                max_face_probability_common_frames=maximum * weight,
                face_probability_bound_any_sharp_destination_frames=str(Fraction(maximum, 384)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    result = run()
    with args.output.open('x') as f:
        json.dump(result, f, indent=2)
        f.write('\n')
