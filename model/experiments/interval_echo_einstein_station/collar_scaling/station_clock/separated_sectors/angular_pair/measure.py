"""Compare addressed entropy with radial and angular energy records on one source."""

import argparse
from collections import Counter
from fractions import Fraction
from itertools import combinations
import json
from pathlib import Path

import numpy as np
from scipy.sparse.csgraph import shortest_path

from ... import measure as clock
from ... import word_entropy
from .. import measure as sectors

HERE = Path(__file__).resolve().parent
RADII = tuple(range(2, 7))
REFERENCE = Fraction(1, 4)
TRIALS = 2**46


def signature(distance, block, far):
    return Counter(
        (int(distance[0, x]), int(distance[0, y]), int(distance[x, y]))
        for x in range(block)
        for y in range(block * far, block * (far + 1))
    )


def native_regions():
    board = clock.federation(2)
    graph = clock.graph(np.asarray(board.external) // 12, board.cells)
    dist = shortest_path(graph, directed=False, unweighted=True).astype(int)
    block = board.cells // 20
    left, right = (signature(dist, block, face) for face in (2, 3))
    if left != right:
        raise ValueError("a three-range bin resolves the reflected faces")
    total = board.cells * (board.cells - 1)
    rows = []
    for radius in RADII:
        pairs = sum(n for (x, y, _), n in left.items() if max(x, y) <= radius)
        populations = [
            int(np.sum(dist[0, i * block : (i + 1) * block] <= radius))
            for i in (0, 2, 3)
        ]
        if pairs != populations[0] * populations[1] or populations[1] != populations[2]:
            raise ValueError("range signatures disagree with regional admissions")
        rows.append((radius, populations, Fraction(2 * pairs, total)))
    return int(dist[0].max()), len(left), rows


def atlas_reflection():
    vertices, edges = clock.icosahedron_from_coordinates()
    graph = clock.graph(np.asarray(edges), len(vertices))
    distance = shortest_path(graph, directed=False, unweighted=True).astype(int)
    signature_map = {tuple(row[:3]): i for i, row in enumerate(distance)}
    perm = [signature_map[a, c, b] for a, b, c in distance[:, :3]]
    if not np.array_equal(distance[np.ix_(perm, perm)], distance):
        raise ValueError("vertex reflection fails adjacency")
    edges = set(edges)
    faces = [f for f in combinations(range(12), 3) if set(combinations(f, 2)) <= edges]
    mirror = [faces.index(tuple(sorted(perm[i] for i in face))) for face in faces]
    if (perm[0], mirror[0], mirror[3], mirror[4]) != (0, 0, 4, 3):
        raise ValueError("reflection misses the corner or oriented sectors")
    return perm


def exact_derivative(groups):
    # At t=1/4, grouped logarithms have integer exponents in units of 1/144.
    base = Fraction(1)
    for p, d, count in groups:
        exponent, rem = divmod(d * count, 18)
        if rem:
            raise ArithmeticError("unexpected log exponent")
        base *= ((p - REFERENCE * d) / (p + REFERENCE * d)) ** exponent
    return base, word_entropy.IV.ln(word_entropy.interval(base)) / 144


def table(source, groups, angular):
    if angular:
        panel = sectors.word_panel()
        orientations = np.asarray(panel["orientation_on_energy_labels"]).reshape(9)
        score = np.outer(orientations, orientations)
        _, _, _, native = sectors.volume.radar_source_and_hand(clock.federation(0))
        raw, _ = sectors.shared_variance.source_words(native)
        reverse, fresh, _ = np.rint(2592 * raw).astype(int)
        law = fresh + float(source) * (fresh - reverse)
        return [float(np.sum(law[score == i]) / 2592) for i in (-1, 0, 1)]
    return [float(c * (p + source * d) / 2592) for p, d, c in groups]


def draw_bank(rng, rows, source, groups, angular=False):
    # One stream of attempts feeds all five nested ball reads. Each admitted
    # pair is assigned the first registered ball in which both stations enter.
    eligible = [row[2] for row in rows]
    shells = [p - prev for p, prev in zip(eligible, [0] + eligible[:-1], strict=True)]
    location = [float(x) for x in shells for _face in (1, -1)]
    admissions = rng.multinomial(TRIALS, [*location, float(1 - 2 * eligible[-1])])
    laws = [table(sgn * source, groups, angular) for sgn in (1, -1)]
    records = [
        rng.multinomial(int(n), laws[i % 2]).tolist()
        for i, n in enumerate(admissions[:-1])
    ]
    return {"first_admitted_shell_AB_AC": records, "other": int(admissions[-1])}


def aggregate(bank, index, face):
    return np.sum(bank["first_admitted_shell_AB_AC"][face : 2 * index + 2 : 2], axis=0)


def estimates(counts):
    total = int(sum(counts))
    radius = word_entropy.IV.sqrt(word_entropy.interval(Fraction(8, total)))
    error = word_entropy.IV.mpf([-1, 1]) * radius
    return [word_entropy.interval(Fraction(int(k), total)) + error for k in counts]


def sampled_slopes(banks, index, norm):
    # Two disjoint entropy-direction arms, in addition to the reference and
    # angular banks. The class multiplicities cancel inside the log ratio.
    ref_a, ref_c = (estimates(aggregate(banks["reference"], index, f)) for f in (0, 1))
    plus, minus = (
        estimates(aggregate(banks[arm], index, 0)) for arm in ("plus", "minus")
    )
    if any(x.a <= 0 for x in (*ref_a, *ref_c)):
        raise ArithmeticError("unresolved positive reference probabilities")
    slope = sum(
        (hi - lo) * word_entropy.IV.ln(c / a)
        for hi, lo, a, c in zip(plus, minus, ref_a, ref_c, strict=True)
    )
    oriented = []
    for f in (0, 1):
        neg, _, pos = estimates(aggregate(banks["angular"], index, f))
        oriented.append(pos - neg)
    gain = (oriented[0] - oriented[1]) / word_entropy.interval(REFERENCE)
    return {
        "addressed_entropy_99": word_entropy.endpoints(norm * slope),
        "angular_energy_score_99": word_entropy.endpoints(norm * gain),
    }


def run():
    height, bins, rows = native_regions()
    groups = word_entropy.probabilities()
    base, ds = exact_derivative(groups)
    rng = np.random.Generator(np.random.PCG64(134611))
    banks = {
        k: draw_bank(rng, rows, t, groups, k == "angular")
        for k, t in (
            ("reference", REFERENCE),
            ("plus", Fraction(1, 2)),
            ("minus", Fraction(-1, 2)),
            ("angular", REFERENCE),
        )
    }
    nested = []
    for i, (r, pop, prob) in enumerate(rows):
        scale = prob / Fraction(r, height) ** 4
        norm = word_entropy.interval(scale)
        nested.append(
            {
                "R": r,
                "populations_ABC": pop,
                "eligible_per_sector_pair": str(prob),
                "entropy_slope_per_fourth": word_entropy.endpoints(norm * ds),
                "angular_energy_slope_per_fourth": str(scale / 18),
                **sampled_slopes(banks, i, norm),
            }
        )
    surveyed = json.loads((sectors.HERE / "results.json").read_text())
    midpoint = []
    for row in surveyed["surveys"]:
        if row["region"] == "ball":
            scale = Fraction(row["normalized_eligible"])
            midpoint.append(
                {
                    "level": row["level"],
                    "radius": row["radius"],
                    "entropy_slope_per_fourth": word_entropy.endpoints(
                        word_entropy.interval(scale) * ds
                    ),
                    "angular_energy_slope_per_fourth": str(scale / 18),
                }
            )
    limit = Fraction(surveyed["limits"]["ball"]["eligible_per_fourth"])
    return {
        "bundled_cells": 320,
        "source_reference_t": str(REFERENCE),
        "clock_appointments_per_hop": clock.CLOCK_TICKS,
        "bundled_corner_eccentricity": height,
        "equal_full_range_triple_bins_per_sector_pair": bins,
        "midpoint_vertex_reflection": atlas_reflection(),
        "conditional_entropy_derivative_log_base_per_144": str(base),
        "conditional_angular_derivative": "1/18",
        "trials_per_independent_bank": TRIALS,
        "terminal_outcomes": banks,
        "fixed_board_nested_readings": nested,
        "midpoint_readings": midpoint,
        "entropy_slope_per_fourth_limit": word_entropy.endpoints(
            word_entropy.interval(limit) * ds
        ),
        "angular_energy_slope_per_fourth_limit": str(limit / 18),
    }


if __name__ == "__main__":
    arg = argparse.ArgumentParser(description=__doc__)
    arg.add_argument("--write", action="store_true")
    opts = arg.parse_args()
    report = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if opts.write:
        (HERE / "results.json").write_text(report)
    else:
        print(report, end="")
