"""Pair admission with the established independently tagged macroscopic faces.

The home station is the corner descendant of the first original face. A ball
with radius R <= 2m-2 meets only the five faces at its distinguished vertex.
For their cyclic offsets 0,1,1,2,2, let a_j=F(R-j)-F(R/2-j). A common tag is
available only when both distinct source cells are in the same original face.
For R=4k those five annular counts are 3k²+2k, 3k²+k twice, and 3k² twice.
Thus sum a_j(a_j-1) = 45k⁴+24k³-9k²-4k. Since N=5R⁴, ℓ=2/(bR) and b=H/m→6,
the ratio of the ordered within-face probability to ℓ⁴ tends to 729/1280.
The finite coefficients use the actual measured H and the original source's
ancestor-conditional reversed reference, with zero charge across ancestors.
"""

import argparse
from collections import deque
from fractions import Fraction
from itertools import combinations
import json
from pathlib import Path

import numpy as np

from . import measure


HERE = Path(__file__).parent


def offsets(faces):
    distinguished = faces[0][0]
    incident = [j for j, f in enumerate(faces) if distinguished in f]
    assert len(incident) == 5
    hops = {0: 0}
    q = deque([0])
    while q:
        i = q.popleft()
        for j in incident:
            if j not in hops and len(set(faces[i]) & set(faces[j])) == 2:
                hops[j] = hops[i] + 1
                q.append(j)
    assert sorted(hops.values()) == [0, 1, 1, 2, 2]
    return hops


def ancestor_counts(dist, parents, radius):
    # Count each original macroface from station enrollment, without using
    # the predicted sector populations or the conditional charge.
    return np.bincount(
        parents[(dist > radius // 2) & (dist <= radius)], minlength=20
    ).tolist()


def polynomial(k):
    return 45 * k**4 + 24 * k**3 - 9 * k**2 - 4 * k


def run():
    ledger = json.loads((HERE / "results.json").read_text())["midpoint_surveys"]
    saved = {row["level"]: row for row in ledger}
    vertices, seams = measure.icosahedron_from_coordinates()
    seams = set(seams)
    faces = [f for f in combinations(range(12), 3) if set(combinations(f, 2)) <= seams]
    corner = faces[0][0]
    turns = offsets(faces)
    parents = np.arange(20, dtype=int)
    panels = []
    native = {}
    pair_limit = Fraction(729, 1280)
    for level in range(max(saved) + 1):
        m = 2**level
        # Subdivision emits four children per parent, keeping its first
        # vertex in the first child. Follow the physical home to this corner.
        assert faces[0][0] == corner and len(parents) == len(faces) == 20 * m * m
        if level in saved:
            row = saved[level]
            radius = row["radius"]
            assert radius <= 2 * m - 2
            d = measure.dual_distances(faces)
            counts = ancestor_counts(d, parents, radius)
            h = int(d.max())
            size = len(d)
            assert (h, size, sum(counts)) == (
                row["eccentricity"],
                row["cells"],
                row["annulus_stations"],
            )
            predicted = [
                measure.corner_sector(radius - turns[j])
                - measure.corner_sector(radius // 2 - turns[j])
                if j in turns
                else 0
                for j in range(20)
            ]
            assert counts == predicted
            pairs = sum(a * (a - 1) for a in counts)
            if radius % 4 == 0:
                assert pairs == polynomial(radius // 4)
            ell = Fraction(radius, h)
            norm = Fraction(pairs, size * (size - 1)) / ell**4
            if radius % 4 == 0:
                b, t = Fraction(h, m), Fraction(1, radius)
                expected = (
                    b**4
                    * (
                        Fraction(45, 256)
                        + Fraction(3, 8) * t
                        - Fraction(9, 16) * t * t
                        - t**3
                    )
                    / (400 * (1 - t**4 / 5))
                )
                assert expected == norm
            panels.append(
                dict(
                    level=level,
                    frequency=m,
                    radius=radius,
                    eccentricity=h,
                    ordered_pair_count=pairs,
                    ancestor_populations=counts,
                    normalized_pair=str(norm),
                )
            )
        if level in (1, 2):
            # Test the ancestor counts separately on the original native
            # external contacts up through the local bundled 320-cell board.
            board = measure.federation(level)
            nd = measure.scan(
                measure.graph(np.asarray(board.external) // 12, board.cells)
            )
            test_radius = 2 * m - 2
            nc = ancestor_counts(nd, np.arange(board.cells) // (m * m), test_radius)
            # The bundled parent names may be permuted relative to the atlas.
            predicted = sorted(
                measure.corner_sector(test_radius - j)
                - measure.corner_sector(test_radius // 2 - j)
                for j in turns.values()
            )
            assert sorted([n for n in nc if n]) == [x for x in predicted if x]
            native[str(level)] = dict(radius=test_radius, counts_by_native_ancestor=nc)
        if level < max(saved):
            vertices, faces = measure._subdivide(vertices, faces)
            parents = np.repeat(parents, 4)
    charge, lower = measure.clock_balls.source_contrast()
    return dict(
        ancestor_tags="independent fair tag per original face; ancestor-conditional reversed reference",
        corner_vertex=corner,
        sector_offsets={str(j): v for j, v in sorted(turns.items())},
        midpoint_surveys=panels,
        bundled_native_surveys=native,
        single_region_pair_limit=str(pair_limit),
        fraction_of_global_pair_limit=str(pair_limit / Fraction(729, 256)),
        collar_contrast_per_fourth_power_limit_nats=charge * float(pair_limit),
        strict_rational_lower_nats=str(lower * pair_limit),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    opts = parser.parse_args()
    result = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if opts.write:
        (HERE / "ancestor_limit.json").write_text(result)
    else:
        print(result, end="")
