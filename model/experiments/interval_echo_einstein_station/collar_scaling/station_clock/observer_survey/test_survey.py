"""Independent counts on the tiling and on the native federations."""

from collections import Counter
from fractions import Fraction
import json

from .. import measure as station
from . import measure


def test_hexagonal_counts_and_local_state():
    # Direct two-sublattice graph, independent of the atlas and count formulas.
    frontier = {(0, 0, 0)}
    seen = set(frontier)
    for s in range(1, 36):
        nxt = set()
        for i, j, side in frontier:
            d = 1 if side else -1
            nxt.update(((i, j, 1 - side), (i + d, j, 1 - side), (i, j + d, 1 - side)))
        frontier = nxt - seen
        seen.update(frontier)
        assert len(frontier) == 3 * s
        assert len(seen) == measure.flat_ball(s)
    data = json.loads((measure.HERE / "results.json").read_text())
    assert data == measure.run()
    assert data["regular_pair_limit"] == str(Fraction(81, 40) ** 2)
    assert data["mean_regional_entropy_slope_nats"][0] > 0
    for level, hist in enumerate(data["native_radius_two_annuli"]):
        adjacency = measure.atlas(level)
        counted = Counter(
            len(measure.annulus(adjacency, i, 2)) for i in range(len(adjacency))
        )
        assert hist == {str(k): v for k, v in counted.items()}


def test_label_comparison_and_corner_leading_coefficients():
    triples = measure.gate_triples()
    assert set(triples) == {1, 2}
    assert all(4 not in t for t in triples.values())
    for c, expected in ((1, Fraction(243, 64)), (2, Fraction(27, 16))):
        w = 1 + c
        t = 4000
        ratio = Fraction(station.corner_ball(t) - station.corner_ball(t // 2), t**2)
        assert abs(ratio * Fraction(18**2, 20 * w**2) - expected) < Fraction(1, 400)
    assert all(
        row["changed_station_ranges"] > 0
        for row in json.loads((measure.HERE / "results.json").read_text())[
            "comparison_labels"
        ]
    )
