"""Independent port contacts, corner layers and the finite normalized source."""

from fractions import Fraction
import json

import pytest
from scipy.sparse.csgraph import shortest_path

from model.executable_model import federation
from . import measure


def test_material_station_services_and_independent_native_ball_counts():
    for level in range(3):
        board = federation(level)
        result = measure.hardware(level)
        local = measure.cell_paths()
        # Independent graph distances on the actual protected template.
        template = measure.cell().adjacency().toarray()
        template[4, :] = template[:, 4] = 0
        distances = shortest_path(template, unweighted=True)[0]
        assert all(distances[p] == len(path) - 1 for p, path in local.items())
        assert int(max(len(path) for path in local.values())) == 4
        # Direct frontier expansion from the native external list.
        neighbors = [set() for _ in range(board.cells)]
        for a, b in board.external:
            neighbors[a // 12].add(b // 12)
            neighbors[b // 12].add(a // 12)
        layer = visited = {0}
        histogram = []
        while layer:
            histogram.append(len(layer))
            layer = set.union(*(neighbors[c] for c in layer)) - visited
            visited = visited | layer
        assert sum(histogram) == board.cells
        assert histogram == result["distance_histogram"]
        assert sum(result["spare_ticks"].values()) == 3 * board.cells
        m = 2**level
        for radius in range(2 * m - 1):
            assert sum(histogram[: radius + 1]) == measure.corner_ball(radius)
    standard, changed = measure.transfer_tables()
    assert standard == pytest.approx([1, 1], abs=1e-12)
    assert changed == pytest.approx([1, 0], abs=1e-12)


def test_exact_survey_and_limit_normalization():
    stored = json.loads((measure.HERE / "results.json").read_text())
    rebuilt = json.loads(json.dumps(measure.run()))
    for key in (
        "hardware",
        "midpoint_surveys",
        "single_limit",
        "pair_limit",
        "strict_rational_lower_nats",
    ):
        assert rebuilt[key] == stored[key]
    assert stored["conditional_collar_contrast_nats"] == pytest.approx(
        rebuilt["conditional_collar_contrast_nats"], rel=1e-12
    )
    assert Fraction(stored["single_limit"]) ** 2 == Fraction(729, 256)
    for row in stored["midpoint_surveys"]:
        n, size, r, h = (
            row[k] for k in ("annulus_stations", "cells", "radius", "eccentricity")
        )
        m = 2 ** row["level"]
        assert 6 * m - 3 <= h <= 6 * m
        if r % 4 == 0:
            assert n == 15 * r * r // 16 + r
            t, b = Fraction(1, r), Fraction(h, m)
            expected = (
                (Fraction(15, 16) + t)
                * (Fraction(15, 16) + t - t * t)
                * b**4
                / (400 * (1 - t**4 / 5))
            )
            assert expected == Fraction(row["normalized_pair"])
        ell = Fraction(r, h)
        assert Fraction(n * (n - 1), size * (size - 1)) / ell**4 == Fraction(
            row["normalized_pair"]
        )
    lower = Fraction(101, 33696) * Fraction(729, 256)
    assert Fraction(stored["strict_rational_lower_nats"]) == lower
    assert (
        float(lower) < stored["collar_contrast_per_fourth_power_limit_nats"] < 0.00942
    )
