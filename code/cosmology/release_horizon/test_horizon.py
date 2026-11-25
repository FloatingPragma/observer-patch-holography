"""Native gate checks, causal cuts and verification of the deposited samples."""

import json
import math
import sys
from fractions import Fraction as Q
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from contrast import select_edges
from horizon import (
    HERE,
    centered_power,
    exact_table,
    moments,
    past_cones,
    selected_pairs,
)
from relay import STATES, gates, propagate, topology

PROTOCOL = json.loads((HERE / "protocol.json").read_text())
RESULT = json.loads((HERE / "results.json").read_text())


def test_native_registers_and_forbidden_writes():
    menu = gates()
    assert len(menu) == len(STATES)
    for x, d in STATES:
        options = menu[str((x, d))]
        for row in options:
            y, f = row["output"]
            assert f <= d and row["shared"] + row["private"] <= 1
        assert (1 - x, 1 - d) not in [r["output"] for r in options]
        if d == 0:
            assert (x, 1) not in [r["output"] for r in options]


def test_distance_deadline_and_trigger_tail():
    # A time-2 release at a neighboring receiver can be seen but cannot affect
    # its time-2 energy. Missing and receiving both have positive probability.
    tbl = exact_table("release", 2, 1, 1)
    assert sum(map(sum, tbl[2:])) == Q(1, 4)
    assert moments(tbl)[-1] == 0
    assert moments(exact_table("release", 3, 1, 1))[-1] > 0
    for h in (2, 4, 8):
        assert centered_power(np.array([2, 2, 2]), h) == 0
        assert moments(exact_table("shared_bit", h, 1, 4))[-1] == Q(1, 225 * 4**h)
        assert moments(exact_table("independent_bits", h, 1, 4))[-1] == 0


def test_finite_deadline_bounded_by_unlimited_covariance():
    for i, j in ((1, 1), (1, 4), (3, 4), (6, 17)):
        limit = Q(4, 24375 * 2 ** (i + j))
        seq = [moments(exact_table("release", h, i, j))[-1] for h in range(1, 36)]
        assert seq == sorted(seq) and seq[-1] <= limit
        assert limit - seq[-1] < limit / 100000


def test_graphs_and_pair_archive():
    for level, record in zip(PROTOCOL["levels"], RESULT["boards"], strict=True):
        board, groups, dist = topology(level)
        root = record["root"]
        np.testing.assert_array_equal(
            propagate(groups, board.cells, root, record["eccentricity"]), dist[root]
        )
        assert record["distance_by_cell"] == dist[root].tolist()
        for dl in record["deadlines"]:
            h = dl["sweeps"]
            cones = past_cones(groups, board.cells, h)
            disjoint, pairs = selected_pairs(dist, cones, root)
            assert dl["disjoint_pairs"] == disjoint
            # Every cell whose energy can depend on this post-reset trigger
            # has the root in its circuit past, on every measured board.
            for i, radius in enumerate(dist[root]):
                if radius < h - 1:
                    assert cones[i] & (1 << root)
            for sample in dl["pairs"]:
                i, j = pairs[sample["selection"]]
                assert sample["cells"] == [i, j]
                assert (
                    sample["common_initial_cells"] == (cones[i] & cones[j]).bit_count()
                )
                for bank, stats in sample["banks"].items():
                    age = h + record["eccentricity"] if bank == "earlier_release" else h
                    table = exact_table(bank, age, *sample["radii"])
                    assert stats["probabilities"] == [
                        [str(x) for x in row] for row in table
                    ]
                    counts = np.array(stats["counts"])
                    assert (
                        np.all(counts >= 0)
                        and counts.sum() == PROTOCOL["precision_originals_per_pair"]
                    )
                    measured = moments(counts / counts.sum())[-1]
                    assert math.isclose(
                        measured, stats["sampled_covariance"], abs_tol=1e-15
                    )
                    assert (
                        abs(float(moments(table)[-1]) - measured)
                        <= stats["covariance_half_width"]
                    )


def test_compensated_signal_against_four_pair_expansion():
    for level, rec in zip(PROTOCOL["levels"], RESULT["boards"], strict=True):
        board, groups, dist = topology(level)
        for dl in rec["deadlines"]:
            h = dl["sweeps"]
            archive = dl["local_differences"]
            total, cells = select_edges(
                groups, dist[0], past_cones(groups, board.cells, h), 0
            )
            assert total == archive["eligible_edge_pairs"]
            if cells is None:
                assert "banks" not in archive
                continue
            assert list(cells) == archive["cells"]
            r = dist[0, list(cells)].tolist()
            for mode, stats in archive["banks"].items():
                age = h + rec["eccentricity"] if mode == "earlier_release" else h
                cross = sum(
                    (-1) ** (i + j)
                    * moments(exact_table(mode, age, r[i], r[j + 2]))[-1]
                    for i in range(2)
                    for j in range(2)
                )
                assert cross == Q(stats["covariance"])
                counts = np.array(stats["counts"])
                assert (
                    len(counts) == 16
                    and np.all(counts >= 0)
                    and counts.sum() == PROTOCOL["contrast_trials"]
                )
                x = np.arange(16)
                a, b = (x & 1) - ((x >> 1) & 1), ((x >> 2) & 1) - ((x >> 3) & 1)
                p = counts / counts.sum()
                estimate = p @ (a * b) - (p @ a) * (p @ b)
                assert math.isclose(
                    estimate, stats["sampled_covariance"], abs_tol=1e-15
                )
                assert abs(estimate - float(cross)) <= stats["covariance_half_width"]
                if h == 2 and mode == "earlier_release":
                    assert estimate > stats["covariance_half_width"]
