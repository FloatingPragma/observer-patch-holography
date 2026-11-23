"""Independent geometry, regional modular charge, completeness and locality checks."""

from itertools import permutations
import json

import numpy as np
import pytest
from scipy.sparse.csgraph import shortest_path

from model.executable_model import federation
from model.executable_model.refinement import Refinement
from . import ancestor, measure
from .. import assert_same_replay


def test_protected_reduction_matches_every_reached_station_on_direct_board():
    board = federation(1)
    compressed = measure.ProtectedRanges(
        Refinement(1, board.cells, "", np.asarray(board.external))
    )
    stations = np.arange(0, board.n, 12)
    kept = np.arange(board.n)[np.arange(board.n) % 12 != 4]
    graph = board.adjacency()[kept][:, kept]
    for origin in (0, 37, 79):
        distances = shortest_path(
            graph,
            directed=False,
            unweighted=True,
            indices=np.searchsorted(kept, 12 * origin),
        )
        np.testing.assert_array_equal(
            compressed.survey(origin), distances[np.searchsorted(kept, stations)]
        )


def test_tagged_admission_equals_direct_ordered_pair_census():
    table = [[2, 1, 1], [1, 0, 2], [0, 1, 1]]
    cells = [
        (g, j)
        for g, counts in enumerate(table)
        for j, count in enumerate(counts)
        for _ in range(count)
    ]
    exact = np.zeros((2, 3, 3))
    for (group, a), (other, b) in permutations(cells, 2):
        exact[int(group != other), a, b] += 1
    np.testing.assert_allclose(
        ancestor.tagged_slots(table), exact / exact.sum(), atol=2e-16, rtol=0
    )


def test_charge_difference_equals_full_reference_modular_expectation():
    laws, _, _ = measure.source_reader()
    fractions = measure.slot_membership([4, 7, 11], False)
    score = measure.modular_score(laws[0])
    # For an observer reading just A, assemble the 100 complete outcomes in the
    # same algebra (0, 1st only, 2nd only, both). Probability of admitted flags is
    # independent of history. The shared and fresh one-word marginals coincide.
    p = fractions[0, 0]
    left = fractions[0, 1:].sum()
    right = fractions[1:, 0].sum()

    def region(law):
        return np.r_[
            1 - p - left - right,
            left * law.sum(axis=1),
            right * law.sum(axis=0),
            p * law.ravel(),
        ]

    direct = -np.log(region(laws[0])) @ (region(laws[2]) - region(laws[1]))
    np.testing.assert_allclose(
        direct, p * np.sum((laws[2] - laws[1]) * score), rtol=1e-12
    )
    # Annihilation of additive single-word energies across the complete word law.
    sums = np.add.outer([0, 32, 33], [0, 32, 33]).ravel()
    for a in np.unique(sums):
        for b in np.unique(sums):
            mask = (sums == a)[:, None] & (sums == b)[None, :]
            np.testing.assert_allclose(
                laws[2][mask].sum(), laws[1][mask].sum(), atol=2e-16
            )


def test_raw_sample_completeness_and_frozen_controls():
    primary = json.loads((measure.HERE / "results.json").read_text())
    measure.verify_replay(measure.experiment([r["geometry"] for r in primary["panels"]]), primary, name="reads")
    for panel in primary["panels"]:
        for row in panel["settings"]:
            for study in row["arms"].values():
                for arm in study["reads"]:
                    assert sum(sum(bins) for bins in arm.values()) == measure.SHOTS
                for q in study["contrasts"].values():
                    expected = q["population_shared_change"]
                    radius = q["simultaneous_radius"]
                    assert abs(q["measured_shared_change"] - expected) <= radius
                    assert abs(q["measured_reverse_change"] + expected) <= radius
    local = json.loads((measure.HERE / "ancestor.json").read_text())
    ancestor.verify_replay(ancestor.experiment([r["geometry"] for r in local["panels"]]), local, name="readings")
    for panel in local["panels"]:
        for setting, data in zip(
            panel["geometry"]["settings"], panel["arms"], strict=True
        ):
            for arm in data["readings"]:
                assert sum(sum(bins) for bins in arm.values()) == measure.SHOTS
            cross = data["contrasts"]["cross"]
            if setting["partner"] == "far":
                assert cross["shared"] == cross["expected"] == 0
            else:
                assert cross["shared"] > cross["simultaneous_radius"]
            for q in data["contrasts"].values():
                assert abs(q["shared"] - q["expected"]) <= q["simultaneous_radius"]
                assert abs(q["reversed"] + q["expected"]) <= q["simultaneous_radius"]


def test_quartic_transfer_and_colocated_rejection():
    data = json.loads((measure.HERE / "results.json").read_text())["panels"]

    def read(index, scale, arm):
        return data[index]["settings"][scale]["arms"][arm]["contrasts"]["home"][
            "normalized_population"
        ]

    assert abs(read(-1, 1, "separated") / read(-1, 0, "separated") - 1) < 0.2
    assert abs(read(-1, 0, "separated") / read(-2, 0, "separated") - 1) < 0.2
    assert abs(read(-1, 1, "colocated") / read(-1, 0, "colocated") - 1) > 0.2
    assert abs(read(-1, 0, "colocated") / read(-2, 0, "colocated") - 1) > 0.2
    a = json.loads((measure.HERE / "ancestor.json").read_text())["panels"]

    def near(index, scale):
        return a[index]["arms"][scale]["contrasts"]["home"]["normalized"]

    assert (
        abs(near(-1, 1) / near(-1, 0) - 1) < 0.2
        and abs(near(-1, 0) / near(-2, 0) - 1) < 0.2
    )


def test_same_total_single_trial_reassignment_is_rejected():
    raw = (measure.HERE / "results.json").read_bytes()
    original, changed = json.loads(raw), json.loads(raw)
    counts = changed["panels"][0]["settings"][0]["arms"]["separated"]["reads"][0]["00"]
    counts[0] -= 1
    counts[1] += 1
    for fresh, archive in ((changed, original), (original, changed), (changed, changed)):
        with pytest.raises(ValueError, match="count stream"):
            measure.verify_replay(fresh, archive, name="reads")


@pytest.mark.parametrize("value,stored", [(True, 1), (1.0, 1), (True, 1.0), (float("inf"), float("inf"))])
def test_numeric_type_and_finiteness_are_part_of_the_receipt(value, stored):
    with pytest.raises(AssertionError):
        assert_same_replay({"value": value}, {"value": stored})


def test_sharp_energy_controls_and_radial_mark_density():
    from . import source_scale

    record = source_scale.evaluate()
    assert_same_replay(record, json.loads((measure.HERE / "source_scale.json").read_text()))
    low, high = (record["energy_conditioned"][i] for i in (0, 2))
    assert sum(low["word_energy_sums"]) == 64 and sum(high["word_energy_sums"]) == 130
    assert [r["population"] for r in low["results"]] == [
        r["population"] for r in high["results"]
    ]
    assert (
        low["results"][0]["population"]
        < low["results"][1]["population"]
        < low["results"][2]["population"]
    )
    for row in (low, high):
        for arm in row["results"]:
            assert (
                abs(arm["measured"] - arm["population_conditional_on_reference"])
                < 5 * arm["standard_error"]
            )
    inner, outer = record["mark_density"][-1]["measurements"]
    assert (
        inner["distinct_occupied_shells_per_routed_length_bounds"][1]
        < outer["distinct_occupied_shells_per_routed_length_bounds"][0]
    )
