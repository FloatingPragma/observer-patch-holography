"""Independent identities for conserved source counts and native regional energy."""

from fractions import Fraction
from itertools import product
import json

import numpy as np

from . import exact, measure, modular, reservoir


def test_exact_native_projection_against_flight_sectors():
    saved = json.loads((measure.HERE / "exact.json").read_text())
    assert exact.run() == saved
    for row in saved:
        panel = measure.run_case(1, row["stiffness"])
        profile = panel["average"]["mean"]
        for q, s in zip(row["shells"], profile, strict=True):
            for key in ("energy", "load"):
                np.testing.assert_allclose(s[key], float(Fraction(q[key])), atol=3e-12)
        assert Fraction(row["shells"][0]["energy"]) < 0
        assert Fraction(row["shells"][2]["energy"]) > 0
        np.testing.assert_allclose(
            panel["ordinary_single_load"]["mean"][-1]["energy"],
            5 * (row["stiffness"] + 0.5),
        )


def test_counts_against_enumeration_and_first_resonance():
    histogram = {}
    for word in product(((0, 0), (8, 1), (9, 1), (0, 2)), repeat=5):
        e, q = map(sum, zip(*word, strict=True))
        histogram[e, q] = histogram.get((e, q), 0) + 1
    for e, q in product(range(46), range(11)):
        assert reservoir.counts(5, e, q) == histogram.get((e, q), 0)
    assert reservoir.counts(18, 144, 18) == 1 + 18 * 17


def test_energy_matched_source_and_independent_calibration():
    data = reservoir.run()
    for a in data["held_equal_energy_challenge"]:
        assert a["modular_difference"] > 0.22
        np.testing.assert_allclose(
            a["energy_difference_inferred_from_masked_log"], 8 / 9
        )
    assert all(
        r["physical_difference"] == 0 and r["predicted_difference_from_masked_log"] < 0
        for r in data["material_endpoints"]
    )
    highest = data["contact_bank"][-1]
    assert highest["fractional_error_billionths"] == [44487, 44488]
    for s in highest["source"]:
        assert 44487 < 1e9 * s["fractional_error"] < 44488
    for s in modular.run():
        for c in s["free_contacts"]:
            for label in ("formed", "independent"):
                v = c[label]
                assert not v["singular"]
                assert abs(v["beta_even_to_odd"] - v["beta_equal_to_even"]) > 0.016


def test_physical_reference_full_matrix_and_rejected_preparations():
    from . import herald

    receipt = herald.run()
    for row in receipt["panels"]:
        for starts, stops, kept in zip(
            row["ordinary_inputs"], row["stops_by_gate"], row["survivors"], strict=True
        ):
            assert starts == sum(stops) + kept
        if row["design"] == "ordinary energy":
            assert row["full_generator_residual"] < 1e-14
            assert abs(row["matched_banks_modular"]) < 1e-14
        else:
            assert abs(row["matched_banks_modular"]) > 0.1
    assert (
        6.6e-6
        < receipt["microcanonical_selection"][-1][
            "probability_totals_given_accepted_contacts"
        ]
        < 6.7e-6
    )
