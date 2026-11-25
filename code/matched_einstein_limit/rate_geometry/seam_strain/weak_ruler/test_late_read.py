"""Rank, held-preparation transfer and confidence checks on the nine-clock experiment."""

from fractions import Fraction
from math import factorial
import json

import numpy as np
import pytest

from . import late_read as read


RECEIPT = read.flight.HERE / "late_read.json"


@pytest.fixture(scope="module")
def archive():
    return json.loads(RECEIPT.read_text())


@pytest.fixture(scope="module")
def rebuilt():
    return read.run()


def test_rebuilds_both_rank_certificates_and_the_minimal_front_supplement(
    archive, rebuilt
):
    # Re-run discovery independently of the frozen integer proof. Equivalent
    # pivots can exchange rows across floating-point implementations.
    _, calibration = read.menu()
    menu = read.matrix(calibration)
    proposed = read.select(menu)
    assert len(set(proposed)) == 9
    assert np.linalg.svd(menu[proposed], compute_uv=False)[-1] > 8e-7
    assert (
        rebuilt["settings"]
        == archive["settings"]
        == [read.setting(s) for s in read.SELECTED]
    )
    assert len({tuple(sorted(s.items())) for s in archive["settings"]}) == 9
    for out in (archive, rebuilt):
        primary = out["certified_least_singular_value_lower"]
        continuous = out["continuous_monitor"]
        assert out["rank"] == continuous["rank"] == 9
        assert 0 < primary < out["selected_singular_values"][-1]
        assert (
            0
            < continuous["least_singular_value_lower"]
            < continuous["singular_values"][-1]
        )
        assert (
            out["independent_first_hit_versus_survival_max_error"]
            < out["entry_error_upper"]
        )
        assert 8e4 < out["selected_condition"] < 9e4
        extra = out["three_bin_supplement"]
        assert extra["mixed_rank"] == 9 and len(extra["added_settings"]) == 3
        assert extra["mixed_least_singular_value_lower"] > 0
        assert len(extra["front_links"]) == 6
    # Every midpoint and independence margin below is derived from checked
    # integer propagation/rational inversion, independent of local BLAS.
    for k in (
        "response_matrix",
        "certified_least_singular_value_lower",
        "entry_error_upper",
    ):
        assert rebuilt[k] == archive[k]
    assert (
        rebuilt["continuous_monitor"]["response_matrix"]
        == archive["continuous_monitor"]["response_matrix"]
    )
    assert (
        rebuilt["three_bin_supplement"]["added_settings"]
        == archive["three_bin_supplement"]["added_settings"]
    )
    assert (
        rebuilt["three_bin_supplement"]["mixed_least_singular_value_lower"]
        == archive["three_bin_supplement"]["mixed_least_singular_value_lower"]
    )


def test_held_preparations_transfer_from_separate_positive_clock_simulation(
    archive, rebuilt
):
    for result in (archive, rebuilt):
        held = result["held"]
        full = json.loads((read.flight.HERE / "flight.json").read_text())
        panels = {p["charge"]: p for p in full["panels"] if p["divisor"] == 16}
        # Compare each positive reference probability with the independent
        # eight-station flight archive, which also counts every failed trial.
        for i, (q, a, b, t) in enumerate(read.SELECTED):
            expected = panels[q]["states"]["reference"]["first_hist"][a][b][t]
            assert abs(held["probabilities"]["reference"][i] - expected) < 1e-13
        for name, r in held["positive_preparations"].items():
            assert r["input_minimum_eigenvalue"] > 0
            assert r["maximum_marginal_change"] < 4e-13
            assert r["maximum_coordinate_recovery_error"] < 1e-8
            direct = r["direct_retained_modular_change"]
            assert abs(r["reconstructed_retained_modular_change"] - direct) < 1e-8
            assert (
                abs(r["reconstructed_entropy_change"] - r["actual_entropy_change"])
                < 1e-8
            )
            lo, hi = r["counted_modular_interval"]
            assert lo < direct < hi
            truth = archive["held"]["positive_preparations"][name]
            assert abs(direct - truth["direct_retained_modular_change"]) < 3e-11


def test_frozen_counts_use_every_launch_and_reproduce_the_modular_interval(archive):
    held = archive["held"]
    n = held["shots_per_bin"]
    width = held["count_mean_half_width"]
    assert n == 2**54 and held["simultaneous_confidence"] == 0.99
    for name, rows in held["counts"].items():
        p = held["probabilities"][name]
        for count, prob in zip(rows, p, strict=True):
            assert type(count) is int and 0 < count < n
            assert abs(count / n - prob) < width
    # A strict Bernstein check at 99 percent of the displayed width leaves
    # room for binary64 evaluation of the nominal Hoeffding radius. All
    # Bernoulli variances are <= p. For 36 settings the required exp >7200.
    pmax = max(Fraction(p) for row in held["probabilities"].values() for p in row)
    u = Fraction(width) * Fraction(99, 100)
    assert n * u**2 / (2 * (pmax + u / 3)) > 10
    assert sum(Fraction(10**k, factorial(k)) for k in range(15)) > 7200
    weights = np.array(archive["log_weights"])
    for name, r in held["positive_preparations"].items():
        observed = (
            np.array(held["counts"][name]) - np.array(held["counts"]["reference"])
        ) / n
        pred = weights @ observed
        lo, hi = r["counted_modular_interval"]
        assert abs((lo + hi) / 2 - pred) < 1e-12
        assert abs((hi - lo) / 2 - held["modular_half_width"]) < 1e-12
    ratio = np.linalg.norm(weights, 1)
    assert abs(ratio - archive["weight_l1"]) < 1e-8
    assert np.allclose(
        weights @ np.array(archive["response_matrix"]),
        archive["nominal_modular_covector"],
        rtol=1e-9,
        atol=1e-12,
    )


def test_larger_count_budget_resolves_quarter_turns_without_using_their_known_amplitude(
    archive, rebuilt
):
    for result in (archive, rebuilt):
        fine = result["fine_count_check"]
        name = fine["held"]
        held = result["held"]["positive_preparations"][name]
        assert (
            fine["shots_per_bin"] == 2**64 and fine["simultaneous_confidence"] == 0.99
        )
        assert (
            held["counted_modular_interval"][0]
            < 0
            < held["counted_modular_interval"][1]
        )
        lo, hi = fine["counted_modular_interval"]
        assert 0 < lo < held["direct_retained_modular_change"] < hi
        assert (
            np.linalg.norm(held["true_input_coordinates"], 1)
            < fine["coefficient_l1_upper_from_counts"]
            < 0.04
        )
        assert fine["modular_bias_upper_from_counts"] < 2e-7


def test_native_front_blind_lotteries_have_a_signed_independent_late_event(archive):
    readout = read.front_blind_read(archive)
    assert readout == archive["front_blind_read"]
    lo, hi = readout["blind_input_probability_slope"]
    assert -1.58e-6 < lo < hi < -1.56e-6
    lo, hi = readout["exact_native_lottery_probability_difference"]
    assert -6.09e-11 < lo < hi < -6.05e-11
