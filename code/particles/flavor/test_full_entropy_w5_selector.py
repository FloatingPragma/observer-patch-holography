"""Reject false minima, altered geometry and physical promotion of the selector."""
from __future__ import annotations

import copy
from decimal import Decimal, localcontext
import json
from pathlib import Path
import sys

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import full_entropy_w5_selector as producer
import verify_full_entropy_w5_selector as independent


@pytest.fixture(scope="module")
def receipt():
    return json.loads(producer.OUT_PATH.read_text())


def test_rebuild_and_independent_verification(receipt):
    assert producer.build_payload() == receipt
    checked = independent.verify_payload(receipt)
    assert checked["stationary_branches_checked"] == 16
    assert checked["independent_axis_quadrupoles"] == 6


@pytest.mark.parametrize("path,replacement", [
    (("selector", "minimizer", "high"), "sqrt(15)*r/18 + 1/6"),
    (("selector", "minimizer", "high_multiplicity"), 2),
    (("selector", "purity"), "r**2/36 + 1/6"),
    (("selector", "amplitude_domain"), "all real r"),
    (("selector", "zero_amplitude_boundary"), "all normalized directions are axial at r=0"),
    (("selector", "stationary_branches", 2, "high_amplitude"), "1/2"),
    (("selector", "stationary_branches", 1, "minimum_classification"), "unique_global_minimum_up_to_permutation"),
    (("selector", "boundary_certificate", "equal_support_rows", 0, "delta_squared"), "epsilon/2"),
    (("geometry", "winning_axis_rows", 0, "exact_Q_entries", 0, 0), "0"),
    (("geometry", "winning_axis_rows", 0, "discriminant"), "1"),
    (("geometry", "winning_axis_rows", 0, "characteristic_polynomial"), "lambda**3 - lambda"),
    (("geometry", "eigenvalues", 1), "-sqrt(15)/15"),
    (("geometry", "minimum_adjacent_gap"), "1"),
    (("scope", "all_OPH_completions_excluded"), True),
    (("scope", "core_axiom_derivation_claimed"), True),
    (("scope", "source_selected_physical_readout"), True),
    (("scope", "empirical_target_consumed"), True),
    (("scope", "extra_continuous_fitted_coefficients"), 1),
    (("physical_boundary", "quark_ratios_emitted"), True),
    (("physical_boundary", "absolute_masses_emitted"), True),
    (("physical_boundary", "positive_physical_evidence_claimed"), True),
    (("physical_boundary", "prediction_freeze_claimed"), True),
    (("source_pins", 0, "path"), "../outside-source.json"),
    (("source_pins", 1, "sha256"), "0" * 64),
    (("proof_grade",), "formal_proof_of_all_OPH"),
])
def test_independent_verifier_rejects_false_green_mutations(receipt, path, replacement):
    mutated = copy.deepcopy(receipt)
    cursor = mutated
    for key in path[:-1]:
        cursor = cursor[key]
    cursor[path[-1]] = replacement
    with pytest.raises((ValueError, KeyError)):
        independent.verify_payload(mutated)


@pytest.mark.parametrize("path,suffix", [
    # These interpolants vanished at EVERY sample used by the old verifier.
    (("selector", "minimizer", "high"), "+r*(r-1)*(r-3)*(r-7)*(r**2-60)"),
    (("selector", "stationary_branches", 0, "high_probability"), "+r*(r-1)"),
    (("selector", "boundary_certificate", "equal_support_rows", 0, "delta_squared"),
     "+epsilon*(1000*1000*epsilon-1)"),
    (("geometry", "characteristic_polynomial"), "+lambda*(lambda-1)*(lambda-2)*(lambda-3)"),
    (("geometry", "winning_axis_rows", 0, "characteristic_polynomial"),
     "+lambda*(lambda-1)*(lambda-2)*(lambda-3)"),
    # Increasing sample counts or precision does not authenticate exact formulas.
    (("selector", "minimizer", "high"), "+r/(100000**6*100000**6*100000**6)"),
    (("selector", "purity"), "+r/(100000**6*100000**6*100000**6)"),
    (("geometry", "winning_axis_rows", 0, "exact_Q_entries", 0, 0),
     "+1/(100000**6*100000**6*100000**6)"),
])
def test_prescribed_normal_forms_reject_interpolation_and_tiny_corrections(receipt, path, suffix):
    mutated = copy.deepcopy(receipt)
    cursor = mutated
    for key in path[:-1]:
        cursor = cursor[key]
    cursor[path[-1]] += suffix
    with pytest.raises(ValueError, match="noncanonical"):
        independent.verify_payload(mutated)


@pytest.mark.parametrize("path,replacement", [
    (("selector", "differential_certificate", "repeated_high_second_variation"), "1"),
    (("selector", "differential_certificate", "log_secant_gap_derivative"), "-1"),
    (("selector", "differential_certificate", "strict_sign_domain"), "all a, b"),
    (("selector", "boundary_certificate", "unequal_support_purity_residual"), "0"),
    (("selector", "boundary_certificate", "unequal_support_implicit_derivative_at_origin"), "0"),
    (("selector", "amplitude_map"), "P already selects the physical amplitude"),
    (("geometry", "winning_axis_scope"), "all normalized directions are axial at r=0"),
    (("geometry", "single_spectral_readout_boundary"), "all OPH readouts are excluded"),
    (("geometry", "forward_map_readback"), "n_i^T Q n_i=a_i"),
    (("geometry", "unoriented_axis_stabilizer"), "C5 of order 5"),
    (("scope", "source_law"), "source-derived universal quark mass law"),
])
def test_displayed_proof_and_scope_cannot_drift_while_flags_stay_true(receipt, path, replacement):
    mutated = copy.deepcopy(receipt)
    cursor = mutated
    for key in path[:-1]:
        cursor = cursor[key]
    cursor[path[-1]] = replacement
    with pytest.raises(ValueError):
        independent.verify_payload(mutated)


def test_missing_branch_is_not_exhaustive(receipt):
    mutated = copy.deepcopy(receipt)
    mutated["selector"]["stationary_branches"].pop()
    with pytest.raises(ValueError, match="branch enumeration"):
        independent.verify_payload(mutated)


def test_duplicate_axis_cannot_replace_an_unchecked_axis(receipt):
    mutated = copy.deepcopy(receipt)
    row = copy.deepcopy(mutated["geometry"]["winning_axis_rows"][0])
    row["axis_index"] = 1
    mutated["geometry"]["winning_axis_rows"][1] = row
    with pytest.raises(ValueError, match="duplicate"):
        independent.verify_payload(mutated)


def test_quartic_support_two_boundary_does_not_minimize_full_entropy():
    # This is the explicit escape left by the old quartic truncation at r=5.
    # It obeys the same moments but loses under the actual logarithmic entropy.
    with localcontext() as context:
        context.prec = 80
        r = Decimal(5)
        purity = Decimal(1) / 6 + r * r / 72
        quartic = [Decimal(7) / 12, Decimal(5) / 12] + [Decimal(0)] * 4
        full = independent.stationary_profile(purity, 1)
        assert abs(independent.moments(quartic)[1] - purity) < Decimal("1e-70")
        assert independent.entropy(full) < independent.entropy(quartic)


def test_zero_amplitude_does_not_select_an_angular_direction(receipt):
    with localcontext() as context:
        context.prec = 80
        for count in range(1, 6):
            assert independent.stationary_profile(Decimal(1) / 6, count) == [Decimal(1) / 6] * 6
    assert "no direction is selected" in receipt["selector"]["zero_amplitude_boundary"]
    assert receipt["physical_boundary"]["actual_record_quadrupole"] == "r*Q"


@pytest.mark.parametrize("bad", ["__import__('os').getcwd()", "sqrt(5, 6)", "1e100000", "10**999", "P"])
def test_receipt_expressions_are_restricted_algebra(bad):
    with pytest.raises(ValueError):
        independent.expression(bad)
