"""Tests for the exact entropic cubic global-extremum receipt."""

from __future__ import annotations

import derive_charged_entropic_branch_no_go as lane


def test_exact_cubic_global_extremum_is_narrowly_certified():
    artifact = lane.build()
    assert artifact["status"] == "ENTROPIC_CUBIC_AND_QUARTIC_EXACT_GLOBAL_CLASSIFICATION"
    assert artifact["checks_pass"] is True
    assert artifact["promotion_allowed"] is False
    assert artifact["cubic_truncation_theorem_established"] is True
    assert artifact["epistemic_gates"]["cubic_global_extremum_proved"] is True
    assert artifact["epistemic_gates"]["cubic_no_go_certified"] is True
    assert artifact["epistemic_gates"]["finite_seed_search_can_close_route"] is False
    assert artifact["quartic_global_minimizer_theorem_established"] is True
    assert artifact["epistemic_gates"]["quartic_global_optimality_proved"] is True
    assert artifact["epistemic_gates"][
        "strict_full_support_quartic_no_go_conditional"
    ] is True
    assert artifact["quartic_packet_globally_excluded"] is False
    assert artifact["epistemic_gates"]["full_entropic_mechanism_no_go_certified"] is False
    assert artifact["epistemic_gates"][
        "quadrupole_to_physical_log_mass_attachment_established"
    ] is False
    assert artifact["numerical_replay"]["role"] == "redundant_diagnostic_only_not_used_in_proof"
    assert artifact["scope"]["universal_impossibility_claimed"] is False
    assert "separate conditional attachment" in artifact["scope"]["physical_mass_attachment"]


def test_exact_stationary_enumeration_proves_multiplicity_one_maximum():
    certificate = lane.exact_cubic_certificate()
    rows = certificate["stationary_enumeration"]
    assert [row["positive_root_multiplicity"] for row in rows] == [1, 2, 3, 4, 5]
    assert [row["objective_sign"] for row in rows] == [
        "positive",
        "positive",
        "zero",
        "negative",
        "negative",
    ]
    assert rows[0]["objective_squared"] == {"numerator": 4, "denominator": 15}
    assert rows[1]["objective_squared"] == {"numerator": 1, "denominator": 24}
    assert certificate["global_maximum"]["multiplicity"] == 1
    assert certificate["quadrupole_degeneracy"]["double_eigenvalue_exact"] is True
    assert certificate["quadrupole_degeneracy"]["axis_frame_check"][
        "sum_over_six_axes_pipt_equals_2I"
    ] is True
    assert all(
        all(row["exact_checks"].values())
        for row in rows
    )
    assert all(certificate["checks"].values())
