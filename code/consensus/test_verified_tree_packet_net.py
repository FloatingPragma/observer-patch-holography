#!/usr/bin/env python3
"""Tests for the verified tree packet-net repair domain."""

from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import export_verified_tree_packet_net as domain  # noqa: E402


def test_tree_packet_domain_checks_close() -> None:
    payload = domain.build_payload()
    checks = payload["theorem_checks"]
    assert checks["total_states_checked"] == 3**4 * 2**4
    assert checks["repair_completeness"] is True
    assert checks["strict_lyapunov_descent"] is True
    assert checks["unique_terminal_normal_form"] is True
    assert checks["max_terminal_count_seen"] == 1


def test_petz_recovery_is_computed_from_explicit_classical_channel() -> None:
    witness = domain.compute_petz_recovery_witness()

    assert witness["witness_kind"] == "independently_computed_classical_recovery_channel"
    assert witness["audit_issue"] == 518
    assert witness["recovery_matrix_shape"] == [9, 3]
    assert witness["recovery_matrix"][0] == ["1/2", "0", "0"]
    assert witness["recovery_matrix"][4] == ["0", "1/2", "0"]
    assert witness["recovery_matrix"][8] == ["0", "0", "1/2"]
    assert len(witness["classical_choi_diagonal"]) == 27
    assert witness["classical_choi_diagonal_nonnegative"] is True
    assert witness["column_sums_exact"] == {"0": "1", "1": "1", "2": "1"}
    assert witness["column_residuals_exact"] == {"0": "0", "1": "0", "2": "0"}
    assert witness["trace_preserving"] is True
    assert witness["basis_marginal_recovery_matrix"] == [
        ["1", "0", "0"],
        ["0", "1", "0"],
        ["0", "0", "1"],
    ]
    assert witness["basis_marginal_recovery_identity"] is True
    assert witness["induced_l1_norm"] == "1"
    assert witness["cptp"] is True
    assert witness["trace_norm_contractive"] is True
    assert witness["reference_state_normalized"] is True
    assert witness["support_gap_exact"] == "1/3"
    assert witness["full_support"] is True
    assert witness["failures"] == []
    assert witness["pass"] is True
    assert "remaining receipt-class" in witness["claim_boundary"]["not_closed_here"][2]


def _conditional_fixture() -> dict[int, dict[int, Fraction]]:
    return {
        b: dict(row)
        for b, row in domain.DEFAULT_PETZ_CONDITIONAL.items()
    }


def test_petz_negative_control_rejects_broken_column_normalization() -> None:
    conditional = _conditional_fixture()
    conditional[0][0] = Fraction(3, 5)

    witness = domain.compute_petz_recovery_witness(conditional=conditional)

    assert witness["classical_choi_diagonal_nonnegative"] is True
    assert witness["column_sums_exact"]["0"] == "11/10"
    assert witness["trace_preserving"] is False
    assert witness["basis_marginal_recovery_identity"] is False
    assert witness["induced_l1_norm"] == "11/10"
    assert witness["cptp"] is False
    assert witness["trace_norm_contractive"] is False
    assert witness["failures"] == [
        "recovery_column_normalization_failed",
        "basis_marginal_recovery_failed",
        "induced_l1_norm_exceeds_one",
    ]
    assert witness["pass"] is False


def test_petz_negative_control_rejects_negative_choi_entry() -> None:
    conditional = _conditional_fixture()
    conditional[0] = {
        0: Fraction(1),
        1: Fraction(-1, 4),
        2: Fraction(1, 4),
    }

    witness = domain.compute_petz_recovery_witness(conditional=conditional)

    assert witness["column_sums_exact"]["0"] == "1"
    assert witness["trace_preserving"] is True
    assert witness["basis_marginal_recovery_identity"] is True
    assert witness["classical_choi_diagonal_nonnegative"] is False
    assert witness["induced_l1_norm"] == "3/2"
    assert witness["cptp"] is False
    assert witness["trace_norm_contractive"] is False
    assert witness["failures"] == [
        "classical_choi_diagonal_has_negative_entry",
        "induced_l1_norm_exceeds_one",
    ]
    assert witness["pass"] is False


def test_petz_negative_control_rejects_zero_reference_support() -> None:
    reference = {
        0: Fraction(0),
        1: Fraction(1, 2),
        2: Fraction(1, 2),
    }

    witness = domain.compute_petz_recovery_witness(reference_sigma_b=reference)

    assert witness["cptp"] is True
    assert witness["trace_norm_contractive"] is True
    assert witness["basis_marginal_recovery_identity"] is True
    assert witness["reference_state_normalized"] is True
    assert witness["support_gap_exact"] == "0"
    assert witness["full_support"] is False
    assert witness["inverse_sqrt_bound"] is None
    assert witness["failures"] == ["reference_state_not_full_support"]
    assert witness["pass"] is False


def test_tracked_petz_receipt_matches_independent_recomputation() -> None:
    tracked = json.loads(
        (HERE / "runs" / "verified_tree_packet_net_domain.json").read_text(
            encoding="utf-8"
        )
    )

    assert tracked["petz_domain"] == domain.compute_petz_recovery_witness()


def test_exported_payload_roundtrips(tmp_path: Path) -> None:
    out = tmp_path / "verified_tree_packet_net_domain.json"
    assert domain.main.__module__ == "export_verified_tree_packet_net"
    payload = domain.build_payload()
    out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    loaded = json.loads(out.read_text(encoding="utf-8"))
    assert loaded["artifact"] == "oph_verified_tree_packet_net_domain"
    assert loaded["issue"] == 238
    assert loaded["petz_domain"]["cptp"] is True
    assert loaded["petz_domain"]["support_gap_gamma_sigma"] > 0.0
    assert loaded["petz_domain"]["trace_norm_contractive"] is True
    assert loaded["petz_domain"]["basis_marginal_recovery_identity"] is True
    assert loaded["petz_domain"]["pass"] is True
    assert loaded["quotient_compatibility"]["descends_to_quotient"] is True
    assert "quotient normal form" in loaded["quotient_compatibility"]["physical_law_use"]
