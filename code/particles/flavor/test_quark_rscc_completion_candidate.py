"""Regression and fail-closed tests for the retrospective RSCC candidate."""

from __future__ import annotations

import ast
from pathlib import Path
import subprocess
import sys

import pytest


HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from audit_quark_rscc_completion_candidate import build_audit  # noqa: E402
from quark_rscc_completion_candidate import evaluate  # noqa: E402
from verify_quark_flavor_source_closure import build_artifact as build_source_closure  # noqa: E402
from verify_quark_rscc_module_arithmetic import build_artifact  # noqa: E402


def test_rscc_coordinates_reproduce_submitted_dimensionless_packet() -> None:
    payload = evaluate()
    coordinates = {
        key: float(value)
        for key, value in payload["dimensionless_output_coordinates_over_E_star"].items()
    }
    assert coordinates == pytest.approx(
        {
            "u": 1.7681821156784641e-22,
            "c": 1.0419911054963275e-19,
            "t": 1.4092333531091165e-17,
            "d": 3.8583169244260653e-22,
            "s": 7.622440077369575e-21,
            "b": 3.431198578910962e-19,
        },
        rel=1e-14,
    )


def test_default_candidate_has_no_gev_display_and_cannot_promote() -> None:
    payload = evaluate()
    assert payload["conditional_D10_gev_display"] is None
    assert payload["promotion_allowed"] is False
    assert payload["provenance"]["public_prediction_allowed"] is False


def test_runtime_removes_old_decimals_but_ancestry_remains_unclean() -> None:
    payload = evaluate()
    provenance = payload["provenance"]
    assert provenance["old_candidate_flavor_decimals_consumed"] is False
    assert provenance["pixel_source_uses_internal_stage5_quark_model"] is True
    assert provenance["pixel_branch_eligible_for_live_particle_prediction"] is False
    assert provenance["d10_mixed_sources_detected"] is True
    assert provenance["rscc_P_matches_d10_P"] is False
    assert provenance["rscc_alpha_U_matches_d10_alpha_U"] is False


def test_exact_verifier_scopes_arithmetic_away_from_physics() -> None:
    artifact = build_artifact()
    assert artifact["proof_status"] == (
        "exact_arithmetic_passed_physical_module_selection_not_proved"
    )
    assert artifact["exact_checks"]["dimension_arithmetic_passes"] is True
    assert artifact["exact_checks"]["F_contains_sign_sector_or_laplacian_6_mode"] is False
    assert artifact["exact_checks"]["dimension_End_S3_of_F"] == 5
    assert artifact["exact_checks"][
        "S3_invariance_selects_unique_scalar_response_on_F"
    ] is False
    assert artifact["exact_checks"][
        "negative_w2_terms_require_external_signed_response_law"
    ] is True
    assert all(value is False for value in artifact["non_theorems"].values())
    assert artifact["existing_count_only_register_boundary"][
        "rscc_requires_new_family_non_singlet_attachment"
    ] is True


def test_compare_only_audit_reproduces_arithmetic_and_keeps_F1_to_F6_open() -> None:
    audit = build_audit()
    comparison = audit["descriptive_mixed_chart_comparison"]
    assert comparison["max_abs_relative_error_percent"] == pytest.approx(
        0.29435865035596365
    )
    assert comparison["raw_diagonal_residual_sum"] == pytest.approx(
        1.161396614166569
    )
    assert comparison["statistical_interpretation_allowed"] is False
    assert audit["negative_control"]["max_abs_relative_error_percent"] == pytest.approx(
        0.21422928629970528
    )
    assert audit["negative_control"]["raw_diagonal_residual_sum"] == pytest.approx(
        0.6694007350278313
    )
    assert audit["negative_control"]["beats_full_rscc_maximum_error"] is True
    assert audit["negative_control"]["beats_full_rscc_raw_residual_sum"] is True
    assert all(
        audit["physical_nonclosure"][f"F{index}_{name}"] is False
        for index, name in (
            (1, "clean_source_root"),
            (2, "refinement_natural_family_carrier"),
            (3, "physical_channel_functor"),
            (4, "affine_mean_law"),
            (5, "common_scale_physical_readout"),
            (6, "RG_threshold_scheme_packet"),
        )
    )


def test_predictor_contains_no_reference_table_literal() -> None:
    source = (HERE / "quark_rscc_completion_candidate.py").read_text(
        encoding="utf-8"
    )
    tree = ast.parse(source)
    names = {node.id for node in ast.walk(tree) if isinstance(node, ast.Name)}
    assert "REFERENCE" not in names
    for literal in ("0.002160", "0.004700", "0.092900", "1.272900", "4.186000", "172.100000"):
        assert literal not in source


def test_cli_requires_explicit_retrospective_acknowledgement(
    tmp_path: Path,
) -> None:
    script = HERE / "quark_rscc_completion_candidate.py"
    output = tmp_path / "rscc.json"
    result = subprocess.run(
        [sys.executable, str(script), "--output", str(output)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode != 0
    assert "--allow-retrospective-rscc" in result.stderr
    assert not output.exists()


def test_cli_dimensionless_default_after_acknowledgement(tmp_path: Path) -> None:
    script = HERE / "quark_rscc_completion_candidate.py"
    output = tmp_path / "rscc.json"
    subprocess.run(
        [
            sys.executable,
            str(script),
            "--allow-retrospective-rscc",
            "--output",
            str(output),
        ],
        check=True,
    )
    import json

    payload = json.loads(output.read_text(encoding="utf-8"))
    assert payload["conditional_D10_gev_display"] is None


def test_rscc_discharges_no_flavor_source_closure_receipt() -> None:
    closure = build_source_closure()
    rscc = closure["rscc_formula_level_candidate"]
    assert rscc["promotion_allowed"] is False
    assert rscc["all_F1_to_F6_receipts_remain_open"] is True
    assert rscc["negative_control_beats_full_rscc"] is True
    assert closure["all_six_receipts_closed"] is False
    assert all(
        receipt["closed"] is False
        for receipt in closure["flavor_source_closure_receipts"].values()
    )
