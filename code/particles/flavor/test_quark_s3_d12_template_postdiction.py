"""Regression and fail-closed policy tests for the S3/D12 quark diagnostic."""

from __future__ import annotations

import math
from pathlib import Path
import subprocess
import sys

import pytest


HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from audit_quark_s3_d12_template_postdiction import build_audit  # noqa: E402
from derive_quark_kernel_normalization_acceptance_harness import evaluate as evaluate_legacy_acceptance  # noqa: E402
from derive_quark_kernel_three_scalar_interface_theorem import build as build_ray_interface  # noqa: E402
from quark_s3_d12_template_postdiction import evaluate, load_repository_inputs  # noqa: E402
from verify_quark_flavor_source_closure import build_artifact as build_source_closure  # noqa: E402
from verify_s3_transposition_heat_shape import build_artifact as build_s3_artifact  # noqa: E402


def test_frozen_template_formula_reproduces_submitted_coordinates() -> None:
    payload = evaluate()
    coordinates = {
        key: float(value) for key, value in payload["dimensionless_output_coordinates"].items()
    }
    assert coordinates == pytest.approx(
        {
            "u": 0.0021602387376279098,
            "d": 0.004713845068652681,
            "s": 0.09312611565085318,
            "c": 1.2730350538241826,
            "b": 4.192020570504227,
            "t": 172.17119584025592,
        },
        rel=1e-14,
        abs=1e-16,
    )


def test_repository_ancestry_is_explicit_and_non_promotable() -> None:
    inputs, provenance = load_repository_inputs()
    assert float(inputs["g_ch_template"]) == pytest.approx(0.9231656602589082)
    assert provenance["family_transport_kernel_status"] == "template"
    assert provenance["numerical_flavor_template_consumed"] is True
    assert provenance["p_source_uses_stage5_quark_model"] is True
    assert provenance["dimensionful_scale_emitted"] is False
    assert provenance["public_prediction_allowed"] is False

    payload = evaluate()
    assert payload["promotion_allowed"] is False
    assert payload["unit_status"]["emitted_coordinates"] == "dimensionless"
    assert payload["unit_status"]["physical_mass_or_yukawa_claim_allowed"] is False


def test_compare_only_audit_reproduces_arithmetic_without_promoting_it() -> None:
    audit = build_audit()
    comparison = audit["descriptive_mixed_chart_comparison"]
    assert comparison["max_abs_relative_error_percent"] == pytest.approx(
        0.29457592878043837, rel=1e-12
    )
    assert comparison["raw_diagonal_residual_sum"] == pytest.approx(
        1.1653240915049539, rel=1e-12
    )
    assert comparison["statistical_interpretation_allowed"] is False
    assert audit["ancestry_audit"]["source_only_ancestry_passes"] is False
    assert audit["runtime_separation_audit"]["runtime_comparison_separation_passes"] is True
    assert audit["local_grammar_audit_correction"]["models_tied_for_best_maximum_error"] == 8


def test_lq_basis_is_centered_and_nonsingular_at_candidate_x() -> None:
    payload = evaluate()
    linear = [float(value) for value in payload["basis"]["L=ctr(-1,x,1)"]]
    quadratic = [float(value) for value in payload["basis"]["Q=ctr(1,x^2,1)"]]
    x = float(payload["derived_scalars"]["x"])
    ll = sum(value * value for value in linear)
    lq = sum(a * b for a, b in zip(linear, quadratic))
    qq = sum(value * value for value in quadratic)
    assert sum(linear) == pytest.approx(0.0, abs=1e-15)
    assert sum(quadratic) == pytest.approx(0.0, abs=1e-15)
    assert ll * qq - lq * lq == pytest.approx(4 * (1 - x * x) ** 2 / 3, rel=1e-14)
    assert not math.isclose(abs(x), 1.0)


def test_s3_heat_identity_is_exact_but_physically_scoped() -> None:
    artifact = build_s3_artifact()
    assert artifact["proof_status"] == "exact_symbolic_verification_passed"
    assert artifact["laplacian_eigenvalues_with_multiplicity"] == {"0": 1, "3": 4, "6": 1}
    assert "does not derive" in str(artifact["scope_boundary"])


def test_legacy_three_scalar_theorem_is_scoped_to_ray_subfamily() -> None:
    artifact = build_ray_interface()
    assert artifact["proof_status"] == "closed_ray_subfamily_interface_theorem"
    assert artifact["guards"]["general_quark_interface_claim_allowed"] is False
    assert artifact["general_interface_boundary"]["common_scale_eigenvalue_coordinates"] == 6


def test_legacy_acceptance_harness_rejects_template_ancestry_at_g1() -> None:
    result = evaluate_legacy_acceptance(
        {
            "kernel": {"refinements": [{}]},
            "ancestry": {
                "attestations": {"numerical_flavor_template_consumed": True},
                "artifacts": [],
            },
        }
    )
    assert result["acceptance_passed"] is False
    assert result["first_failed_gate"] == "G1_ancestry"
    assert result["gates"]["G1_ancestry"]["template_ancestry_declared"] is True


def test_flavor_source_closure_is_exact_algebra_plus_open_physics_contract() -> None:
    artifact = build_source_closure()
    exact = artifact["exact_algebra"]
    assert exact["proof_status"] == "exact_symbolic_verification_passed"
    assert exact["normalized_trace_theorem"]["rank_one_or_slot_coefficients"] == {
        "isotropic_heat_slot": "1/25",
        "up_odd_rank_one": "1/2",
        "down_odd_rank_one": "1/5",
        "up_even_rank_one": "1/10",
        "down_even_rank_one": "1/4",
    }
    assert exact["normalized_trace_theorem"]["conditional_denominator_tuple"] == [
        5,
        2,
        5,
        10,
        4,
    ]
    assert artifact["all_six_receipts_closed"] is False
    assert artifact["conditional_mass_evaluator_unique_if_all_receipts_close"] is True
    assert artifact["current_repository_emits_physical_quark_sextet"] is False
    assert artifact["promotion_allowed"] is False
    assert artifact["submitted_proof_bundle_receipt"]["expected_sha256"] == (
        "86e476d4b9c7910aefa0d259395e63e40b0b71c1534cf192b4f578cb21320759"
    )
    assert all(
        receipt["closed"] is False
        for receipt in artifact["flavor_source_closure_receipts"].values()
    )


def test_cli_requires_explicit_template_acknowledgement(tmp_path: Path) -> None:
    script = HERE / "quark_s3_d12_template_postdiction.py"
    result = subprocess.run(
        [sys.executable, str(script), "--output", str(tmp_path / "candidate.json")],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode != 0
    assert "--allow-template-ancestry" in result.stderr
    assert not (tmp_path / "candidate.json").exists()
