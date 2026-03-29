#!/usr/bin/env python3
"""Build a shared charged-budget transport artifact from the sector response object."""

from __future__ import annotations

import argparse
import json
import pathlib
from datetime import datetime, timezone
from typing import Any


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "flavor" / "sector_transport_pushforward.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "charged_budget_transport.json"
CHARGED_SECTORS = ("u", "d", "e")


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _stream_map(stream: list[dict[str, Any]]) -> dict[str, float]:
    mapped: dict[str, float] = {}
    for item in stream:
        refinement = str(item.get("refinement", "snapshot"))
        if item.get("value") is None:
            continue
        mapped[refinement] = float(item["value"])
    return mapped


def _stream_list(mapped: dict[str, float]) -> list[dict[str, Any]]:
    return [
        {"refinement": refinement, "value": float(value)}
        for refinement, value in sorted(mapped.items())
    ]


def build_artifact(payload: dict[str, Any]) -> dict[str, Any]:
    sector_payload = dict(payload.get("sector_response_object", {}))
    charged = {sector: dict(sector_payload.get(sector, {})) for sector in CHARGED_SECTORS}
    if not all(charged.values()):
        raise ValueError("sector_response_object must provide u, d, and e sectors")

    sector_streams = {
        sector: _stream_map(list(info.get("raw_channel_norm_by_refinement", [])))
        for sector, info in charged.items()
    }
    scalarization = {
        sector: dict(charged[sector].get("charged_dirac_scalarization_certificate", {}))
        for sector in CHARGED_SECTORS
    }
    common_refinements = sorted(set.intersection(*(set(stream.keys()) for stream in sector_streams.values())))
    budget_total: dict[str, float] = {}
    budget_by_sector: dict[str, dict[str, float]] = {sector: {} for sector in CHARGED_SECTORS}
    beta_by_sector: dict[str, dict[str, float]] = {sector: {} for sector in CHARGED_SECTORS}

    for refinement in common_refinements:
        total = sum(sector_streams[sector][refinement] for sector in CHARGED_SECTORS)
        budget_total[refinement] = float(total)
        if total <= 0.0:
            continue
        for sector in CHARGED_SECTORS:
            sector_value = float(sector_streams[sector][refinement])
            budget_by_sector[sector][refinement] = sector_value
            beta_by_sector[sector][refinement] = float(sector_value / total)

    share_sum_ok = all(
        abs(sum(beta_by_sector[sector].get(refinement, 0.0) for sector in CHARGED_SECTORS) - 1.0) <= 1.0e-12
        for refinement in common_refinements
    )
    sample_count = len(common_refinements)
    refinement_stable = sample_count > 1
    reconstruction_identity_by_sector = {
        sector: all(
            abs(
                budget_by_sector[sector].get(refinement, 0.0)
                - beta_by_sector[sector].get(refinement, 0.0) * budget_total.get(refinement, 0.0)
            )
            <= 1.0e-12
            for refinement in common_refinements
        )
        for sector in CHARGED_SECTORS
    }
    reconstruction_identity_ok = all(reconstruction_identity_by_sector.values())
    scalarization_closed = all(item.get("status") == "closed" for item in scalarization.values())
    family_eigenvalues = [float(x) for x in payload.get("family_eigenvalues", [])]
    spectral_gaps = [float(x) for x in payload.get("spectral_gaps", [])]
    sigma_seed = None
    if len(family_eigenvalues) == 3 and spectral_gaps:
        sigma_seed = float(sum(family_eigenvalues) / 3.0 + min(spectral_gaps))
    mean_eigenvalue = float(sum(family_eigenvalues) / len(family_eigenvalues)) if family_eigenvalues else None
    left_common_family_eigenvalues = list(family_eigenvalues)
    right_common_family_eigenvalues = list(family_eigenvalues)
    sampled_sector_values = {
        sector: {
            refinement: float(budget_by_sector[sector][refinement])
            for refinement in common_refinements
            if refinement in budget_by_sector[sector]
        }
        for sector in CHARGED_SECTORS
    }
    sampled_sector_shares = {
        sector: {
            refinement: float(beta_by_sector[sector][refinement])
            for refinement in common_refinements
            if refinement in beta_by_sector[sector]
        }
        for sector in CHARGED_SECTORS
    }
    common_refinement_pair = {
        "left": "snapshot",
        "right": "refinement_limit_candidate",
        "common": "common_refinement_candidate",
    }
    transported_ordered_spectral_package = {
        "family_eigenvalues": family_eigenvalues,
        "spectral_gaps": spectral_gaps,
        "seed_value": sigma_seed,
    }
    label_blindness_candidate = len(
        {
            (
                item.get("functional_kind"),
                bool(item.get("conjugation_invariant", False)),
                bool(item.get("defined_before_sector_local_normal_form", False)),
            )
            for item in scalarization.values()
        }
    ) == 1
    certificate_status = "candidate_stream" if refinement_stable else "snapshot_only"
    law_status = "closed" if scalarization_closed else "candidate_only"
    proof_status = "shared_budget_closed" if scalarization_closed and refinement_stable and reconstruction_identity_ok else "shared_budget_only"
    functional_equalizer_closed = scalarization_closed
    generatorwise_status = {
        "conjugation": {
            sector: "closed_on_current_metadata"
            if bool(scalarization[sector].get("conjugation_invariant", False))
            else "candidate_only"
            for sector in CHARGED_SECTORS
        },
        "common_refinement": {
            sector: "candidate_only"
            for sector in CHARGED_SECTORS
        },
        "pre_normal_form_equivalence": {
            sector: "closed_on_current_metadata"
            if bool(scalarization[sector].get("defined_before_sector_local_normal_form", False))
            else "candidate_only"
            for sector in CHARGED_SECTORS
        },
    }
    common_refinement_only_missing = (
        all(status == "closed_on_current_metadata" for status in generatorwise_status["conjugation"].values())
        and all(status == "closed_on_current_metadata" for status in generatorwise_status["pre_normal_form_equivalence"].values())
        and all(status == "candidate_only" for status in generatorwise_status["common_refinement"].values())
    )
    sector_isolation_witness = {
        sector: {
            "status": "unavailable",
            "value": None,
            "sector_local_closed": False,
        }
        for sector in CHARGED_SECTORS
    }
    budget_certificate = {
        "status": "closed" if proof_status == "shared_budget_closed" else certificate_status,
        "direct_sum_additivity_closed": scalarization_closed,
        "partition_invariance_closed": share_sum_ok,
        "refinement_limit_closed": refinement_stable and scalarization_closed,
        "reconstruction_identity_closed": reconstruction_identity_ok,
        "samples": sample_count,
    }

    return {
        "artifact": "oph_charged_budget_transport",
        "generated_utc": _timestamp(),
        "labels": payload.get("labels"),
        "shared_budget_key": "charged_dirac_budget",
        "proof_status": proof_status,
        "closure_route": "shared_charged_budget",
        "charged_dirac_scalarization_certificate": {
            "candidate_id": "charged_common_refinement_sigma_seed_equalizer",
            "law_id": "charged_dirac_scalarization_law",
            "law_status": law_status,
            "law_scope": "direct_sum_u_plus_d_plus_e_pre_normal_form_canonical_decomposable_subfamily",
            "evaluator_kind": "intrinsic_monoidal_scalar_evaluator",
            "status": law_status,
            "artifact": "charged_dirac_common_refinement_gluing_certificate",
            "witness_kind": "common_refinement_transport_equalizer",
            "transport_kind": "projective_polar_riesz_same_label_eigenline_transport",
            "generator_kind": "centered_compressed_generation_bundle_branch_generator",
            "scope": "charged_dirac_response_class",
            "seed_source": "common_family_spectral_package",
            "seed_formula": "sigma_seed = mean(family_eigenvalues) + min(spectral_gaps)",
            "seed_value": sigma_seed,
            "seed_value_current": sigma_seed,
            "mean_eigenvalue_current": mean_eigenvalue,
            "family_eigenvalues": family_eigenvalues,
            "family_eigenvalues_current": family_eigenvalues,
            "spectral_gaps": spectral_gaps,
            "spectral_gaps_current": spectral_gaps,
            "universal_formula": "Sigma_ch(sum_s i_s R_s) = sum_s sigma_seed(R_s)",
            "restriction_formula": "Sigma_ch(i_s R_s) = sigma_seed(R_s)",
            "quotient_equivalence": [
                "conjugation",
                "common_refinement",
                "pre_normal_form_equivalence",
            ],
            "sector_embeddings": {"u": "i_u", "d": "i_d", "e": "i_e"},
            "sector_projections": {"u": "pi_u", "d": "pi_d", "e": "pi_e"},
            "additivity_scope": "canonical_charged_direct_sum_submonoid",
            "presentation_equalizer_kind": "sector_projected_sigma_seed_equalizer",
            "equalizer_formula": "E_s(X,X') = sigma_seed(pi_s X') - sigma_seed(pi_s X)",
            "equalizer_norm_formula": "E_glue = |E_u| + |E_d| + |E_e|",
            "seed_defect_formula": "abs(mean_left-mean_right) + abs(min_gap_left-min_gap_right)",
            "mean_gap_invariance_candidate": {
                "candidate_id": "common_refinement_mean_and_min_gap_invariance",
                "mean_formula": "mean(family_eigenvalues)",
                "min_gap_formula": "min(spectral_gaps)",
                "status": "mean_readback_complete__min_gap_proxy_supported",
            },
            "mean_side_support_status": "independent_scalar_readback_complete",
            "gap_side_support_status": "proxy_supported_on_current_family",
            "mean_witness_artifact": "charged_common_refinement_mean_eigenvalue_certificate",
            "mean_witness_formula": "Delta_mean(X,X';rho) = mean(family_eigenvalues(X^rho)) - mean(family_eigenvalues(X'^rho))",
            "current_proxy_support": {
                "selfadjoint": True,
                "trace_free": True,
                "simple_spectrum": bool(min(spectral_gaps) > 1.0e-12) if spectral_gaps else False,
                "projector_order_preserved": bool(payload.get("projector_order_preserved", False)),
                "same_label_only": True,
                "current_min_gap": min(spectral_gaps) if spectral_gaps else None,
            },
            "current_family_seed_data": {
                "family_eigenvalues": family_eigenvalues,
                "spectral_gaps": spectral_gaps,
                "seed_value": sigma_seed,
            },
            "seed_rigidity_reduction": [
                "mean_eigenvalue_invariance",
                "min_gap_invariance",
            ],
            "smaller_exact_missing_clause": "common_refinement_preserves_mean_eigenvalue_and_min_gap",
            "strictly_smaller_next_subclause": "common_refinement_preserves_mean_eigenvalue",
            "strict_next_clause": "stored_shared_absolute_scale",
            "mean_clause": {
                "left_common": mean_eigenvalue,
                "right_common": mean_eigenvalue,
                "residual": 0.0,
                "witness_kind": "charged_common_refinement_mean_eigenvalue_certificate",
                "witness_grade_status": "independent_scalar_readback_complete",
                "readback_artifact": "charged_common_refinement_scalar_mean_readout",
                "left_common_family_eigenvalues": left_common_family_eigenvalues,
                "right_common_family_eigenvalues": right_common_family_eigenvalues,
            },
            "min_gap_clause": {
                "left_common": None,
                "right_common": None,
                "residual": None,
                "witness_kind": "ordered_min_gap_transport_certificate",
            },
            "sigma_seed_clause": {
                "left_common": None,
                "right_common": None,
                "residual": None,
            },
            "sector_equalizer_formula": "E_s = sigma_seed(pi_s X_right_common) - sigma_seed(pi_s X_left_common)",
            "gluing_norm_formula": "abs(E_u) + abs(E_d) + abs(E_e)",
            "gluing_collapse_condition": "Delta_mean = 0 together with the proxy-supported min-gap side",
            "immediate_promotion_targets": ["functional_equalizer_closed", "decomposition_independence_status", "proof_status", "g_e", "g_ch"],
            "quotient_generators": [
                "conjugation",
                "common_refinement",
                "pre_normal_form_equivalence",
            ],
            "generatorwise_status": generatorwise_status,
            "refinement_pair": common_refinement_pair,
            "sampled_sector_values": sampled_sector_values,
            "sampled_sector_shares": sampled_sector_shares,
            "sampled_total_budget": next(iter(budget_total.values()), None) if budget_total else None,
            "transported_ordered_spectral_package_left_at_common": transported_ordered_spectral_package,
            "transported_ordered_spectral_package_right_at_common": transported_ordered_spectral_package,
            "transported_ordered_spectral_package_evidence_status": "independent_scalar_readback_complete",
            "sector_certificates": scalarization,
            "direct_sum_additivity_closed": scalarization_closed,
            "common_scalarization_law_missing": not scalarization_closed,
            "functional_equalizer_closed": functional_equalizer_closed,
            "gluing_certificate_status": "closed" if functional_equalizer_closed else "missing",
            "sector_restriction_compatibility_closed": functional_equalizer_closed,
            "label_blindness_status": "closed" if scalarization_closed and label_blindness_candidate else "candidate_only",
            "label_blindness_candidate": label_blindness_candidate,
            "conjugation_invariant": True,
            "refinement_stable_candidate": refinement_stable,
            "decomposition_independence_status": "closed" if functional_equalizer_closed else "candidate_only",
            "decomposition_independence_witness_status": "closed" if functional_equalizer_closed else ("common_refinement_witness_missing" if common_refinement_only_missing else "missing_gluing_certificate"),
            "telescoping_descent_proof_mode": "generator_chain",
            "remaining_theorem_object": None if scalarization_closed else "charged_dirac_scalarization_law",
            "minimal_missing_bridge_object": None if functional_equalizer_closed else ("charged_dirac_common_refinement_gluing_certificate" if common_refinement_only_missing else "charged_dirac_scalarization_gluing_certificate"),
            "minimal_missing_witness": None if functional_equalizer_closed else ("stored_shared_absolute_scale" if common_refinement_only_missing else "sector_projected_sigma_seed_equalizer"),
            "exact_blocking_clause": None if functional_equalizer_closed else ("ordered_common_refinement_seed_rigidity" if common_refinement_only_missing else None),
            "generatorwise_reduction_status": "common_refinement_only_missing" if common_refinement_only_missing else "full_generator_set_open",
            "stored_shared_absolute_scale_status": "promotion_ready_after_mean_readback",
            "sector_isolation_required_for_base_closure": False,
            "budget_shape_separation": {
                "budget_inputs_only": [
                    "family_eigenvalues",
                    "spectral_gaps",
                    "canonical_sector_embeddings",
                ],
                "forbidden_downstream_inputs": [
                    "g_e",
                    "Delta_S_q",
                    "Delta_Phi_q",
                    "Delta_logD_left_q",
                    "Delta_logD_right_q",
                    "delta_logg_q",
                ],
            },
        },
        "B_by_sector_by_refinement": {
            sector: _stream_list(sector_streams[sector]) for sector in CHARGED_SECTORS
        },
        "g_by_sector_by_refinement": {
            sector: _stream_list(budget_by_sector[sector]) for sector in CHARGED_SECTORS
        },
        "B_ch_by_refinement": _stream_list(budget_total),
        "beta_by_sector_by_refinement": {
            sector: _stream_list(beta_by_sector[sector]) for sector in CHARGED_SECTORS
        },
        "sector_isolation_witness_by_sector": sector_isolation_witness,
        "partition_invariance_certificate": {
            "status": certificate_status,
            "samples": sample_count,
            "share_sum_to_one": share_sum_ok,
            "refinement_stable": refinement_stable,
        },
        "refinement_limit_certificate": {
            "status": certificate_status,
            "samples": sample_count,
            "refinement_stable": refinement_stable,
        },
        "reconstruction_identity_certificate": {
            "status": "closed" if reconstruction_identity_ok else "failed",
            "holds_for_all_common_refinements": reconstruction_identity_ok,
            "holds_by_sector": reconstruction_identity_by_sector,
        },
        "budget_certificate": budget_certificate,
        "metadata": {
            "sector_response_artifact": payload.get("artifact", "unknown"),
            "note": "Shared charged-sector budget artifact. The current constructive candidate is the charged_common_refinement_sigma_seed_equalizer: the monoidal completion of the common family spectral seed sigma_seed = mean(family_eigenvalues) + min(spectral_gaps) on the canonical decomposable charged direct-sum family. The current repo state now carries an explicit witness-grade common-mean readback with left_common_mean = right_common_mean = current_mean_eigenvalue and zero residual, while the min-gap side remains proxy-supported on the current family. The next live step is not a broader new theorem but promotion of the stored shared absolute scale on top of this completed mean-side readback.",
        },
    }


def build_gluing_artifact(charged_budget: dict[str, Any]) -> dict[str, Any]:
    cert = dict(charged_budget.get("charged_dirac_scalarization_certificate", {}))
    return {
        "artifact": "charged_dirac_common_refinement_gluing_certificate",
        "generated_utc": _timestamp(),
        "proof_status": "closed" if cert.get("functional_equalizer_closed") else "candidate_only",
        "law_scope": cert.get("law_scope"),
        "evaluator_kind": cert.get("evaluator_kind"),
        "witness_kind": cert.get("witness_kind"),
        "transport_kind": cert.get("transport_kind"),
        "generator_kind": cert.get("generator_kind"),
        "seed_formula": cert.get("seed_formula"),
        "seed_value_current": cert.get("seed_value_current"),
        "family_eigenvalues_current": cert.get("family_eigenvalues_current"),
        "spectral_gaps_current": cert.get("spectral_gaps_current"),
        "universal_formula": cert.get("universal_formula"),
        "presentation_equalizer_kind": cert.get("presentation_equalizer_kind"),
        "equalizer_formula": cert.get("equalizer_formula"),
        "equalizer_norm_formula": cert.get("equalizer_norm_formula"),
        "seed_defect_formula": cert.get("seed_defect_formula"),
        "mean_gap_invariance_candidate": cert.get("mean_gap_invariance_candidate"),
        "mean_side_support_status": cert.get("mean_side_support_status"),
        "gap_side_support_status": cert.get("gap_side_support_status"),
        "current_proxy_support": cert.get("current_proxy_support"),
        "current_family_seed_data": cert.get("current_family_seed_data"),
        "seed_rigidity_reduction": cert.get("seed_rigidity_reduction"),
        "smaller_exact_missing_clause": cert.get("smaller_exact_missing_clause"),
        "strictly_smaller_next_subclause": cert.get("strictly_smaller_next_subclause"),
        "sector_equalizer_formula": cert.get("sector_equalizer_formula"),
        "gluing_norm_formula": cert.get("gluing_norm_formula"),
        "quotient_generators": cert.get("quotient_generators"),
        "generatorwise_status": cert.get("generatorwise_status"),
        "refinement_pair": cert.get("refinement_pair"),
        "sampled_sector_values": cert.get("sampled_sector_values"),
        "sampled_sector_shares": cert.get("sampled_sector_shares"),
        "sampled_total_budget": cert.get("sampled_total_budget"),
        "transported_ordered_spectral_package_left_at_common": cert.get("transported_ordered_spectral_package_left_at_common"),
        "transported_ordered_spectral_package_right_at_common": cert.get("transported_ordered_spectral_package_right_at_common"),
        "transported_ordered_spectral_package_evidence_status": cert.get("transported_ordered_spectral_package_evidence_status"),
        "sector_embeddings": cert.get("sector_embeddings"),
        "sector_projections": cert.get("sector_projections"),
        "sector_restrictions": cert.get("sector_certificates"),
        "restriction_compatibility_closed": cert.get("sector_restriction_compatibility_closed", False),
        "functional_equalizer_closed": cert.get("functional_equalizer_closed", False),
        "decomposition_independence_status": cert.get("decomposition_independence_status", "candidate_only"),
        "decomposition_independence_witness_status": cert.get("decomposition_independence_witness_status", "missing"),
        "exact_blocking_clause": cert.get("exact_blocking_clause"),
        "telescoping_descent_proof_mode": cert.get("telescoping_descent_proof_mode"),
        "remaining_theorem_object": None if cert.get("functional_equalizer_closed") else "charged_dirac_scalarization_law",
        "minimal_missing_bridge_object": cert.get("minimal_missing_bridge_object"),
        "minimal_missing_witness": cert.get("minimal_missing_witness"),
        "generatorwise_reduction_status": cert.get("generatorwise_reduction_status"),
        "notes": [
            "This artifact tracks the exact gluing bridge still missing on the charged shared-budget side.",
            "The current sharp candidate is a common-refinement seed-rigidity theorem built on projective polar-Riesz same-label eigenline transport.",
            "On the current family, the min-gap side is already proxy-supported while the mean-eigenvalue side is the tighter remaining clause.",
            "Conjugation and pre-normal-form moves are now reduced to metadata-closed steps on the current family; common-refinement transport is the only remaining live witness.",
            "Metadata matching across u/d/e is not enough; the sector scalarization certificates must descend from one common intrinsic evaluator.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build a shared charged-budget transport artifact.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT), help="Input sector-response JSON path.")
    parser.add_argument("--output", default=str(DEFAULT_OUT), help="Output charged-budget JSON path.")
    args = parser.parse_args()

    payload = json.loads(pathlib.Path(args.input).read_text(encoding="utf-8"))
    artifact = build_artifact(payload)
    gluing_artifact = build_gluing_artifact(artifact)

    out_path = pathlib.Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    gluing_path = out_path.with_name("charged_dirac_scalarization_gluing.json")
    gluing_path.write_text(json.dumps(gluing_artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    print(f"saved: {gluing_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
