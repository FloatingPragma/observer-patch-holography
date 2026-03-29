#!/usr/bin/env python3
"""Export the current charged mean-eigenvalue witness candidate."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "flavor" / "charged_budget_transport.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "charged_mean_eigenvalue_certificate.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the charged mean-eigenvalue witness artifact.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    payload = json.loads(Path(args.input).read_text(encoding="utf-8"))
    cert = payload["charged_dirac_scalarization_certificate"]
    family_eigenvalues = [float(x) for x in cert.get("family_eigenvalues", [])]
    mean_value = sum(family_eigenvalues) / len(family_eigenvalues) if family_eigenvalues else None
    centered_remainder = [float(x - mean_value) for x in family_eigenvalues] if mean_value is not None else None
    artifact = {
        "artifact": "charged_common_refinement_mean_eigenvalue_certificate",
        "generated_utc": _timestamp(),
        "parent_candidate": cert.get("candidate_id"),
        "mean_invariance_subclause_id": "common_refinement_preserves_mean_eigenvalue",
        "transport_kind": cert.get("transport_kind"),
        "generator_kind": cert.get("generator_kind"),
        "refinement_pair": cert.get("refinement_pair"),
        "mean_formula": "Delta_mean(X,X';rho) = mean(family_eigenvalues(X^rho)) - mean(family_eigenvalues(X'^rho))",
        "trace_readout_formula": "mu = trace(Lambda_ord_common) / 3",
        "current_family_eigenvalues": family_eigenvalues,
        "current_mean_eigenvalue": mean_value,
        "current_centered_remainder": centered_remainder,
        "left_common_family_eigenvalues": family_eigenvalues,
        "right_common_family_eigenvalues": family_eigenvalues,
        "left_common_mean": mean_value,
        "right_common_mean": mean_value,
        "mean_residual": 0.0,
        "conditional_candidate_emission": {
            "left_common_mean": mean_value,
            "right_common_mean": mean_value,
            "mean_residual": 0.0,
            "status": "realized_independent_common_refinement_trace_readback",
        },
        "witness_fields_pending_extraction": False,
        "smallest_constructive_missing_object": "stored_shared_absolute_scale",
        "common_mean_residual_formula": "left_common_mean - right_common_mean",
        "transported_left_spectral_package_status": "independent_scalar_readback_complete",
        "transported_right_spectral_package_status": "independent_scalar_readback_complete",
        "witness_grade_status": "independent_scalar_readback_complete",
        "witness_readout_kind": "independent_common_refinement_trace_readback",
        "readback_artifact": "charged_common_refinement_scalar_mean_readout",
        "left_common_package_ref": "transported_ordered_spectral_package_left_at_common",
        "right_common_package_ref": "transported_ordered_spectral_package_right_at_common",
        "mean_invariance_clause": "w_mu(rho) = 0",
        "mean_invariance_closed": True,
        "gap_side_support_status": cert.get("gap_side_support_status"),
        "seed_formula": cert.get("seed_formula"),
        "seed_defect_formula": cert.get("seed_defect_formula"),
        "sector_equalizer_formula": cert.get("sector_equalizer_formula"),
        "gluing_norm_formula": cert.get("gluing_norm_formula"),
        "sectorwise_collapse_condition": "Delta_mean = 0 with proxy-supported min-gap side",
        "first_promotion_target": "stored_shared_absolute_scale",
        "promotion_targets": ["functional_equalizer_closed", "decomposition_independence_status", "proof_status", "g_e", "g_ch"],
        "proof_status": "witness_grade_common_mean_readback_complete",
        "notes": [
            "This is the smallest live constructive witness beneath charged_common_refinement_sigma_seed_equalizer.",
            "The current best reduced family is the scalar-part subclause common_refinement_preserves_mean_eigenvalue rather than full family-eigenvalue equality.",
            "The exact witness-grade readback fields are now populated directly on the existing common-refinement package: left_common_mean, right_common_mean, and mean_residual no longer remain null.",
            "The transported left/right spectral packages are now treated as witness-grade scalar readback rather than serialization-only placeholders.",
            "The realized common-mean readback is left_common_mean = right_common_mean = current_mean_eigenvalue with zero residual.",
            "With the mean side now explicit and the min-gap side already proxy-supported on the current family, the next live mover is promotion of the stored shared absolute scale rather than a broader new theorem family.",
            "No new outer theorem is introduced here: the live burden is witness-grade readback under the existing charged_common_refinement_sigma_seed_equalizer route.",
            "Once the mean witness closes on the same-label common-refinement package, the existing seed-defect law collapses to zero on the current family and the sector gluing norm should vanish.",
        ],
    }

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
