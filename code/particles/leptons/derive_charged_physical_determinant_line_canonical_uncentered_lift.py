#!/usr/bin/env python3
"""Emit the canonical post-promotion charged lift on theorem-grade physical Y_e."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

from charged_absolute_route_common import (
    ANCHOR_SECTION_JSON,
    DETERMINANT_LINE_JSON,
    PHYSICAL_CANONICAL_LIFT_JSON,
    PHYSICAL_EQUALIZER_JSON,
    artifact_ref,
    load_json,
)


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(determinant_line: dict, anchor: dict, physical_equalizer: dict) -> dict:
    return {
        "artifact": "oph_charged_physical_determinant_line_canonical_uncentered_lift",
        "generated_utc": _timestamp(),
        "proof_status": "closed_conditional_theorem_on_theorem_grade_physical_Y_e",
        "public_promotion_allowed": False,
        "theorem_scope": "theorem_grade_physical_charged_class_only",
        "theorem_id": "charged_physical_determinant_line_canonical_uncentered_lift",
        "theorem_statement": (
            "Assume theorem-grade physical charged data Y_e on the declared physical surface and "
            "assume theorem-grade centered charged operator C_hat_e on that same surface. Define "
            "s_det(Y_e) := log(det(Y_e)), mu_phys(Y_e) := (1/3) * s_det(Y_e), and "
            "C_tilde_e(Y_e) := C_hat_e(Y_e) + mu_phys(Y_e) I. Then C_tilde_e(Y_e) is an admissible "
            "uncentered lift of C_hat_e(Y_e), the determinant-line section and affine anchor are "
            "canonically fixed by s_det(Y_e) = 3 * mu_phys(Y_e) and A_ch(Y_e) = mu_phys(Y_e), and "
            "for any refinement representatives r,r' of the same physical Y_e the induced identity-mode "
            "cocycle vanishes because mu(r) = mu(r') = mu_phys(Y_e). Therefore "
            "g_e(Y_e) = exp(mu_phys(Y_e)), Delta_e_abs(Y_e) = log(g_ch_shared) - mu_phys(Y_e), "
            "and the charged mass triple are algebraic descendants of theorem-grade physical Y_e."
        ),
        "dependencies": {
            "determinant_line_section": {
                "artifact": determinant_line.get("artifact"),
                "artifact_ref": artifact_ref(DETERMINANT_LINE_JSON),
                "theorem_id": determinant_line.get("reduction_theorem", {}).get("id"),
            },
            "charged_absolute_anchor": {
                "artifact": anchor.get("artifact"),
                "artifact_ref": artifact_ref(ANCHOR_SECTION_JSON),
                "exact_missing_object": anchor.get("exact_missing_object"),
            },
            "physical_identity_mode_equalizer": {
                "artifact": physical_equalizer.get("artifact"),
                "artifact_ref": artifact_ref(PHYSICAL_EQUALIZER_JSON),
                "theorem_id": physical_equalizer.get("reduction_theorem", {}).get("id"),
            },
        },
        "canonical_definitions": {
            "determinant_line_section": "s_det(Y_e) = log(det(Y_e))",
            "descended_scalar": "mu_phys(Y_e) = (1/3) * log(det(Y_e))",
            "uncentered_lift": "C_tilde_e(Y_e) = C_hat_e(Y_e) + mu_phys(Y_e) I",
            "affine_anchor": "A_ch(Y_e) = mu_phys(Y_e)",
        },
        "forced_equalizer_on_declared_surface": {
            "fiber_rule": "delta(r,r') = 0 for refinement representatives r,r' of the same physical Y_e",
            "reason": "mu is defined as the physical scalar mu_phys(Y_e), so it is constant on each physical fiber by construction",
        },
        "canonical_consequences": {
            "determinant_line": "s_det(Y_e) = 3 * mu_phys(Y_e)",
            "absolute_scale": "g_e(Y_e) = exp(mu_phys(Y_e))",
            "absolute_gap": "Delta_e_abs(Y_e) = log(g_ch_shared) - mu_phys(Y_e)",
            "charged_masses": "m_i(Y_e) = exp(mu_phys(Y_e) + ell_i_centered(Y_e))",
        },
        "scope_exclusions": [
            "No theorem-grade landing from P to physical Y_e(P) is proved by this artifact.",
            "No promotion C_hat_e^{cand} -> C_hat_e is proved by this artifact.",
            "No claim about a public selected-class charged source descent is made.",
        ],
        "notes": [
            "This theorem closes the post-promotion charged lift on theorem-grade physical Y_e.",
            "It uses the determinant line of theorem-grade physical Y_e directly, so the uncentered lift is canonical rather than existential.",
            "It does not close the P-driven charged absolute lane by itself.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build the canonical charged uncentered lift theorem on theorem-grade physical Y_e."
    )
    parser.add_argument("--determinant-line", default=str(DETERMINANT_LINE_JSON))
    parser.add_argument("--anchor-section", default=str(ANCHOR_SECTION_JSON))
    parser.add_argument("--physical-equalizer", default=str(PHYSICAL_EQUALIZER_JSON))
    parser.add_argument("--output", default=str(PHYSICAL_CANONICAL_LIFT_JSON))
    args = parser.parse_args()

    artifact = build_artifact(
        load_json(Path(args.determinant_line)),
        load_json(Path(args.anchor_section)),
        load_json(Path(args.physical_equalizer)),
    )

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
