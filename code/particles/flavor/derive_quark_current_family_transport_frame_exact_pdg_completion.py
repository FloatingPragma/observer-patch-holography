#!/usr/bin/env python3
"""Emit the merged exact-PDG completion theorem on the current-family transport-frame surface."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
STRENGTHENED_THEOREM_JSON = (
    ROOT / "particles" / "runs" / "flavor" / "quark_current_family_transport_frame_strengthened_physical_sigma_lift_theorem.json"
)
ABS_COLLAPSE_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_absolute_readout_algebraic_collapse.json"
PDG_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_exact_pdg_theorem.json"
DEFAULT_OUT = (
    ROOT / "particles" / "runs" / "flavor" / "quark_current_family_transport_frame_exact_pdg_completion.json"
)


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def build_artifact(strengthened_theorem: dict, collapse: dict, pdg: dict) -> dict:
    sigma_data = strengthened_theorem["theorem_grade_physical_sigma_datum"]
    exact_sigma = collapse["exact_current_family_sigma_specialization"]
    if abs(float(sigma_data["sigma_u"]) - float(exact_sigma["sigma_u_target"])) > 1.0e-12:
        raise ValueError("lift sigma_u target does not match algebraic-collapse specialization")
    if abs(float(sigma_data["sigma_d"]) - float(exact_sigma["sigma_d_target"])) > 1.0e-12:
        raise ValueError("lift sigma_d target does not match algebraic-collapse specialization")

    return {
        "artifact": "oph_quark_current_family_transport_frame_exact_pdg_completion",
        "generated_utc": _timestamp(),
        "proof_status": "closed_target_anchored_mixed_convention_coordinate_audit",
        "theorem_scope": "current_family_common_refinement_transport_frame_only",
        "public_promotion_allowed": False,
        "supporting_artifacts": {
            "strengthened_physical_sigma_lift_theorem": strengthened_theorem["artifact"],
            "absolute_readout_algebraic_collapse": collapse["artifact"],
            "exact_pdg_reconstruction": pdg["artifact"],
        },
        "theorem_statement": (
            "On the declared current-family common-refinement target-audit surface, the target-attached sigma lift "
            "and affine mean law reconstruct the chosen quark coordinates algebraically. The light, heavy, and top "
            "rows use different comparison conventions, so this is not one source-emitted running sextet."
        ),
        "target_anchored": True,
        "source_only_prediction": False,
        "single_running_quark_sextet_claim_allowed": False,
        "strengthened_physical_sigma_lift": {
            "artifact": strengthened_theorem["artifact"],
            "sigma_id": strengthened_theorem["restricted_sigma_ud_phys_element"]["sigma_id"],
            "canonical_token": strengthened_theorem["restricted_sigma_ud_phys_element"]["canonical_token"],
            "ckm_invariants": strengthened_theorem["restricted_sigma_ud_phys_element"]["ckm_invariants"],
            "theorem_grade_physical_sigma_datum": strengthened_theorem["theorem_grade_physical_sigma_datum"],
        },
        "forced_absolute_sector_scales": {
            "artifact": collapse["artifact"],
            "g_u_forced": float(exact_sigma["g_u_forced"]),
            "g_d_forced": float(exact_sigma["g_d_forced"]),
        },
        "exact_running_values_gev": pdg["reconstructed_current_family_running_values_gev"],
        "notes": [
            "This is a merged exact completion theorem only on the explicit current-family/common-refinement surface.",
            "It does not by itself promote the public target-free physical-sheet quark lane.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the merged current-family transport-frame exact PDG completion artifact.")
    parser.add_argument("--strengthened-theorem", default=str(STRENGTHENED_THEOREM_JSON))
    parser.add_argument("--absolute-collapse", default=str(ABS_COLLAPSE_JSON))
    parser.add_argument("--pdg", default=str(PDG_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    payload = build_artifact(
        _load_json(Path(args.strengthened_theorem)),
        _load_json(Path(args.absolute_collapse)),
        _load_json(Path(args.pdg)),
    )
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
