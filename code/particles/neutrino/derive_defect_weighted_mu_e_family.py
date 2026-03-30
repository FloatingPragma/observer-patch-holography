#!/usr/bin/env python3
"""Export the defect-weighted Majorana edge-weight family candidate."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SCALAR = ROOT / "particles" / "runs" / "neutrino" / "majorana_overlap_defect_scalar_evaluator.json"
DEFAULT_FORWARD = ROOT / "particles" / "runs" / "neutrino" / "forward_neutrino_closure_bundle.json"
DEFAULT_READBACK = ROOT / "particles" / "runs" / "neutrino" / "realized_same_label_gap_defect_readback.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "neutrino" / "defect_weighted_mu_e_family.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the defect-weighted mu_e family artifact.")
    parser.add_argument("--scalar", default=str(DEFAULT_SCALAR))
    parser.add_argument("--forward", default=str(DEFAULT_FORWARD))
    parser.add_argument("--readback", default=str(DEFAULT_READBACK))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    scalar = json.loads(Path(args.scalar).read_text(encoding="utf-8"))
    forward = json.loads(Path(args.forward).read_text(encoding="utf-8"))
    readback_path = Path(args.readback)
    readback = json.loads(readback_path.read_text(encoding="utf-8")) if readback_path.exists() else {}
    masses = [float(x) for x in forward.get("masses_gev_sorted", [])]
    m0 = sum(masses[:2]) / 2.0 if len(masses) >= 2 else None
    heavy_light_gap = (masses[2] - m0) if len(masses) >= 3 and m0 is not None else None
    mu_nu = float(scalar.get("mu_nu", 0.0))
    readback_complete = readback.get("payload_status") == "complete_from_live_flavor_artifacts"
    smallest_constructive_missing = None if readback_complete else "oph_realized_same_label_gap_defect_readback"
    strict_repo_missing_object = None if readback_complete else "oph_same_label_overlap_defect_log_source"

    artifact = {
        "artifact": "oph_defect_weighted_majorana_edge_weight_family",
        "generated_utc": _timestamp(),
        "parent_theorem_id": scalar.get("theorem_candidate_id"),
        "upstream_exact_clause": scalar.get("required_overlap_certificate"),
        "normalizer_artifact": "oph_same_label_overlap_defect_weight_normalizer",
        "realized_same_label_gap_defect_readback_status": (
            readback.get("payload_status") if readback else "missing_readback_artifact"
        ),
        "smallest_constructive_missing_object": smallest_constructive_missing,
        "strict_repo_missing_object": strict_repo_missing_object,
        "raw_edge_score_emitter": "oph_same_label_overlap_defect_log_source",
        "selector_center": scalar.get("selector_center"),
        "selector_point_absolute": scalar.get("selector_point_absolute"),
        "kernel_choice": "1-cos",
        "base_mu_nu": mu_nu,
        "raw_defect_source": scalar.get("overlap_nonvanishing_witness_hint"),
        "realized_same_label_arrows": ["psi12", "psi23", "psi31"],
        "neutrino_only_isotropy_obstruction": True,
        "same_label_readback_origin": "realized_arrow_pullback_from_flavor_gap_and_defect_certificates",
        "majorana_matrix_real": forward.get("majorana_matrix_real"),
        "majorana_matrix_imag": forward.get("majorana_matrix_imag"),
        "raw_edge_score_symbol": "q_e > 0",
        "raw_edge_score_family": "q_e(alpha) = gap_e^alpha * defect_e^(1-alpha)",
        "canonical_alpha": 0.5,
        "raw_edge_score_rule": "q_e = sqrt(gap_e * defect_e)",
        "same_label_overlap_sq": {"psi12": None, "psi23": None, "psi31": None},
        "same_label_gap_witness": {"psi12": None, "psi23": None, "psi31": None},
        "same_label_defect_rule": "d_e = 0.5 * ||P_head(e) - L_e P_tail(e) L_e^*||_HS^2 = 1 - overlap_sq_e",
        "same_label_defect_witness": {"psi12": None, "psi23": None, "psi31": None},
        "raw_edge_score": {"psi12": None, "psi23": None, "psi31": None},
        "centered_log_rule": "eta_e = 0.5 * ((log g_e + log d_e) - mean_f(log g_f + log d_f))",
        "defect_log_centered": {"psi12": None, "psi23": None, "psi31": None},
        "weight_rule": "mu_e = base_mu_nu * exp(delta_e) / mean_f(exp(delta_f))",
        "edge_weights": {"psi12": None, "psi23": None, "psi31": None},
        "residual_hessian_formula_2x2": [
            ["mu12 + mu31", "mu31"],
            ["mu31", "mu23 + mu31"],
        ],
        "mean_preserved": True,
        "positive_weights": True,
        "isotropic_limit_recovered": True,
        "current_doublet_center_gev": m0,
        "current_heavy_light_gap_gev": heavy_light_gap,
        "current_delta_m21_sq_gev2": forward.get("delta_m21_sq_gev2"),
        "current_delta_m31_sq_gev2": forward.get("delta_m31_sq_gev2"),
        "delta_m21_sq_expected_status": "strictly_positive_for_generic_nonconstant_delta",
        "first_order_solar_response_formula": "delta_m21_sq ~= 4 * m0 * sigma_q",
        "first_order_solar_response_coefficient_gev": 4.0 * m0 if m0 is not None else None,
        "ordering_expected_status": "normal_like_stable_for_small_defect_anisotropy",
        "proof_status": "candidate_only",
        "notes": [
            "This is the first local mass-moving object that can lift the current 1-2 near-degeneracy without changing the centered selector, kernel choice, or edge-character origin candidate.",
            "The current forward neutrino bundle is S_3-isotropic, so any same-label scalar readback built only from the neutrino payload is forced to stay edge-constant.",
            "The best reduced family is a realized-arrow readback of same-label gap and defect witnesses, followed by the canonical raw score q_e = sqrt(gap_e * defect_e), centered-log lift, and mean-preserving mu_e family.",
            "The current canonical no-new-parameter point is q_e = sqrt(gap_e * defect_e), but the actual next mover is the realized flavor-side same-label gap/defect readback that feeds that rule.",
            "The exact theorem blocker remains on same-label overlap / edge-bundle normalization.",
            (
                "On the live corpus the realized same-label gap/defect readback is already complete from flavor-side certificates, "
                "so this family no longer treats that readback as the next missing constructive object."
                if readback_complete
                else "The smallest spectrum-moving local object is still the realized same-label gap/defect readback beneath the repo-facing log-source wrapper."
            ),
        ],
    }

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
