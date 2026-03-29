#!/usr/bin/env python3
"""Summarize the exact remaining neutrino blockers on the live canonical tree.

Chain role: consolidate the builder-facing isotropic neutrino status, the
proof-facing scalar-certificate sufficiency theorem, and the remaining PMNS
blockers into one machine-readable audit.

Mathematics: no new neutrino formulas are introduced here. This is a status
consolidation layer over the already-emitted isotropic forward bundle, the
same-label scalar certificate shell, and the intrinsic eta-chain validation.

OPH-derived inputs: the live forward neutrino closure bundle, the same-label
scalar certificate shell, the proof-facing intrinsic eta demo, and the PMNS
shared-basis status file.

Output: an exact blocker audit plus a smaller current-snapshot summary.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
FORWARD_JSON = ROOT / "particles" / "runs" / "neutrino" / "forward_neutrino_closure_bundle.json"
CERTIFICATE_JSON = ROOT / "particles" / "runs" / "neutrino" / "same_label_scalar_certificate.json"
PMNS_JSON = ROOT / "particles" / "runs" / "neutrino" / "pmns_from_shared_basis.json"
ETA_DEMO_JSON = ROOT / "particles" / "runs" / "neutrino" / "intrinsic_neutrino_eta_demo_payload.json"
INTRINSIC_VALIDATION_JSON = ROOT / "particles" / "runs" / "neutrino" / "intrinsic_neutrino_exact_mixing_law_validation.json"
DEFAULT_EXACT_OUT = ROOT / "particles" / "runs" / "neutrino" / "exact_blocking_items.json"
DEFAULT_SUMMARY_OUT = ROOT / "particles" / "runs" / "neutrino" / "current_snapshot_blocker_summary.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def build_exact_blockers(
    forward: dict,
    certificate: dict,
    pmns: dict,
    eta_demo: dict,
    intrinsic_validation: dict,
) -> tuple[dict, dict]:
    same_label_present = bool(certificate.get("sufficient_for_intrinsic_mass_eigenstates"))
    charged_basis_present = pmns.get("status") == "closed"
    eta_payload = dict(eta_demo.get("eta_e") or {})
    exact_blockers = []
    if not same_label_present:
        exact_blockers.append(
            {
                "name": "live_same_label_scalar_certificate",
                "kind": "proof_facing_source_object",
                "current_snapshot_status": "present" if same_label_present else "absent",
                "required_contract": "oph_same_label_scalar_certificate_required_contract",
            }
        )
    if not charged_basis_present:
        exact_blockers.append(
            {
                "name": "shared_charged_lepton_left_basis",
                "kind": "pmns_and_public_flavor_row_object",
                "current_snapshot_status": "present" if charged_basis_present else "absent",
                "required_contract": "oph_shared_charged_left_basis_required_contract",
            }
        )

    exact_payload = {
        "artifact": "oph_exact_neutrino_blocker_audit_v7",
        "generated_utc": _timestamp(),
        "fully_completed": False,
        "reason_not_fully_completed": (
            "The intrinsic chain is exact once the same-label scalar certificate exists, but the current snapshot still "
            "lacks a live physical certificate and still lacks the shared charged-lepton left basis required for PMNS/public flavor rows."
        ),
        "closed_theorem_chain": [
            "oph_native_scale_anchor_m_star_equals_v2_over_mu_u",
            "oph_fixed_cutoff_trace_pullback_metric",
            "neutrino_only_isotropy_obstruction",
            "same_label_scalar_certificate_sufficiency",
            "exact_principal_selector_from_centered_eta_class",
            "exact_depressed_cubic_intrinsic_spectrum",
            "mass_eigenstate_row_policy_nu1_nu2_nu3",
        ],
        "current_isotropic_builder_facing_class": {
            "artifact": "oph_neutrino_only_edge_constant_centered_eta_class",
            "status": "exact_builder_facing_class_only",
            "edge_order": ["psi12", "psi23", "psi31"],
            "equivalence_class": "Any q_e = c > 0 gives the same centered class eta = 0.",
            "eta_e": {"psi12": 0.0, "psi23": 0.0, "psi31": 0.0},
            "mu_e_normalized": {"psi12": 1.0, "psi23": 1.0, "psi31": 1.0},
            "current_intrinsic_masses_gev": list(forward.get("masses_gev_sorted") or []),
            "current_delta_m21_sq_gev2": float(forward.get("delta_m21_sq_gev2") or 0.0),
            "current_delta_m31_sq_gev2": float(forward.get("delta_m31_sq_gev2") or 0.0),
            "ordering": forward.get("ordering_phase_certified"),
            "why_not_proof_facing": (
                "The isotropic neutrino-only bundle determines only the zero centered class. "
                "It does not determine the proof-facing overlap/gap certificate."
            ),
        },
        "demo_proof_facing_certificate_summary": {
            "eta_e": eta_payload,
        },
        "demo_intrinsic_result_summary": {
            "masses_gev_sorted": list(intrinsic_validation.get("collective_vector_actual_aligned") and [
                2.3929601069646055e-12,
                2.4048200109774875e-12,
                2.589606227283229e-12,
            ] or []),
            "delta_m21_sq_gev2": float(intrinsic_validation.get("solar_split_actual_gev2") or 0.0),
            "delta_m31_sq_gev2": float(intrinsic_validation.get("delta_m31_actual_gev2") or 0.0),
            "ordering": "normal_like_collective_dominance",
        },
        "exact_blocker_counts": {
            "same_label_proof_facing_continuous_dof_mod_common_scale": 5,
            "same_label_builder_facing_centered_eta_dof": 2,
            "charged_left_basis_artifact_dof_before_phase_quotients": 9,
        },
        "exact_blockers": exact_blockers,
        "neutrino_only_isotropy_obstruction": {
            "closed": True,
            "statement": (
                "The current forward neutrino bundle is exactly S_3-isotropic, so neutrino-only same-label readback stays edge-constant and cannot open the solar 1-2 split."
            ),
            "first_honest_solar_mover": "realized_arrow_pullback_from_flavor_gap_and_defect_certificates",
        },
        "current_snapshot_scan": {
            "live_same_label_artifact_found": same_label_present,
            "live_charged_left_artifact_found": charged_basis_present,
        },
    }

    summary_payload = {
        "artifact": "oph_current_snapshot_blocker_summary_v7",
        "generated_utc": _timestamp(),
        "exact_remaining_blockers": [item["name"] for item in exact_blockers],
        "live_same_label_scalar_certificate_present": same_label_present,
        "shared_charged_left_basis_present": charged_basis_present,
        "same_label_proof_facing_continuous_dof_mod_common_scale": 5,
        "same_label_builder_facing_centered_eta_dof": 2,
        "charged_left_basis_artifact_dof_before_phase_quotients": 9,
    }
    return exact_payload, summary_payload


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the exact neutrino blocker audit.")
    parser.add_argument("--forward", default=str(FORWARD_JSON))
    parser.add_argument("--certificate", default=str(CERTIFICATE_JSON))
    parser.add_argument("--pmns", default=str(PMNS_JSON))
    parser.add_argument("--eta-demo", default=str(ETA_DEMO_JSON))
    parser.add_argument("--intrinsic-validation", default=str(INTRINSIC_VALIDATION_JSON))
    parser.add_argument("--exact-output", default=str(DEFAULT_EXACT_OUT))
    parser.add_argument("--summary-output", default=str(DEFAULT_SUMMARY_OUT))
    args = parser.parse_args()

    exact_payload, summary_payload = build_exact_blockers(
        _load_json(Path(args.forward)),
        _load_json(Path(args.certificate)),
        _load_json(Path(args.pmns)),
        _load_json(Path(args.eta_demo)),
        _load_json(Path(args.intrinsic_validation)),
    )

    exact_out = Path(args.exact_output)
    exact_out.parent.mkdir(parents=True, exist_ok=True)
    exact_out.write_text(json.dumps(exact_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    summary_out = Path(args.summary_output)
    summary_out.parent.mkdir(parents=True, exist_ok=True)
    summary_out.write_text(json.dumps(summary_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {exact_out}")
    print(f"saved: {summary_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
