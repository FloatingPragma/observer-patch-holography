#!/usr/bin/env python3
"""Build the particle derivation-chain closure matrix."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
PARTICLES_ROOT = ROOT / "particles"
FINAL_PREDICTIONS = PARTICLES_ROOT / "runs" / "status" / "final_end_to_end_predictions.json"
BLIND_PROVENANCE = PARTICLES_ROOT / "runs" / "status" / "blind_prediction_provenance.json"
CARRIER_ACCEPTANCE = PARTICLES_ROOT / "runs" / "status" / "carrier_mode_acceptance.json"
CHARGED_NONCLOSURE = PARTICLES_ROOT / "runs" / "leptons" / "charged_end_to_end_impossibility_theorem.json"
QUARK_GLOBAL_OBSTRUCTION = (
    PARTICLES_ROOT / "runs" / "flavor" / "quark_class_uniform_public_frame_descent_obstruction.json"
)
DIRECT_TOP_CONTRACT = PARTICLES_ROOT / "runs" / "calibration" / "direct_top_bridge_contract.json"
HIERARCHY_RESONANCE = PARTICLES_ROOT / "hierarchy" / "certificates" / "R_local_global_hierarchy_resonance_closeout_335.json"
HIERARCHY_EW_CAPACITY = PARTICLES_ROOT / "hierarchy" / "certificates" / "R_EW_global_capacity_certificate.json"
HIERARCHY_NATURALITY = PARTICLES_ROOT / "hierarchy" / "issue_332_rg_naturality_certificate.json"
DEFAULT_JSON_OUT = PARTICLES_ROOT / "runs" / "status" / "derivation_chain_closure_matrix.json"
DEFAULT_MD_OUT = PARTICLES_ROOT / "DERIVATION_CHAIN_CLOSURE_MATRIX.md"


def _now_utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _prediction_map(payload: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {entry["particle_id"]: entry for entry in payload["predictions"]}


def _display_status(status: str) -> str:
    return status.replace("current_corpus", "corpus_limited").replace("not_live", "not_certified")


def build_payload() -> dict[str, Any]:
    final_predictions = _load_json(FINAL_PREDICTIONS)
    provenance = _load_json(BLIND_PROVENANCE)
    carrier_acceptance = _load_json(CARRIER_ACCEPTANCE)
    charged_nonclosure = _load_json(CHARGED_NONCLOSURE)
    quark_global = _load_json(QUARK_GLOBAL_OBSTRUCTION)
    direct_top = _load_json(DIRECT_TOP_CONTRACT)
    hierarchy_resonance = _load_json(HIERARCHY_RESONANCE)
    hierarchy_capacity = _load_json(HIERARCHY_EW_CAPACITY)
    hierarchy_naturality = _load_json(HIERARCHY_NATURALITY)
    predictions = _prediction_map(final_predictions)
    conditional_candidates = {
        entry["particle_id"]: entry
        for entry in final_predictions.get("conditional_nonpromotable_rows", [])
    }
    higgs_output = predictions.get("higgs") or conditional_candidates.get("higgs")
    if higgs_output is None:
        raise ValueError("missing Higgs prediction or conditional candidate row")
    hierarchy_capacity_root = hierarchy_capacity["exact_capacity_fixed_point"]
    withheld_by_id = {
        entry["particle_id"]: entry for entry in final_predictions.get("withheld_non_prediction_rows", [])
    }
    quark_withheld = [pid for pid in ("up_quark", "down_quark", "strange_quark", "charm_quark", "bottom_quark", "top_quark") if pid in withheld_by_id]

    rows = [
        {
            "chain": "p_closure_root",
            "status": "candidate_nonpromoting_root",
            "claim_level": final_predictions["p_closure"]["claim_status"],
            "outputs": {
                "P": final_predictions["p_closure"]["P"],
                "alpha_inv": final_predictions["p_closure"]["alpha_inv"],
            },
            "promotable": False,
            "required_receipts": [
                "source_emitted_same_scheme_Ward_projected_R_Q",
                "full_interval_composition_certificate",
            ],
            "scientific_boundary": "The fixed-point root is a guarded candidate and cannot feed promoted particle rows.",
            "evidence_artifacts": [
                "code/P_derivation/runtime/r_q_residual_contract_current.json",
                "code/P_derivation/runtime/fine_structure_interval_certificate_current.json",
            ],
        },
        {
            "chain": "conditional_classical_carrier_modes",
            "status": "conditional_classical_modes__quantum_particle_receipts_absent",
            "claim_level": "conditional_classical_or_perturbative_mode",
            "outputs": {
                row["carrier_id"]: {
                    "hard_quadratic_mass_parameter_squared": row[
                        "hard_quadratic_mass_parameter_squared"
                    ],
                    "classical_carrier_gate": row["classical_carrier_gate"]["status"],
                    "quantum_particle_gate": row["quantum_particle_gate"]["status"],
                }
                for row in carrier_acceptance["carriers"]
            },
            "promotable": False,
            "required_receipts": list(carrier_acceptance["quantum_required_receipts"]),
            "scientific_boundary": (
                "The conditional classical-mode theorem does not satisfy the quantum-particle "
                "gate; no 0 GeV particle row is emitted."
            ),
            "evidence_artifacts": [
                "code/particles/runs/status/carrier_mode_acceptance.json"
            ],
        },
        {
            "chain": "electroweak_massive_bosons",
            "status": "no_public_prediction__target_free_source_law_absent",
            "claim_level": "gap_frontier",
            "outputs": {},
            "promotable": False,
            "required_receipts": [
                "target_free_D10_source_law",
                "certified_P_root_endpoint_stack",
            ],
            "scientific_boundary": "No source-promoted W/Z mass row is available.",
            "evidence_artifacts": [
                "code/P_derivation/runtime/r_q_residual_contract_current.json"
            ],
        },
        {
            "chain": "hierarchy_naturality_bridge",
            "status": "conditional_hierarchy_bridge__declared_naturality_identity_not_promoted",
            "claim_level": "conditional_bridge_plus_declared_map_packaging_identity",
            "outputs": {
                "N_CRC_EW": hierarchy_capacity_root["N_CRC_EW"],
                "bridge_residual": hierarchy_capacity_root["bridge_residual"],
                "declared_map_epsilon_H": hierarchy_naturality["epsilon_H"],
            },
            "promotable": False,
            "required_receipts": hierarchy_resonance["work_in_progress_receipts"]
            + [
                "source-derived RG/Higgs comparison maps",
                "independent evaluation of both naturality residuals",
            ],
            "evidence_artifacts": [
                "code/particles/hierarchy/certificates/R_local_global_hierarchy_resonance_closeout_335.json",
                "code/particles/hierarchy/certificates/R_EW_global_capacity_certificate.json",
                "code/particles/hierarchy/issue_332_rg_naturality_certificate.json",
            ],
            "resonance_status": hierarchy_resonance["status"],
            "full_theorem_grade_resonance_promoted": hierarchy_resonance[
                "full_theorem_grade_resonance_promoted"
            ],
            "work_in_progress_receipts": hierarchy_resonance["work_in_progress_receipts"],
            "claim_boundary": (
                "The hierarchy implication is exact under the declared screen premises. The stored "
                "zero naturality defect follows only for separately declared comparison maps; the "
                "source antecedents admit nonzero-defect completions, so source-derived maps and an "
                "independent residual evaluation remain open. The public Thomson endpoint, W/Z "
                "mass promotion, charged-lepton absolute masses, source-only hadrons, Strong CP, "
                "and full SI G are separate surfaces."
            ),
        },
        {
            "chain": "higgs_top_declared_surface",
            "status": (
                "declared_d10_d11_surface__direct_top_no_go"
                if higgs_output["promotable"]
                else "conditional_declared_surface_higgs_top_candidate"
            ),
            "claim_level": higgs_output["exact_kind"],
            "outputs": {
                "higgs": higgs_output["value"],
            },
            "promotable": higgs_output["promotable"],
            "required_receipts": []
            if higgs_output["promotable"]
            else ["EWTargetFreeRepairValueLaw_D10"],
            "evidence_artifacts": [
                "code/particles/runs/calibration/direct_top_bridge_contract.json"
            ],
            "codomain_obstruction": direct_top.get("status"),
        },
        {
            "chain": "charged_leptons",
            "status": "corpus_limited_charged_end_to_end_no_go",
            "claim_level": "target_anchored_witness_withheld_no_public_charged_mass_output",
            "outputs": {},
            "withheld_non_prediction_rows": [
                withheld_by_id[pid] for pid in ("electron", "muon", "tau") if pid in withheld_by_id
            ],
            "promotable": False,
            "required_receipts": [],
            "evidence_artifacts": [
                "code/particles/runs/leptons/charged_end_to_end_impossibility_theorem.json"
            ],
            "nonclosure_theorem": charged_nonclosure.get("artifact"),
        },
        {
            "chain": "selected_class_quarks",
            "status": "quark_source_nonidentifiability__numeric_rows_withheld",
            "claim_level": (
                withheld_by_id.get("top_quark", {}).get("claim_tier")
                or predictions.get("top_quark", {}).get("exact_kind")
                or "source_spread_nonidentifiability_obstruction"
            ),
            "outputs": {},
            "withheld_non_prediction_rows": [withheld_by_id[pid] for pid in quark_withheld],
            "promotable": bool(predictions.get("top_quark", {}).get("promotable", False)),
            "required_receipts": [
                "QUARK_SOURCE_SPREAD_PAIR_ACTION_BREAKING_THEOREM",
                "QUARK_RG_COVARIANT_TRAJECTORY_OR_INVARIANT",
                "QUARK_OPERATIONAL_SCHEME_AND_SCALE_SECTION",
                "QUARK_THRESHOLD_AND_TOP_CONVERSION",
                "QUARK_COMMON_SCALE_DIMENSIONLESS_YUKAWA_CERTIFICATE",
            ],
            "evidence_artifacts": [
                "code/particles/runs/flavor/quark_sigma_source_nonidentifiability_obstruction.json",
                "code/particles/runs/flavor/quark_running_mass_scheme_convention_obstruction.json",
                "code/particles/runs/flavor/quark_class_uniform_public_frame_descent_obstruction.json",
            ],
            "global_classification_obstruction": quark_global.get("proof_status"),
            "source_spread_obstruction": final_predictions["quark_sigma_source_boundary"]["obstruction_artifact"],
            "scheme_and_yukawa_obstruction": final_predictions[
                "quark_scheme_and_yukawa_boundary"
            ]["artifact"],
        },
        {
            "chain": "neutrino_absolute_attachment",
            "status": (
                "closed_weighted_cycle_absolute_attachment_with_comparison_tension_visible"
                if predictions.get("electron_neutrino", {}).get("promotable")
                else "rejected_target_informed_weighted_cycle_candidate"
            ),
            "claim_level": (
                predictions.get("electron_neutrino", {}).get("exact_kind")
                or withheld_by_id.get("electron_neutrino", {}).get("exact_kind")
                or "rejected_target_informed_weighted_cycle_candidate"
            ),
            "outputs": {
                particle: predictions[particle]["value"]
                for particle in ("electron_neutrino", "muon_neutrino", "tau_neutrino")
                if particle in predictions
            },
            "withheld_non_prediction_rows": [
                withheld_by_id[pid]
                for pid in ("electron_neutrino", "muon_neutrino", "tau_neutrino")
                if pid in withheld_by_id
            ],
            "unit": "eV",
            "promotable": bool(predictions.get("electron_neutrino", {}).get("promotable", False)),
            "required_receipts": []
            if predictions.get("electron_neutrino", {}).get("promotable")
            else [
                "source_emitted_family_transport_kernel",
                "source_derived_weight_exponent_and_cycle_matrix_law",
                "source_derived_basis_permutation_and_holonomy_orientation",
                "pre_reference_hash_lock",
            ],
            "evidence_artifacts": [],
        },
        {
            "chain": "hadrons",
            "status": "source_backend_absent__empirical_output_class_separate",
            "claim_level": "source_only_absent_empirical_rows_separate",
            "outputs": {},
            "promotable": False,
            "required_receipts": ["source_only_hadron_backend"],
            "scientific_boundary": (
                "Source-only hadron prediction requires a working OPH hadron backend. Empirical "
                "hadron closure uses a separate e+e- payload class."
            ),
            "evidence_artifacts": [
                "docs/HADRON.md",
                "code/particles/hadron/empirical_ee_hadrons_sources.yaml",
                "code/particles/hadron/empirical_ee_hadronic_spectral_measure.schema.json",
            ],
        },
    ]

    return {
        "artifact": "oph_particle_derivation_chain_classification_matrix",
        "generated_utc": _now_utc(),
        "classification": "scientific_chain_boundaries_and_required_receipts",
        "classification_summary": {
            "promotable_chains": [row["chain"] for row in rows if row["promotable"]],
            "nonpromoting_chains": [row["chain"] for row in rows if not row["promotable"]],
            "source_backend_absent_chains": ["hadrons"],
            "empirical_output_class_chains": ["hadrons"],
            "policy": (
                "Do not promote candidate, compare-only, witness-only, corpus-limited no-go, or "
                "source-backend-absent chains as theorem predictions."
            ),
        },
        "source_artifacts": {
            "final_predictions": "code/particles/runs/status/final_end_to_end_predictions.json",
            "blind_provenance": "code/particles/runs/status/blind_prediction_provenance.json",
            "hierarchy_local_global_resonance": (
                "code/particles/hierarchy/certificates/R_local_global_hierarchy_resonance_closeout_335.json"
            ),
            "hierarchy_ew_capacity": (
                "code/particles/hierarchy/certificates/R_EW_global_capacity_certificate.json"
            ),
            "hierarchy_higgs_naturality": "code/particles/hierarchy/issue_332_rg_naturality_certificate.json",
        },
        "provenance": {
            "artifact": provenance["artifact"],
            "promotion_allowed": provenance["promotion_allowed"],
        },
        "rows": rows,
    }


def render_markdown(payload: dict[str, Any]) -> str:
    lines = [
        "# Particle Derivation Chain Classification Matrix",
        "",
        f"Generated: `{payload['generated_utc']}`",
        "",
        f"Classification: `{payload['classification']}`",
        f"Promotable chains: `{', '.join(payload['classification_summary']['promotable_chains']) or 'none'}`",
        f"Non-promoting chains: `{', '.join(payload['classification_summary']['nonpromoting_chains'])}`",
        "",
        "QFT status is orthogonal to this numeric closure matrix. The namespaced",
        "`SM_QFT_*` oracle proves conditional implications and exact strict W/Z algebra,",
        "but QFT-Q2 and QFT-Q3 remain parallel and QFT-Q4 consumes a separate",
        "nonperturbative tower. None of those theorem checks supplies the source action,",
        "physical-current amplitude, numerical freeze, or resonance sheet required to",
        "promote `electroweak_massive_bosons`.",
        "",
        "| Chain | Classification | Promotable | Required receipts | Outputs | Evidence artifacts |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for row in payload["rows"]:
        outputs = ", ".join(f"`{key}={value}`" for key, value in row["outputs"].items()) or "n/a"
        receipts = ", ".join(f"`{item}`" for item in row["required_receipts"]) or "none"
        evidence = ", ".join(f"`{item}`" for item in row["evidence_artifacts"]) or "none"
        lines.append(
            f"| `{row['chain']}` | `{_display_status(row['status'])}` | `{row['promotable']}` | {receipts} | {outputs} | {evidence} |"
        )
    return "\n".join(lines)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build particle derivation-chain closure matrix.")
    parser.add_argument("--json-out", default=str(DEFAULT_JSON_OUT))
    parser.add_argument("--markdown-out", default=str(DEFAULT_MD_OUT))
    parser.add_argument("--print-json", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    payload = build_payload()
    json_text = json.dumps(payload, indent=2, sort_keys=True) + "\n"

    json_out = Path(args.json_out)
    json_out.parent.mkdir(parents=True, exist_ok=True)
    json_out.write_text(json_text, encoding="utf-8")

    markdown_out = Path(args.markdown_out)
    markdown_out.write_text(render_markdown(payload) + "\n", encoding="utf-8")

    if args.print_json:
        print(json_text, end="")
    else:
        print(f"saved: {json_out}")
        print(f"saved: {markdown_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
