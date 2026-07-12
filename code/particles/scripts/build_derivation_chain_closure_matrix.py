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
PIPELINE_STATUS = PARTICLES_ROOT / "runs" / "status" / "particle_pipeline_closure_status.json"
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


def _issue_map(payload: dict[str, Any]) -> dict[int, dict[str, Any]]:
    return {int(gate["issue"]): gate for gate in payload["issue_gates"]}


def _display_status(status: str) -> str:
    return status.replace("current_corpus", "corpus_limited").replace("not_live", "not_certified")


def build_payload() -> dict[str, Any]:
    final_predictions = _load_json(FINAL_PREDICTIONS)
    pipeline = _load_json(PIPELINE_STATUS)
    provenance = _load_json(BLIND_PROVENANCE)
    carrier_acceptance = _load_json(CARRIER_ACCEPTANCE)
    charged_nonclosure = _load_json(CHARGED_NONCLOSURE)
    quark_global = _load_json(QUARK_GLOBAL_OBSTRUCTION)
    direct_top = _load_json(DIRECT_TOP_CONTRACT)
    hierarchy_resonance = _load_json(HIERARCHY_RESONANCE)
    hierarchy_capacity = _load_json(HIERARCHY_EW_CAPACITY)
    hierarchy_naturality = _load_json(HIERARCHY_NATURALITY)
    predictions = _prediction_map(final_predictions)
    gates = _issue_map(pipeline)
    hierarchy_capacity_root = hierarchy_capacity["exact_capacity_fixed_point"]
    withheld_by_id = {
        entry["particle_id"]: entry for entry in final_predictions.get("withheld_non_prediction_rows", [])
    }
    quark_withheld = [pid for pid in ("up_quark", "down_quark", "strange_quark", "charm_quark", "bottom_quark", "top_quark") if pid in withheld_by_id]

    rows = [
        {
            "chain": "p_closure_root",
            "status": "candidate_not_live_root",
            "claim_level": final_predictions["p_closure"]["claim_status"],
            "outputs": {
                "P": final_predictions["p_closure"]["P"],
                "alpha_inv": final_predictions["p_closure"]["alpha_inv"],
            },
            "promotable": False,
            "open_gates": [],
            "closed_issue_refs": [224],
            "stage_gate": "populated source spectral measure payload + full interval certificate",
            "next_artifact": "code/P_derivation/runtime/r_q_residual_contract_current.json + code/P_derivation/runtime/fine_structure_interval_certificate_current.json",
        },
        {
            "chain": "conditional_classical_carrier_modes",
            "status": "closed_scoped_classical_modes_quantum_particle_gate_open",
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
            "open_gates": list(carrier_acceptance["quantum_required_receipts"]),
            "closure_reason": (
                "Issue #536 is closed by separating the conditional classical-mode theorem from the "
                "unmet quantum-particle gate; no 0 GeV particle row is emitted."
            ),
            "next_artifact": "code/particles/runs/status/carrier_mode_acceptance.json",
        },
        {
            "chain": "electroweak_massive_bosons",
            "status": "no_public_prediction_row_target_free_d10_repair_not_promoted",
            "claim_level": "gap_frontier",
            "outputs": {},
            "promotable": False,
            "open_gates": [],
            "closed_issue_refs": [224],
            "stage_gate": "promotable target-free D10 repair plus certified P/root endpoint stack",
            "next_artifact": "code/P_derivation/runtime/r_q_residual_contract_current.json",
        },
        {
            "chain": "hierarchy_naturality_bridge",
            "status": "closed_selected_branch_local_global_hierarchy_naturality",
            "claim_level": "selected_branch_source_side_theorem",
            "outputs": {
                "N_CRC_EW": hierarchy_capacity_root["N_CRC_EW"],
                "bridge_residual": hierarchy_capacity_root["bridge_residual"],
                "epsilon_H": hierarchy_naturality["epsilon_H"],
            },
            "promotable": True,
            "open_gates": [],
            "closed_issue_refs": [332, 335, 337],
            "next_artifact": None,
            "resonance_status": hierarchy_resonance["status"],
            "full_theorem_grade_resonance_promoted": hierarchy_resonance[
                "full_theorem_grade_resonance_promoted"
            ],
            "remaining_promotion_gates": hierarchy_resonance["remaining_promotion_gates"],
            "claim_boundary": (
                "Promotes the hierarchy and Higgs naturality rows only. The public Thomson endpoint, "
                "W/Z mass promotion, charged-lepton absolute masses, source-only hadrons, Strong CP, "
                "and full SI G remain separate surfaces."
            ),
        },
        {
            "chain": "higgs_top_declared_surface",
            "status": (
                "closed_on_declared_d10_d11_surface_direct_top_no_go"
                if predictions["higgs"]["promotable"]
                else "conditional_declared_surface_higgs_top_candidate"
            ),
            "claim_level": predictions["higgs"]["exact_kind"],
            "outputs": {
                "higgs": predictions["higgs"]["value"],
            },
            "promotable": predictions["higgs"]["promotable"],
            "open_gates": []
            if predictions["higgs"]["promotable"]
            else ["closed_promotable_EWTargetFreeRepairValueLaw_D10"],
            "closed_issue_refs": [207],
            "next_artifact": "code/particles/runs/calibration/direct_top_bridge_contract.json",
            "codomain_obstruction": direct_top.get("status"),
        },
        {
            "chain": "charged_leptons",
            "status": "closed_current_corpus_charged_end_to_end_no_go",
            "claim_level": "target_anchored_witness_withheld_no_public_charged_mass_output",
            "outputs": {},
            "withheld_non_prediction_rows": [
                withheld_by_id[pid] for pid in ("electron", "muon", "tau") if pid in withheld_by_id
            ],
            "promotable": False,
            "open_gates": [],
            "closed_issue_refs": [201],
            "next_artifact": "code/particles/runs/leptons/charged_end_to_end_impossibility_theorem.json",
            "nonclosure_theorem": charged_nonclosure.get("artifact"),
        },
        {
            "chain": "selected_class_quarks",
            "status": "closed_sharper_quark_source_obstructions_numeric_rows_withheld",
            "claim_level": (
                withheld_by_id.get("top_quark", {}).get("claim_tier")
                or predictions.get("top_quark", {}).get("exact_kind")
                or "source_spread_nonidentifiability_obstruction"
            ),
            "outputs": {},
            "withheld_non_prediction_rows": [withheld_by_id[pid] for pid in quark_withheld],
            "promotable": bool(predictions.get("top_quark", {}).get("promotable", False)),
            "open_gates": [
                "QUARK_SOURCE_SPREAD_PAIR_ACTION_BREAKING_THEOREM",
                "QUARK_RG_COVARIANT_TRAJECTORY_OR_INVARIANT",
                "QUARK_OPERATIONAL_SCHEME_AND_SCALE_SECTION",
                "QUARK_THRESHOLD_AND_TOP_CONVERSION",
                "QUARK_COMMON_SCALE_DIMENSIONLESS_YUKAWA_CERTIFICATE",
            ],
            "closed_issue_refs": [199, 207, 212, 377, 379, 380, 381, 382],
            "next_artifact": "code/particles/runs/flavor/quark_sigma_source_nonidentifiability_obstruction.json",
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
            "open_gates": []
            if predictions.get("electron_neutrino", {}).get("promotable")
            else [
                "source_emitted_family_transport_kernel",
                "source_derived_weight_exponent_and_cycle_matrix_law",
                "source_derived_basis_permutation_and_holonomy_orientation",
                "pre_reference_hash_lock",
            ],
            "next_artifact": None,
        },
        {
            "chain": "hadrons",
            "status": "source_backend_absent_empirical_closure_policy_emitted",
            "claim_level": "source_only_absent_empirical_rows_separate",
            "outputs": {},
            "promotable": False,
            "open_gates": [],
            "closed_issue_refs": [153, 157],
            "closure_reason": (
                "Source-only hadron prediction requires a working OPH hadron backend. Empirical "
                "hadron closure uses a separate e+e- payload class."
            ),
            "next_artifact": (
                "HADRON.md + code/particles/hadron/empirical_ee_hadrons_sources.yaml + "
                "code/particles/hadron/empirical_ee_hadronic_spectral_measure.schema.json"
            ),
        },
    ]
    closed_or_scoped = [
        row["chain"]
        for row in rows
        if row["status"]
        in {
            "closed_scoped_classical_modes_quantum_particle_gate_open",
            "closed_on_declared_d10_d11_surface_direct_top_no_go",
            "closed_current_corpus_charged_end_to_end_no_go",
            "closed_selected_public_class_global_classification_no_go",
            "closed_selected_branch_local_global_hierarchy_naturality",
            "closed_weighted_cycle_absolute_attachment_with_comparison_tension_visible",
        }
    ]
    remaining_nonclosed = [
        row["chain"]
        for row in rows
        if row["chain"] not in closed_or_scoped
        and row["status"] != "source_backend_absent_empirical_closure_policy_emitted"
    ]

    return {
        "artifact": "oph_particle_derivation_chain_closure_matrix",
        "generated_utc": _now_utc(),
        "status": "executable_nonhadron_chain_matrix_emitted",
        "closure_summary": {
            "all_derivation_chains_claimed_closed": False,
            "closed_or_scoped_chains": closed_or_scoped,
            "remaining_nonclosed_chains": remaining_nonclosed,
            "source_backend_absent_chains": ["hadrons"],
            "empirical_closure_policy_chains": ["hadrons"],
            "policy": (
                "Do not promote candidate, compare-only, witness-only, corpus-limited no-go, or "
                "source-backend-absent chains as closed theorem predictions."
            ),
        },
        "source_artifacts": {
            "final_predictions": "code/particles/runs/status/final_end_to_end_predictions.json",
            "pipeline_status": "code/particles/runs/status/particle_pipeline_closure_status.json",
            "blind_provenance": "code/particles/runs/status/blind_prediction_provenance.json",
            "hierarchy_local_global_resonance": (
                "code/particles/hierarchy/certificates/R_local_global_hierarchy_resonance_closeout_335.json"
            ),
            "hierarchy_ew_capacity": (
                "code/particles/hierarchy/certificates/R_EW_global_capacity_certificate.json"
            ),
            "hierarchy_higgs_naturality": "code/particles/hierarchy/issue_332_rg_naturality_certificate.json",
        },
        "worker_policy": {
            "chrome_pro_workers_needed_now": False,
            "reason": (
                "Hadron issues #153/#157 have a source-backend boundary and an empirical closure policy surface; "
                "the P/electromagnetic source-only endpoint needs a populated source spectral measure payload and interval certificate. "
                "The hierarchy/naturality bridge is closed on the selected branch. "
                "The charged absolute-anchor, quark global-classification, and direct-top auxiliary-codomain lanes "
                "have corpus-limited no-go boundaries."
            ),
        },
        "particle_five_gates": {
            str(issue): gates[issue]
            for issue in (
                32,
                153,
                199,
                201,
                207,
                223,
                224,
                225,
                234,
                235,
                332,
                335,
                337,
                377,
                378,
                379,
                380,
                381,
                382,
            )
            if issue in gates
        },
        "provenance_status": provenance["status"],
        "rows": rows,
    }


def render_markdown(payload: dict[str, Any]) -> str:
    lines = [
        "# Particle Derivation Chain Closure Matrix",
        "",
        f"Generated: `{payload['generated_utc']}`",
        "",
        f"Status: `{payload['status']}`",
        f"All derivation chains claimed closed: `{payload['closure_summary']['all_derivation_chains_claimed_closed']}`",
        f"Remaining nonclosed chains: `{', '.join(payload['closure_summary']['remaining_nonclosed_chains'])}`",
        f"Chrome Pro workers needed: `{payload['worker_policy']['chrome_pro_workers_needed_now']}`",
        f"Reason: {payload['worker_policy']['reason']}",
        "",
        "| Chain | Status | Promotable | Open gates | Outputs | Next artifact |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for row in payload["rows"]:
        outputs = ", ".join(f"`{key}={value}`" for key, value in row["outputs"].items()) or "n/a"
        gates = ", ".join(f"#{issue}" for issue in row["open_gates"]) or "none"
        next_artifact = row["next_artifact"] or "none"
        lines.append(
            f"| `{row['chain']}` | `{_display_status(row['status'])}` | `{row['promotable']}` | {gates} | {outputs} | {next_artifact} |"
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
