#!/usr/bin/env python3
"""Build the final end-to-end particle prediction bundle."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
PARTICLES_ROOT = ROOT / "particles"
P_ROOT = ROOT / "P_derivation"
P_TRUNK = P_ROOT / "runtime" / "p_closure_trunk_current.json"
MEASURED_ENDPOINT = P_ROOT / "runtime" / "measured_endpoint_calibration_current.json"
PIPELINE_STATUS = PARTICLES_ROOT / "runs" / "status" / "particle_pipeline_closure_status.json"
EXACT_NONHADRON = PARTICLES_ROOT / "exact_nonhadron_masses.json"
RESULTS_STATUS = PARTICLES_ROOT / "results_status.json"
DIRECT_TOP = PARTICLES_ROOT / "runs" / "calibration" / "direct_top_bridge_contract.json"
EMPIRICAL_EE_REGISTRY = PARTICLES_ROOT / "hadron" / "empirical_ee_hadrons_sources.yaml"
EMPIRICAL_EE_SCHEMA = PARTICLES_ROOT / "hadron" / "empirical_ee_hadronic_spectral_measure.schema.json"
HIERARCHY_ROOT = PARTICLES_ROOT / "hierarchy"
HIERARCHY_RESONANCE = HIERARCHY_ROOT / "certificates" / "R_local_global_hierarchy_resonance_closeout_335.json"
HIERARCHY_EW_CAPACITY = HIERARCHY_ROOT / "certificates" / "R_EW_global_capacity_certificate.json"
HIERARCHY_READBACK = HIERARCHY_ROOT / "certificates" / "R_readback_resolution_certificate.json"
HIERARCHY_M_REP = HIERARCHY_ROOT / "certificates" / "R_m_rep_24_certificate.json"
HIERARCHY_SCREEN_SIEVE = HIERARCHY_ROOT / "certificates" / "R_screen_sieve_icosahedral_certificate.json"
HIERARCHY_NATURALITY = HIERARCHY_ROOT / "issue_332_rg_naturality_certificate.json"
DEFAULT_JSON_OUT = PARTICLES_ROOT / "runs" / "status" / "final_end_to_end_predictions.json"
DEFAULT_MD_OUT = PARTICLES_ROOT / "FINAL_END_TO_END_PREDICTIONS.md"


PARTICLE_ORDER = [
    "photon",
    "gluon",
    "graviton",
    "w_boson",
    "z_boson",
    "higgs",
    "electron",
    "muon",
    "tau",
    "up_quark",
    "down_quark",
    "strange_quark",
    "charm_quark",
    "bottom_quark",
    "top_quark",
    "electron_neutrino",
    "muon_neutrino",
    "tau_neutrino",
]


def _now_utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _load_optional_json(path: Path) -> dict[str, Any] | None:
    if not path.exists():
        return None
    return _load_json(path)


def _rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def _display_status(status: str) -> str:
    if status is None:
        return "missing"
    return status.replace("current_corpus", "corpus_limited")


def _prediction_entry(entry: dict[str, Any]) -> dict[str, Any]:
    if entry.get("mass_gev") is not None:
        value = float(entry["mass_gev"])
        unit = "GeV"
    else:
        value = float(entry["mass_eV"])
        unit = "eV"
    return {
        "particle_id": entry["particle_id"],
        "label": entry.get("label"),
        "value": value,
        "unit": unit,
        "exact_kind": entry.get("exact_kind"),
        "scope": entry.get("scope"),
        "promotable": entry.get("promotable"),
        "source_artifact": entry.get("source_artifact"),
        "supporting_scope_closure_artifact": entry.get("supporting_scope_closure_artifact"),
    }


def _fine_structure_surface(measured_endpoint: dict[str, Any]) -> dict[str, Any]:
    calibrated = measured_endpoint["calibrated_values"]
    source_candidate = measured_endpoint["current_source_candidate"]
    endpoint_requirement = measured_endpoint["codata_mapped_endpoint_requirement"]
    return {
        "source_only_oph": {
            "row_class": "source_only_oph",
            "status": "source_audit_without_hadronic_spectral_payload",
            "alpha_inv": source_candidate["alpha_inv"],
            "alpha": source_candidate["alpha"],
            "P": source_candidate["P"],
            "missing_inverse_alpha_units": source_candidate["missing_inverse_alpha_units"],
            "p_gap_implemented_minus_empirical": source_candidate["p_gap_implemented_minus_calibrated"],
            "promotable_as_exact_source_theorem": False,
        },
        "oph_plus_empirical_hadron_closure": {
            "row_class": "oph_plus_empirical_hadron_closure",
            "status": measured_endpoint["status"],
            "alpha_inv": calibrated["alpha_inv_0"],
            "alpha_inv_standard_uncertainty": calibrated["alpha_inv_0_standard_uncertainty"],
            "alpha": calibrated["alpha_0"],
            "P": calibrated["P_from_outer_equation"],
            "P_standard_uncertainty": calibrated["P_standard_uncertainty"],
            "required_transport_delta_alpha_inv": endpoint_requirement["required_transport_delta_alpha_inv"],
            "implemented_transport_delta_alpha_inv": endpoint_requirement["implemented_transport_delta_alpha_inv"],
            "missing_source_transport_delta_alpha_inv": endpoint_requirement[
                "missing_source_transport_delta_alpha_inv"
            ],
            "required_screening_factor": endpoint_requirement["required_screening_factor"],
            "residual_second_order_coefficient": endpoint_requirement[
                "residual_second_order_coefficient"
            ],
            "promotable_as_exact_source_theorem": False,
        },
        "empirical_payload_policy": measured_endpoint["empirical_hadron_closure"],
    }


def _hierarchy_surface() -> dict[str, Any]:
    resonance = _load_optional_json(HIERARCHY_RESONANCE)
    ew_capacity = _load_optional_json(HIERARCHY_EW_CAPACITY)
    readback = _load_optional_json(HIERARCHY_READBACK)
    m_rep = _load_optional_json(HIERARCHY_M_REP)
    screen_sieve = _load_optional_json(HIERARCHY_SCREEN_SIEVE)
    naturality = _load_optional_json(HIERARCHY_NATURALITY)

    exact_capacity = dict((ew_capacity or {}).get("exact_capacity_fixed_point") or {})
    source_values = dict((ew_capacity or {}).get("source_values") or {})
    exact_statement = dict((resonance or {}).get("exact_surviving_statement") or {})
    target_relation = dict((resonance or {}).get("target_relation") or {})
    readback_resolution = dict((resonance or {}).get("finite_readback_resolution_certificate") or {})
    round_count = dict((resonance or {}).get("round_count_certificate") or {})
    screen_sieve_summary = dict((resonance or {}).get("screen_sieve_certificate") or {})

    return {
        "status": {
            "resonance_status": (resonance or {}).get("status"),
            "accepted": bool((resonance or {}).get("accepted", False)),
            "full_theorem_grade_resonance_promoted": bool(
                (resonance or {}).get("full_theorem_grade_resonance_promoted", False)
            ),
            "remaining_promotion_gates": list((resonance or {}).get("remaining_promotion_gates") or []),
            "ew_capacity_status": (ew_capacity or {}).get("status"),
            "readback_status": (readback or {}).get("status") or readback_resolution.get("status"),
            "m_rep_status": (m_rep or {}).get("status") or round_count.get("status"),
            "screen_sieve_status": (screen_sieve or {}).get("status") or screen_sieve_summary.get("status"),
            "higgs_naturality_status": (naturality or {}).get("status") or (
                "closed_exact_selected_branch" if (naturality or {}).get("accepted") else None
            ),
        },
        "source_values": {
            "P_star": source_values.get("P_star"),
            "alpha_U": source_values.get("alpha_U"),
            "alpha_U_interval": source_values.get("alpha_U_interval"),
        },
        "local_global_bridge": {
            "target_relation": target_relation,
            "projection_map": exact_statement.get("projection_map"),
            "bridge_residual_formula": exact_statement.get("bridge_residual"),
            "exact_bridge_target": exact_statement.get("exact_bridge_target"),
            "N_CRC_EW": exact_capacity.get("N_CRC_EW") or exact_statement.get("N_EW_public_endpoint"),
            "bridge_residual": exact_capacity.get("bridge_residual"),
            "fixed_point_residual_x": exact_capacity.get("fixed_point_residual_x"),
            "v_over_E_cell_source": exact_capacity.get("v_over_E_cell_source"),
            "v_identity_error": exact_capacity.get("v_identity_error"),
        },
        "factor_origins": {
            "screen_ports": screen_sieve_summary.get("orbit_size", 12),
            "curvature_charge": screen_sieve_summary.get("total_curvature_charge", 12),
            "m_rep": round_count.get("m_rep", 24),
            "specialized_exponent": round_count.get("specialized_exponent", "-1/48"),
            "higgs_naturality_defect": (naturality or {}).get("epsilon_H"),
        },
        "claim_boundary": {
            "improves": (
                "The hierarchy and Higgs naturality rows are promoted as selected-branch "
                "source-side results with zero bridge residual and epsilon_H=0."
            ),
            "does_not_promote": [
                "public Thomson endpoint without the missing hadronic spectral payload",
                "W/Z compare-only mass rows",
                "charged-lepton absolute masses",
                "source-only hadron masses",
                "strong CP",
                "SI G without the full no-G clock stack",
            ],
        },
        "artifacts": {
            "local_global_resonance": _rel(HIERARCHY_RESONANCE),
            "ew_capacity": _rel(HIERARCHY_EW_CAPACITY),
            "finite_readback": _rel(HIERARCHY_READBACK),
            "m_rep": _rel(HIERARCHY_M_REP),
            "screen_sieve": _rel(HIERARCHY_SCREEN_SIEVE),
            "higgs_naturality": _rel(HIERARCHY_NATURALITY),
        },
    }


def build_payload() -> dict[str, Any]:
    p_trunk = _load_json(P_TRUNK)
    measured_endpoint = _load_json(MEASURED_ENDPOINT)
    pipeline = _load_json(PIPELINE_STATUS)
    exact = _load_json(EXACT_NONHADRON)
    results = _load_json(RESULTS_STATUS)
    direct_top = _load_json(DIRECT_TOP)
    by_id = {entry["particle_id"]: _prediction_entry(entry) for entry in exact["entries"]}
    predictions = [by_id[particle_id] for particle_id in PARTICLE_ORDER if particle_id in by_id]
    particle_five_gates = [
        gate
        for gate in pipeline["issue_gates"]
        if gate["issue"] in {223, 224, 225, 234, 235, 32, 153, 199, 201, 207}
    ]

    return {
        "artifact": "oph_final_current_end_to_end_particle_predictions",
        "generated_utc": _now_utc(),
        "scope": "nonhadron_particle_pipeline_with_empirical_hadron_closure_policy",
        "claim_status": "final_nonhadron_predictions_with_separate_empirical_hadron_closure_surface",
        "source_surfaces": {
            "p_trunk": "code/P_derivation/runtime/p_closure_trunk_current.json",
            "measured_endpoint_calibration": (
                "code/P_derivation/runtime/measured_endpoint_calibration_current.json"
            ),
            "thomson_endpoint_package": "code/P_derivation/runtime/thomson_endpoint_package_current.json",
            "pipeline_status": "code/particles/runs/status/particle_pipeline_closure_status.json",
            "exact_nonhadron": "code/particles/exact_nonhadron_masses.json",
            "results_status": "code/particles/results_status.json",
            "direct_top_bridge": "code/particles/runs/calibration/direct_top_bridge_contract.json",
            "hadron_policy": "HADRON.md",
            "empirical_ee_hadrons_source_registry": (
                "code/particles/hadron/empirical_ee_hadrons_sources.yaml"
            ),
            "empirical_ee_hadronic_spectral_measure_schema": (
                "code/particles/hadron/empirical_ee_hadronic_spectral_measure.schema.json"
            ),
            "hierarchy_local_global_resonance": (
                "code/particles/hierarchy/certificates/R_local_global_hierarchy_resonance_closeout_335.json"
            ),
            "hierarchy_ew_capacity": (
                "code/particles/hierarchy/certificates/R_EW_global_capacity_certificate.json"
            ),
            "hierarchy_higgs_naturality": (
                "code/particles/hierarchy/issue_332_rg_naturality_certificate.json"
            ),
        },
        "p_closure": {
            "P": p_trunk["fixed_point_candidate"]["P"],
            "alpha_inv": p_trunk["fixed_point_candidate"]["alpha_inv"],
            "claim_status": p_trunk["claim_status"],
            "may_feed_live_particle_predictions": p_trunk["consumer_policy"]["may_feed_live_particle_predictions"],
        },
        "runtime_inputs": results.get("inputs", {}),
        "output_classes": [
            "source_only_oph",
            "oph_plus_empirical_hadron_closure",
            "compare_only",
            "work_in_progress",
        ],
        "fine_structure": _fine_structure_surface(measured_endpoint),
        "hierarchy_and_naturality": _hierarchy_surface(),
        "finalization_gates": pipeline["finalization_gates"],
        "particle_five_issue_gates": particle_five_gates,
        "companion_open_branches": list(pipeline.get("companion_status_branches", [])),
        "predictions": predictions,
        "hadron_policy": {
            "source_only_hadron_predictions_emitted": False,
            "empirical_hadron_closure_allowed_for_display": True,
            "policy_artifact": "HADRON.md",
            "source_registry": str(EMPIRICAL_EE_REGISTRY.relative_to(ROOT)),
            "empirical_payload_schema": str(EMPIRICAL_EE_SCHEMA.relative_to(ROOT)),
            "reason": (
                "Source-only hadron outputs require a working OPH hadron backend. Empirical "
                "hadron closure values stay in a separate output class; the e+e- spectral payload "
                "has a source registry and schema."
            ),
            "github_issues": [153, 157],
        },
        "direct_top_auxiliary_comparison": {
            "current_top_coordinate_gev": direct_top["current_theorem_coordinate"]["value_gev"],
            "current_top_codomain": direct_top["current_theorem_coordinate"]["pdg_summary_id"],
            "auxiliary_direct_top_gev": direct_top["auxiliary_direct_top_coordinate"]["value_gev"],
            "auxiliary_direct_top_codomain": direct_top["auxiliary_direct_top_coordinate"]["pdg_summary_id"],
            "direct_minus_current_coordinate_gev": direct_top["comparison_only_readout"]["direct_minus_current_coordinate_gev"],
            "pull_in_combined_sigma": direct_top["comparison_only_readout"]["pull_in_combined_sigma"],
            "bridge_status": direct_top["status"],
        },
    }


def render_markdown(payload: dict[str, Any]) -> str:
    lines = [
        "# Final End-to-End Particle Predictions",
        "",
        f"Generated: `{payload['generated_utc']}`",
        "",
        f"Scope: `{payload['scope']}`",
        f"Claim label: `{payload['claim_status']}`",
        "",
        "## P Closure",
        "",
        f"- Candidate `P`: `{payload['p_closure']['P']}`",
        f"- Candidate `alpha^-1`: `{payload['p_closure']['alpha_inv']}`",
        f"- Claim label: `{payload['p_closure']['claim_status']}`",
        f"- May feed promoted particle predictions: `{payload['p_closure']['may_feed_live_particle_predictions']}`",
        "",
        "## Particle-Five Receipts",
        "",
        "| Receipt label | Closable | Local artifact | Worker policy |",
        "| --- | --- | --- | --- |",
    ]
    for gate in payload["particle_five_issue_gates"]:
        lines.append(
            f"| `{_display_status(gate['state'])}` | `{gate['closable_now']}` | "
            f"`{gate['local_next_artifact']}` | {gate['chrome_workers']} |"
        )
    companion_open_branches = payload.get("companion_open_branches") or []
    if companion_open_branches:
        lines.extend(
            [
                "",
                "## Companion Claim Boundaries",
                "",
                "| Topic | Claim label | Boundary | Gate |",
                "| --- | --- | --- | --- |",
            ]
        )
        for branch in companion_open_branches:
            lines.append(
                f"| {branch['label']} | `{_display_status(branch['state'])}` | {branch['summary']} | {branch['next_action']} |"
            )
    lines.extend(
        [
            "",
            "## Predictions",
            "",
            "| Particle | Prediction | Claim label | Scope | Promotable |",
            "| --- | ---: | --- | --- | --- |",
        ]
    )
    for entry in payload["predictions"]:
        lines.append(
            f"| `{entry['particle_id']}` | `{entry['value']} {entry['unit']}` | "
            f"`{entry['exact_kind']}` | `{entry['scope']}` | `{entry['promotable']}` |"
        )
    direct = payload["direct_top_auxiliary_comparison"]
    hierarchy = payload["hierarchy_and_naturality"]
    hierarchy_status = hierarchy["status"]
    hierarchy_bridge = hierarchy["local_global_bridge"]
    hierarchy_factors = hierarchy["factor_origins"]
    lines.extend(
        [
            "",
            "## Fine Structure",
            "",
            "| Output class | alpha^-1(0) | P | Claim label |",
            "| --- | ---: | ---: | --- |",
        ]
    )
    fine = payload["fine_structure"]
    source_only = fine["source_only_oph"]
    empirical = fine["oph_plus_empirical_hadron_closure"]
    lines.append(
        f"| `source_only_oph` | `{source_only['alpha_inv']}` | `{source_only['P']}` | "
        f"`{source_only['status']}` |"
    )
    lines.append(
        f"| `oph_plus_empirical_hadron_closure` | `{empirical['alpha_inv']}` | "
        f"`{empirical['P']}` | `{empirical['status']}` |"
    )
    lines.extend(
        [
            "",
            f"- Empirical residual at the public endpoint pixel: "
            f"`{empirical['missing_source_transport_delta_alpha_inv']}` inverse-alpha units",
            f"- Empirical payload policy: `{fine['empirical_payload_policy']['dispersion_payload_status']}`",
            "",
            "## Hierarchy And Naturality",
            "",
            f"- Resonance label: `{_display_status(hierarchy_status['resonance_status'])}`",
            f"- Full theorem-grade resonance promoted: `{hierarchy_status['full_theorem_grade_resonance_promoted']}`",
            f"- Remaining promotion gates: `{hierarchy_status['remaining_promotion_gates']}`",
            f"- Exact EW bridge capacity: `{hierarchy_bridge['N_CRC_EW']}`",
            f"- Bridge residual: `{hierarchy_bridge['bridge_residual']}`",
            f"- Source `v/E_cell`: `{hierarchy_bridge['v_over_E_cell_source']}`",
            f"- Factor origins: `ports={hierarchy_factors['screen_ports']}`, "
            f"`m_rep={hierarchy_factors['m_rep']}`, "
            f"`exponent={hierarchy_factors['specialized_exponent']}`",
            f"- Higgs naturality defect: `epsilon_H={hierarchy_factors['higgs_naturality_defect']}`",
            f"- Boundary: {hierarchy['claim_boundary']['improves']}",
            f"- Not promoted by this bridge: {', '.join(hierarchy['claim_boundary']['does_not_promote'])}",
            "",
            "## Direct-Top Auxiliary Comparison",
            "",
            f"- Top theorem coordinate: `{direct['current_top_coordinate_gev']} GeV` on `{direct['current_top_codomain']}`",
            f"- Auxiliary direct-top coordinate: `{direct['auxiliary_direct_top_gev']} GeV` on `{direct['auxiliary_direct_top_codomain']}`",
            f"- Difference: `{direct['direct_minus_current_coordinate_gev']} GeV`",
            f"- Pull: `{direct['pull_in_combined_sigma']}` combined sigma",
            f"- Bridge label: `{_display_status(direct['bridge_status'])}`",
            "",
            "## Hadrons",
            "",
            f"- Source-only hadron predictions emitted: `{payload['hadron_policy']['source_only_hadron_predictions_emitted']}`",
            f"- Empirical hadron closure allowed for display: `{payload['hadron_policy']['empirical_hadron_closure_allowed_for_display']}`",
            f"- Policy artifact: `{payload['hadron_policy']['policy_artifact']}`",
            f"- Reason: {payload['hadron_policy']['reason']}",
        ]
    )
    return "\n".join(lines)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build final current end-to-end particle predictions.")
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
