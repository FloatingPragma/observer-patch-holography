#!/usr/bin/env python3
"""Build the public non-hadron mass-output bundle.

The public `entries` list is deliberately stricter than the exact-fit audit
surface: target-anchored witnesses and compare-only absolute attachments are
withheld from prediction tables. They remain available in `EXACT_FITS_ONLY`
and in the generated `withheld_entries` metadata for audit/debug use.
"""

from __future__ import annotations

import argparse
import json
import pathlib
from datetime import datetime, timezone
from typing import Any


ROOT = pathlib.Path(__file__).resolve().parents[2]
REFERENCE_JSON = ROOT / "particles" / "data" / "particle_reference_values.json"
EW_EXACT_JSON = ROOT / "particles" / "runs" / "calibration" / "d10_ew_w_anchor_neutral_shear_factorization_official_pdg_2025_update.json"
D11_EXACT_JSON = ROOT / "particles" / "runs" / "calibration" / "d11_reference_exact_adapter.json"
D11_EXACT_HIGGS_PROMOTION_JSON = ROOT / "particles" / "runs" / "calibration" / "d11_live_exact_higgs_promotion.json"
D11_EXACT_SPLIT_PAIR_JSON = ROOT / "particles" / "runs" / "calibration" / "d11_live_exact_split_pair_theorem.json"
CHARGED_JSON = ROOT / "particles" / "runs" / "leptons" / "lepton_current_family_exact_readout.json"
QUARK_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_exact_readout.json"
CHARGED_THEOREM_JSON = ROOT / "particles" / "runs" / "leptons" / "lepton_current_family_quadratic_readout_theorem.json"
QUARK_THEOREM_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_quadratic_readout_theorem.json"
CHARGED_AFFINE_JSON = ROOT / "particles" / "runs" / "leptons" / "lepton_current_family_affine_anchor_theorem.json"
CHARGED_TRACE_LIFT_REQUIRED_JSON = (
    ROOT / "particles" / "runs" / "leptons" / "charged_determinant_trace_lift_attachment_required.json"
)
QUARK_CLOSURE_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_selected_sheet_closure.json"
QUARK_TRANSPORT_LIFT_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_transport_frame_sector_attached_lift.json"
QUARK_TRANSPORT_COMPLETION_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_transport_frame_exact_pdg_completion.json"
QUARK_TRANSPORT_FORWARD_YUKAWAS_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_transport_frame_exact_forward_yukawas.json"
QUARK_END_TO_END_CHAIN_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_end_to_end_exact_pdg_derivation_chain.json"
QUARK_PUBLIC_YUKAWA_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_public_exact_yukawa_end_to_end_theorem.json"
QUARK_SIGMA_REQUIRED_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_sigma_source_datum_no_target_leak_required.json"
NEUTRINO_JSON = ROOT / "particles" / "runs" / "neutrino" / "neutrino_absolute_attachment_theorem.json"
NEUTRINO_BRIDGE_RIGIDITY_JSON = ROOT / "particles" / "runs" / "neutrino" / "neutrino_bridge_rigidity_theorem.json"
CARRIER_ACCEPTANCE_JSON = ROOT / "particles" / "runs" / "status" / "carrier_mode_acceptance.json"
DEFAULT_MD_OUT = ROOT / "particles" / "EXACT_NONHADRON_MASSES.md"
DEFAULT_JSON_OUT = ROOT / "particles" / "exact_nonhadron_masses.json"
DEFAULT_FORWARD_OUT = ROOT / "particles" / "runs" / "status" / "exact_nonhadron_masses_current.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: pathlib.Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _load_optional_json(path: pathlib.Path) -> dict[str, Any] | None:
    if not path.exists():
        return None
    return _load_json(path)


def _repo_ref(path: pathlib.Path) -> str:
    return str(path.relative_to(ROOT.parent))


def _carrier_summary(entry: dict[str, Any]) -> dict[str, Any]:
    return {
        "carrier_id": entry["carrier_id"],
        "label": entry["label"],
        "claim_kind": entry["claim_kind"],
        "branch": entry["branch"],
        "hard_quadratic_mass_parameter_squared": entry[
            "hard_quadratic_mass_parameter_squared"
        ],
        "classical_carrier_gate": entry["classical_carrier_gate"],
        "quantum_particle_gate": entry["quantum_particle_gate"],
        "particle_promotion_allowed": entry["particle_promotion_allowed"],
        "source_artifact": "code/particles/runs/status/carrier_mode_acceptance.json",
    }


def _non_circularity_promotable(payload: dict[str, Any] | None, default: bool = False) -> bool:
    if not payload:
        return False
    status = payload.get("non_circularity_status")
    if isinstance(status, dict):
        return status.get("promotion_allowed") is True
    return default


def _d11_split_promotable(payload: dict[str, Any] | None) -> bool:
    return (
        bool(payload)
        and payload.get("status") == "closed"
        and payload.get("proof_status") == "closed_source_only_live_exact_split_pair"
        and payload.get("public_surface_candidate_allowed") is True
        and _non_circularity_promotable(payload, default=True)
    )


def _quark_public_promotable(payload: dict[str, Any] | None) -> bool:
    return (
        bool(payload)
        and payload.get("public_promotion_allowed") is True
        and payload.get("proof_status") == "closed_source_only_public_exact_yukawa_end_to_end_theorem"
        and _non_circularity_promotable(payload, default=True)
    )


def _neutrino_absolute_promotable(payload: dict[str, Any] | None) -> bool:
    return (
        bool(payload)
        and payload.get("artifact") == "oph_neutrino_absolute_attachment_theorem"
        and payload.get("status") == "theorem_grade_emitted"
        and payload.get("weighted_cycle_base_eligible") is True
        and payload.get("prediction_promotion_allowed") is True
        and payload.get("public_surface_candidate_allowed") is True
        and _non_circularity_promotable(payload, default=False)
    )


def _is_public_mass_output(entry: dict[str, Any]) -> bool:
    if entry.get("particle_id") in {"photon", "gluon", "graviton"}:
        # A zero hard parameter in a declared quadratic action is not a public
        # quantum-particle rest-mass prediction.  Those modes live in the
        # separate carrier acceptance receipt.
        return False
    if entry.get("legacy_particle_id_slot") is True:
        # Oscillation eigenmasses are not masses of flavor states. A future
        # physical theorem must use dedicated mass-eigenstate identifiers.
        return False
    exact_kind = str(entry.get("exact_kind", ""))
    if "target_anchored" in exact_kind:
        return False
    if "compare_only" in exact_kind:
        return False
    if "target_informed" in exact_kind or "rejected" in exact_kind:
        return False
    return True


def _withheld_entry(entry: dict[str, Any]) -> dict[str, Any]:
    exact_kind = str(entry.get("exact_kind", ""))
    if "target_anchored" in exact_kind:
        reason = "target_anchored_witness_kept_in_exact_fit_audit_not_public_prediction"
    elif "compare_only" in exact_kind:
        reason = "compare_only_absolute_or_adapter_surface_kept_out_of_public_prediction_table"
    elif "target_informed" in exact_kind or "rejected" in exact_kind:
        reason = "target_informed_candidate_rejected_by_correlated_profile"
    else:
        reason = "not_public_prediction_output"
    withheld = {
        "particle_id": entry["particle_id"],
        "label": entry.get("label"),
        "exact_kind": exact_kind,
        "scope": entry.get("scope"),
        "promotable": bool(entry.get("promotable", False)),
        "source_artifact": entry.get("source_artifact"),
        "reason": reason,
    }
    for key in (
        "public_theorem_value",
        "source_only",
        "centered_log",
        "formula_if_anchor_exists",
        "conditional_anchor_symbol",
        "missing_for_promotion",
        "promotion_gate_artifact",
        "claim_tier",
        "source_only_sigma_emitted",
        "downstream_algebra_closed",
        "selected_class",
        "exact_sigma_target",
        "strongest_source_candidate",
        "legacy_particle_id_slot",
        "mass_basis_semantics",
    ):
        if key in entry:
            withheld[key] = entry[key]
    return withheld


CHARGED_MISSING_FOR_PROMOTION = [
    "charged_branch_generator_splitting_promotion",
    "sector_isolated_charged_determinant_exponent_vector_M_ch",
    "source_side_same_label_q_psi_readout_certificate",
    "charged_determinant_trace_lift_attachment",
    "NO_TARGET_LEAK_DAG_CHARGED_A_CH",
]


def _charged_formula(symbol: str, centered_log: float) -> str:
    sign = "+" if centered_log >= 0.0 else "-"
    return f"{symbol}(P)=exp(A_ch(P){sign}{abs(centered_log):.15f})"


def build_all_entries() -> list[dict[str, Any]]:
    references = _load_json(REFERENCE_JSON)["entries"]
    ew_exact = _load_json(EW_EXACT_JSON)
    d11_exact = _load_json(D11_EXACT_JSON)
    d11_exact_higgs = _load_optional_json(D11_EXACT_HIGGS_PROMOTION_JSON)
    d11_exact_pair = _load_optional_json(D11_EXACT_SPLIT_PAIR_JSON)
    charged = _load_json(CHARGED_JSON)
    quark = _load_json(QUARK_JSON)
    charged_affine = _load_optional_json(CHARGED_AFFINE_JSON)
    quark_closure = _load_optional_json(QUARK_CLOSURE_JSON)
    quark_transport_lift = _load_optional_json(QUARK_TRANSPORT_LIFT_JSON)
    quark_transport_completion = _load_optional_json(QUARK_TRANSPORT_COMPLETION_JSON)
    quark_transport_forward_yukawas = _load_optional_json(QUARK_TRANSPORT_FORWARD_YUKAWAS_JSON)
    quark_end_to_end_chain = _load_optional_json(QUARK_END_TO_END_CHAIN_JSON)
    quark_public_exact_yukawa = _load_optional_json(QUARK_PUBLIC_YUKAWA_JSON)
    quark_sigma_required = _load_optional_json(QUARK_SIGMA_REQUIRED_JSON)
    neutrino = _load_json(NEUTRINO_JSON)
    neutrino_bridge_rigidity = _load_json(NEUTRINO_BRIDGE_RIGIDITY_JSON)
    quark_exact_values = (
        dict(quark_public_exact_yukawa["public_exact_outputs"]["exact_running_values_gev"])
        if quark_public_exact_yukawa
        else {
            "u": quark["predicted_singular_values_u"][0],
            "c": quark["predicted_singular_values_u"][1],
            "t": quark["predicted_singular_values_u"][2],
            "d": quark["predicted_singular_values_d"][0],
            "s": quark["predicted_singular_values_d"][1],
            "b": quark["predicted_singular_values_d"][2],
        }
    )
    quark_exact_promotable = _quark_public_promotable(quark_public_exact_yukawa)
    quark_exact_kind = (
        "selected_class_theorem_grade_exact_forward_quark_closure"
        if quark_exact_promotable
        else (
            "selected_class_target_anchored_exact_witness"
            if quark_public_exact_yukawa
            else "exact_target_anchored_current_family_witness"
        )
    )
    quark_exact_scope = (
        (
            "selected_public_physical_quark_frame_class_only_but_sigma_datum_target_derived"
            if quark_public_exact_yukawa and not quark_exact_promotable
            else quark_public_exact_yukawa["theorem_scope"]
        )
        if quark_public_exact_yukawa
        else quark["theorem_scope"]
    )
    quark_exact_source = _repo_ref(QUARK_PUBLIC_YUKAWA_JSON) if quark_public_exact_yukawa else _repo_ref(QUARK_JSON)
    quark_missing_for_promotion = list(
        (quark_sigma_required or {}).get(
            "missing_for_promotion",
            [
                "QUARK_SIGMA_SOURCE_QUOTIENT",
                "QUARK_SIGMA_SOURCE_SELECTOR",
                "QUARK_EDGE_STATISTICS_CORRECTION_THEOREM",
                "QUARK_SIGMA_REFINEMENT_COMPATIBILITY",
                "NO_TARGET_LEAK_DAG_QUARK_SIGMA_SOURCE",
            ],
        )
    )
    quark_common_metadata = {
        "claim_tier": (
            (quark_public_exact_yukawa or {}).get("claim_tier")
            or (quark_sigma_required or {}).get("claim_tier")
            or "selected_class_conditional_on_source_sigma"
        ),
        "source_only_sigma_emitted": bool(
            (quark_public_exact_yukawa or {}).get(
                "source_only_sigma_emitted",
                (quark_sigma_required or {}).get("source_only_sigma_emitted", False),
            )
        ),
        "downstream_algebra_closed": bool(
            (quark_public_exact_yukawa or {}).get(
                "downstream_algebra_closed",
                (quark_sigma_required or {}).get("downstream_algebra_closed", True),
            )
        ),
        "selected_class": (quark_sigma_required or {}).get("selected_class", "f_P"),
        "missing_for_promotion": quark_missing_for_promotion,
        "promotion_gate_artifact": _repo_ref(QUARK_SIGMA_REQUIRED_JSON),
        "exact_sigma_target": (quark_sigma_required or {}).get("target_values_for_future_source_theorem"),
        "strongest_source_candidate": (quark_sigma_required or {}).get("strongest_current_source_candidate"),
    }
    quark_public_note_prefix = (
        "Theorem-grade exact quark output on the selected public physical quark frame class chosen by `P`, "
        if quark_exact_promotable
        else "Selected-class exact quark witness on the public physical quark frame class chosen by `P`, "
    )
    quark_public_nonpromotion_note = (
        ""
        if quark_exact_promotable
        else "Under the strict non-circularity audit it is not promotable because the public sigma descent proves representative independence only; the source-only sigma selector is still missing. "
    )
    quark_exact_note = (
        (
            quark_public_note_prefix
            +
            f"emitted by `{quark_public_exact_yukawa['artifact']}`. The exact sextet matches the official PDG 2025 "
            "API running-quark target surface, and the theorem emits explicit exact forward Yukawas `Y_u` and `Y_d`. "
            + quark_public_nonpromotion_note
        )
        if quark_public_exact_yukawa
        else (
            "Exact current-family quark witness on `current_family_only`. The exact sextet matches the official "
            "PDG 2025 API running-quark target surface on that declared carrier. "
        )
    ) + (
        "The top coordinate uses PDG "
        f"summary `{references['top_quark']['source']['summary_id']}`. The auxiliary direct-top "
        f"entry `{references['top_quark_direct_aux']['source']['summary_id']}` is compare-only; "
        "[#207](https://github.com/FloatingPragma/observer-patch-holography/issues/207) is closed as a corpus-limited no-go by "
        "`code/particles/runs/calibration/direct_top_bridge_contract.json`. The same exact sextet is also "
        "realized on `current_family_only`. A separate restricted theorem "
        "surface emits a sector-attached `Sigma_ud^phys` element on "
        f"`{(quark_transport_lift or {}).get('theorem_scope', 'current_family_common_refinement_transport_frame_only')}`, "
        "the merged transport-frame completion closes the same "
        f"running sextet on `{(quark_transport_completion or {}).get('theorem_scope', 'current_family_common_refinement_transport_frame_only')}`, and the declared transport-frame "
        "chain emits explicit exact forward Yukawas `Y_u` and `Y_d` with certification status "
        f"`{(quark_transport_forward_yukawas or {}).get('certification_status', 'forward_matrix_certified')}`. The full declared-carrier chain is "
        f"recorded in `{(quark_end_to_end_chain or {}).get('artifact', 'oph_quark_current_family_end_to_end_exact_pdg_derivation_chain')}`. The target-free mass bridge closes separately "
        "on the emitted D12 ray, but it does not emit the physical sigma/spread datum. "
        + (
            "This theorem is selected-class closure only. It does not claim a global classification of all quark frame classes."
            if quark_public_exact_yukawa
            else "This entry is carrier-restricted and does not replace the selected-class public theorem."
        )
    )
    d11_split_promotable = _d11_split_promotable(d11_exact_pair)
    higgs_promotable = d11_split_promotable
    if d11_exact_pair:
        higgs_exact_kind = (
            "exact_source_only_higgs_top_split_calibration_theorem"
            if d11_split_promotable
            else "conditional_declared_surface_higgs_top_candidate"
        )
    elif d11_exact_higgs:
        higgs_exact_kind = "exact_target_anchored_higgs_calibration_theorem"
    else:
        higgs_exact_kind = "exact_target_anchored_compare_only_inverse_slice"
    neutrino_promotable = _neutrino_absolute_promotable(neutrino)
    neutrino_exact_kind = (
        "theorem_grade_source_closed_physical_neutrino_eigenstate"
        if neutrino_promotable
        else "rejected_target_informed_weighted_cycle_candidate"
    )
    neutrino_note_prefix = (
        "Source-closed physical neutrino mass eigenvalue with an independent operator, basis, label, and scale contract, "
        if neutrino_promotable
        else "Rejected target-informed weighted-cycle candidate mass eigenvalue, not a flavor-neutrino mass, with compare-only absolute attachment, "
    )
    c_nu_display = neutrino_bridge_rigidity.get("emitted_value")
    if c_nu_display is None:
        c_nu_display = neutrino_bridge_rigidity.get("display_value")

    return [
        {
            "particle_id": "higgs",
            "label": "Higgs Boson",
            "mass_gev": (
                d11_exact_pair["exact_split_pair"]["mH_gev"]
                if d11_exact_pair
                else (
                d11_exact_higgs["mass_readout"]["mH_gev"]
                if d11_exact_higgs
                else d11_exact["predicted_outputs"]["mH_gev"]
                )
            ),
            "exact_kind": higgs_exact_kind,
            "scope": (
                d11_exact_pair["theorem_scope"]
                if d11_exact_pair
                else (
                d11_exact_higgs["theorem_scope"]
                if d11_exact_higgs
                else d11_exact["scope"]
                )
            ),
            "promotable": higgs_promotable,
            "source_artifact": (
                _repo_ref(D11_EXACT_SPLIT_PAIR_JSON)
                if d11_exact_pair
                else (
                _repo_ref(D11_EXACT_HIGGS_PROMOTION_JSON)
                if d11_exact_higgs
                else _repo_ref(D11_EXACT_JSON)
                )
            ),
            "note": (
                (
                "Exact source-only Higgs theorem on the declared D10/D11 surface. "
                "The same theorem emits the exact Higgs row together with a companion top coordinate by direct Jacobian readout from the live D10 repair tuple. "
                "At PDG quoting precision the Higgs row lands on the 2025 Higgs average. "
                "The target-anchored selected-class running-top witness uses the PDG 2025 cross-section entry `Q007TP4`, but it is audit-only under the strict public-output policy. "
                "The auxiliary direct-top average `Q007TP` is compare-only; "
                "[#207](https://github.com/FloatingPragma/observer-patch-holography/issues/207) is closed as a corpus-limited no-go by "
                "`code/particles/runs/calibration/direct_top_bridge_contract.json`."
                if d11_split_promotable
                else
                "Conditional Higgs/top split candidate on the declared D10/D11 surface. It emits the same numeric split pair, but strict promotion is blocked until the upstream D10 target-free repair law is closed and promotable."
                )
                if d11_exact_pair
                else (
                "Exact target-anchored Higgs calibration theorem on the declared D10/D11 surface."
                if d11_exact_higgs
                else "Exact compare-only inverse slice on the D11 Jacobian."
                )
            ),
        },
        {
            "particle_id": "electron",
            "label": "Electron",
            "mass_gev": charged["predicted_singular_values_abs"][0],
            "exact_kind": "exact_target_anchored_current_family_witness",
            "scope": charged["theorem_scope"],
            "promotable": False,
            "source_only": False,
            "public_theorem_value": None,
            "centered_log": charged["centered_log_shape_exact"][0],
            "formula_if_anchor_exists": _charged_formula("m_e", float(charged["centered_log_shape_exact"][0])),
            "conditional_anchor_symbol": "A_ch(P)",
            "missing_for_promotion": CHARGED_MISSING_FOR_PROMOTION,
            "promotion_gate_artifact": _repo_ref(CHARGED_TRACE_LIFT_REQUIRED_JSON),
            "source_artifact": _repo_ref(CHARGED_JSON),
            "supporting_theorem_artifact": _repo_ref(CHARGED_THEOREM_JSON),
            "supporting_scope_closure_artifact": _repo_ref(CHARGED_AFFINE_JSON),
            "note": "Exact `current_family_only` charged-lepton witness. Public promotion requires a source-only charged determinant trace-lift attachment emitting A_ch(P); the same-family A_ch_current_family checksum is forbidden ancestry for that theorem.",
        },
        {
            "particle_id": "muon",
            "label": "Muon",
            "mass_gev": charged["predicted_singular_values_abs"][1],
            "exact_kind": "exact_target_anchored_current_family_witness",
            "scope": charged["theorem_scope"],
            "promotable": False,
            "source_only": False,
            "public_theorem_value": None,
            "centered_log": charged["centered_log_shape_exact"][1],
            "formula_if_anchor_exists": _charged_formula("m_mu", float(charged["centered_log_shape_exact"][1])),
            "conditional_anchor_symbol": "A_ch(P)",
            "missing_for_promotion": CHARGED_MISSING_FOR_PROMOTION,
            "promotion_gate_artifact": _repo_ref(CHARGED_TRACE_LIFT_REQUIRED_JSON),
            "source_artifact": _repo_ref(CHARGED_JSON),
            "supporting_theorem_artifact": _repo_ref(CHARGED_THEOREM_JSON),
            "supporting_scope_closure_artifact": _repo_ref(CHARGED_AFFINE_JSON),
            "note": "Exact `current_family_only` charged-lepton witness. Public promotion requires a source-only charged determinant trace-lift attachment emitting A_ch(P); the same-family A_ch_current_family checksum is forbidden ancestry for that theorem.",
        },
        {
            "particle_id": "tau",
            "label": "Tau",
            "mass_gev": charged["predicted_singular_values_abs"][2],
            "exact_kind": "exact_target_anchored_current_family_witness",
            "scope": charged["theorem_scope"],
            "promotable": False,
            "source_only": False,
            "public_theorem_value": None,
            "centered_log": charged["centered_log_shape_exact"][2],
            "formula_if_anchor_exists": _charged_formula("m_tau", float(charged["centered_log_shape_exact"][2])),
            "conditional_anchor_symbol": "A_ch(P)",
            "missing_for_promotion": CHARGED_MISSING_FOR_PROMOTION,
            "promotion_gate_artifact": _repo_ref(CHARGED_TRACE_LIFT_REQUIRED_JSON),
            "source_artifact": _repo_ref(CHARGED_JSON),
            "supporting_theorem_artifact": _repo_ref(CHARGED_THEOREM_JSON),
            "supporting_scope_closure_artifact": _repo_ref(CHARGED_AFFINE_JSON),
            "note": "Exact `current_family_only` charged-lepton witness. Public promotion requires a source-only charged determinant trace-lift attachment emitting A_ch(P); the same-family A_ch_current_family checksum is forbidden ancestry for that theorem.",
        },
        {
            "particle_id": "up_quark",
            "label": "Up Quark",
            "mass_gev": quark_exact_values["u"],
            "exact_kind": quark_exact_kind,
            "scope": quark_exact_scope,
            "promotable": quark_exact_promotable,
            "source_artifact": quark_exact_source,
            "supporting_theorem_artifact": _repo_ref(QUARK_THEOREM_JSON),
            "supporting_scope_closure_artifact": _repo_ref(QUARK_CLOSURE_JSON),
            "supporting_transport_frame_artifact": _repo_ref(QUARK_TRANSPORT_LIFT_JSON),
            "supporting_transport_frame_completion_artifact": _repo_ref(QUARK_TRANSPORT_COMPLETION_JSON),
            "supporting_transport_frame_forward_yukawas_artifact": _repo_ref(QUARK_TRANSPORT_FORWARD_YUKAWAS_JSON),
            "supporting_end_to_end_chain_artifact": _repo_ref(QUARK_END_TO_END_CHAIN_JSON),
            **quark_common_metadata,
            "note": quark_exact_note,
        },
        {
            "particle_id": "charm_quark",
            "label": "Charm Quark",
            "mass_gev": quark_exact_values["c"],
            "exact_kind": quark_exact_kind,
            "scope": quark_exact_scope,
            "promotable": quark_exact_promotable,
            "source_artifact": quark_exact_source,
            "supporting_theorem_artifact": _repo_ref(QUARK_THEOREM_JSON),
            "supporting_scope_closure_artifact": _repo_ref(QUARK_CLOSURE_JSON),
            "supporting_transport_frame_artifact": _repo_ref(QUARK_TRANSPORT_LIFT_JSON),
            "supporting_transport_frame_completion_artifact": _repo_ref(QUARK_TRANSPORT_COMPLETION_JSON),
            "supporting_transport_frame_forward_yukawas_artifact": _repo_ref(QUARK_TRANSPORT_FORWARD_YUKAWAS_JSON),
            "supporting_end_to_end_chain_artifact": _repo_ref(QUARK_END_TO_END_CHAIN_JSON),
            **quark_common_metadata,
            "note": quark_exact_note,
        },
        {
            "particle_id": "top_quark",
            "label": "Top Quark",
            "mass_gev": quark_exact_values["t"],
            "exact_kind": quark_exact_kind,
            "scope": quark_exact_scope,
            "promotable": quark_exact_promotable,
            "source_artifact": quark_exact_source,
            "supporting_theorem_artifact": _repo_ref(QUARK_THEOREM_JSON),
            "supporting_scope_closure_artifact": _repo_ref(QUARK_CLOSURE_JSON),
            "supporting_transport_frame_artifact": _repo_ref(QUARK_TRANSPORT_LIFT_JSON),
            "supporting_transport_frame_completion_artifact": _repo_ref(QUARK_TRANSPORT_COMPLETION_JSON),
            "supporting_transport_frame_forward_yukawas_artifact": _repo_ref(QUARK_TRANSPORT_FORWARD_YUKAWAS_JSON),
            "supporting_end_to_end_chain_artifact": _repo_ref(QUARK_END_TO_END_CHAIN_JSON),
            **quark_common_metadata,
            "note": quark_exact_note,
        },
        {
            "particle_id": "down_quark",
            "label": "Down Quark",
            "mass_gev": quark_exact_values["d"],
            "exact_kind": quark_exact_kind,
            "scope": quark_exact_scope,
            "promotable": quark_exact_promotable,
            "source_artifact": quark_exact_source,
            "supporting_theorem_artifact": _repo_ref(QUARK_THEOREM_JSON),
            "supporting_scope_closure_artifact": _repo_ref(QUARK_CLOSURE_JSON),
            "supporting_transport_frame_artifact": _repo_ref(QUARK_TRANSPORT_LIFT_JSON),
            "supporting_transport_frame_completion_artifact": _repo_ref(QUARK_TRANSPORT_COMPLETION_JSON),
            "supporting_transport_frame_forward_yukawas_artifact": _repo_ref(QUARK_TRANSPORT_FORWARD_YUKAWAS_JSON),
            "supporting_end_to_end_chain_artifact": _repo_ref(QUARK_END_TO_END_CHAIN_JSON),
            **quark_common_metadata,
            "note": quark_exact_note,
        },
        {
            "particle_id": "strange_quark",
            "label": "Strange Quark",
            "mass_gev": quark_exact_values["s"],
            "exact_kind": quark_exact_kind,
            "scope": quark_exact_scope,
            "promotable": quark_exact_promotable,
            "source_artifact": quark_exact_source,
            "supporting_theorem_artifact": _repo_ref(QUARK_THEOREM_JSON),
            "supporting_scope_closure_artifact": _repo_ref(QUARK_CLOSURE_JSON),
            "supporting_transport_frame_artifact": _repo_ref(QUARK_TRANSPORT_LIFT_JSON),
            "supporting_transport_frame_completion_artifact": _repo_ref(QUARK_TRANSPORT_COMPLETION_JSON),
            "supporting_transport_frame_forward_yukawas_artifact": _repo_ref(QUARK_TRANSPORT_FORWARD_YUKAWAS_JSON),
            "supporting_end_to_end_chain_artifact": _repo_ref(QUARK_END_TO_END_CHAIN_JSON),
            **quark_common_metadata,
            "note": quark_exact_note,
        },
        {
            "particle_id": "bottom_quark",
            "label": "Bottom Quark",
            "mass_gev": quark_exact_values["b"],
            "exact_kind": quark_exact_kind,
            "scope": quark_exact_scope,
            "promotable": quark_exact_promotable,
            "source_artifact": quark_exact_source,
            "supporting_theorem_artifact": _repo_ref(QUARK_THEOREM_JSON),
            "supporting_scope_closure_artifact": _repo_ref(QUARK_CLOSURE_JSON),
            "supporting_transport_frame_artifact": _repo_ref(QUARK_TRANSPORT_LIFT_JSON),
            "supporting_transport_frame_completion_artifact": _repo_ref(QUARK_TRANSPORT_COMPLETION_JSON),
            "supporting_transport_frame_forward_yukawas_artifact": _repo_ref(QUARK_TRANSPORT_FORWARD_YUKAWAS_JSON),
            "supporting_end_to_end_chain_artifact": _repo_ref(QUARK_END_TO_END_CHAIN_JSON),
            **quark_common_metadata,
            "note": quark_exact_note,
        },
        {
            "particle_id": "electron_neutrino",
            "label": "Neutrino Mass Eigenstate s0 (Legacy nu_e Slot)",
            "mass_eV": neutrino["outputs"]["masses_eV"][0],
            "exact_kind": neutrino_exact_kind,
            "scope": (
                "source_closed_weighted_cycle_mass_eigenstate"
                if neutrino_promotable
                else "rejected_weighted_cycle_mass_eigenstate_candidate"
            ),
            "promotable": neutrino_promotable,
            "legacy_particle_id_slot": True,
            "mass_basis_semantics": "ascending_candidate_mass_eigenstate_s0_not_electron_flavor_mass",
            "source_artifact": _repo_ref(NEUTRINO_JSON),
            "supporting_bridge_rigidity_artifact": _repo_ref(NEUTRINO_BRIDGE_RIGIDITY_JSON),
            "note": (
                neutrino_note_prefix
                +
                f"with `C_nu = {float(c_nu_display):.16f}`, "
                f"`P_nu = {neutrino_bridge_rigidity['emitted_proxy']['value']:.15f}`, and "
                f"`B_nu = {neutrino['outputs']['B_nu']:.15f}`."
            ),
        },
        {
            "particle_id": "muon_neutrino",
            "label": "Neutrino Mass Eigenstate s1 (Legacy nu_mu Slot)",
            "mass_eV": neutrino["outputs"]["masses_eV"][1],
            "exact_kind": neutrino_exact_kind,
            "scope": (
                "source_closed_weighted_cycle_mass_eigenstate"
                if neutrino_promotable
                else "rejected_weighted_cycle_mass_eigenstate_candidate"
            ),
            "promotable": neutrino_promotable,
            "legacy_particle_id_slot": True,
            "mass_basis_semantics": "ascending_candidate_mass_eigenstate_s1_not_muon_flavor_mass",
            "source_artifact": _repo_ref(NEUTRINO_JSON),
            "supporting_bridge_rigidity_artifact": _repo_ref(NEUTRINO_BRIDGE_RIGIDITY_JSON),
            "note": (
                neutrino_note_prefix
                +
                f"with `C_nu = {float(c_nu_display):.16f}`, "
                f"`P_nu = {neutrino_bridge_rigidity['emitted_proxy']['value']:.15f}`, and "
                f"`B_nu = {neutrino['outputs']['B_nu']:.15f}`."
            ),
        },
        {
            "particle_id": "tau_neutrino",
            "label": "Neutrino Mass Eigenstate s2 (Legacy nu_tau Slot)",
            "mass_eV": neutrino["outputs"]["masses_eV"][2],
            "exact_kind": neutrino_exact_kind,
            "scope": (
                "source_closed_weighted_cycle_mass_eigenstate"
                if neutrino_promotable
                else "rejected_weighted_cycle_mass_eigenstate_candidate"
            ),
            "promotable": neutrino_promotable,
            "legacy_particle_id_slot": True,
            "mass_basis_semantics": "ascending_candidate_mass_eigenstate_s2_not_tau_flavor_mass",
            "source_artifact": _repo_ref(NEUTRINO_JSON),
            "supporting_bridge_rigidity_artifact": _repo_ref(NEUTRINO_BRIDGE_RIGIDITY_JSON),
            "note": (
                neutrino_note_prefix
                +
                f"with `C_nu = {float(c_nu_display):.16f}`, "
                f"`P_nu = {neutrino_bridge_rigidity['emitted_proxy']['value']:.15f}`, and "
                f"`B_nu = {neutrino['outputs']['B_nu']:.15f}`."
            ),
        },
    ]


def build_entries() -> list[dict[str, Any]]:
    return [entry for entry in build_all_entries() if _is_public_mass_output(entry)]


def build_markdown(
    generated_utc: str,
    entries: list[dict[str, Any]],
    carrier_modes: list[dict[str, Any]],
) -> str:
    lines = [
        "# Public Non-Hadron Mass Outputs",
        "",
        f"Generated: `{generated_utc}`",
        "",
        "This bundle gives numeric non-hadron mass outputs that are not target-anchored witness rows.",
        "Target-anchored exact fits and compare-only absolute attachments are withheld from this prediction table and kept in `EXACT_FITS_ONLY.md` for audit/debug use.",
        "Charged-lepton and quark exact same-family/current-family surfaces remain available as exact-fit audit witnesses, but their numeric values are not public mass outputs under the strict non-circularity policy.",
        "Absolute neutrino masses are likewise withheld here while the absolute attachment remains compare-only; dimensionless PMNS and mass-splitting-ratio comparisons stay in `RESULTS_STATUS.md`.",
        "Photon, gluon, and graviton labels are not emitted as `0 GeV` particle predictions. Their zero hard quadratic parameters are tracked separately as conditional classical/perturbative carrier modes in `CARRIER_MODE_ACCEPTANCE.md`; the quantum-particle gate remains open.",
        "",
        "| Particle | Mass Output | Kind | Scope | Source |",
        "| --- | ---: | --- | --- | --- |",
    ]
    for entry in entries:
        if "mass_gev" in entry:
            value = f"`{entry['mass_gev']} GeV`"
        else:
            value = f"`{entry['mass_eV']} eV`"
        lines.append(
            f"| {entry['label']} | {value} | `{entry['exact_kind']}` | `{entry['scope']}` | `{entry['source_artifact']}` |"
        )
    lines.extend(
        [
            "",
            "## Separated Classical Carrier Modes",
            "",
            "| Carrier | Hard quadratic mass parameter squared | Classical mode gate | Quantum particle gate |",
            "| --- | ---: | --- | --- |",
        ]
    )
    for carrier in carrier_modes:
        lines.append(
            f"| {carrier['label']} | `{carrier['hard_quadratic_mass_parameter_squared']}` | "
            f"`{carrier['classical_carrier_gate']['status']}` | "
            f"`{carrier['quantum_particle_gate']['status']}` |"
        )
    return "\n".join(lines).rstrip()


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the exact non-hadron mass bundle.")
    parser.add_argument("--markdown-out", default=str(DEFAULT_MD_OUT))
    parser.add_argument("--json-out", default=str(DEFAULT_JSON_OUT))
    parser.add_argument("--forward-out", default=str(DEFAULT_FORWARD_OUT))
    args = parser.parse_args()

    generated_utc = _timestamp()
    carrier_acceptance = _load_json(CARRIER_ACCEPTANCE_JSON)
    carrier_modes = [
        _carrier_summary(entry) for entry in carrier_acceptance.get("carriers", [])
    ]
    all_entries = build_all_entries()
    entries = [entry for entry in all_entries if _is_public_mass_output(entry)]
    withheld_entries = [_withheld_entry(entry) for entry in all_entries if not _is_public_mass_output(entry)]
    payload = {
        "artifact": "oph_exact_nonhadron_mass_bundle",
        "generated_utc": generated_utc,
        "status": "public_mass_outputs_with_classical_carriers_separated_and_target_anchored_witnesses_withheld",
        "entries": entries,
        "withheld_entries": withheld_entries,
        "classical_carrier_modes": carrier_modes,
        "carrier_acceptance_artifact": "code/particles/runs/status/carrier_mode_acceptance.json",
        "excluded_lane": "hadrons_compute_bound",
        "audit_surface": "code/particles/EXACT_FITS_ONLY.md",
    }

    markdown_out = pathlib.Path(args.markdown_out)
    markdown_out.parent.mkdir(parents=True, exist_ok=True)
    markdown_out.write_text(
        build_markdown(generated_utc, entries, carrier_modes) + "\n",
        encoding="utf-8",
    )

    json_out = pathlib.Path(args.json_out)
    json_out.parent.mkdir(parents=True, exist_ok=True)
    json_out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    forward_out = pathlib.Path(args.forward_out)
    forward_out.parent.mkdir(parents=True, exist_ok=True)
    forward_out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(f"saved: {markdown_out}")
    print(f"saved: {json_out}")
    print(f"saved: {forward_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
