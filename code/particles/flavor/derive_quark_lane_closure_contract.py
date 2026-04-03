#!/usr/bin/env python3
"""Emit the exact three-step quark closure contract.

Chain role: collect the current theorem boundary and the strongest exact sidecar
surfaces for the quark lane into one machine-readable contract.

Mathematics: the current corpus fixes the exact form of three missing theorem
objects without proving them:
1. a one-scalar D12 mass theorem `quark_d12_t1_value_law`,
2. a sector-attached same-label left-handed lift to the physical CKM shell,
3. a target-free absolute sector readout theorem on that physical sheet.

The same corpus also fixes:
1. the negative local selector closure on `sigma_ref`,
2. the emitted D12 mass ray,
3. the current-family affine sector-mean split,
4. the exact selected-sheet current-family sextet,
5. the continuation-only D12 backread sidecar.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
T1_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_d12_t1_value_law.json"
PHYSICAL_BRANCH_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_physical_branch_repair_theorem.json"
SELECTED_SHEET_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_selected_sheet_closure.json"
EXACT_READOUT_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_exact_readout.json"
BACKREAD_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_d12_internal_backread_continuation_closure.json"
SECTOR_MEAN_SPLIT_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_sector_mean_split.json"
FORWARD_YUKAWAS_JSON = ROOT / "particles" / "runs" / "flavor" / "forward_yukawas.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "quark_lane_closure_contract.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def build_payload(
    t1_value_law: dict[str, Any],
    physical_branch: dict[str, Any],
    selected_sheet: dict[str, Any],
    exact_readout: dict[str, Any],
    backread: dict[str, Any],
    sector_mean_split: dict[str, Any],
    forward_yukawas: dict[str, Any],
) -> dict[str, Any]:
    selected_sigma = selected_sheet["selected_sheet"]["sigma_id"]
    exact_masses = {
        "u": float(exact_readout["predicted_singular_values_u"][0]),
        "d": float(exact_readout["predicted_singular_values_d"][0]),
        "s": float(exact_readout["predicted_singular_values_d"][1]),
        "c": float(exact_readout["predicted_singular_values_u"][1]),
        "b": float(exact_readout["predicted_singular_values_d"][2]),
        "t": float(exact_readout["predicted_singular_values_u"][2]),
    }
    rho_ord = float(sector_mean_split["rho_ord"])
    x2 = float(t1_value_law["sample_same_family_point"]["x2"])
    return {
        "artifact": "oph_quark_lane_closure_contract",
        "base_theorem_emitted_package_artifact": "oph_quark_maximal_theorem_emitted_package",
        "generated_utc": _timestamp(),
        "scope": "quark_lane_theorem_boundary_plus_exact_sidecars",
        "proof_status": "maximal_theorem_emitted_package_identified_three_exact_missing_theorems_recorded_exact_sidecar_masses_recorded",
        "public_promotion_allowed": False,
        "mass_comparison_surface": {
            "kind": "running_mass_comparison_surface",
            "note": "Quark references are running masses, not asymptotic free-particle pole masses.",
        },
        "selected_local_sheet_status": {
            "sigma_id": selected_sigma,
            "proof_status": selected_sheet["proof_status"],
            "theorem_scope": selected_sheet["theorem_scope"],
            "wrong_branch_for_physical_ckm_shell": True,
            "why_not_enough": physical_branch["insufficiency_theorem"]["statement"],
        },
        "exact_sidecar_mass_surface": {
            "artifact": exact_readout["artifact"],
            "scope": exact_readout["theorem_scope"],
            "selected_sheet": selected_sigma,
            "exact_outputs_gev": exact_masses,
            "exact_sector_geometric_means": {
                "g_u": float(exact_readout["g_u"]),
                "g_d": float(exact_readout["g_d"]),
            },
            "closure_statement": selected_sheet["theorem_statement"],
        },
        "continuation_only_mass_sidecar": {
            "artifact": backread["artifact"],
            "scope": backread["scope"],
            "closed_mass_side_package": backread["closed_mass_side_package"],
            "closed_source_side_package": backread["closed_source_side_package"],
            "theorem_boundary_note": backread["theorem_boundary_note"],
        },
        "exact_missing_theorems": [
            {
                "step": 1,
                "id": "quark_d12_t1_value_law",
                "role": "mass_side_value_law_on_emitted_D12_mass_ray",
                "mathematical_name": "Theta_ud^mass",
                "target_free_map": {
                    "domain": "G_ud^light",
                    "codomain": "R_{>= 0}",
                    "formula": "Theta_ud^mass(OPH light data) = t1",
                },
                "must_emit": "t1",
                "equivalently_emits": [
                    "ray_modulus",
                    "Delta_ud_overlap",
                    "eta_Q_centered",
                    "kappa_Q",
                ],
                "must_not_use_target_masses": True,
                "must_not_use_ckm_cp": True,
                "minimal_light_branch": {
                    "epsilon": "1/6",
                    "textures": [
                        "y_u = c_u * epsilon^6",
                        "y_d = c_d * epsilon^6",
                    ],
                    "light_isospin_breaking_scalar": "Delta_ud_overlap = (1/6) * log(c_d / c_u)",
                },
                "same_family_d12_hypercharge_law": {
                    "C1": "C_1(Y) = (3/5) * Y^2",
                    "difference": "C_1(u^c) - C_1(d^c) = 1/5",
                    "ray_identity": "Delta_ud_overlap = t1 / 5",
                    "equivalent_scalar_identities": [
                        "log(c_d / c_u) = (6/5) * t1",
                        "t1 = 5 * Delta_ud_overlap = (5/6) * log(c_d / c_u)",
                    ],
                },
                "induced_formulas": t1_value_law["one_scalar_reduction"]["induced_formulas"],
                "odd_source_pair_forced_by_mass_scalar": {
                    "beta_u_diag_B_source": "t1 / 10",
                    "beta_d_diag_B_source": "-t1 / 10",
                    "source_readback_u_log_per_side": "(-t1/10, 0, +t1/10)",
                    "source_readback_d_log_per_side": "(+t1/10, 0, -t1/10)",
                    "tau_u": "sigma_d * t1 / (10 * (sigma_u + sigma_d))",
                    "tau_d": "sigma_u * t1 / (10 * (sigma_u + sigma_d))",
                },
                "builder_reduction_beneath_public_frontier": {
                    "conditional_on": "pure_B_readback_payload_law",
                    "source_readback_u_log_per_side": "(-t1/10, 0, +t1/10)",
                    "source_readback_d_log_per_side": "(+t1/10, 0, -t1/10)",
                },
            },
            {
                "step": 2,
                "id": "quark_physical_sigma_ud_lift",
                "role": "sector_attached_same_label_left_handed_lift_to_physical_ckm_shell",
                "mathematical_name": "Theta_ud^phys",
                "physical_carrier": {
                    "name": "Sigma_ud^phys",
                    "definition": (
                        "{(sigma_id, tau, U_uL, U_dL, V_CKM, I_CKM) : "
                        "V_CKM = U_uL^dagger U_dL} / ~"
                    ),
                    "equivalence_relation": (
                        "(U_uL, U_dL, V) ~ (U_uL D_u, U_dL D_d, D_u^dagger V D_d) "
                        "for diagonal D_u, D_d in U(1)^3"
                    ),
                },
                "must_emit": "{sigma_id, canonical_token, U_u_left, U_d_left, V_CKM, ckm_invariants}",
                "must_not_use_compare_fit_masses": True,
                "quotient_section_form": {
                    "projection": (
                        "pi: Sigma_ud^phys -> F_ud / (U(1)^3 x U(1)^3), "
                        "sending a same-label left-handed datum to its frame-overlap class [V_CKM]"
                    ),
                    "missing_theorem": (
                        "a sector-attached section s_phys with pi o s_phys = id on the physical frame class"
                    ),
                    "compare_only_frame_signal": "[F_0^dagger F_1]",
                },
                "current_exact_obstruction": (
                    "No sector-attached lift or equivalence currently identifies the common-refinement "
                    "self-overlap matrix with an emitted same-label left-handed Sigma_ud element."
                ),
                "current_negative_closure": {
                    "selected_sigma": selected_sigma,
                    "statement": physical_branch["insufficiency_theorem"]["statement"],
                },
            },
            {
                "step": 3,
                "id": "quark_absolute_sector_readout_theorem",
                "role": "target_free_absolute_sector_scales_on_physical_sheet",
                "mathematical_name": "Theta_ud^abs",
                "must_emit": ["g_u", "g_d"],
                "must_be_target_free": True,
                "must_be_on_physical_sheet": True,
                "current_candidate_surface": {
                    "artifact": sector_mean_split["artifact"],
                    "theorem_scope": sector_mean_split["theorem_scope"],
                    "g_u_candidate": float(sector_mean_split["g_u_candidate"]),
                    "g_d_candidate": float(sector_mean_split["g_d_candidate"]),
                    "candidate_formulas": {
                        "g_u": "g_ch * exp(-(A_ud * sigma_seed_ud - B_ud * eta_ud))",
                        "g_d": "g_ch * exp(-(A_ud * sigma_seed_ud + B_ud * eta_ud))",
                    },
                },
                "affine_mean_split_law": {
                    "sigma_seed_ud": "(sigma_u + sigma_d) / 2",
                    "eta_ud": "(sigma_u - sigma_d) / 2",
                    "A_ud": "1 / (2 * (1 + rho_ord - x2^2))",
                    "B_ud": "1 / (2 * (1 - x2^2 - x2^2 / (1 + rho_ord)))",
                    "rho_ord": rho_ord,
                    "x2": x2,
                    "g_u": "g_ch * exp(-(A_ud * sigma_seed_ud - B_ud * eta_ud))",
                    "g_d": "g_ch * exp(-(A_ud * sigma_seed_ud + B_ud * eta_ud))",
                },
                "ordered_three_point_readout": {
                    "u_sector": "(m_u, m_c, m_t) = g_u * (exp(E_u,1^log), exp(E_u,2^log), exp(E_u,3^log))",
                    "d_sector": "(m_d, m_s, m_b) = g_d * (exp(E_d,1^log), exp(E_d,2^log), exp(E_d,3^log))",
                },
                "current_exact_obstruction": {
                    "forward_artifact": forward_yukawas["artifact"],
                    "certification_status": forward_yukawas["certification_status"],
                    "promotion_blockers": forward_yukawas["promotion_blockers"],
                    "reason": (
                        "The current absolute sector scales live only on a current-family candidate surface; "
                        "the forward Yukawa artifact remains non-promotable until the spread emitter is theorem-fed."
                    ),
                },
                "selected_sheet_alignment": {
                    "g_u_exact_selected_sheet": float(exact_readout["g_u"]),
                    "g_d_exact_selected_sheet": float(exact_readout["g_d"]),
                    "g_u_candidate": float(sector_mean_split["g_u_candidate"]),
                    "g_d_candidate": float(sector_mean_split["g_d_candidate"]),
                },
            },
        ],
        "closure_chain": [
            "(axioms + light-data transport) => Theta_ud^mass = quark_d12_t1_value_law => t1 => (Delta_ud_overlap, eta_Q_centered, kappa_Q, tau_u, tau_d)",
            "(axioms + same-label left-handed sector data) => Theta_ud^phys = quark_physical_sigma_ud_lift => Sigma_ud^phys = (sigma_id, tau, U_uL, U_dL, V_CKM, I_CKM)",
            "(axioms + Sigma_ud^phys) => Theta_ud^abs = quark_absolute_sector_readout_theorem => (g_u, g_d) => (m_u, m_d, m_s, m_c, m_b, m_t)",
        ],
        "notes": [
            "The maximal theorem-emitted quark package on the present ledger is the D12 mass ray, the negative selector closure sigma_ref, and the restricted-scope affine mean package; the three exact missing theorem objects recorded here sit above that package.",
            "The exact current-family witness and the D12 internal backread sidecar exhibit the mass data on sidecar surfaces, but they do not promote the public quark theorem lane.",
            "The selected local same-label left-handed sheet closes negatively to sigma_ref; the remaining physical quark burden is not another local selector search.",
            "Full end-to-end physical quark closure waits on three exact theorem objects: the D12 one-scalar value law, a sector-attached physical CKM lift, and a target-free absolute sector readout on that physical sheet.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the exact three-step quark closure contract.")
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    payload = build_payload(
        _load_json(T1_JSON),
        _load_json(PHYSICAL_BRANCH_JSON),
        _load_json(SELECTED_SHEET_JSON),
        _load_json(EXACT_READOUT_JSON),
        _load_json(BACKREAD_JSON),
        _load_json(SECTOR_MEAN_SPLIT_JSON),
        _load_json(FORWARD_YUKAWAS_JSON),
    )

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
