#!/usr/bin/env python3
"""Build an exact-fits-only surface from particle artifacts.

This surface lists exact target matches on their declared OPH carriers. It
includes theorem-grade selected-class quark closure together with exact
carrier-restricted and compare-only supporting surfaces.
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
CHARGED_JSON = ROOT / "particles" / "runs" / "leptons" / "lepton_current_family_exact_readout.json"
QUARK_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_exact_readout.json"
CHARGED_AFFINE_JSON = ROOT / "particles" / "runs" / "leptons" / "lepton_current_family_affine_anchor_theorem.json"
QUARK_CLOSURE_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_selected_sheet_closure.json"
QUARK_TRANSPORT_LIFT_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_transport_frame_sector_attached_lift.json"
QUARK_TRANSPORT_COMPLETION_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_transport_frame_exact_pdg_completion.json"
QUARK_TRANSPORT_FORWARD_YUKAWAS_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_transport_frame_exact_forward_yukawas.json"
QUARK_END_TO_END_CHAIN_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_end_to_end_exact_pdg_derivation_chain.json"
QUARK_PUBLIC_SIGMA_DESCENT_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_public_physical_sigma_datum_descent.json"
QUARK_PUBLIC_EXACT_YUKAWA_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_public_exact_yukawa_end_to_end_theorem.json"
NEUTRINO_JSON = ROOT / "particles" / "runs" / "neutrino" / "neutrino_compare_only_scale_fit.json"
NEUTRINO_TWO_PARAMETER_JSON = ROOT / "particles" / "runs" / "neutrino" / "neutrino_two_parameter_exact_adapter.json"
NEUTRINO_BRIDGE_COORDINATE_JSON = ROOT / "particles" / "runs" / "neutrino" / "neutrino_exact_adapter_bridge_coordinate.json"
DEFAULT_MD_OUT = ROOT / "particles" / "EXACT_FITS_ONLY.md"
DEFAULT_JSON_OUT = ROOT / "particles" / "exact_fits_only.json"
DEFAULT_FORWARD_OUT = ROOT / "particles" / "runs" / "status" / "exact_fits_only_current.json"


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


def _max_abs(values: list[float]) -> float:
    return max(abs(value) for value in values) if values else 0.0


def build_entries() -> list[dict[str, Any]]:
    references = _load_json(REFERENCE_JSON)["entries"]
    ew_exact = _load_json(EW_EXACT_JSON)
    d11_exact = _load_json(D11_EXACT_JSON)
    charged = _load_json(CHARGED_JSON)
    quark = _load_json(QUARK_JSON)
    charged_affine = _load_optional_json(CHARGED_AFFINE_JSON)
    quark_closure = _load_optional_json(QUARK_CLOSURE_JSON)
    quark_transport_lift = _load_optional_json(QUARK_TRANSPORT_LIFT_JSON)
    quark_transport_completion = _load_optional_json(QUARK_TRANSPORT_COMPLETION_JSON)
    quark_transport_forward_yukawas = _load_optional_json(QUARK_TRANSPORT_FORWARD_YUKAWAS_JSON)
    quark_end_to_end_chain = _load_optional_json(QUARK_END_TO_END_CHAIN_JSON)
    quark_public_sigma_descent = _load_optional_json(QUARK_PUBLIC_SIGMA_DESCENT_JSON)
    quark_public_exact_yukawa = _load_optional_json(QUARK_PUBLIC_EXACT_YUKAWA_JSON)
    neutrino = _load_json(NEUTRINO_JSON)
    neutrino_two_parameter = _load_json(NEUTRINO_TWO_PARAMETER_JSON)
    neutrino_bridge_coordinate = _load_json(NEUTRINO_BRIDGE_COORDINATE_JSON)

    entries: list[dict[str, Any]] = [
        {
            "id": "electroweak_frozen_target_exact_pair",
            "label": "Electroweak Frozen-Target Exact Pair",
            "fit_kind": "exact_frozen_target_compare_only_adapter",
            "scope": "frozen_authoritative_target_surface",
            "promotable": False,
            "matched_observables": ["m_W", "m_Z"],
            "units": "GeV",
            "values": {
                "m_W": ew_exact["coherent_repaired_quintet"]["MW_pole"],
                "m_Z": ew_exact["coherent_repaired_quintet"]["MZ_pole"],
            },
            "references": {
                "m_W": references["w_boson"]["value_gev"],
                "m_Z": references["z_boson"]["value_gev"],
            },
            "max_abs_residual": max(
                abs(ew_exact["coherent_repaired_quintet"]["MW_pole"] - references["w_boson"]["value_gev"]),
                abs(ew_exact["coherent_repaired_quintet"]["MZ_pole"] - references["z_boson"]["value_gev"]),
            ),
            "source_artifact": _repo_ref(EW_EXACT_JSON),
            "note": (
                "Exact on the frozen-authoritative D10 repair surface. The public theorem surface is the "
                "target-free source-only emission, which is separate and differs only at the `1e-8 GeV` scale."
            ),
        },
        {
            "id": "higgs_top_reference_exact_adapter",
            "label": "Higgs/Top Reference Exact Adapter",
            "fit_kind": "exact_target_anchored_compare_only_inverse_slice",
            "scope": d11_exact["scope"],
            "promotable": False,
            "matched_observables": ["m_H", "m_t"],
            "units": "GeV",
            "values": {
                "m_H": d11_exact["predicted_outputs"]["mH_gev"],
                "m_t": d11_exact["predicted_outputs"]["mt_pole_gev"],
            },
            "references": {
                "m_H": d11_exact["exact_reference_targets"]["mH_gev"],
                "m_t": d11_exact["exact_reference_targets"]["mt_pole_gev"],
            },
            "max_abs_residual": max(
                abs(d11_exact["exact_fit_residuals_gev"]["mH_gev"]),
                abs(d11_exact["exact_fit_residuals_gev"]["mt_pole_gev"]),
            ),
            "source_artifact": _repo_ref(D11_EXACT_JSON),
            "note": (
                "Exact only as a compare-only inverse slice on the D11 Jacobian. The live D11 theorem lane uses "
                "the source-only split theorem on the declared D10/D11 surface, not this adapter. "
                "That theorem emits `m_H = 125.1995304097179 GeV` and a companion top coordinate "
                "`m_t = 172.3523553288312 GeV`. The Higgs row lands on the PDG 2025 Higgs average at quoted precision. "
                "The exact public running-top row uses the PDG 2025 cross-section entry `Q007TP4`. "
                "The bridge to the auxiliary direct-top average `Q007TP` is open and is tracked in "
                "[#207](https://github.com/FloatingPragma/observer-patch-holography/issues/207)."
            ),
        },
        {
            "id": "charged_current_family_exact_witness",
            "label": "Charged Current-Family Exact Witness",
            "fit_kind": "exact_target_anchored_current_family_witness",
            "scope": charged.get("theorem_scope", "current_family_only"),
            "promotable": False,
            "matched_observables": ["m_e", "m_mu", "m_tau"],
            "units": "GeV",
            "values": {
                "m_e": charged["predicted_singular_values_abs"][0],
                "m_mu": charged["predicted_singular_values_abs"][1],
                "m_tau": charged["predicted_singular_values_abs"][2],
            },
            "references": {
                "m_e": charged["reference_targets"][0],
                "m_mu": charged["reference_targets"][1],
                "m_tau": charged["reference_targets"][2],
            },
            "max_abs_residual": _max_abs(charged["exact_fit_residuals_abs"]),
            "source_artifact": _repo_ref(CHARGED_JSON),
            "supporting_scope_closure_artifact": _repo_ref(CHARGED_AFFINE_JSON),
            "note": (
                "Exact on the ordered charged eigenvalue triple, with a closed ordered-three-point "
                "readout theorem inside `current_family_only`, and with the scoped affine coordinate "
                "`A_ch_current_family` closed on that same exact family. The charged theorem lane does not "
                "emit a theorem-grade absolute anchor."
            ),
        },
        {
            "id": "quark_current_family_exact_witness",
            "label": "Quark Current-Family Exact Witness",
            "fit_kind": "exact_target_anchored_current_family_witness",
            "scope": quark.get("theorem_scope", "current_family_only"),
            "promotable": False,
            "matched_observables": ["m_u", "m_c", "m_t", "m_d", "m_s", "m_b"],
            "units": "GeV",
            "values": {
                "m_u": quark["predicted_singular_values_u"][0],
                "m_c": quark["predicted_singular_values_u"][1],
                "m_t": quark["predicted_singular_values_u"][2],
                "m_d": quark["predicted_singular_values_d"][0],
                "m_s": quark["predicted_singular_values_d"][1],
                "m_b": quark["predicted_singular_values_d"][2],
            },
            "references": {
                "m_u": quark["reference_targets_u"][0],
                "m_c": quark["reference_targets_u"][1],
                "m_t": quark["reference_targets_u"][2],
                "m_d": quark["reference_targets_d"][0],
                "m_s": quark["reference_targets_d"][1],
                "m_b": quark["reference_targets_d"][2],
            },
            "max_abs_residual": max(
                _max_abs(quark["exact_fit_residuals_u"]),
                _max_abs(quark["exact_fit_residuals_d"]),
            ),
            "source_artifact": _repo_ref(QUARK_JSON),
            "supporting_scope_closure_artifact": _repo_ref(QUARK_CLOSURE_JSON),
            "supporting_transport_frame_artifact": _repo_ref(QUARK_TRANSPORT_LIFT_JSON),
            "supporting_transport_frame_completion_artifact": _repo_ref(QUARK_TRANSPORT_COMPLETION_JSON),
            "supporting_transport_frame_forward_yukawas_artifact": _repo_ref(QUARK_TRANSPORT_FORWARD_YUKAWAS_JSON),
            "supporting_end_to_end_chain_artifact": _repo_ref(QUARK_END_TO_END_CHAIN_JSON),
            "supporting_public_sigma_descent_artifact": _repo_ref(QUARK_PUBLIC_SIGMA_DESCENT_JSON),
            "supporting_public_exact_yukawa_artifact": _repo_ref(QUARK_PUBLIC_EXACT_YUKAWA_JSON),
            "note": (
                "Exact on the official PDG 2025 API running-quark target surface on the ordered three-point "
                "quark family witness, with the internal same-family quadratic readout closed on the fixed "
                "carrier and the selected-sheet exact closure packaged on `sigma_ref`. The top coordinate "
                f"uses PDG summary `{references['top_quark']['source']['summary_id']}`. "
                "The bridge to the auxiliary direct-top entry "
                f"`{references['top_quark_direct_aux']['source']['summary_id']}` is open and is tracked in "
                "[#207](https://github.com/FloatingPragma/observer-patch-holography/issues/207). "
                "The declared scope is `current_family_only`. A separate restricted theorem chain emits a "
                "sector-attached `Sigma_ud^phys` element on the explicit "
                f"`{(quark_transport_lift or {}).get('theorem_scope', 'current_family_common_refinement_transport_frame_only')}` carrier, and the merged transport-frame theorem "
                f"reconstructs the same running sextet exactly on `{(quark_transport_completion or {}).get('theorem_scope', 'current_family_common_refinement_transport_frame_only')}`. "
                "The declared transport-frame chain also closes explicit exact forward Yukawas `Y_u` and "
                f"`Y_d` with certification status `{(quark_transport_forward_yukawas or {}).get('certification_status', 'forward_matrix_certified')}`, "
                f"and the full declared-carrier chain is recorded in `{(quark_end_to_end_chain or {}).get('artifact', 'oph_quark_current_family_end_to_end_exact_pdg_derivation_chain')}`. "
                "A separate target-free mass bridge closes `Delta_ud_overlap = (1/6) * log(c_d / c_u)`, "
                "equivalently `quark_d12_t1_value_law`, on the emitted D12 ray. "
                + (
                    "A separate public theorem closes on the selected public physical quark frame class chosen by `P`: "
                    f"`{quark_public_sigma_descent['artifact']}` makes the exact physical sigma datum target-free "
                    "public on that selected class, and "
                    f"`{quark_public_exact_yukawa['artifact']}` emits the same exact sextet together with explicit "
                    "exact forward Yukawas `Y_u` and `Y_d`. This entry remains an exact-fit surface rather than that "
                    "public selected-class theorem."
                    if quark_public_sigma_descent and quark_public_exact_yukawa
                    else "This entry is the strongest exact carrier-restricted quark surface present in the local artifact set."
                )
            ),
        },
        {
            "id": "neutrino_two_parameter_exact_adapter",
            "label": "Neutrino Two-Parameter Exact Adapter",
            "fit_kind": "exact_two_observable_compare_only_segment_adapter",
            "scope": neutrino_two_parameter["scope"],
            "promotable": False,
            "matched_observables": ["Delta m21^2", "Delta m32^2"],
            "units": "eV / eV^2",
            "values": {
                "m1_eV": neutrino_two_parameter["exact_outputs"]["masses_eV"][0],
                "m2_eV": neutrino_two_parameter["exact_outputs"]["masses_eV"][1],
                "m3_eV": neutrino_two_parameter["exact_outputs"]["masses_eV"][2],
                "delta_m21_sq_eV2": neutrino_two_parameter["exact_outputs"]["delta_m_sq_eV2"]["21"],
                "delta_m31_sq_eV2": neutrino_two_parameter["exact_outputs"]["delta_m_sq_eV2"]["31"],
                "delta_m32_sq_eV2": neutrino_two_parameter["exact_outputs"]["delta_m_sq_eV2"]["32"],
            },
            "references": {
                "delta_m21_sq_eV2": neutrino_two_parameter["reference_central_values"]["delta_m21_sq_eV2"],
                "delta_m31_sq_eV2": neutrino_two_parameter["reference_central_values"]["delta_m31_sq_eV2"],
                "delta_m32_sq_eV2": neutrino_two_parameter["reference_central_values"]["delta_m32_sq_eV2"],
            },
            "max_abs_residual": max(
                abs(neutrino_two_parameter["exact_fit_residuals_eV2"]["21"]),
                abs(neutrino_two_parameter["exact_fit_residuals_eV2"]["31"]),
                abs(neutrino_two_parameter["exact_fit_residuals_eV2"]["32"]),
            ),
            "source_artifact": _repo_ref(NEUTRINO_TWO_PARAMETER_JSON),
            "supporting_bridge_coordinate_artifact": _repo_ref(NEUTRINO_BRIDGE_COORDINATE_JSON),
            "note": (
                "Exact compare-only fit to both representative PDG central splittings by moving along the explicit "
                "positive selector segment and then rescaling with one positive `lambda_nu`. It is diagnostic-only "
                "after the emitted weighted-cycle bridge-rigidity and absolute-attachment theorems. On that same "
                "exact compare-only branch, the explicit bridge coordinates are "
                f"`B_nu = {neutrino_bridge_coordinate['bridge_coordinates']['paper_facing_amplitude']['value']:.8f}` and "
                f"`C_nu = {neutrino_bridge_coordinate['bridge_coordinates']['reduced_correction_invariant']['value']:.8f}`, "
                "but they remain sidecars and must not feed back into theorem state."
            ),
        },
    ]

    if quark_public_sigma_descent and quark_public_exact_yukawa:
        entries.insert(
            3,
            {
                "id": "quark_selected_class_exact_theorem",
                "label": "Quark Selected-Class Exact Theorem",
                "fit_kind": "selected_class_theorem_grade_exact_forward_quark_closure",
                "scope": quark_public_exact_yukawa["theorem_scope"],
                "promotable": True,
                "matched_observables": ["m_u", "m_c", "m_t", "m_d", "m_s", "m_b"],
                "units": "GeV",
                "values": {
                    "m_u": quark_public_exact_yukawa["public_exact_outputs"]["exact_running_values_gev"]["u"],
                    "m_c": quark_public_exact_yukawa["public_exact_outputs"]["exact_running_values_gev"]["c"],
                    "m_t": quark_public_exact_yukawa["public_exact_outputs"]["exact_running_values_gev"]["t"],
                    "m_d": quark_public_exact_yukawa["public_exact_outputs"]["exact_running_values_gev"]["d"],
                    "m_s": quark_public_exact_yukawa["public_exact_outputs"]["exact_running_values_gev"]["s"],
                    "m_b": quark_public_exact_yukawa["public_exact_outputs"]["exact_running_values_gev"]["b"],
                },
                "references": {
                    "m_u": quark["reference_targets_u"][0],
                    "m_c": quark["reference_targets_u"][1],
                    "m_t": quark["reference_targets_u"][2],
                    "m_d": quark["reference_targets_d"][0],
                    "m_s": quark["reference_targets_d"][1],
                    "m_b": quark["reference_targets_d"][2],
                },
                "max_abs_residual": max(
                    abs(
                        quark_public_exact_yukawa["public_exact_outputs"]["exact_running_values_gev"]["u"]
                        - quark["reference_targets_u"][0]
                    ),
                    abs(
                        quark_public_exact_yukawa["public_exact_outputs"]["exact_running_values_gev"]["c"]
                        - quark["reference_targets_u"][1]
                    ),
                    abs(
                        quark_public_exact_yukawa["public_exact_outputs"]["exact_running_values_gev"]["t"]
                        - quark["reference_targets_u"][2]
                    ),
                    abs(
                        quark_public_exact_yukawa["public_exact_outputs"]["exact_running_values_gev"]["d"]
                        - quark["reference_targets_d"][0]
                    ),
                    abs(
                        quark_public_exact_yukawa["public_exact_outputs"]["exact_running_values_gev"]["s"]
                        - quark["reference_targets_d"][1]
                    ),
                    abs(
                        quark_public_exact_yukawa["public_exact_outputs"]["exact_running_values_gev"]["b"]
                        - quark["reference_targets_d"][2]
                    ),
                ),
                "source_artifact": _repo_ref(QUARK_PUBLIC_EXACT_YUKAWA_JSON),
                "supporting_public_sigma_descent_artifact": _repo_ref(QUARK_PUBLIC_SIGMA_DESCENT_JSON),
                "supporting_scope_closure_artifact": _repo_ref(QUARK_CLOSURE_JSON),
                "supporting_transport_frame_artifact": _repo_ref(QUARK_TRANSPORT_LIFT_JSON),
                "supporting_transport_frame_completion_artifact": _repo_ref(QUARK_TRANSPORT_COMPLETION_JSON),
                "supporting_transport_frame_forward_yukawas_artifact": _repo_ref(QUARK_TRANSPORT_FORWARD_YUKAWAS_JSON),
                "supporting_end_to_end_chain_artifact": _repo_ref(QUARK_END_TO_END_CHAIN_JSON),
                "note": (
                    "Exact theorem on the selected public physical quark frame class chosen by `P`. "
                    f"`{quark_public_sigma_descent['artifact']}` makes the exact physical sigma datum target-free public "
                    "on that selected class, and "
                    f"`{quark_public_exact_yukawa['artifact']}` emits the exact PDG 2025 running-quark sextet together "
                    "with explicit exact forward Yukawas `Y_u` and `Y_d`. The top coordinate uses PDG summary "
                    f"`{references['top_quark']['source']['summary_id']}`. The bridge to the auxiliary direct-top entry "
                    f"`{references['top_quark_direct_aux']['source']['summary_id']}` is open and is tracked in "
                    "[#207](https://github.com/FloatingPragma/observer-patch-holography/issues/207). This is selected-class closure "
                    "only. It does not claim a global classification of all quark frame classes."
                ),
            },
        )

    fits = neutrino.get("fits", {})
    for fit_name, matched_key, key_label in (
        ("atmospheric_only", "32", "Delta m32^2"),
        ("solar_only", "21", "Delta m21^2"),
    ):
        fit = dict(fits.get(fit_name, {}))
        if not fit:
            continue
        entries.append(
            {
                "id": f"neutrino_{fit_name}_exact_adapter",
                "label": f"Neutrino {fit_name.replace('_', ' ').title()} Exact Adapter",
                "fit_kind": "exact_single_observable_compare_only_adapter",
                "scope": "compare_only",
                "promotable": False,
                "matched_observables": [key_label],
                "units": "eV / eV^2",
                "values": {
                    "m1_eV": fit["masses_eV"][0],
                    "m2_eV": fit["masses_eV"][1],
                    "m3_eV": fit["masses_eV"][2],
                    "delta_m21_sq_eV2": fit["delta_m_sq_eV2"]["21"],
                    "delta_m31_sq_eV2": fit["delta_m_sq_eV2"]["31"],
                    "delta_m32_sq_eV2": fit["delta_m_sq_eV2"]["32"],
                },
                "references": {
                    "delta_m21_sq_eV2": neutrino["reference_central_values"]["delta_m21_sq_eV2"],
                    "delta_m32_sq_eV2": neutrino["reference_central_values"]["delta_m32_sq_eV2"],
                },
                "exact_match_observable": key_label,
                "source_artifact": _repo_ref(NEUTRINO_JSON),
                "note": (
                    "Exact only for one splitting observable on the repaired weighted-cycle family. "
                    "The same artifact states that no single `lambda_nu` hits both central splittings exactly."
                ),
            }
        )

    return entries


def build_markdown(generated_utc: str, entries: list[dict[str, Any]]) -> str:
    has_quark_selected_class_theorem = any(entry["id"] == "quark_selected_class_exact_theorem" for entry in entries)
    lines = [
        "# Exact Fits Only",
        "",
        f"Generated: `{generated_utc}`",
        "",
        "This surface lists exact target matches on declared OPH carriers. "
        "It separates theorem-grade selected-class outputs from compare-only and carrier-restricted exact surfaces.",
        (
            "For quarks, the selected-class theorem and its supporting exact carriers coincide with the official PDG 2025 API running-quark target surface."
            if has_quark_selected_class_theorem
            else "For quarks, the exact carrier-restricted witnesses coincide with the official PDG 2025 API running-quark target surface on their declared carriers."
        ),
        "",
    ]
    for entry in entries:
        lines.extend(
            [
                f"## {entry['label']}",
                "",
                f"- Fit kind: `{entry['fit_kind']}`",
                f"- Scope: `{entry['scope']}`",
                f"- Promotable: `{str(entry['promotable']).lower()}`",
                f"- Source artifact: `{entry['source_artifact']}`",
            ]
        )
        if "max_abs_residual" in entry:
            lines.append(f"- Max absolute residual: `{entry['max_abs_residual']}`")
        if "exact_match_observable" in entry:
            lines.append(f"- Exact matched observable: `{entry['exact_match_observable']}`")
        lines.append(f"- Note: {entry['note']}")
        lines.append("")
        lines.append("| Observable | Value | Reference |")
        lines.append("| --- | ---: | ---: |")
        for key, value in entry["values"].items():
            reference = entry["references"].get(key, "n/a")
            lines.append(f"| `{key}` | `{value}` | `{reference}` |")
        lines.append("")
    return "\n".join(lines).rstrip()


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the exact-fits-only particle surface.")
    parser.add_argument("--markdown-out", default=str(DEFAULT_MD_OUT))
    parser.add_argument("--json-out", default=str(DEFAULT_JSON_OUT))
    parser.add_argument("--forward-out", default=str(DEFAULT_FORWARD_OUT))
    args = parser.parse_args()

    generated_utc = _timestamp()
    entries = build_entries()
    payload = {
        "artifact": "oph_exact_fits_only_surface",
        "generated_utc": generated_utc,
        "status": "exact_target_anchored_or_single_observable_diagnostic_surface",
        "entries": entries,
    }

    markdown_out = pathlib.Path(args.markdown_out)
    markdown_out.parent.mkdir(parents=True, exist_ok=True)
    markdown_out.write_text(build_markdown(generated_utc, entries) + "\n", encoding="utf-8")

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
