#!/usr/bin/env python3
"""Export a blind flavor-dictionary artifact for downstream flavor lanes."""

from __future__ import annotations

import argparse
import json
import pathlib
from datetime import datetime, timezone


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_OBSERVABLE = ROOT / "particles" / "runs" / "flavor" / "flavor_observable_artifact.json"
DEFAULT_COCYCLE = ROOT / "particles" / "runs" / "flavor" / "overlap_edge_transport_cocycle.json"
DEFAULT_LINE_LIFT = ROOT / "particles" / "runs" / "flavor" / "overlap_edge_line_lift.json"
DEFAULT_PUSHFORWARD = ROOT / "particles" / "runs" / "flavor" / "sector_transport_pushforward.json"
DEFAULT_CHARGED_BUDGET = ROOT / "particles" / "runs" / "flavor" / "charged_budget_transport.json"
DEFAULT_TENSORS = ROOT / "particles" / "runs" / "flavor" / "suppression_phase_tensors.json"
DEFAULT_QUARK_ODD_FORM = ROOT / "particles" / "runs" / "flavor" / "charged_dirac_odd_deformation_form.json"
DEFAULT_QUARK_RESPONSE_LAW = ROOT / "particles" / "runs" / "flavor" / "quark_odd_response_law.json"
DEFAULT_QUARK_DESCENT = ROOT / "particles" / "runs" / "flavor" / "quark_sector_descent.json"
DEFAULT_YUKAWAS = ROOT / "particles" / "runs" / "flavor" / "forward_yukawas.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "flavor_dictionary_artifact.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def main() -> int:
    parser = argparse.ArgumentParser(description="Export a blind flavor-dictionary artifact.")
    parser.add_argument("--observable", default=str(DEFAULT_OBSERVABLE), help="Input flavor-observable JSON path.")
    parser.add_argument("--cocycle", default=str(DEFAULT_COCYCLE), help="Input overlap-edge transport cocycle JSON path.")
    parser.add_argument("--line-lift", default=str(DEFAULT_LINE_LIFT), help="Input overlap-edge line-lift JSON path.")
    parser.add_argument("--pushforward", default=str(DEFAULT_PUSHFORWARD), help="Input sector-response JSON path.")
    parser.add_argument("--charged-budget", default=str(DEFAULT_CHARGED_BUDGET), help="Input charged-budget JSON path.")
    parser.add_argument("--tensors", default=str(DEFAULT_TENSORS), help="Input suppression/phase tensor JSON path.")
    parser.add_argument("--quark-odd-form", default=str(DEFAULT_QUARK_ODD_FORM), help="Input charged Dirac odd deformation-form JSON path.")
    parser.add_argument("--quark-response-law", default=str(DEFAULT_QUARK_RESPONSE_LAW), help="Input quark odd-response-law JSON path.")
    parser.add_argument("--quark-descent", default=str(DEFAULT_QUARK_DESCENT), help="Input quark-sector-descent JSON path.")
    parser.add_argument("--yukawas", default=str(DEFAULT_YUKAWAS), help="Input forward Yukawa JSON path.")
    parser.add_argument("--output", default=str(DEFAULT_OUT), help="Output JSON path.")
    args = parser.parse_args()

    observable = json.loads(pathlib.Path(args.observable).read_text(encoding="utf-8"))
    cocycle_path = pathlib.Path(args.cocycle)
    cocycle = json.loads(cocycle_path.read_text(encoding="utf-8")) if cocycle_path.exists() else {}
    line_lift_path = pathlib.Path(args.line_lift)
    line_lift = json.loads(line_lift_path.read_text(encoding="utf-8")) if line_lift_path.exists() else {}
    pushforward = json.loads(pathlib.Path(args.pushforward).read_text(encoding="utf-8"))
    charged_budget_path = pathlib.Path(args.charged_budget)
    charged_budget = json.loads(charged_budget_path.read_text(encoding="utf-8")) if charged_budget_path.exists() else {}
    tensors = json.loads(pathlib.Path(args.tensors).read_text(encoding="utf-8"))
    quark_odd_form_path = pathlib.Path(args.quark_odd_form)
    quark_odd_form = json.loads(quark_odd_form_path.read_text(encoding="utf-8")) if quark_odd_form_path.exists() else {}
    quark_response_path = pathlib.Path(args.quark_response_law)
    quark_response = json.loads(quark_response_path.read_text(encoding="utf-8")) if quark_response_path.exists() else {}
    quark_descent_path = pathlib.Path(args.quark_descent)
    quark_descent = json.loads(quark_descent_path.read_text(encoding="utf-8")) if quark_descent_path.exists() else {}
    yukawa_path = pathlib.Path(args.yukawas)
    yukawas = json.loads(yukawa_path.read_text(encoding="utf-8")) if yukawa_path.exists() else {}
    artifact = {
        "artifact": "oph_flavor_dictionary",
        "generated_utc": _timestamp(),
        "observable_kind": observable.get("observable_kind"),
        "observable_proof_status": observable.get("proof_status"),
        "labels": observable.get("labels"),
        "intrinsic_label_order": observable.get("intrinsic_label_order"),
        "family_projectors": observable.get("family_projectors"),
        "family_eigenvalues": observable.get("family_eigenvalues"),
        "spectral_gaps": observable.get("spectral_gaps"),
        "pairwise_suppression": observable.get("pairwise_suppression"),
        "cycle_phases": observable.get("cycle_phases"),
        "persistent_projector_certificate": observable.get("persistent_projector_certificate"),
        "persistent_spectral_triple": observable.get("persistent_spectral_triple"),
        "overlap_edge_transport_cocycle": cocycle,
        "overlap_edge_line_lift": line_lift,
        "pairwise_suppression_directed": pushforward.get("pairwise_suppression_directed"),
        "cycle_holonomy_class": pushforward.get("cycle_holonomy_class"),
        "pushforward_proof_status": pushforward.get("proof_status"),
        "refinement_stability": observable.get("refinement_stability"),
        "gauge_invariance_checks": observable.get("gauge_invariance_checks"),
        "sector_response_object": pushforward.get("sector_response_object"),
        "charged_budget_transport": charged_budget,
        "charged_dirac_scalarization_certificate": {
            sector: (pushforward.get("sector_response_object") or {}).get(sector, {}).get("charged_dirac_scalarization_certificate")
            for sector in ("u", "d", "e")
        },
        "S_u": tensors.get("S_u"),
        "S_d": tensors.get("S_d"),
        "S_e": tensors.get("S_e"),
        "S_nu": tensors.get("S_nu"),
        "Phi_u": tensors.get("Phi_u"),
        "Phi_d": tensors.get("Phi_d"),
        "Phi_e": tensors.get("Phi_e"),
        "phi_nu_unclosed": tensors.get("phi_nu_unclosed"),
        "majorana_lift_rule": tensors.get("majorana_lift_rule"),
        "u_normalization_class": tensors.get("u_normalization_class"),
        "d_normalization_class": tensors.get("d_normalization_class"),
        "e_normalization_class": tensors.get("e_normalization_class"),
        "nu_normalization_class": tensors.get("nu_normalization_class"),
        "u_raw_channel_norm_candidate": tensors.get("u_raw_channel_norm_candidate"),
        "d_raw_channel_norm_candidate": tensors.get("d_raw_channel_norm_candidate"),
        "e_raw_channel_norm_candidate": tensors.get("e_raw_channel_norm_candidate"),
        "nu_raw_channel_norm_candidate": tensors.get("nu_raw_channel_norm_candidate"),
        "charged_dirac_odd_deformation_form": quark_odd_form,
        "quark_odd_response_law": quark_response,
        "quark_sector_descent": quark_descent,
        "Y_u": yukawas.get("Y_u"),
        "Y_d": yukawas.get("Y_d"),
        "quark_certification_status": yukawas.get("certification_status"),
        "quark_promotion_blockers": yukawas.get("promotion_blockers"),
        "template_amplitude_fallback_used": yukawas.get("template_amplitude_fallback_used"),
        "up_down_sector_distinct": yukawas.get("up_down_sector_distinct"),
        "singular_values_u": yukawas.get("singular_values_u"),
        "singular_values_d": yukawas.get("singular_values_d"),
        "U_u_left": yukawas.get("U_u_left"),
        "U_d_left": yukawas.get("U_d_left"),
        "V_CKM": yukawas.get("V_CKM"),
        "jarlskog": yukawas.get("jarlskog"),
        "commutator_invariant": yukawas.get("commutator_invariant"),
        "metadata": {
            "observable_artifact": observable.get("artifact"),
            "cocycle_artifact": cocycle.get("artifact"),
            "line_lift_artifact": line_lift.get("artifact"),
            "pushforward_artifact": pushforward.get("artifact"),
            "charged_budget_artifact": charged_budget.get("artifact"),
            "tensor_artifact": tensors.get("artifact"),
            "quark_odd_form_artifact": quark_odd_form.get("artifact"),
            "quark_response_artifact": quark_response.get("artifact"),
            "quark_descent_artifact": quark_descent.get("artifact"),
            "yukawa_artifact": yukawas.get("artifact"),
            "note": "Blind forward dictionary artifact for flavor-theorem development.",
        },
    }

    out_path = pathlib.Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
