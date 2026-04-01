#!/usr/bin/env python3
"""Build a canonical exact non-hadron mass bundle.

This bundle consolidates the strongest exact non-hadron mass outputs currently
on disk into one deduplicated surface: structural zeros, exact electroweak
sidecar masses, exact Higgs readout, exact charged-lepton current-family
witness, exact quark current-family witness, and exact neutrino compare-only
adapter masses.
"""

from __future__ import annotations

import argparse
import json
import pathlib
from datetime import datetime, timezone
from typing import Any


ROOT = pathlib.Path(__file__).resolve().parents[2]
EW_EXACT_JSON = ROOT / "particles" / "runs" / "calibration" / "d10_ew_w_anchor_neutral_shear_factorization_official_pdg_2025_update.json"
D11_EXACT_JSON = ROOT / "particles" / "runs" / "calibration" / "d11_reference_exact_adapter.json"
CHARGED_JSON = ROOT / "particles" / "runs" / "leptons" / "lepton_current_family_exact_readout.json"
QUARK_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_exact_readout.json"
CHARGED_THEOREM_JSON = ROOT / "particles" / "runs" / "leptons" / "lepton_current_family_quadratic_readout_theorem.json"
QUARK_THEOREM_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_quadratic_readout_theorem.json"
CHARGED_AFFINE_JSON = ROOT / "particles" / "runs" / "leptons" / "lepton_current_family_affine_anchor_theorem.json"
QUARK_CLOSURE_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_selected_sheet_closure.json"
NEUTRINO_JSON = ROOT / "particles" / "runs" / "neutrino" / "neutrino_two_parameter_exact_adapter.json"
DEFAULT_MD_OUT = ROOT / "particles" / "EXACT_NONHADRON_MASSES.md"
DEFAULT_JSON_OUT = ROOT / "particles" / "exact_nonhadron_masses.json"
DEFAULT_FORWARD_OUT = ROOT / "particles" / "runs" / "status" / "exact_nonhadron_masses_current.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: pathlib.Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _repo_ref(path: pathlib.Path) -> str:
    return str(path.relative_to(ROOT.parent))


def build_entries() -> list[dict[str, Any]]:
    ew_exact = _load_json(EW_EXACT_JSON)
    d11_exact = _load_json(D11_EXACT_JSON)
    charged = _load_json(CHARGED_JSON)
    quark = _load_json(QUARK_JSON)
    charged_affine = _load_json(CHARGED_AFFINE_JSON)
    quark_closure = _load_json(QUARK_CLOSURE_JSON)
    neutrino = _load_json(NEUTRINO_JSON)

    return [
        {
            "particle_id": "photon",
            "label": "Photon",
            "mass_gev": 0.0,
            "exact_kind": "structural_zero",
            "scope": "structural",
            "promotable": True,
            "source_artifact": "structural_gauge_redundancy_surface",
            "note": "Exact structural zero from the nonbroken U(1) overlap-gluing redundancy.",
        },
        {
            "particle_id": "gluon",
            "label": "Gluon",
            "mass_gev": 0.0,
            "exact_kind": "structural_zero",
            "scope": "structural",
            "promotable": True,
            "source_artifact": "structural_color_gauge_surface",
            "note": "Exact structural zero for the color gauge sector.",
        },
        {
            "particle_id": "graviton",
            "label": "Graviton",
            "mass_gev": 0.0,
            "exact_kind": "structural_zero",
            "scope": "structural",
            "promotable": True,
            "source_artifact": "structural_diffeomorphism_redundancy_surface",
            "note": "Exact structural zero from bulk diffeomorphism redundancy.",
        },
        {
            "particle_id": "w_boson",
            "label": "W Boson",
            "mass_gev": ew_exact["coherent_repaired_quintet"]["MW_pole"],
            "exact_kind": "exact_frozen_target_compare_only_adapter",
            "scope": "frozen_authoritative_target_surface",
            "promotable": False,
            "source_artifact": _repo_ref(EW_EXACT_JSON),
            "note": "Exact on the frozen D10 repair surface.",
        },
        {
            "particle_id": "z_boson",
            "label": "Z Boson",
            "mass_gev": ew_exact["coherent_repaired_quintet"]["MZ_pole"],
            "exact_kind": "exact_frozen_target_compare_only_adapter",
            "scope": "frozen_authoritative_target_surface",
            "promotable": False,
            "source_artifact": _repo_ref(EW_EXACT_JSON),
            "note": "Exact on the frozen D10 repair surface.",
        },
        {
            "particle_id": "higgs",
            "label": "Higgs Boson",
            "mass_gev": d11_exact["predicted_outputs"]["mH_gev"],
            "exact_kind": "exact_target_anchored_compare_only_inverse_slice",
            "scope": d11_exact["scope"],
            "promotable": False,
            "source_artifact": _repo_ref(D11_EXACT_JSON),
            "note": "Exact compare-only inverse slice on the D11 Jacobian.",
        },
        {
            "particle_id": "electron",
            "label": "Electron",
            "mass_gev": charged["predicted_singular_values_abs"][0],
            "exact_kind": "exact_target_anchored_current_family_witness",
            "scope": charged["theorem_scope"],
            "promotable": False,
            "source_artifact": _repo_ref(CHARGED_JSON),
            "supporting_theorem_artifact": _repo_ref(CHARGED_THEOREM_JSON),
            "supporting_scope_closure_artifact": _repo_ref(CHARGED_AFFINE_JSON),
            "note": "Exact current-family charged-lepton witness on a closed ordered-three-point readout chain, with the scoped affine coordinate A_ch_current_family closed on the same exact family.",
        },
        {
            "particle_id": "muon",
            "label": "Muon",
            "mass_gev": charged["predicted_singular_values_abs"][1],
            "exact_kind": "exact_target_anchored_current_family_witness",
            "scope": charged["theorem_scope"],
            "promotable": False,
            "source_artifact": _repo_ref(CHARGED_JSON),
            "supporting_theorem_artifact": _repo_ref(CHARGED_THEOREM_JSON),
            "supporting_scope_closure_artifact": _repo_ref(CHARGED_AFFINE_JSON),
            "note": "Exact current-family charged-lepton witness on a closed ordered-three-point readout chain, with the scoped affine coordinate A_ch_current_family closed on the same exact family.",
        },
        {
            "particle_id": "tau",
            "label": "Tau",
            "mass_gev": charged["predicted_singular_values_abs"][2],
            "exact_kind": "exact_target_anchored_current_family_witness",
            "scope": charged["theorem_scope"],
            "promotable": False,
            "source_artifact": _repo_ref(CHARGED_JSON),
            "supporting_theorem_artifact": _repo_ref(CHARGED_THEOREM_JSON),
            "supporting_scope_closure_artifact": _repo_ref(CHARGED_AFFINE_JSON),
            "note": "Exact current-family charged-lepton witness on a closed ordered-three-point readout chain, with the scoped affine coordinate A_ch_current_family closed on the same exact family.",
        },
        {
            "particle_id": "up_quark",
            "label": "Up Quark",
            "mass_gev": quark["predicted_singular_values_u"][0],
            "exact_kind": "exact_target_anchored_current_family_witness",
            "scope": quark["theorem_scope"],
            "promotable": False,
            "source_artifact": _repo_ref(QUARK_JSON),
            "supporting_theorem_artifact": _repo_ref(QUARK_THEOREM_JSON),
            "supporting_scope_closure_artifact": _repo_ref(QUARK_CLOSURE_JSON),
            "note": "Exact current-family quark witness on the selected sigma_ref sheet, with the selected-sheet exact readout chain closed on current_family_only.",
        },
        {
            "particle_id": "charm_quark",
            "label": "Charm Quark",
            "mass_gev": quark["predicted_singular_values_u"][1],
            "exact_kind": "exact_target_anchored_current_family_witness",
            "scope": quark["theorem_scope"],
            "promotable": False,
            "source_artifact": _repo_ref(QUARK_JSON),
            "supporting_theorem_artifact": _repo_ref(QUARK_THEOREM_JSON),
            "supporting_scope_closure_artifact": _repo_ref(QUARK_CLOSURE_JSON),
            "note": "Exact current-family quark witness on the selected sigma_ref sheet, with the selected-sheet exact readout chain closed on current_family_only.",
        },
        {
            "particle_id": "top_quark",
            "label": "Top Quark",
            "mass_gev": quark["predicted_singular_values_u"][2],
            "exact_kind": "exact_target_anchored_current_family_witness",
            "scope": quark["theorem_scope"],
            "promotable": False,
            "source_artifact": _repo_ref(QUARK_JSON),
            "supporting_theorem_artifact": _repo_ref(QUARK_THEOREM_JSON),
            "supporting_scope_closure_artifact": _repo_ref(QUARK_CLOSURE_JSON),
            "note": "Exact current-family quark witness on the selected sigma_ref sheet, with the selected-sheet exact readout chain closed on current_family_only.",
        },
        {
            "particle_id": "down_quark",
            "label": "Down Quark",
            "mass_gev": quark["predicted_singular_values_d"][0],
            "exact_kind": "exact_target_anchored_current_family_witness",
            "scope": quark["theorem_scope"],
            "promotable": False,
            "source_artifact": _repo_ref(QUARK_JSON),
            "supporting_theorem_artifact": _repo_ref(QUARK_THEOREM_JSON),
            "supporting_scope_closure_artifact": _repo_ref(QUARK_CLOSURE_JSON),
            "note": "Exact current-family quark witness on the selected sigma_ref sheet, with the selected-sheet exact readout chain closed on current_family_only.",
        },
        {
            "particle_id": "strange_quark",
            "label": "Strange Quark",
            "mass_gev": quark["predicted_singular_values_d"][1],
            "exact_kind": "exact_target_anchored_current_family_witness",
            "scope": quark["theorem_scope"],
            "promotable": False,
            "source_artifact": _repo_ref(QUARK_JSON),
            "supporting_theorem_artifact": _repo_ref(QUARK_THEOREM_JSON),
            "supporting_scope_closure_artifact": _repo_ref(QUARK_CLOSURE_JSON),
            "note": "Exact current-family quark witness on the selected sigma_ref sheet, with the selected-sheet exact readout chain closed on current_family_only.",
        },
        {
            "particle_id": "bottom_quark",
            "label": "Bottom Quark",
            "mass_gev": quark["predicted_singular_values_d"][2],
            "exact_kind": "exact_target_anchored_current_family_witness",
            "scope": quark["theorem_scope"],
            "promotable": False,
            "source_artifact": _repo_ref(QUARK_JSON),
            "supporting_theorem_artifact": _repo_ref(QUARK_THEOREM_JSON),
            "supporting_scope_closure_artifact": _repo_ref(QUARK_CLOSURE_JSON),
            "note": "Exact current-family quark witness on the selected sigma_ref sheet, with the selected-sheet exact readout chain closed on current_family_only.",
        },
        {
            "particle_id": "electron_neutrino",
            "label": "Electron Neutrino",
            "mass_eV": neutrino["exact_outputs"]["masses_eV"][0],
            "exact_kind": "exact_two_observable_compare_only_segment_adapter",
            "scope": neutrino["scope"],
            "promotable": False,
            "source_artifact": _repo_ref(NEUTRINO_JSON),
            "note": "Exact compare-only neutrino mass from the two-parameter positive-segment adapter.",
        },
        {
            "particle_id": "muon_neutrino",
            "label": "Muon Neutrino",
            "mass_eV": neutrino["exact_outputs"]["masses_eV"][1],
            "exact_kind": "exact_two_observable_compare_only_segment_adapter",
            "scope": neutrino["scope"],
            "promotable": False,
            "source_artifact": _repo_ref(NEUTRINO_JSON),
            "note": "Exact compare-only neutrino mass from the two-parameter positive-segment adapter.",
        },
        {
            "particle_id": "tau_neutrino",
            "label": "Tau Neutrino",
            "mass_eV": neutrino["exact_outputs"]["masses_eV"][2],
            "exact_kind": "exact_two_observable_compare_only_segment_adapter",
            "scope": neutrino["scope"],
            "promotable": False,
            "source_artifact": _repo_ref(NEUTRINO_JSON),
            "note": "Exact compare-only neutrino mass from the two-parameter positive-segment adapter.",
        },
    ]


def build_markdown(generated_utc: str, entries: list[dict[str, Any]]) -> str:
    lines = [
        "# Exact Non-Hadron Masses",
        "",
        f"Generated: `{generated_utc}`",
        "",
        "This bundle gives one exact mass output for every non-hadron particle currently covered by the OPH particle stack.",
        "It closes the exact-output lane, not the theorem-grade derivation lane.",
        "",
        "| Particle | Exact Mass | Kind | Scope | Source |",
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
    return "\n".join(lines).rstrip()


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the exact non-hadron mass bundle.")
    parser.add_argument("--markdown-out", default=str(DEFAULT_MD_OUT))
    parser.add_argument("--json-out", default=str(DEFAULT_JSON_OUT))
    parser.add_argument("--forward-out", default=str(DEFAULT_FORWARD_OUT))
    args = parser.parse_args()

    generated_utc = _timestamp()
    entries = build_entries()
    payload = {
        "artifact": "oph_exact_nonhadron_mass_bundle",
        "generated_utc": generated_utc,
        "status": "exact_output_lane_closed_nonhadron_only",
        "entries": entries,
        "excluded_lane": "hadrons_compute_bound",
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
