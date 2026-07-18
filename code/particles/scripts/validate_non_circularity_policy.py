#!/usr/bin/env python3
"""Validate strict non-circularity promotion boundaries for particle artifacts."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
PARTICLES = ROOT / "particles"

ARTIFACTS = {
    "d10_wz_factorization": PARTICLES
    / "runs"
    / "calibration"
    / "d10_ew_w_anchor_neutral_shear_factorization_official_pdg_2025_update.json",
    "d10_readout": PARTICLES / "runs" / "calibration" / "d10_ew_source_transport_readout.json",
    "d11_split": PARTICLES / "runs" / "calibration" / "d11_live_exact_split_pair_theorem.json",
    "quark_sigma": PARTICLES / "runs" / "flavor" / "quark_public_physical_sigma_datum_descent.json",
    "quark_sigma_required": PARTICLES
    / "runs"
    / "flavor"
    / "quark_sigma_source_datum_no_target_leak_required.json",
    "quark_yukawa": PARTICLES / "runs" / "flavor" / "quark_public_exact_yukawa_end_to_end_theorem.json",
    "neutrino_bridge": PARTICLES / "runs" / "neutrino" / "neutrino_bridge_rigidity_theorem.json",
    "neutrino_absolute": PARTICLES / "runs" / "neutrino" / "neutrino_absolute_attachment_theorem.json",
    "exact_nonhadron": PARTICLES / "exact_nonhadron_masses.json",
    "carrier_acceptance": PARTICLES / "runs" / "status" / "carrier_mode_acceptance.json",
}


def _load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _require(condition: bool, failures: list[str], message: str) -> None:
    if not condition:
        failures.append(message)


def validate() -> dict[str, Any]:
    payloads = {name: _load(path) for name, path in ARTIFACTS.items()}
    failures: list[str] = []

    wz = payloads["d10_wz_factorization"]
    _require(wz.get("public_surface_candidate_allowed") is False, failures, "W/Z reference-fitted adapter is public")
    _require(wz.get("prediction_promotion_allowed") is False, failures, "W/Z reference-fitted adapter is promotable")
    _require(wz.get("display_allowed_as_compare_only") is True, failures, "W/Z reference-fitted adapter lacks compare-only display flag")

    readout = payloads["d10_readout"]
    if readout.get("target_free_repair_value_law_status") != "closed":
        _require(readout.get("public_surface_candidate_allowed") is False, failures, "D10 readout promotes without closed target-free repair")
        _require(readout.get("display_allowed_as_compare_only") is True, failures, "D10 readout lacks compare-only display flag")

    d11 = payloads["d11_split"]
    gate = dict(d11.get("upstream_promotion_gate") or {})
    if gate.get("passed") is not True:
        _require(d11.get("status") == "candidate_only", failures, "D11 split is closed despite failed D10 repair gate")
        _require(d11.get("public_surface_candidate_allowed") is False, failures, "D11 split promotes despite failed D10 repair gate")
        _require(d11.get("prediction_promotion_allowed") is False, failures, "D11 split prediction promotion leaked")

    quark_sigma = payloads["quark_sigma"]
    sigma_nc = dict(quark_sigma.get("non_circularity_status") or {})
    if sigma_nc.get("target_derived_sigma_datum_used") is True:
        _require(quark_sigma.get("public_promotion_allowed") is False, failures, "target-derived quark sigma datum promotes")
        _require(sigma_nc.get("promotion_allowed") is False, failures, "target-derived quark sigma non-circularity gate promotes")
        _require(
            quark_sigma.get("source_only_sigma_emitted") is False,
            failures,
            "target-derived quark sigma is marked source-only",
        )
        _require(
            "QUARK_SOURCE_SPREAD_PAIR_ACTION_BREAKING_THEOREM"
            in quark_sigma.get("missing_for_promotion", []),
            failures,
            "quark source-spread action-breaking gate is not recorded",
        )

    quark_sigma_required = payloads["quark_sigma_required"]
    _require(
        quark_sigma_required.get("status") == "closed_current_corpus_nonidentifiability_obstruction",
        failures,
        "quark sigma required artifact does not record the closed current-corpus obstruction",
    )
    _require(
        quark_sigma_required.get("source_only_sigma_emitted") is False,
        failures,
        "quark sigma required artifact emits source sigma",
    )
    _require(
        "QUARK_SOURCE_SPREAD_PAIR_ACTION_BREAKING_THEOREM"
        in quark_sigma_required.get("missing_for_promotion", []),
        failures,
        "quark sigma required artifact lacks the action-breaking gate",
    )

    quark_yukawa = payloads["quark_yukawa"]
    yukawa_nc = dict(quark_yukawa.get("non_circularity_status") or {})
    if yukawa_nc.get("target_derived_sigma_datum_used") is True:
        _require(quark_yukawa.get("public_promotion_allowed") is False, failures, "target-derived quark Yukawa theorem promotes")
        _require(
            quark_yukawa.get("proof_status")
            in {
                "blocked_target_derived_sigma_source_missing",
                "blocked_inconsistent_exact_mass_yukawa_surfaces",
                "blocked_mixed_scheme_dimensionful_mass_textures_not_physical_yukawas",
            },
            failures,
            "target-derived quark Yukawa proof status is not fail-closed",
        )

    bridge = payloads["neutrino_bridge"]
    bridge_nc = dict(bridge.get("non_circularity_status") or {})
    if bridge_nc.get("compare_only_correction_audit_used") is True:
        _require(bridge.get("public_surface_candidate_allowed") is False, failures, "compare-only neutrino C_nu promotes")
        _require(bridge.get("emitted_value") is None, failures, "compare-only neutrino C_nu is emitted as theorem value")
        _require(bridge.get("display_value") is not None, failures, "compare-only neutrino C_nu lost display value")

    absolute = payloads["neutrino_absolute"]
    absolute_nc = dict(absolute.get("non_circularity_status") or {})
    if absolute_nc.get("compare_only_C_nu_used") is True:
        _require(absolute.get("public_surface_candidate_allowed") is False, failures, "compare-only neutrino absolute attachment promotes")
        _require(
            absolute.get("status") == "conditional_absolute_family_blocked_by_compare_only_C_nu",
            failures,
            "compare-only neutrino absolute attachment status is not blocked",
        )

    entries = {entry["particle_id"]: entry for entry in payloads["exact_nonhadron"]["entries"]}
    withheld_entries = {
        entry["particle_id"]: entry for entry in payloads["exact_nonhadron"].get("withheld_entries", [])
    }
    carrier_modes = {
        entry["carrier_id"]: entry for entry in payloads["carrier_acceptance"].get("carriers", [])
    }
    for particle_id in ("photon", "gluon", "graviton"):
        _require(particle_id not in entries, failures, f"{particle_id} leaked into public particle mass entries")
        _require(particle_id in carrier_modes, failures, f"{particle_id} conditional carrier mode is absent")
        _require(
            carrier_modes.get(particle_id, {}).get("hard_quadratic_mass_parameter_squared") == 0,
            failures,
            f"{particle_id} hard quadratic zero is absent",
        )
        _require(
            carrier_modes.get(particle_id, {}).get("abstract_symmetry_group_alone_sufficient") is False,
            failures,
            f"{particle_id} treats the abstract group as sufficient",
        )
        _require(
            carrier_modes.get(particle_id, {}).get("quantum_particle_gate", {}).get("passed") is False,
            failures,
            f"{particle_id} passes the quantum-particle gate without receipts",
        )
        _require(
            carrier_modes.get(particle_id, {}).get("particle_promotion_allowed") is False,
            failures,
            f"{particle_id} conditional carrier mode promotes as a particle",
        )
    for particle_id in ("up_quark", "down_quark", "strange_quark", "charm_quark", "bottom_quark", "top_quark"):
        _require(particle_id not in entries, failures, f"{particle_id} target-anchored quark witness is a public entry")
        _require(
            withheld_entries.get(particle_id, {}).get("exact_kind")
            == "selected_class_target_anchored_mixed_convention_mass_texture_audit",
            failures,
            f"{particle_id} target-anchored mixed-convention audit is not withheld",
        )
        _require(
            withheld_entries.get(particle_id, {}).get("promotable") is False,
            failures,
            f"{particle_id} target-anchored quark witness promotes",
        )
        _require(
            withheld_entries.get(particle_id, {}).get("claim_tier")
            == "mixed_scheme_dimensionful_mass_texture_audit",
            failures,
            f"{particle_id} quark witness has wrong claim tier",
        )
        _require(
            withheld_entries.get(particle_id, {}).get("source_only_sigma_emitted") is False,
            failures,
            f"{particle_id} quark witness is marked source-only sigma",
        )
        _require(
            "QUARK_SOURCE_SPREAD_PAIR_ACTION_BREAKING_THEOREM"
            in withheld_entries.get(particle_id, {}).get("missing_for_promotion", []),
            failures,
            f"{particle_id} missing quark source-spread action-breaking gate is not recorded",
        )
    for particle_id in ("electron", "muon", "tau"):
        _require(particle_id not in entries, failures, f"{particle_id} target-anchored charged witness is a public entry")
        _require(
            withheld_entries.get(particle_id, {}).get("exact_kind") == "exact_target_anchored_current_family_witness",
            failures,
            f"{particle_id} target-anchored charged witness is not withheld",
        )
        _require(
            withheld_entries.get(particle_id, {}).get("source_only") is False,
            failures,
            f"{particle_id} target-anchored charged witness is marked source-only",
        )
        _require(
            "charged_determinant_trace_lift_attachment"
            in withheld_entries.get(particle_id, {}).get("missing_for_promotion", []),
            failures,
            f"{particle_id} missing determinant trace-lift gate is not recorded",
        )
    for particle_id in ("electron_neutrino", "muon_neutrino", "tau_neutrino"):
        _require(particle_id not in entries, failures, f"{particle_id} compare-only absolute attachment is a public entry")
        _require(
            withheld_entries.get(particle_id, {}).get("exact_kind")
            == "rejected_target_informed_weighted_cycle_candidate",
            failures,
            f"{particle_id} rejected target-informed candidate is not withheld",
        )
        _require(
            withheld_entries.get(particle_id, {}).get("promotable") is False,
            failures,
            f"{particle_id} compare-only absolute attachment promotes",
        )

    return {
        "artifact": "oph_particle_non_circularity_policy_validation",
        "pass": not failures,
        "failures": failures,
        "checked_artifacts": {name: path.relative_to(ROOT).as_posix() for name, path in ARTIFACTS.items()},
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate particle non-circularity promotion policy.")
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON.")
    args = parser.parse_args()
    payload = validate()
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0 if payload["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
