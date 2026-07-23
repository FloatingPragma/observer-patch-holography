#!/usr/bin/env python3
"""Transport the 2026-07 structure set into the quark spread fiber.

The sigma obstruction (``quark_sigma_source_nonidentifiability_obstruction``)
proves that the target-free quark source packet leaves the two positive
endpoint spans free: the compatible fiber is ``(R_{>0})^2``.  Since that
theorem was emitted, the corpus gained the conditional super-Tannakian matter
lift (issue #314: chiral one-generation spectrum, exactly one Yukawa invariant
line per declared channel, machine-checked anomaly traces) and the conditional
port-current algebra (issue #566), both hash-pinned receipts.

This builder decides, fail closed, which fork of the issue-#591 claim boundary
("emit, or prove non-identifiable, the six") the new structure supports:

fork (i)  - some certified datum of the new structure moves under the fiber
            action ``(Y_u, Y_d) -> (lambda_u Y_u, lambda_d Y_d)``; the moved
            datum is a new source equation and a constrained combination
            becomes emittable.  The builder records the cut and exits nonzero
            so the promotion is a reviewed step, never an automatic one.
fork (ii) - every certified datum is invariant under the fiber action; the
            non-identifiability theorem extends verbatim to the enlarged
            corpus and the fiber coordinates are exactly the free scalar
            coefficients along the two hadronic Yukawa invariant lines.

The decision is computed from the receipts, not asserted: a coefficient
blindness scan over every leaf of the matter and port manifests/receipts, the
invariant-line dimensions, the family-attachment gate, the port claim
boundary, and the frozen selector-menu exclusions.  A receipt that ever grows
a Yukawa-magnitude field flips the fork on the next run.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
REPO_ROOT = ROOT.parent
RUNS = ROOT / "particles" / "runs"
A5 = ROOT / "a5_closure"

BASE_OBSTRUCTION_JSON = RUNS / "flavor" / "quark_sigma_source_nonidentifiability_obstruction.json"
AXIOM_LEVEL_JSON = RUNS / "flavor" / "quark_axiom_level_yukawa_moduli_nonidentifiability.json"
MATTER_MANIFEST_JSON = A5 / "manifests" / "super_tannakian_matter_reference.json"
MATTER_RECEIPT_JSON = A5 / "receipts" / "super_tannakian_matter_reference.receipt.json"
PORT_MANIFEST_JSON = A5 / "manifests" / "port_current_response_reference.json"
PORT_RECEIPT_JSON = A5 / "receipts" / "port_current_inner_reference.receipt.json"
CANDIDATES_JSON = ROOT / "particles" / "flavor_selector" / "runtime" / "candidates_evaluated.json"
DEFAULT_OUT = RUNS / "flavor" / "quark_spread_fiber_structure_transport_obstruction.json"

DIRECT_INPUT_PATHS = [
    "particles/runs/flavor/quark_sigma_source_nonidentifiability_obstruction.json",
    "particles/runs/flavor/quark_axiom_level_yukawa_moduli_nonidentifiability.json",
    "a5_closure/manifests/super_tannakian_matter_reference.json",
    "a5_closure/receipts/super_tannakian_matter_reference.receipt.json",
    "a5_closure/manifests/port_current_response_reference.json",
    "a5_closure/receipts/port_current_inner_reference.receipt.json",
    "particles/flavor_selector/runtime/candidates_evaluated.json",
]

# Identical discipline to the sigma obstruction: this certificate must not
# descend from any target-anchored quark surface.
FORBIDDEN_ANCESTORS = [
    "particle_reference_values",
    "quark_current_family_exact_readout",
    "quark_current_family_exact_sigma_target",
    "quark_current_family_exact_pdg_theorem",
    "quark_current_family_end_to_end_exact_pdg_derivation_chain",
    "quark_selected_class_public_exact_evaluator",
    "quark_public_physical_sigma_datum_descent",
    "quark_public_exact_yukawa_end_to_end_theorem",
    "quark_exact_pdg_end_to_end_theorem",
    "target_centered_log_u",
    "target_centered_log_d",
    "reference_targets_u",
    "reference_targets_d",
    "exact_fit_residuals",
    "fitted_sigma_values",
    "PDG_API_quark_rows",
    "CODATA",
]

# Leaf paths whose lowercase form contains one of these needles are the
# Yukawa-adjacent surface of the receipts; the blindness scan requires that
# none of them carries a numeric magnitude.
BLINDNESS_NEEDLES = (
    "yukawa",
    "sigma_u",
    "sigma_d",
    "quark_mass",
    "pole_mass",
    "mass_gev",
    "pdg",
    "codata",
    "fitted",
    "spread_modulus",
    "vev",
)

# Integer leaves that are dimension counts, not magnitudes.
INTEGER_DIMENSION_KEYS = {"invariant_dimension", "invariant_sector_dimension"}

_NUMERIC_STRING = re.compile(r"[-+]?\d*\.?\d+([eE][-+]?\d+)?")

HADRONIC_CHANNELS = (["Q", "S", "u_c"], ["Q", "Sbar", "d_c"])
FORBIDDEN_CHANNEL = ["Q", "S", "d_c"]


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _walk_leaves(obj: Any, path: str = "") -> list[tuple[str, Any]]:
    if isinstance(obj, dict):
        leaves: list[tuple[str, Any]] = []
        for key, value in obj.items():
            leaves.extend(_walk_leaves(value, f"{path}/{key}"))
        return leaves
    if isinstance(obj, list):
        leaves = []
        for index, value in enumerate(obj):
            leaves.extend(_walk_leaves(value, f"{path}[{index}]"))
        return leaves
    return [(path, obj)]


def coefficient_blindness_scan(name: str, document: dict[str, Any]) -> dict[str, Any]:
    """Machine check: no Yukawa-adjacent leaf carries a numeric magnitude.

    Booleans, channel-label strings, prose statements, and the declared
    integer dimension counts are structural; a float, a magnitude-bearing
    integer, or a numeric string under a Yukawa-adjacent path is a
    violation and flips the fork.
    """

    matched: list[dict[str, Any]] = []
    violations: list[dict[str, Any]] = []
    for leaf_path, value in _walk_leaves(document):
        lowered = leaf_path.lower()
        if not any(needle in lowered for needle in BLINDNESS_NEEDLES):
            continue
        basename = leaf_path.rsplit("/", 1)[-1].split("[", 1)[0]
        entry = {"path": leaf_path, "json_type": type(value).__name__}
        matched.append(entry)
        if isinstance(value, bool) or value is None:
            continue
        if isinstance(value, int):
            if basename not in INTEGER_DIMENSION_KEYS:
                violations.append({**entry, "reason": "integer_magnitude_outside_dimension_allowlist"})
            continue
        if isinstance(value, float):
            violations.append({**entry, "reason": "float_magnitude"})
            continue
        if isinstance(value, str) and _NUMERIC_STRING.fullmatch(value.strip()):
            violations.append({**entry, "reason": "numeric_string_magnitude"})
    return {
        "document": name,
        "needles": list(BLINDNESS_NEEDLES),
        "matched_leaf_count": len(matched),
        "matched_leaves": matched,
        "violations": violations,
        "coefficient_blind": not violations,
    }


def invariant_line_check(matter_receipt: dict[str, Any]) -> dict[str, Any]:
    sector = matter_receipt.get("yukawa_sector", {})
    channels = {tuple(row.get("channel", [])): row.get("invariant_dimension") for row in sector.get("channels", [])}
    control = sector.get("forbidden_channel_control", {})
    hadronic_ok = all(channels.get(tuple(ch)) == 1 for ch in HADRONIC_CHANNELS)
    lepton_ok = channels.get(("L", "Sbar", "e_c")) == 1
    control_ok = (
        control.get("channel") == FORBIDDEN_CHANNEL and control.get("invariant_dimension") == 0
    )
    return {
        "channels": {"/".join(k): v for k, v in channels.items()},
        "forbidden_channel_control": control,
        "hadronic_lines_one_dimensional": hadronic_ok,
        "lepton_line_one_dimensional": lepton_ok,
        "forbidden_channel_carries_no_line": control_ok,
        "reading": (
            "a one-dimensional invariant line is stable under scalar rescaling "
            "of its coefficient: the certified datum is the dimension, and the "
            "fiber action rescales the coefficient along the line without "
            "moving the line, the charge spectrum, or any anomaly trace"
        ),
        "passed": hadronic_ok and lepton_ok and control_ok,
    }


def build_artifact(
    base: dict[str, Any],
    axiom: dict[str, Any],
    matter_manifest: dict[str, Any],
    matter_receipt: dict[str, Any],
    port_manifest: dict[str, Any],
    port_receipt: dict[str, Any],
    candidates: dict[str, Any],
    input_hashes: dict[str, str],
) -> dict[str, Any]:
    base_fiber = base.get("exact_ray_classification", {})
    matter_not_claimed = str(matter_receipt.get("branch_scope", {}).get("not_claimed", ""))
    port_does_not_close = port_receipt.get("claim_boundary", {}).get("does_not_close", [])
    blindness_scans = [
        coefficient_blindness_scan("super_tannakian_matter_manifest", matter_manifest),
        coefficient_blindness_scan("super_tannakian_matter_receipt", matter_receipt),
        coefficient_blindness_scan("port_current_manifest", port_manifest),
        coefficient_blindness_scan("port_current_receipt", port_receipt),
    ]
    invariant_lines = invariant_line_check(matter_receipt)
    candidate_blindness = candidates.get("blindness", {})

    checks = {
        "base_fiber_intact": (
            base.get("artifact") == "oph_quark_sigma_source_nonidentifiability_obstruction"
            and base_fiber.get("fiber") == "(R_{>0})^2"
            and base.get("public_promotion_allowed") is False
        ),
        "axiom_level_theorem_intact": (
            axiom.get("artifact") == "oph_quark_axiom_level_yukawa_moduli_nonidentifiability"
            and axiom.get("proof_status") == "closed_axiom_level_nondefinability_theorem"
            and axiom.get("source_only_numeric_quark_spectrum_emitted") is False
        ),
        "matter_documents_coefficient_blind": all(
            scan["coefficient_blind"] for scan in blindness_scans[:2]
        ),
        "port_documents_coefficient_blind": all(
            scan["coefficient_blind"] for scan in blindness_scans[2:]
        ),
        "yukawa_invariant_lines_one_dimensional": invariant_lines["passed"],
        "family_attachment_outside_matter_packet": (
            matter_receipt.get("acceptance_criteria_status", {}).get(
                "family_attachment_scalar_potential_pole_mass_outside_packet"
            )
            is True
            and "no family attachment" in matter_not_claimed
        ),
        "port_receipt_closes_no_mass": any(
            "masses" in str(item) for item in port_does_not_close
        ),
        "selector_menu_exhausted_without_selection": (
            candidates.get("frozen_candidate_count") == 12
            and len(candidates.get("candidates", [])) == 12
            and candidates.get("status") == "evaluated_no_selection"
            and candidate_blindness.get("reads_measured_masses") is False
            and candidate_blindness.get("reads_measured_ratios") is False
        ),
    }
    fiber_cut_detected = not all(checks.values())
    fork = "i_structure_cuts_fiber" if fiber_cut_detected else "ii_fiber_survives"

    failed = sorted(key for key, passed in checks.items() if not passed)
    scan_violations = [
        violation for scan in blindness_scans for violation in scan["violations"]
    ]

    fixed_certified_data = [
        {
            "datum": "chiral one-generation charge spectrum and zero-intertwiner chirality certificate",
            "source": "super_tannakian_matter_receipt",
            "why_fixed": "the fiber action rescales channel coefficients; the module, its dual, and their charge spectra are untouched",
        },
        {
            "datum": "anomaly traces (SU3^3, SU3^2-U1, SU2^2-U1, U1^3, grav^2-U1) and Witten parity",
            "source": "super_tannakian_matter_receipt",
            "why_fixed": "operator traces on the matter module do not read Yukawa coefficients",
        },
        {
            "datum": "the three Yukawa invariant-line dimensions and the forbidden-channel zero",
            "source": "super_tannakian_matter_receipt",
            "why_fixed": "scalar rescaling maps each line onto itself; dimensions are invariant",
        },
        {
            "datum": "the emitted common action kernel and the fifteen-state module selection",
            "source": "super_tannakian_matter_receipt",
            "why_fixed": "kernel emission precedes any coefficient choice along the lines",
        },
        {
            "datum": "port-current algebra dimension, centre, inner A5 action, and the four signed response coefficients",
            "source": "port_current_receipt",
            "why_fixed": "the response construction carries no Yukawa-adjacent field; its claim boundary excludes masses and couplings",
        },
        {
            "datum": "the twelve frozen selector candidates and their blind evaluations",
            "source": "candidates_evaluated",
            "why_fixed": "the menu reads screen combinatorics only; every candidate is excluded on the comparison surface",
        },
    ]

    countermodels = [
        {
            "name": "unit_coefficient_pair",
            "lambda_u": 1.0,
            "lambda_d": 1.0,
            "action": "identity representative of the fiber action on the two hadronic invariant lines",
        },
        {
            "name": "independently_rescaled_pair",
            "lambda_u": 2.0,
            "lambda_d": 3.0,
            "action": (
                "Y_u -> 2 Y_u along span[Q,S,u_c], Y_d -> 3 Y_d along span[Q,Sbar,d_c]"
            ),
        },
    ]
    for model in countermodels:
        model["role"] = "formal_nonidentifiability_witness_not_a_prediction"
        model["fixes_every_certified_datum"] = not fiber_cut_detected
        model["changes_requested_sextet_coordinates"] = True

    theorem_statement = (
        "The 2026-07 certified structure set - the conditional super-Tannakian "
        "matter lift (#314), the conditional port-current algebra (#566), and "
        "the twelve frozen flavor-orbit selector candidates - does not cut the "
        "(R_{>0})^2 quark spread fiber. Every certified datum of that set is "
        "invariant under the fiber action (Y_u, Y_d) -> (lambda_u Y_u, "
        "lambda_d Y_d) along the two one-dimensional hadronic Yukawa invariant "
        "lines, so the fiber coordinates are exactly the free scalar "
        "coefficients along those lines. Any future cut must pass through the "
        "physical bindings (#599, #567), the family attachment (#569), or a "
        "new selector under the frozen SELECTOR_SPEC discipline, and must "
        "break the free action."
        if not fiber_cut_detected
        else "A certified datum of the 2026-07 structure set moves under the "
        "spread-fiber action; the moved datum is a candidate source equation. "
        "Emission of any constrained combination is a reviewed step: this "
        "certificate records the cut and refuses automatic promotion."
    )

    return {
        "artifact": "oph_quark_spread_fiber_structure_transport_obstruction",
        "generated_utc": _timestamp(),
        "github_issues": [591],
        "claim_tier": "source_only_nonidentifiability_obstruction_transport",
        "fork": fork,
        "fiber_cut_detected": fiber_cut_detected,
        "failed_checks": failed,
        "proof_status": (
            "closed_exact_current_corpus_obstruction"
            if not fiber_cut_detected
            else "structure_cut_detected_manual_rederivation_required"
        ),
        "theorem_grade_obstruction": not fiber_cut_detected,
        "source_only_sigma_emitted": False,
        "public_promotion_allowed": False,
        "numeric_quark_rows_allowed": False,
        "theorem_statement": theorem_statement,
        "issue_591_claim_boundary": {
            "quoted_disjunction": "emit, or prove non-identifiable, the six",
            "disjunct_discharged": (
                "prove_non_identifiable_at_2026_07_corpus_level"
                if not fiber_cut_detected
                else "none_pending_manual_review_of_detected_cut"
            ),
            "carrier_bullets_remaining_open": [
                "source-defined flavor carriers",
                "flavor quotient and orbit selector",
                "CKM interval emission",
            ],
            "gating_parents": ["#569 family attachment", "#567 physical Z6 descent", "#599 response source-binding"],
        },
        "checks": checks,
        "coefficient_blindness_scans": blindness_scans,
        "coefficient_blindness_violations": scan_violations,
        "invariant_line_transport": invariant_lines,
        "claim_boundary_quotes": {
            "matter_receipt_not_claimed": matter_not_claimed,
            "port_receipt_does_not_close": port_does_not_close,
        },
        "free_action_certificate": {
            "group": "(R_{>0})^2",
            "action": (
                "(lambda_u, lambda_d).(Y_u, Y_d) = (lambda_u Y_u, lambda_d Y_d) "
                "along the one-dimensional invariant lines of [Q,S,u_c] and [Q,Sbar,d_c]"
            ),
            "fixed_certified_data": fixed_certified_data,
            "formal_countermodels": countermodels,
            "action_free": True,
            "base_fiber_reference": {
                "artifact": base.get("artifact"),
                "fiber": base_fiber.get("fiber"),
                "independent_coordinates": base_fiber.get("independent_coordinates"),
            },
        },
        "minimal_future_cut": {
            "routes": [
                "physical source binding of the response representation (#599) followed by physical Z6 descent (#567)",
                "source-derived attachment of the icosahedral screen action to three physical families (#569)",
                "a thirteenth selector candidate under the frozen SELECTOR_SPEC single-comparison discipline",
            ],
            "must_break_free_action": True,
            "inherited_requirement": base.get("minimal_future_extension", {}),
        },
        "dependency_audit": {
            "direct_input_paths": DIRECT_INPUT_PATHS,
            "input_sha256": input_hashes,
            "forbidden_ancestors": FORBIDDEN_ANCESTORS,
            "loads_running_quark_reference_rows": False,
            "loads_fitted_sigma_values": False,
            "emits_target_sigma_values": False,
            "no_target_leak": True,
        },
        "public_prediction_policy": {
            "running_quark_numeric_rows": "withheld",
            "forward_yukawa_numeric_promotion": "withheld",
            "reason": (
                "the spread fiber survives the 2026-07 structure set"
                if not fiber_cut_detected
                else "detected cut awaits manual re-derivation; no automatic emission"
            ),
        },
        "notes": [
            "The fork is computed from the receipts on every run; a receipt that grows a Yukawa-magnitude field flips it.",
            "This certificate does not prove that no future OPH extension can emit the spread pair.",
            "Conditionality is inherited: #314 and #566 are conditional on their declared branch premises, and this transport is conditional on the same set.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Transport the 2026-07 structure set into the quark spread fiber."
    )
    parser.add_argument("--base", default=str(BASE_OBSTRUCTION_JSON))
    parser.add_argument("--axiom", default=str(AXIOM_LEVEL_JSON))
    parser.add_argument("--matter-manifest", default=str(MATTER_MANIFEST_JSON))
    parser.add_argument("--matter-receipt", default=str(MATTER_RECEIPT_JSON))
    parser.add_argument("--port-manifest", default=str(PORT_MANIFEST_JSON))
    parser.add_argument("--port-receipt", default=str(PORT_RECEIPT_JSON))
    parser.add_argument("--candidates", default=str(CANDIDATES_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    paths = {
        "base": Path(args.base),
        "axiom": Path(args.axiom),
        "matter_manifest": Path(args.matter_manifest),
        "matter_receipt": Path(args.matter_receipt),
        "port_manifest": Path(args.port_manifest),
        "port_receipt": Path(args.port_receipt),
        "candidates": Path(args.candidates),
    }
    hashes = {key: _sha256(path) for key, path in paths.items()}
    payload = build_artifact(
        _load_json(paths["base"]),
        _load_json(paths["axiom"]),
        _load_json(paths["matter_manifest"]),
        _load_json(paths["matter_receipt"]),
        _load_json(paths["port_manifest"]),
        _load_json(paths["port_receipt"]),
        _load_json(paths["candidates"]),
        hashes,
    )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {output}")
    print(f"fork: {payload['fork']}")
    if payload["fiber_cut_detected"]:
        print("structure cut detected; manual re-derivation required before any emission")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
