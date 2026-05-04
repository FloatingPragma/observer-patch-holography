#!/usr/bin/env python3
"""Build the particle-derivation gap ledger.

The ledger separates the new compressed P-trunk from the remaining theorem,
codepath, and execution gaps in the particle program. It is intentionally
static and explicit: changing a status should be a conscious edit, not an
accidental side effect of a numeric rebuild.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
PARTICLES_ROOT = ROOT / "particles"
P_TRUNK = ROOT / "P_derivation" / "runtime" / "p_closure_trunk_current.json"
DEFAULT_JSON_OUT = PARTICLES_ROOT / "runs" / "status" / "particle_derivation_gap_ledger.json"
DEFAULT_MD_OUT = PARTICLES_ROOT / "DERIVATION_GAP_LEDGER.md"


def _now_utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_p_trunk_summary() -> dict[str, Any]:
    if not P_TRUNK.exists():
        return {
            "artifact_path": str(P_TRUNK.relative_to(ROOT)),
            "exists": False,
            "claim_status": "not_emitted",
            "may_feed_live_particle_predictions": False,
        }
    payload = json.loads(P_TRUNK.read_text(encoding="utf-8"))
    return {
        "artifact_path": str(P_TRUNK.relative_to(ROOT)),
        "exists": True,
        "claim_status": payload.get("claim_status"),
        "P": payload.get("fixed_point_candidate", {}).get("P"),
        "alpha_inv": payload.get("fixed_point_candidate", {}).get("alpha_inv"),
        "source_report_mode": payload.get("source_report_mode"),
        "may_feed_live_particle_predictions": payload.get("consumer_policy", {}).get(
            "may_feed_live_particle_predictions",
            False,
        ),
    }


def build_gap_rows() -> list[dict[str, Any]]:
    return [
        {
            "id": "pclosure.compressed-trunk-artifact",
            "lane": "P closure",
            "status": "candidate_artifact",
            "github_issue": 224,
            "title": "Make the five-equation P trunk the canonical audit artifact",
            "current_boundary": (
                "The compressed trunk can organize the code path, but it is not a certified live "
                "particle root until the endpoint and interval-certificate gates close."
            ),
            "next_action": "Keep emitting p_closure_trunk_current.json and use it only as audit metadata for now.",
            "target_surfaces": ["code/P_derivation", "code/particles"],
        },
        {
            "id": "d10.ward-projected-thomson-endpoint",
            "lane": "D10 electromagnetic endpoint",
            "status": "open_theorem_gap",
            "github_issue": 223,
            "title": "Close the Ward-projected U(1)_Q Thomson endpoint",
            "current_boundary": (
                "The current structured-running law is a continuation: lepton one-loop transport is "
                "implemented, while hadronic spectral transport, scheme matching, and interval error "
                "control are not theorem-grade."
            ),
            "next_action": (
                "Derive Delta_Th(P) from the same Ward-projected source family as a0(P), including "
                "rho_had(s;P), matching remainder, and certified quadrature bounds."
            ),
            "target_surfaces": ["code/P_derivation/THOMSON_TRANSPORT_THEOREMS.md", "code/particles/calibration"],
        },
        {
            "id": "d10.rg-matching-threshold-scheme",
            "lane": "D10 running and matching",
            "status": "open_theorem_gap",
            "github_issue": 32,
            "title": "Internalize RG matching, threshold placement, and scheme conversion",
            "current_boundary": (
                "The D10 branch still uses declared running/matching conventions. The compressed P "
                "trunk makes this easier to see but does not prove those conventions from OPH."
            ),
            "next_action": (
                "Turn the running/matching package into an OPH edge-sector theorem or explicitly "
                "keep it as a declared convention in every prediction surface."
            ),
            "target_surfaces": ["paper compact D10 section", "code/P_derivation", "code/particles/calibration"],
        },
        {
            "id": "pclosure.live-codepath-adoption",
            "lane": "P closure",
            "status": "blocked_pending_certified_root",
            "github_issue": 224,
            "title": "Replace live particle P consumers with the certified P root",
            "current_boundary": (
                "Particle status surfaces still expose legacy metadata P=1.63094. This is acceptable "
                "only while the compressed trunk remains candidate/audit metadata."
            ),
            "next_action": (
                "After issues 223 and 32 close, switch live particle builders to the certified trunk "
                "artifact and make compare-only or historical P paths non-default."
            ),
            "target_surfaces": ["code/particles/scripts", "code/particles/runs/status", "WebProjects OPH summaries"],
        },
        {
            "id": "charged.determinant-normalization-transport",
            "lane": "Charged leptons",
            "status": "open_source_theorem",
            "github_issue": None,
            "title": "Close the P-to-charged affine anchor bridge",
            "current_boundary": (
                "No public charged-lepton values are emitted on the theorem lane. The open object is "
                "the determinant-normalization / sector-isolated trace-lift identity beneath A_ch(P)."
            ),
            "next_action": (
                "Prove 3 mu(r) = sum_e M_e^ch log q_e(r), equivalently zero normalization defect "
                "N_det(P), on the physical charged branch."
            ),
            "target_surfaces": ["code/particles/leptons", "code/particles/runs/leptons"],
        },
        {
            "id": "quark.selected-class-vs-global-classification",
            "lane": "Quarks",
            "status": "selected_class_closed_global_classification_open",
            "github_issue": 198,
            "title": "Keep exact quark rows scoped to the selected public frame class",
            "current_boundary": (
                "The selected-class exact Yukawa theorem is live, but it is not a global "
                "classification of all quark frame classes."
            ),
            "next_action": (
                "Either prove the global frame-class classification or keep every public exact-quark "
                "claim explicitly selected-class."
            ),
            "target_surfaces": ["code/particles/flavor", "particle paper quark section"],
        },
        {
            "id": "neutrino.pmns-status-and-absolute-rows",
            "lane": "Neutrinos",
            "status": "theorem_rows_with_visible_comparison_tension",
            "github_issue": 117,
            "title": "Keep neutrino theorem rows and comparison tension separated",
            "current_boundary": (
                "The weighted-cycle absolute splitting lane is emitted, while PMNS angle/phase "
                "comparisons still show visible displacement from current central values."
            ),
            "next_action": (
                "Do not hide PMNS residuals behind the exact absolute-splitting rows; either prove a "
                "better branch or leave the comparison tension explicit."
            ),
            "target_surfaces": ["code/particles/neutrino", "RESULTS_STATUS.md"],
        },
        {
            "id": "hadron.production-backend-systematics",
            "lane": "Hadrons",
            "status": "deferred_execution_contract",
            "github_issue": 153,
            "title": "Execute the production hadron backend and publish systematics",
            "current_boundary": (
                "The hadron lane is execution-contract-frozen. Symbolic simplification of P does not "
                "replace production QCD backend output and declared systematics."
            ),
            "next_action": (
                "Implement/run the production backend dump path, validate writeback, and publish "
                "continuum/volume/chiral/statistical budgets."
            ),
            "target_surfaces": ["code/particles/hadron", "code/particles/qcd"],
        },
    ]


def build_bundles() -> list[dict[str, Any]]:
    return [
        {
            "id": "electroweak-root-closure-bundle",
            "status": "constructive_endpoint_contract_emitted",
            "gap_ids": [
                "pclosure.compressed-trunk-artifact",
                "d10.ward-projected-thomson-endpoint",
                "d10.rg-matching-threshold-scheme",
                "pclosure.live-codepath-adoption",
            ],
            "promotion_question": (
                "Can one source-emitted map Delta_Th(P), with declared matching and interval bounds, "
                "certify the compressed P trunk as the live particle root without importing alpha(0)?"
            ),
            "result": (
                "Constructive result. The admissible endpoint object is now explicit: Delta_Th(P) must split into "
                "source lepton transport, a Ward-projected hadronic spectral density rho_had(s;P), "
                "a certified electroweak/scheme remainder, RG/matching certificates, quadrature bounds, "
                "and an interval-level fixed-point certificate. The local implementation target is "
                "P_derivation/runtime/thomson_endpoint_contract_current.json."
            ),
        },
        {
            "id": "spectrum-source-bundle",
            "status": "constructive_trace_lift_schema_emitted",
            "gap_ids": [
                "charged.determinant-normalization-transport",
                "quark.selected-class-vs-global-classification",
                "neutrino.pmns-status-and-absolute-rows",
            ],
            "promotion_question": (
                "Is there one OPH excitation dictionary and sector-isolated trace-lift theorem that "
                "explains the charged affine anchor, quark selected-class boundary, and neutrino PMNS "
                "comparison surface without hidden target fitting?"
            ),
            "result": (
                "Constructive result. SourceNormalizedTraceLiftDescent is the reusable schema to implement next. Charged leptons "
                "still need N_det(P)=0, quarks remain selected-class on f_P, and neutrino PMNS rows "
                "remain visible comparison-tension rows."
            ),
        },
        {
            "id": "qcd-thomson-backend-bundle",
            "status": "constructive_spectral_measure_contract_emitted",
            "gap_ids": [
                "d10.ward-projected-thomson-endpoint",
                "hadron.production-backend-systematics",
            ],
            "promotion_question": (
                "Can the hadron production backend emit the rho_had(s;P) object and uncertainty budget "
                "needed by the Ward-projected Thomson endpoint, rather than leaving hadrons and alpha(0) "
                "as separate deferred gaps?"
            ),
            "result": (
                "Constructive result. The current stable-channel backend is not the endpoint object. The smallest "
                "missing primitive is production_ward_projected_hadronic_spectral_measure_export with "
                "continuum, volume, chiral, statistical, matching, quadrature, and endpoint budgets. The "
                "local implementation target is particles/hadron/ward_projected_spectral_measure.schema.json."
            ),
        },
        {
            "id": "particle-root-integration-gate",
            "status": "keep_candidate_with_constructive_next_artifacts",
            "gap_ids": [
                "pclosure.compressed-trunk-artifact",
                "d10.ward-projected-thomson-endpoint",
                "d10.rg-matching-threshold-scheme",
                "pclosure.live-codepath-adoption",
                "charged.determinant-normalization-transport",
                "quark.selected-class-vs-global-classification",
                "neutrino.pmns-status-and-absolute-rows",
                "hadron.production-backend-systematics",
            ],
            "promotion_question": (
                "Do the returned packets jointly close the endpoint, matching, interval, and source-object "
                "requirements strongly enough to promote the compressed trunk into live particle builders?"
            ),
            "result": (
                "No promotion. The first wave now emits constructive next artifacts, so the compressed P trunk "
                "remains candidate/audit metadata until those artifacts are populated and certified."
            ),
        },
    ]


def build_ledger() -> dict[str, Any]:
    rows = build_gap_rows()
    return {
        "artifact": "oph_particle_derivation_gap_ledger",
        "generated_utc": _now_utc(),
        "purpose": "Systematic claim-safe queue after the five-equation P-trunk simplification.",
        "p_trunk": _load_p_trunk_summary(),
        "bundles": build_bundles(),
        "rows": rows,
        "promotion_policy": {
            "compressed_p_trunk_is_live_prediction_root": False,
            "reason": "The endpoint, RG/matching, and interval-certificate gates remain open.",
            "torus_mode_language_allowed_in_pipeline": False,
            "address_remaining_blockers_one_by_one": False,
            "obstruction_only_worker_result_allowed": False,
        },
    }


def render_markdown(ledger: dict[str, Any]) -> str:
    p_trunk = ledger["p_trunk"]
    lines = [
        "# Particle Derivation Gap Ledger",
        "",
        f"Generated: `{ledger['generated_utc']}`",
        "",
        ledger["purpose"],
        "",
        "## P-Trunk Status",
        "",
        f"- Artifact: `{p_trunk['artifact_path']}`",
        f"- Exists: `{p_trunk['exists']}`",
        f"- Claim status: `{p_trunk['claim_status']}`",
        f"- May feed live particle predictions: `{p_trunk['may_feed_live_particle_predictions']}`",
    ]
    if p_trunk.get("P") is not None:
        lines.extend(
            [
                f"- Candidate P: `{p_trunk['P']}`",
                f"- Candidate alpha^-1: `{p_trunk['alpha_inv']}`",
                f"- Source report mode: `{p_trunk['source_report_mode']}`",
            ]
        )
    lines.extend(
        [
            "",
            "## Bundle Execution Plan",
            "",
            "The remaining work is grouped into coupled closure packets rather than a one-blocker-at-a-time queue.",
            "",
            "| Bundle | Status | Gaps | Promotion question |",
            "| --- | --- | --- | --- |",
        ]
    )
    for bundle in ledger["bundles"]:
        gaps = ", ".join(f"`{gap}`" for gap in bundle["gap_ids"])
        lines.append(
            f"| `{bundle['id']}` | `{bundle['status']}` | {gaps} | {bundle['promotion_question']} |"
        )
    lines.extend(
        [
            "",
            "## Bundle Packet Results",
            "",
        ]
    )
    for bundle in ledger["bundles"]:
        lines.append(f"- `{bundle['id']}`: `{bundle['status']}`. {bundle['result']}")
    lines.extend(
        [
            "",
            "## Remaining Gaps",
            "",
            "| ID | Lane | Status | Issue | Next action |",
            "| --- | --- | --- | --- | --- |",
        ]
    )
    for row in ledger["rows"]:
        issue = f"#{row['github_issue']}" if row.get("github_issue") is not None else "n/a"
        lines.append(
            f"| `{row['id']}` | {row['lane']} | `{row['status']}` | {issue} | {row['next_action']} |"
        )
    lines.extend(
        [
            "",
            "## Claim Policy",
            "",
            "- The compressed P trunk is an audit/candidate artifact until the endpoint and certificate gates close.",
            "- The remaining blockers should be worked as coupled bundles, not as isolated one-off fixes.",
            "- The particle pipeline must keep compare-only, continuation, selected-class, and theorem-grade rows mechanically distinct.",
            "- Golden-ratio torus or resonance language is not a live derivation input unless a separate representation-to-spectrum theorem is supplied.",
        ]
    )
    return "\n".join(lines)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build the particle derivation gap ledger.")
    parser.add_argument("--json-out", default=str(DEFAULT_JSON_OUT))
    parser.add_argument("--markdown-out", default=str(DEFAULT_MD_OUT))
    parser.add_argument("--print-json", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    ledger = build_ledger()
    json_text = json.dumps(ledger, indent=2, sort_keys=True) + "\n"

    json_out = Path(args.json_out)
    json_out.parent.mkdir(parents=True, exist_ok=True)
    json_out.write_text(json_text, encoding="utf-8")

    markdown_out = Path(args.markdown_out)
    markdown_out.write_text(render_markdown(ledger) + "\n", encoding="utf-8")

    if args.print_json:
        print(json_text, end="")
    else:
        print(f"saved: {json_out}")
        print(f"saved: {markdown_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
