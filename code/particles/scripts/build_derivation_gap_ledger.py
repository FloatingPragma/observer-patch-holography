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


def _display_status(status: str) -> str:
    return status.replace("current_corpus", "corpus_limited")


def build_gap_rows() -> list[dict[str, Any]]:
    return [
        {
            "id": "pclosure.compressed-trunk-artifact",
            "lane": "P closure",
            "status": "candidate_artifact",
            "github_issue": 224,
            "title": "Make the five-equation P trunk the canonical audit artifact",
            "current_boundary": (
                "The compressed trunk organizes the code path. It remains outside the certified "
                "particle-root surface until the endpoint and interval-certificate gates close."
            ),
            "next_action": "Keep emitting p_closure_trunk_current.json and use it as audit metadata.",
            "target_surfaces": ["code/P_derivation", "code/particles"],
        },
        {
            "id": "d10.ward-projected-thomson-endpoint",
            "lane": "D10 electromagnetic endpoint",
            "status": "closed_blocker_isolated_endpoint_package",
            "github_issue": 223,
            "successor_github_issue": 235,
            "title": "Close the Ward-projected U(1)_Q endpoint package",
            "current_boundary": (
                "The structured-running law is a continuation. The endpoint package computes the "
                "residual inverse-alpha transport packet and isolates the first non-internalized "
                "object. The source-emitted QCD screening map, scheme remainder, and interval error "
                "control are carried by issue #235."
            ),
            "next_action": (
                "Keep the package as the closed blocker-isolation artifact for issue #223."
            ),
            "target_surfaces": [
                "code/P_derivation/THOMSON_TRANSPORT_THEOREMS.md",
                "code/P_derivation/runtime/thomson_endpoint_package_current.json",
                "code/particles/calibration",
            ],
        },
        {
            "id": "d10.source-residual-map-and-interval-certificate",
            "lane": "D10 electromagnetic endpoint",
            "status": "open_source_residual_map_and_interval_certificate",
            "github_issue": 235,
            "closed_issue_refs": [223],
            "title": "Emit the source-only Thomson residual map and interval certificate",
            "current_boundary": (
                "The endpoint package fixes the target residual, but no OPH source theorem emits "
                "the Ward-projected hadronic spectral measure, electroweak scheme remainder, or "
                "interval certificate for the full map."
            ),
            "next_action": (
                "Emit a source-only residual map that matches the endpoint package, including "
                "rho_had(s;P) or an equivalent Ward-projected spectral primitive, matching remainder, "
                "and certified quadrature bounds."
            ),
            "target_surfaces": [
                "code/P_derivation/THOMSON_TRANSPORT_THEOREMS.md",
                "code/P_derivation/runtime/thomson_endpoint_contract_current.json",
                "code/particles/hadron/ward_projected_spectral_measure.schema.json",
            ],
        },
        {
            "id": "d10.rg-matching-threshold-scheme",
            "lane": "D10 running and matching",
            "status": "open_theorem_gap",
            "github_issue": 32,
            "title": "Internalize RG matching, threshold placement, and scheme conversion",
            "current_boundary": (
                "The D10 branch uses declared running/matching conventions. The compressed P "
                "trunk exposes that boundary without proving those conventions from OPH."
            ),
            "next_action": (
                "Turn the running/matching package into an OPH edge-sector theorem or explicitly "
                "keep it as a declared convention in every prediction surface."
            ),
            "target_surfaces": ["paper compact D10 section", "code/P_derivation", "code/particles/calibration"],
        },
        {
            "id": "pclosure.certified-codepath-adoption",
            "lane": "P closure",
            "status": "blocked_pending_certified_root",
            "github_issue": 224,
            "title": "Replace particle P consumers with the certified P root",
            "current_boundary": (
                "Particle status surfaces use the public display value P=1.630968209403959 "
                "where the value is reader-facing. Certified particle-root adoption remains "
                "blocked while the compressed trunk remains candidate/audit metadata."
            ),
            "next_action": (
                "When issues 235 and 32 close, switch particle builders to the certified trunk "
                "artifact and make compare-only or historical P paths non-default."
            ),
            "target_surfaces": ["code/particles/scripts", "code/particles/runs/status", "WebProjects OPH summaries"],
        },
        {
            "id": "charged.determinant-normalization-transport",
            "lane": "Charged leptons",
            "status": "closed_current_corpus_charged_end_to_end_no_go",
            "github_issue": 201,
            "title": "Keep the P-to-charged affine anchor bridge scoped as a no-go",
            "current_boundary": (
                "No public charged-lepton values are emitted on the theorem lane. The available corpus "
                "does not prove the determinant-normalization / sector-isolated trace-lift identity "
                "beneath A_ch(P), and the impossibility packet rules out end-to-end charged closure "
                "from the present centered-data surface."
            ),
            "next_action": (
                "Keep charged masses suppressed on the public theorem lane. Reopen only for a "
                "theorem-grade uncentered trace lift proving 3 mu(r) = sum_e M_e^ch log q_e(r), "
                "equivalently zero normalization defect N_det(P), on the physical charged branch."
            ),
            "target_surfaces": ["code/particles/leptons", "code/particles/runs/leptons"],
        },
        {
            "id": "quark.selected-class-vs-global-classification",
            "lane": "Quarks",
            "status": "selected_class_closed_global_classification_no_go",
            "github_issue": 199,
            "title": "Keep exact quark rows scoped to the selected public frame class",
            "current_boundary": (
                "The selected-class exact Yukawa theorem is promoted. The stronger class-uniform/global "
                "classification lane is closed as a corpus-limited no-go because no source-emitted "
                "ambient public-frame classifier or quotient-intrinsic sigma law exists."
            ),
            "next_action": (
                "Keep every public exact-quark claim explicitly selected-class. Reopen only for a new "
                "source-emitted global public-frame classifier artifact."
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
            "id": "calibration.direct-top-bridge",
            "lane": "D11/top codomain",
            "status": "closed_current_corpus_codomain_no_go",
            "github_issue": 207,
            "title": "Bridge the exact top coordinate to the auxiliary direct-top PDG row",
            "current_boundary": (
                "The exact top coordinate uses the PDG cross-section codomain Q007TP4. The auxiliary "
                "direct-top entry Q007TP is a separate extraction codomain and remains compare-only; "
                "the available corpus emits no extraction-response map or uncertainty-propagation certificate."
            ),
            "next_action": (
                "Keep Q007TP compare-only while the theorem row remains anchored on Q007TP4. Reopen only "
                "for a concrete source-side extraction-response kernel."
            ),
            "target_surfaces": ["code/particles/calibration", "code/particles/runs/status"],
        },
        {
            "id": "hadron.production-backend-systematics",
            "lane": "Hadrons",
            "status": "closed_out_of_scope_computationally_blocked",
            "github_issue": 153,
            "related_github_issues": [153, 157],
            "title": "Execute the production hadron backend and publish systematics",
            "current_boundary": (
                "The hadron lane is explicitly outside the local pipeline. Issues #153 and "
                "#157 are closed as out-of-scope/computationally blocked because there "
                "is no working production hadron backend here. A credible backend is hardware-gated "
                "on OPH hardware such as GLORB/Echosahedron, outside local Python and Chrome workers."
            ),
            "next_action": (
                "Keep hadron rows suppressed. Reopen only after a working OPH hadron backend emits "
                "production hadron output, Ward-projected spectral data, and systematics."
            ),
            "target_surfaces": ["code/particles/hadron", "code/particles/qcd"],
        },
    ]


def build_bundles() -> list[dict[str, Any]]:
    return [
        {
            "id": "electroweak-root-closure-bundle",
            "status": "endpoint_package_closed_source_residual_open",
            "gap_ids": [
                "pclosure.compressed-trunk-artifact",
                "d10.ward-projected-thomson-endpoint",
                "d10.source-residual-map-and-interval-certificate",
                "d10.rg-matching-threshold-scheme",
                "pclosure.certified-codepath-adoption",
            ],
            "promotion_question": (
                "Can one source-emitted map Delta_Th(P), with declared matching and interval bounds, "
                "certify the compressed P trunk as the particle root without importing alpha(0)?"
            ),
            "result": (
                "Constructive result. The admissible endpoint object is explicit and the endpoint "
                "package computes the residual inverse-alpha packet. Delta_Th(P) must split into "
                "source lepton transport, a Ward-projected hadronic spectral density rho_had(s;P), "
                "a certified electroweak/scheme remainder, RG/matching certificates, quadrature bounds, "
                "and an interval-level fixed-point certificate. The local implementation targets are "
                "P_derivation/runtime/thomson_endpoint_contract_current.json and "
                "P_derivation/runtime/thomson_endpoint_package_current.json."
            ),
        },
        {
            "id": "spectrum-source-bundle",
            "status": "closed_current_corpus_source_boundaries_emitted",
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
                "No promotion. Charged leptons are closed as a corpus-limited no-go by the end-to-end "
                "impossibility theorem: the same-family witness and conditional algebraic readout remain, "
                "but no theorem-grade A_ch(P) is emitted. Quarks remain selected-class on f_P with global "
                "classification closed as a corpus-limited no-go, and neutrino PMNS rows remain visible "
                "comparison-tension rows."
            ),
        },
        {
            "id": "qcd-thomson-backend-bundle",
            "status": "closed_out_of_scope_scope_lock_emitted",
            "gap_ids": [
                "d10.ward-projected-thomson-endpoint",
                "hadron.production-backend-systematics",
            ],
            "promotion_question": (
                "Can the hadron production backend emit the rho_had(s;P) object and uncertainty budget "
                "needed by the Ward-projected Thomson endpoint with hadrons and alpha(0) handled in one "
                "closed packet?"
            ),
            "result": (
                "Constructive result with a scope lock. The stable-channel backend is outside the endpoint "
                "object, and the real hadron backend is outside local execution. The missing primitive remains "
                "production_ward_projected_hadronic_spectral_measure_export. Issues #153/#157 are closed "
                "out-of-scope/computationally blocked pending OPH hardware and outside Chrome workers or "
                "local surrogate code."
            ),
        },
        {
            "id": "top-codomain-bridge-bundle",
            "status": "closed_current_corpus_codomain_no_go",
            "gap_ids": [
                "calibration.direct-top-bridge",
            ],
            "promotion_question": (
                "Can the exact top coordinate be mapped into the auxiliary direct-top extraction codomain "
                "without using Q007TP as a calibration input?"
            ),
            "result": (
                "No-go result. The exact top theorem row remains on Q007TP4. The auxiliary direct-top row "
                "Q007TP is compare-only because the available corpus emits no source-side extraction-response "
                "kernel into that codomain."
            ),
        },
        {
            "id": "particle-root-integration-gate",
            "status": "keep_candidate_with_constructive_next_artifacts",
            "gap_ids": [
                "pclosure.compressed-trunk-artifact",
                "d10.ward-projected-thomson-endpoint",
                "d10.rg-matching-threshold-scheme",
                "pclosure.certified-codepath-adoption",
                "charged.determinant-normalization-transport",
                "quark.selected-class-vs-global-classification",
                "neutrino.pmns-status-and-absolute-rows",
                "calibration.direct-top-bridge",
                "hadron.production-backend-systematics",
            ],
            "promotion_question": (
                "Do the returned packets jointly close the endpoint, matching, interval, and source-object "
                "requirements strongly enough to promote the compressed trunk into particle builders?"
            ),
            "result": (
                "No promotion. The first wave emits constructive next artifacts, so the compressed P trunk "
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
            "compressed_p_trunk_is_certified_prediction_root": False,
            "reason": "The endpoint, RG/matching, and interval-certificate gates remain open.",
            "hadron_backend_in_current_local_scope": False,
            "hadron_backend_scope_reason": (
                "Production hadrons require a real OPH hardware backend. Issues #153/#157 are closed "
                "as out-of-scope/computationally blocked; local surrogate output is non-promoting."
            ),
            "torus_mode_language_allowed_in_pipeline": False,
            "address_remaining_blockers_one_by_one": False,
            "obstruction_only_worker_result_allowed": True,
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
        f"- May feed promoted particle predictions: `{p_trunk['may_feed_live_particle_predictions']}`",
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
            f"| `{bundle['id']}` | `{_display_status(bundle['status'])}` | {gaps} | {bundle['promotion_question']} |"
        )
    lines.extend(
        [
            "",
            "## Bundle Packet Results",
            "",
        ]
    )
    for bundle in ledger["bundles"]:
        lines.append(f"- `{bundle['id']}`: `{_display_status(bundle['status'])}`. {bundle['result']}")
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
            f"| `{row['id']}` | {row['lane']} | `{_display_status(row['status'])}` | {issue} | {row['next_action']} |"
        )
    lines.extend(
        [
            "",
            "## Claim Policy",
            "",
            "- The compressed P trunk is an audit/candidate artifact until the endpoint and certificate gates close.",
            "- The remaining blockers should be worked as coupled bundles, not as isolated one-off fixes.",
            "- The particle pipeline must keep compare-only, continuation, selected-class, and theorem-grade rows mechanically distinct.",
            "- Golden-ratio torus or resonance language is not a derivation input unless a separate representation-to-spectrum theorem is supplied.",
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
