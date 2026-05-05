#!/usr/bin/env python3
"""Build the public OPH open-problem ledger from live GitHub issues."""

from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
REPO = "FloatingPragma/observer-patch-holography"
DEFAULT_JSON_OUT = ROOT / "tracking" / "open_issues" / "open_problem_ledger.json"
DEFAULT_MD_OUT = ROOT / "OPEN_PROBLEMS.md"


ISSUE_POLICY: dict[int, dict[str, str]] = {
    28: {
        "phase": "compact-structure",
        "claim_level": "theorem gap",
        "blocker": "Hypercharge lattice and exact Z6 quotient proof packet.",
        "closure": "Compact paper and DAG surfaces cite the same theorem-grade lattice/quotient derivation.",
        "falsification": "A required SM charge assignment cannot be represented on the proposed OPH lattice.",
        "chrome_policy": "Use only for proof audit after a local lattice packet exists.",
    },
    32: {
        "phase": "particle-rg",
        "claim_level": "constructive contract",
        "blocker": "Populate scheme lock, threshold map, beta provenance table, and matching interval certificate.",
        "closure": "Every running/matching ingredient is OPH-derived, explicitly borrowed, or declared with interval impact.",
        "falsification": "Hidden threshold or scheme freedom is needed to recover the quantitative rows.",
        "chrome_policy": "Use only after the local RG packet exists.",
    },
    49: {
        "phase": "observer-monograph",
        "claim_level": "claim-hygiene gap",
        "blocker": "Proton stability and proton-spin claims need theorem-grade status or downgrade wording.",
        "closure": "Observer paper and summaries agree on exact theorem, conditional result, or open-continuation status.",
        "falsification": "A retained proton claim depends on unsupported QCD or baryon-structure assumptions.",
        "chrome_policy": "Use for independent claim audit if local wording remains ambiguous.",
    },
    55: {
        "phase": "observer-monograph",
        "claim_level": "decision gate",
        "blocker": "Critical-superstring lift needs either a theorem packet or an explicit non-claim decision.",
        "closure": "Paper and book consistently mark the lift as proved, conditional, or not part of the OPH core.",
        "falsification": "The lift is required for a core OPH claim but cannot be derived or cleanly removed.",
        "chrome_policy": "Use for proof/claim audit after local classification.",
    },
    58: {
        "phase": "observer-monograph",
        "claim_level": "theorem package gap",
        "blocker": "Observer continuation/backup theorem package.",
        "closure": "Continuation, recovery, and backup statements are packaged with hypotheses and proof boundary.",
        "falsification": "Observer continuation requires nonlocal state access outside OPH overlap rules.",
        "chrome_policy": "Use for proof audit after local theorem packet.",
    },
    59: {
        "phase": "observer-monograph",
        "claim_level": "audit gap",
        "blocker": "Line-by-line audit of the Additional Problem Closures table.",
        "closure": "Every table row links to a theorem, conditional claim, or open issue.",
        "falsification": "A table row asserts closure without a traceable theorem or artifact.",
        "chrome_policy": "Use for independent table audit after local row mapping.",
    },
    60: {
        "phase": "observer-monograph",
        "claim_level": "final audit gate",
        "blocker": "Final proof, citation, and reproducibility audit.",
        "closure": "The monograph passes local build, citation, theorem-status, and reproducibility checks.",
        "falsification": "A core claim cannot be linked to a proof or reproducible artifact.",
        "chrome_policy": "Use for final independent audit after local checks are green.",
    },
    66: {
        "phase": "reality-paper",
        "claim_level": "theorem gap",
        "blocker": "Extension from abelian cycle holonomy to the full OPH defect hierarchy.",
        "closure": "Reality paper states the full hierarchy construction with hypotheses and proof boundary.",
        "falsification": "Nonabelian or higher-defect holonomy is incompatible with the proposed reconciliation law.",
        "chrome_policy": "Use for proof search only after local hierarchy packet.",
    },
    70: {
        "phase": "reality-paper",
        "claim_level": "theorem gap",
        "blocker": "Proof that coarse-graining commutes with reconciliation.",
        "closure": "Coarse-graining/reconciliation commutation theorem is stated and synchronized across reality surfaces.",
        "falsification": "A counterexample shows reconciliation depends on refinement order.",
        "chrome_policy": "Use for proof audit after local counterexample search.",
    },
    112: {
        "phase": "archive-speculative",
        "claim_level": "personal/open hypothesis",
        "blocker": "No public theorem package or falsifiable closure artifact is defined.",
        "closure": "Either define a falsifiable OPH theorem target or archive as non-core speculation.",
        "falsification": "The hypothesis cannot be connected to OPH overlap consistency without extra metaphysics.",
        "chrome_policy": "Do not spend Chrome workers until a precise theorem target exists.",
    },
    113: {
        "phase": "observer-monograph",
        "claim_level": "theorem construction gap",
        "blocker": "Closure map and invariant sector construction.",
        "closure": "Observer paper gives a reproducible closure-map/invariant-sector theorem package.",
        "falsification": "The closure map is not well-defined or is not invariant under allowed patch refinements.",
        "chrome_policy": "Use for proof audit after local construction.",
    },
    153: {
        "phase": "hardware-gated-hadrons",
        "claim_level": "out of local scope",
        "blocker": "Working OPH hadron backend on suitable hardware such as GLORB/Echosahedron.",
        "closure": "Production backend output and continuum/volume/chiral/statistical systematics are published.",
        "falsification": "Surrogate local hadron artifacts are required as if they were production QCD outputs.",
        "chrome_policy": "Do not use Chrome workers for backend execution.",
    },
    155: {
        "phase": "continuation",
        "claim_level": "open branch",
        "blocker": "Strong-CP branch theorem or explicit continuation boundary.",
        "closure": "Strong-CP mechanism is derived, falsified, or downgraded consistently across public surfaces.",
        "falsification": "The proposed OPH branch leaves theta_QCD unconstrained while claiming closure.",
        "chrome_policy": "Use only after a concrete branch packet exists.",
    },
    157: {
        "phase": "hardware-gated-hadrons",
        "claim_level": "out of local scope",
        "blocker": "Nonperturbative QCD/hadron backend and systematics.",
        "closure": "Compact paper states the backend output and budgets or marks hadrons out-of-scope.",
        "falsification": "Compact paper promotes hadron masses without production backend evidence.",
        "chrome_policy": "Do not use Chrome workers for backend execution.",
    },
    199: {
        "phase": "quark-global-classification",
        "claim_level": "selected-class theorem closed, global classification open",
        "blocker": "Class-uniform public quark-frame descent and sigma classification.",
        "closure": "Global quark frame classes are classified or every claim remains selected-class.",
        "falsification": "Another admissible public quark class breaks the selected-class exact theorem language.",
        "chrome_policy": "Use for proof audit after local classification packet.",
    },
    201: {
        "phase": "charged-lepton-source",
        "claim_level": "source theorem gap",
        "blocker": "Sector-isolated charged determinant trace-lift attachment and normalization.",
        "closure": "The determinant character lands on the physical charged determinant line without target readback.",
        "falsification": "The additive determinant normalization remains underdetermined on the OPH source data.",
        "chrome_policy": "Use for proof search after local trace-lift packet.",
    },
    207: {
        "phase": "top-codomain-bridge",
        "claim_level": "constructive conversion contract",
        "blocker": "Source-side extraction response kernel from Q007TP4 coordinate to Q007TP codomain.",
        "closure": "Converted top coordinate and uncertainty budget close without using Q007TP as calibration input.",
        "falsification": "The direct-top codomain cannot be mapped from the current theorem coordinate without a free shift.",
        "chrome_policy": "Use only to audit a proposed response kernel.",
    },
    212: {
        "phase": "quark-off-canonical",
        "claim_level": "candidate-only",
        "blocker": "Target-free off-canonical flavor transport law and P-to-sigma evaluator.",
        "closure": "Off-canonical quark masses move by theorem-grade P-driven transport, not default-universe anchoring.",
        "falsification": "Off-canonical motion requires anchoring to the canonical PDG target surface.",
        "chrome_policy": "Use for proof audit after local off-canonical transport packet.",
    },
    223: {
        "phase": "thomson-endpoint",
        "claim_level": "constructive contract",
        "blocker": "Ward-projected source endpoint including hadronic spectral measure and interval certificate.",
        "closure": "alpha_Th(P) is emitted from source objects with certified transport/error bounds.",
        "falsification": "Measured alpha(0) or a free screened ansatz is required to close the endpoint.",
        "chrome_policy": "Use only after the source endpoint packet exists.",
    },
    224: {
        "phase": "p-root-adoption",
        "claim_level": "blocked on certified root",
        "blocker": "Certified P root after #223 and #32 close.",
        "closure": "Live particle consumers read the certified trunk root and legacy P paths are non-default.",
        "falsification": "Derived P cannot replace legacy/candidate paths without weakening prediction rows.",
        "chrome_policy": "Do not use until endpoint and RG interval gates close.",
    },
    225: {
        "phase": "publication-sync",
        "claim_level": "blocked on #223/#224",
        "blocker": "Certified derived P closure values and live consumer adoption.",
        "closure": "Observers and particle papers quote the same certified P, alpha, and dependent quantities.",
        "falsification": "Paper values must be published before the code root is certified.",
        "chrome_policy": "Use for claim-hygiene audit after code root certification.",
    },
    231: {
        "phase": "axiom-status",
        "claim_level": "independence decision gate",
        "blocker": "Axiom 3/Axiom 4 independence proof or downgrade.",
        "closure": "Five-axiom basis language matches the proved independence/dependence status.",
        "falsification": "Axiom independence is asserted while a dependency proof or countermodel is missing.",
        "chrome_policy": "Use for proof/countermodel audit after local packet.",
    },
    232: {
        "phase": "bw-cap-pair",
        "claim_level": "theorem gap",
        "blocker": "Transported BW/geometric cap-pair extraction and ordered cut-pair rigidity.",
        "closure": "BW cap-pair extraction and ordered cut-pair rigidity are emitted as theorem-grade artifacts.",
        "falsification": "The extracted cap/cut structure depends on arbitrary collar choices.",
        "chrome_policy": "Use for proof audit after local BW packet.",
    },
    233: {
        "phase": "mar-formalization",
        "claim_level": "theorem gap",
        "blocker": "MAR realization space, well-founded order, minimality, and uniqueness theorems.",
        "closure": "MAR claims have formal definitions, order, minimality theorem, and uniqueness boundary.",
        "falsification": "MAR has no well-founded order compatible with the stated minimality claims.",
        "chrome_policy": "Use for proof audit after local formalization.",
    },
    234: {
        "phase": "particle-provenance",
        "claim_level": "audit gap",
        "blocker": "Blind-prediction provenance and convention-sensitivity ledger.",
        "closure": "Every quantitative row records input use, blind/compare status, provenance, and convention sensitivity.",
        "falsification": "A public quantitative row hides target leakage or unquantified convention dependence.",
        "chrome_policy": "Use for audit only after local provenance artifact exists.",
    },
    235: {
        "phase": "p-closure-theorem",
        "claim_level": "residual theorem gap",
        "blocker": "Root monotonicity/uniqueness and exact endpoint boundary.",
        "closure": "P-closure has a scalar root theorem and endpoint status boundary synchronized across surfaces.",
        "falsification": "The claimed root is nonunique or endpoint-dependent in a way the papers do not state.",
        "chrome_policy": "Use for proof audit after local monotonicity packet.",
    },
    236: {
        "phase": "publication-ledger",
        "claim_level": "ledger artifact",
        "blocker": "Keep this ledger synchronized across README, papers, book, and public summaries.",
        "closure": "Public OPH surfaces point to the same open-problem ledger and no longer imply hidden closure.",
        "falsification": "A public surface claims a branch is closed while this ledger or its issue says open.",
        "chrome_policy": "Not needed unless a downstream claim-hygiene audit is ambiguous.",
    },
    237: {
        "phase": "screen-microphysics",
        "claim_level": "benchmark suite gap",
        "blocker": "Explicit existence-program benchmark suite for screen microphysics reference architecture.",
        "closure": "Benchmarks specify model class, expected artifacts, pass/fail criteria, and reproducibility path.",
        "falsification": "The reference architecture cannot be exercised by any finite benchmark suite.",
        "chrome_policy": "Use for architecture audit after local benchmark suite skeleton.",
    },
}


CLOSED_OUT_OF_SCOPE_ISSUES: dict[int, dict[str, str]] = {
    153: {
        "title": "hadron backend/systematics",
        "url": "https://github.com/FloatingPragma/observer-patch-holography/issues/153",
        "closed_status": "closed_not_planned_out_of_scope_computationally_blocked",
        "closure_note": (
            "Closed as out-of-scope, not solved. The current pipeline emits no hadron predictions "
            "until a working OPH hadron backend on suitable hardware publishes production output and systematics."
        ),
    },
    157: {
        "title": "nonperturbative QCD hadron branch",
        "url": "https://github.com/FloatingPragma/observer-patch-holography/issues/157",
        "closed_status": "closed_not_planned_out_of_scope_computationally_blocked",
        "closure_note": (
            "Closed as out-of-scope, not solved. Compact-paper hadron claims remain suppressed until "
            "a working OPH backend emits the required nonperturbative QCD/hadron data."
        ),
    },
}


def _now_utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _ascii(text: str) -> str:
    return text.encode("ascii", "ignore").decode("ascii")


def _run_gh() -> list[dict[str, Any]]:
    completed = subprocess.run(
        [
            "gh",
            "issue",
            "list",
            "-R",
            REPO,
            "--state",
            "open",
            "--limit",
            "500",
            "--json",
            "number,title,labels,updatedAt,url",
        ],
        check=True,
        text=True,
        capture_output=True,
    )
    return json.loads(completed.stdout)


def _fallback_policy(issue: dict[str, Any]) -> dict[str, str]:
    return {
        "phase": "unclassified",
        "claim_level": "open",
        "blocker": "Classify blocker from the live issue body.",
        "closure": "Add exact closure criterion to this ledger.",
        "falsification": "Add exact falsification criterion to this ledger.",
        "chrome_policy": "Do not launch workers until the issue has a concrete local packet.",
    }


def build_ledger(issues: list[dict[str, Any]]) -> dict[str, Any]:
    rows = []
    open_numbers = {int(issue["number"]) for issue in issues}
    for issue in sorted(issues, key=lambda item: item["number"]):
        policy = ISSUE_POLICY.get(int(issue["number"]), _fallback_policy(issue))
        labels = [label["name"] for label in issue.get("labels", [])]
        rows.append(
            {
                "number": int(issue["number"]),
                "title": _ascii(issue["title"]),
                "url": issue["url"],
                "labels": labels,
                "updated_at": issue["updatedAt"],
                **policy,
            }
        )
    closed_out_of_scope_records = []
    for number, record in sorted(CLOSED_OUT_OF_SCOPE_ISSUES.items()):
        if number in open_numbers:
            continue
        policy = ISSUE_POLICY[number]
        closed_out_of_scope_records.append(
            {
                "number": number,
                "title": record["title"],
                "url": record["url"],
                "closed_status": record["closed_status"],
                "closure_note": record["closure_note"],
                **policy,
            }
        )
    return {
        "artifact": "oph_open_problem_ledger",
        "generated_utc": _now_utc(),
        "repo": REPO,
        "open_issue_count": len(rows),
        "closed_out_of_scope_count": len(closed_out_of_scope_records),
        "worker_policy": {
            "chrome_pro_workers_default": "local_first",
            "max_parallel_workers": 6,
            "launch_condition": "only after a concrete theorem, audit, or implementation packet exists",
            "obstruction_only_result_allowed": False,
        },
        "closed_out_of_scope_records": closed_out_of_scope_records,
        "rows": rows,
    }


def render_markdown(ledger: dict[str, Any]) -> str:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in ledger["rows"]:
        grouped[row["phase"]].append(row)

    lines = [
        "# OPH Open Problem Ledger",
        "",
        f"Generated: `{ledger['generated_utc']}` from live GitHub issues in `{ledger['repo']}`.",
        "",
        "This is the public boundary between closed OPH claims, conditional claims, and open work. "
        "Dedicated GitHub issues remain canonical for task state; this ledger records the current "
        "claim level, missing artifact, closure criterion, falsification route, and Chrome Pro worker policy.",
        "",
        "Worker policy: local artifacts first; up to six Chrome Pro workers may be used only after a "
        "concrete theorem, audit, or implementation packet exists. Obstruction-only worker output is not accepted.",
        "",
        f"Open issue count: `{ledger['open_issue_count']}`",
        "",
    ]
    closed_records = ledger.get("closed_out_of_scope_records", [])
    if closed_records:
        lines.extend(
            [
                "## Closed Out-Of-Scope Records",
                "",
                "These are not solved derivations. They are closed as non-current-scope tasks because "
                "the required computation depends on a working OPH hadron backend and suitable hardware.",
                "",
                "| Issue | Status | Claim level | Closure note | Reopen criterion | Chrome policy |",
                "| --- | --- | --- | --- | --- | --- |",
            ]
        )
        for row in closed_records:
            lines.append(
                f"| [#{row['number']}]({row['url']}) {row['title']} | `{row['closed_status']}` | "
                f"`{row['claim_level']}` | {row['closure_note']} | {row['closure']} | {row['chrome_policy']} |"
            )
        lines.append("")
    for phase in sorted(grouped):
        lines.extend([f"## {phase}", ""])
        lines.extend(
            [
                "| Issue | Claim level | Missing artifact / blocker | Closure criterion | Falsification route | Chrome policy |",
                "| --- | --- | --- | --- | --- | --- |",
            ]
        )
        for row in grouped[phase]:
            lines.append(
                f"| [#{row['number']}]({row['url']}) {row['title']} | `{row['claim_level']}` | "
                f"{row['blocker']} | {row['closure']} | {row['falsification']} | {row['chrome_policy']} |"
            )
        lines.append("")
    return "\n".join(lines).rstrip()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build the OPH open-problem ledger from GitHub issues.")
    parser.add_argument("--json-out", default=str(DEFAULT_JSON_OUT))
    parser.add_argument("--markdown-out", default=str(DEFAULT_MD_OUT))
    parser.add_argument("--print-json", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    ledger = build_ledger(_run_gh())
    json_text = json.dumps(ledger, indent=2, sort_keys=True, ensure_ascii=True) + "\n"

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
