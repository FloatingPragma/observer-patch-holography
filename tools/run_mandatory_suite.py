#!/usr/bin/env python3
"""Run the documented mandatory scientific suite from a clean checkout (#507).

This is the one command REPRODUCE.md documents and CI enforces:

    python tools/run_mandatory_suite.py                # mandatory suite
    python tools/run_mandatory_suite.py --certificates # + exact certificate suites
    python tools/run_mandatory_suite.py --certificate-smoke-only

The mandatory suite validates the claim registry against its live gates,
validates external inputs, public quantitative surfaces, null-model controls,
and the paper release manifest, then proves those gates reject isolated
false-green mutations. It also proves the full scientific collection imports
cleanly (which is what keeps the optional cloud/hardware lanes fail-closed)
and executes the fast fixture suites. The exact certificate suites (#566
port-current, #314 matter-lift) run in their own CI workflow on their own
triggers; `--certificates` runs them here with the same commands.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

MANDATORY_STEPS: list[tuple[str, list[str]]] = [
    (
        "Validate the committed open-problem ledger offline",
        [sys.executable, "tools/build_open_problem_ledger.py", "--check"],
    ),
    ("Validate claim registry", [sys.executable, "tools/check_claim_registry.py"]),
    (
        "Validate the selection ledger and its generated surface",
        [sys.executable, "tools/build_selection_ledger.py", "--check"],
    ),
    (
        "Validate the frozen-prediction ladder and its generated surface",
        [sys.executable, "tools/build_fz_registry.py", "--check"],
    ),
    (
        "Validate external-data provenance pins and declared source gaps",
        [sys.executable, "tools/check_external_data_provenance.py"],
    ),
    (
        "Validate clean-clone receipt path portability",
        [sys.executable, "tools/check_receipt_portability.py"],
    ),
    (
        "Validate public quantitative claim surfaces",
        [sys.executable, "tools/check_public_surface_claims.py"],
    ),
    (
        "Validate the null-model scorecard",
        [sys.executable, "tools/check_null_models.py", "--check"],
    ),
    ("Audit issue-518 receipt promotion", [sys.executable, "tools/audit_issue_518_receipts.py"]),
    (
        "Validate source-bound canonical book PDF assets",
        [sys.executable, "tools/book_pdf_assets.py"],
    ),
    (
        "Regression-test canonical book PDF assets",
        [sys.executable, "-m", "pytest", "-q", "tools/test_book_pdf_assets.py"],
    ),
    ("Validate paper release manifest", [sys.executable, "tools/validate_paper_release_manifest.py"]),
    ("Regression-test the manifest validator", [sys.executable, "-m", "pytest", "-q", "tools/test_paper_release_manifest.py"]),
    (
        "Regression-test offline Phase-0 scientific gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "tools/test_public_surface_claims.py",
            "tools/test_null_models.py",
            "tools/test_receipt_portability.py",
            "tools/test_github_release_channel.py",
            "tools/test_gates_actually_fail.py",
        ],
    ),
    (
        "Regression-test the deterministic publication gates",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "paper/tools/test_check_build_warnings.py",
            "tools/test_reproducible_build_env.py",
            "tools/test_open_problem_ledger.py",
        ],
    ),
    (
        "Execute the Phase-0 proof and non-identifiability receipts",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/consensus/test_issue_517_proof_obligations.py",
            "code/particles/hierarchy/test_antecedent_only_nonidentifiability.py",
            "code/particles/hierarchy/test_hierarchy_bundle.py",
        ],
    ),
    ("Check the claims scoreboard is regenerated", [sys.executable, "tools/build_scoreboard.py", "--check"]),
    ("Collect the mandatory scientific suite", [sys.executable, "-m", "pytest", "--collect-only", "-q", "code"]),
    ("Execute the audit fixture suite", [sys.executable, "-m", "pytest", "-q", "code/audit"]),
    ("Audit A5 closure ledgers", [sys.executable, "code/a5_closure/test_audit.py"]),
]

CERTIFICATE_STEPS: list[tuple[str, list[str]]] = [
    ("Execute the conditional port-current certificate suite", [sys.executable, "-m", "pytest", "-q", "code/a5_closure/tests/test_port_current_inner_certificate.py"]),
    ("Execute the source-bound matter-lift certificate suite", [sys.executable, "-m", "pytest", "-q", "code/a5_closure/tests/test_super_tannakian_matter_lift_certificate.py"]),
    ("Execute the axis-center-descent certificate suite", [sys.executable, "-m", "pytest", "-q", "code/a5_closure/tests/test_axis_center_descent_certificate.py"]),
]

CERTIFICATE_SMOKE_STEPS: list[tuple[str, list[str]]] = [
    (
        "Recompute and verify the canonical port-current certificate",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/a5_closure/tests/test_port_current_inner_certificate.py"
            "::PortCurrentInnerCertificateTests::test_reference_receipt_is_exactly_recomputable",
        ],
    ),
    (
        "Recompute and verify the canonical matter-lift certificate",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/a5_closure/tests/test_super_tannakian_matter_lift_certificate.py"
            "::SuperTannakianMatterLiftTests::test_reference_receipt_is_exactly_recomputable",
        ],
    ),
    (
        "Recompute and verify the canonical axis-center-descent certificate",
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "code/a5_closure/tests/test_axis_center_descent_certificate.py"
            "::AxisCenterDescentTests::test_reference_receipt_is_exactly_recomputable",
        ],
    ),
]


def run_steps(steps: list[tuple[str, list[str]]]) -> None:
    for title, command in steps:
        print(f"==> {title}", flush=True)
        result = subprocess.run(command, cwd=ROOT)
        if result.returncode != 0:
            raise SystemExit(f"FAILED: {title} (exit {result.returncode})")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--certificates",
        action="store_true",
        help="also execute the exact certificate suites (adds ~26 minutes)",
    )
    parser.add_argument(
        "--certificates-only",
        action="store_true",
        help="execute only the exact certificate suites",
    )
    parser.add_argument(
        "--certificate-smoke-only",
        action="store_true",
        help="recompute and verify only the two canonical exact certificate receipts",
    )
    args = parser.parse_args()

    if args.certificates_only and args.certificate_smoke_only:
        parser.error("--certificates-only and --certificate-smoke-only are mutually exclusive")

    steps = [] if args.certificates_only or args.certificate_smoke_only else list(MANDATORY_STEPS)
    if args.certificates or args.certificates_only:
        steps += CERTIFICATE_STEPS
    if args.certificate_smoke_only:
        steps += CERTIFICATE_SMOKE_STEPS
    run_steps(steps)
    if args.certificates_only:
        print("certificate suites OK")
    elif args.certificate_smoke_only:
        print("certificate smoke suite OK")
    else:
        print("mandatory suite OK")


if __name__ == "__main__":
    main()
