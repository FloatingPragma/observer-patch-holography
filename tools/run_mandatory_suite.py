#!/usr/bin/env python3
"""Run the documented mandatory scientific suite from a clean checkout (#507).

This is the one command REPRODUCE.md documents and CI enforces:

    python tools/run_mandatory_suite.py                # mandatory suite
    python tools/run_mandatory_suite.py --certificates # + exact certificate suites

The mandatory suite validates the claim registry against its live gates,
validates the paper release manifest, proves the full scientific collection
imports cleanly (which is what keeps the optional cloud/hardware lanes
fail-closed), and executes the fast fixture suites. The exact certificate
suites (#566 port-current, #314 matter-lift) run in their own CI workflow on
their own triggers; `--certificates` runs them here with the same commands.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

MANDATORY_STEPS: list[tuple[str, list[str]]] = [
    ("Validate claim registry", [sys.executable, "tools/check_claim_registry.py"]),
    ("Audit issue-518 receipt promotion", [sys.executable, "tools/audit_issue_518_receipts.py"]),
    ("Validate paper release manifest", [sys.executable, "tools/validate_paper_release_manifest.py"]),
    ("Regression-test the manifest validator", [sys.executable, "-m", "pytest", "-q", "tools/test_paper_release_manifest.py"]),
    ("Check the claims scoreboard is regenerated", [sys.executable, "tools/build_scoreboard.py", "--check"]),
    ("Collect the mandatory scientific suite", [sys.executable, "-m", "pytest", "--collect-only", "-q", "code"]),
    ("Execute the audit fixture suite", [sys.executable, "-m", "pytest", "-q", "code/audit"]),
    ("Audit A5 closure ledgers", [sys.executable, "code/a5_closure/test_audit.py"]),
]

CERTIFICATE_STEPS: list[tuple[str, list[str]]] = [
    ("Execute the conditional port-current certificate suite", [sys.executable, "-m", "pytest", "-q", "code/a5_closure/tests/test_port_current_inner_certificate.py"]),
    ("Execute the conditional matter-lift certificate suite", [sys.executable, "-m", "pytest", "-q", "code/a5_closure/tests/test_super_tannakian_matter_lift_certificate.py"]),
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
    args = parser.parse_args()

    steps = [] if args.certificates_only else list(MANDATORY_STEPS)
    if args.certificates or args.certificates_only:
        steps += CERTIFICATE_STEPS
    run_steps(steps)
    print("mandatory suite OK" if not args.certificates_only else "certificate suites OK")


if __name__ == "__main__":
    main()
