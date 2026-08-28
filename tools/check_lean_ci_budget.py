#!/usr/bin/env python3
"""Fail closed if per-change Lean CI can again exceed thirty minutes."""

from __future__ import annotations

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
WORKFLOW = ROOT / ".github" / "workflows" / "lean-ci.yml"
LAKEFILE = ROOT / "Lean" / "lakefile.lean"


def _job_block(workflow: str, name: str) -> str:
    marker = f"  {name}:\n"
    start = workflow.find(marker)
    if start < 0:
        raise ValueError(f"Lean CI is missing the {name!r} job")
    end_match = re.search(r"(?m)^  [A-Za-z0-9_-]+:\s*$", workflow[start + len(marker) :])
    if end_match is None:
        return workflow[start:]
    return workflow[start : start + len(marker) + end_match.start()]


def validate(workflow: str, lakefile: str) -> list[str]:
    errors: list[str] = []

    default_ophgap = re.search(
        r"@\[default_target\]\s*(?:(?:--[^\n]*)?\n\s*)*lean_lib\s+«OphGap»",
        lakefile,
    )
    if default_ophgap:
        errors.append("OphGap must not be a Lake default target")
    if "lean_lib «OphGap»" not in lakefile:
        errors.append("the OphGap Lake library is missing")

    push_header = re.search(r"(?ms)^  push:\n(?P<body>.*?)(?=^  pull_request:)", workflow)
    if push_header is None or "    branches: [main]" not in push_header.group("body"):
        errors.append("Lean CI branch pushes must be limited to main")
    if "cancel-in-progress: true" not in workflow:
        errors.append("Lean CI must cancel superseded runs")

    try:
        build = _job_block(workflow, "build")
        ophgap = _job_block(workflow, "ophgap")
    except ValueError as exc:
        errors.append(str(exc))
        return errors

    for name, block in (("build", build), ("ophgap", ophgap)):
        match = re.search(r"(?m)^    timeout-minutes:\s*(\d+)\s*$", block)
        if match is None:
            errors.append(f"Lean CI {name!r} job has no hard timeout")
        elif int(match.group(1)) > 30:
            errors.append(f"Lean CI {name!r} job exceeds the 30-minute ceiling")

    if "timeout 18m lake build" not in build:
        errors.append("the default Lean build must retain its 18-minute resumable budget")
    if "timeout 23m lake build OphGap" not in ophgap:
        errors.append("the OphGap build must retain its 23-minute resumable budget")
    if "Detect OphGap-relevant changes" not in ophgap:
        errors.append("the OphGap job must remain change-aware")
    for name, block in (("build", build), ("ophgap", ophgap)):
        if "${{ github.run_attempt }}" not in block:
            errors.append(f"Lean CI {name!r} cache keys must support resumable re-runs")
    for path in (
        "Lean/OphGap Lean/OphGap.lean",
        "Lean/lean-toolchain Lean/lakefile.lean Lean/lake-manifest.json",
        ".github/workflows/lean-ci.yml",
    ):
        if path not in ophgap:
            errors.append(f"the OphGap change detector is missing {path!r}")

    full_lake_cache = re.compile(r"(?m)^\s+Lean/(?:ObservableNormalForms/)?\.lake\s*$")
    if full_lake_cache.search(workflow):
        errors.append("Lean CI must cache build artifacts, not the full Lake dependency tree")

    return errors


def main() -> None:
    errors = validate(
        WORKFLOW.read_text(encoding="utf-8"),
        LAKEFILE.read_text(encoding="utf-8"),
    )
    if errors:
        for error in errors:
            print(f"ERROR: {error}")
        raise SystemExit(1)
    print("Lean CI budget wiring: OK")


if __name__ == "__main__":
    main()
