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
        golden = _job_block(workflow, "golden_sector_psl2f5")
        ophgap = _job_block(workflow, "ophgap")
    except ValueError as exc:
        errors.append(str(exc))
        return errors

    jobs = (("build", build), ("golden_sector_psl2f5", golden), ("ophgap", ophgap))
    for name, block in jobs:
        match = re.search(r"(?m)^    timeout-minutes:\s*(\d+)\s*$", block)
        if match is None:
            errors.append(f"Lean CI {name!r} job has no hard timeout")
        elif int(match.group(1)) > 30:
            errors.append(f"Lean CI {name!r} job exceeds the 30-minute ceiling")

    if "Detect changed Lean modules" not in build:
        errors.append("the core Lean build must remain change-aware")
    changed_step = re.search(
        r"(?ms)^      - name: Build changed Lean modules under time budget\n"
        r"(?P<body>.*?)(?=^      - name:)",
        build,
    )
    changed_body = changed_step.group("body") if changed_step else ""
    if "timeout 18m bash -c" not in changed_body or 'lake build "$target"' not in changed_body:
        errors.append("changed Lean modules must retain their 18-minute build budget")
    if re.search(r"(?m)^\s+lake build\s*$", changed_body):
        errors.append("per-change Lean CI must not run an exhaustive default Lake build")
    for source in ("NativeBinaryThermodynamics", "NativeRecordFourLaws"):
        command = f"lake env lean Thermodynamics/{source}.lean"
        if not any(line.strip() == command for line in changed_body.splitlines()):
            errors.append(f"the native thermodynamics attachment {source} must be checked without suppressing failures")
    if "lake build FourLawAdequacySurface CofinalSpectralTailFamily" not in changed_body:
        errors.append("native thermodynamics dependencies must be built before direct source checks")

    cache_step = re.search(
        r"(?ms)^      - name: Save Lake build cache\n"
        r"(?P<body>.*?)(?=^      - name:)",
        build,
    )
    cache_body = cache_step.group("body") if cache_step else ""
    cache_condition = (
        "        if: steps.budget_build.outputs.status == '0' || "
        "steps.budget_build.outputs.status == '124'"
    )
    if (
        "        id: budget_build" not in changed_body.splitlines()
        or cache_condition not in cache_body.splitlines()
        or "uses: actions/cache/save@" not in cache_body
    ):
        errors.append("the core Lean cache must save successful and timed-out builds using the declared budget_build status")

    for excluded in (
        "Screen/OPHScreen.lean|Screen/GoldenSectorPSL2F5Representations.lean) continue ;;",
    ):
        if excluded not in build:
            errors.append("the core Lean build must leave the import-only OPHScreen umbrella and typed golden bridge to their dedicated/nightly lanes")

    if "Detect Golden-sector PSL2F5-relevant changes" not in golden:
        errors.append("the typed golden-sector job must remain change-aware")
    if "timeout 23m lake build GoldenSectorPSL2F5Representations" not in golden:
        errors.append("the typed golden-sector bridge must retain its 23-minute resumable budget")
    if "Verify Golden-sector umbrella registration" not in golden:
        errors.append("the typed golden-sector lane must verify its OPHScreen/Lake registration")

    if "timeout 23m lake build OphGap" not in ophgap:
        errors.append("the OphGap build must retain its 23-minute resumable budget")
    if "Detect OphGap-relevant changes" not in ophgap:
        errors.append("the OphGap job must remain change-aware")

    for name, block, step in (
        ("golden_sector_psl2f5", golden, "golden_build"),
        ("ophgap", ophgap, "ophgap_build"),
    ):
        condition = (
            f"steps.{step}.outputs.status == '0' || "
            f"steps.{step}.outputs.status == '124'"
        )
        if condition not in block:
            errors.append(
                f"the {name!r} cache must save only completed or timed-out builds "
                f"using the declared {step} status"
            )
        if "timed_out" in block:
            errors.append(f"the {name!r} job reads a timed_out output that no step sets")

    for name, block in jobs:
        if "${{ github.run_attempt }}" not in block:
            errors.append(f"Lean CI {name!r} cache keys must support resumable re-runs")

    # A lakefile edit changes every hashed key; without a hash-free fallback
    # the next run starts cold and exhausts its budget.
    for name, block, prefix in (
        ("build", build, "lean-core-v3"),
        ("golden_sector_psl2f5", golden, "golden-psl-v1"),
        ("ophgap", ophgap, "ophgap-v1"),
    ):
        fallback = re.compile(rf"(?m)^\s+{re.escape(prefix)}-\$\{{\{{ runner\.os \}}\}}-\s*$")
        if not fallback.search(block):
            errors.append(f"Lean CI {name!r} cache must restore across lakefile edits")

    for path in (
        "Lean/Screen/GoldenSectorPSL2F5Representations.lean",
        "Lean/Screen/OPHScreen.lean",
        "Lean/lean-toolchain Lean/lakefile.lean Lean/lake-manifest.json",
        ".github/workflows/lean-ci.yml",
    ):
        if path not in golden:
            errors.append(f"the typed golden-sector change detector is missing {path!r}")

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
