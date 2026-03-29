#!/usr/bin/env python3
"""Fail if the flavor sandbox starts choosing labels by experimental identity."""

from __future__ import annotations

import argparse
import json
import pathlib
import re
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_TARGETS = [
    ROOT / "particles" / "flavor" / "derive_overlap_flavor_observable.py",
    ROOT / "particles" / "flavor" / "derive_family_transport_kernel.py",
    ROOT / "particles" / "flavor" / "derive_suppression_phase_tensors.py",
    ROOT / "particles" / "flavor" / "build_forward_yukawas.py",
    ROOT / "particles" / "flavor" / "export_flavor_dictionary_artifact.py",
]
DEFAULT_ARTIFACT = ROOT / "particles" / "runs" / "flavor" / "flavor_observable_artifact.json"
FORBIDDEN_PATTERNS = [
    re.compile(r"\bclosest\b", re.IGNORECASE),
    re.compile(r"\bbest\s*fit\b", re.IGNORECASE),
    re.compile(r"\bmeasured\b", re.IGNORECASE),
    re.compile(r"\bexperimental\b", re.IGNORECASE),
    re.compile(r"\bPDG\b"),
    re.compile(r"\bKoide\b"),
    re.compile(r"\bPMNS\b"),
    re.compile(r"\bWolfenstein\b"),
    re.compile(r"\belectron\b", re.IGNORECASE),
    re.compile(r"\bmuon\b", re.IGNORECASE),
    re.compile(r"\btau\b", re.IGNORECASE),
    re.compile(r"\bcharm\b", re.IGNORECASE),
    re.compile(r"\bstrange\b", re.IGNORECASE),
    re.compile(r"\bbottom\b", re.IGNORECASE),
]


def _validate_artifact(path: pathlib.Path) -> list[str]:
    if not path.exists():
        return []
    payload = json.loads(path.read_text(encoding="utf-8"))
    failures: list[str] = []
    labels = payload.get("labels", [])
    for label in labels:
        for pattern in FORBIDDEN_PATTERNS:
            if pattern.search(str(label)):
                failures.append(f"{path}: label {label!r} matched forbidden pattern {pattern.pattern!r}")
    spectral_gaps = [float(item) for item in payload.get("spectral_gaps", [])]
    if spectral_gaps and any(gap <= 0.0 for gap in spectral_gaps):
        failures.append(f"{path}: non-positive spectral gap detected")
    if "family_projectors" in payload and len(payload.get("family_projectors", [])) not in (0, 3):
        failures.append(f"{path}: expected either zero or three family projectors")
    return failures


def main() -> int:
    parser = argparse.ArgumentParser(description="Check flavor scripts for label-selection leakage.")
    parser.add_argument("targets", nargs="*", help="Files to scan. Defaults to the /particles flavor lane.")
    parser.add_argument("--artifact", default=str(DEFAULT_ARTIFACT), help="Optional flavor-observable artifact to validate.")
    args = parser.parse_args()

    targets = [pathlib.Path(item) for item in args.targets] if args.targets else DEFAULT_TARGETS
    failures: list[str] = []

    for target in targets:
        text = target.read_text(encoding="utf-8")
        for pattern in FORBIDDEN_PATTERNS:
            if pattern.search(text):
                failures.append(f"{target}: matched forbidden pattern {pattern.pattern!r}")

    artifact_path = pathlib.Path(args.artifact)
    failures.extend(_validate_artifact(artifact_path))

    if failures:
        for failure in failures:
            print(failure, file=sys.stderr)
        return 1

    print("no flavor-dictionary disambiguation leaks found")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
