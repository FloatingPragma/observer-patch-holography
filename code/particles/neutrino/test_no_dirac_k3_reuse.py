#!/usr/bin/env python3
"""Fail if the neutrino lane reuses the Dirac K3 phase lift."""

from __future__ import annotations

import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
TARGETS = [
    ROOT / "particles" / "neutrino" / "derive_majorana_holonomy_lift.py",
    ROOT / "particles" / "neutrino" / "derive_family_response_tensor.py",
    ROOT / "particles" / "neutrino" / "build_forward_majorana_matrix.py",
]
FORBIDDEN = [
    "_dirac_phase_lift",
    "K3 lift",
    "antisymmetric",
]


def main() -> int:
    failures = []
    for path in TARGETS:
        text = path.read_text(encoding="utf-8")
        for pattern in FORBIDDEN:
            if pattern in text:
                failures.append(f"{path}: forbidden pattern {pattern!r}")
    if failures:
        print("\n".join(failures), file=sys.stderr)
        return 1
    print("no Dirac K3 reuse in neutrino lane")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
