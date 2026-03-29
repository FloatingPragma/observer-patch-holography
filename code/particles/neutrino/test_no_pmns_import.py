#!/usr/bin/env python3
"""Fail if the neutrino sandbox hardcodes PMNS-style imports."""

from __future__ import annotations

import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
TARGETS = [
    ROOT / "particles" / "neutrino" / "derive_neutrino_scale_anchor.py",
    ROOT / "particles" / "neutrino" / "derive_family_response_tensor.py",
    ROOT / "particles" / "neutrino" / "derive_majorana_holonomy_lift.py",
    ROOT / "particles" / "neutrino" / "build_forward_majorana_matrix.py",
    ROOT / "particles" / "neutrino" / "build_forward_splittings.py",
    ROOT / "particles" / "neutrino" / "build_pmns_from_shared_flavor_basis.py",
]
FORBIDDEN = [
    "pmns_best_fit",
    "theta12",
    "theta13",
    "theta23",
    "U_PMNS = [[",
]


def main() -> int:
    failures = []
    for path in TARGETS:
        text = path.read_text(encoding="utf-8")
        for pattern in FORBIDDEN:
            if pattern in text:
                failures.append(f"{path}: forbidden pattern {pattern!r}")
    if failures:
        print("\n".join(failures))
        return 1
    print("no PMNS-import patterns found")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
