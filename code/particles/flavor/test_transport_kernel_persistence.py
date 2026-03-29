#!/usr/bin/env python3
"""Check the transport-kernel and observable artifacts for persistence basics."""

from __future__ import annotations

import argparse
import json
import math
import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_KERNEL = ROOT / "particles" / "runs" / "flavor" / "family_transport_kernel.json"
DEFAULT_OBSERVABLE = ROOT / "particles" / "runs" / "flavor" / "flavor_observable_artifact.json"


def _cycle_consistent(cycle_phases: dict[str, float], tol: float = 1e-6) -> bool:
    if "012" not in cycle_phases or "021" not in cycle_phases:
        return True
    lhs = float(cycle_phases["021"]) + float(cycle_phases["012"])
    wrapped = math.atan2(math.sin(lhs), math.cos(lhs))
    return abs(wrapped) <= tol


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate flavor transport persistence basics.")
    parser.add_argument("--kernel", default=str(DEFAULT_KERNEL), help="Transport-kernel artifact path.")
    parser.add_argument("--observable", default=str(DEFAULT_OBSERVABLE), help="Flavor-observable artifact path.")
    args = parser.parse_args()

    kernel = json.loads(pathlib.Path(args.kernel).read_text(encoding="utf-8"))
    observable = json.loads(pathlib.Path(args.observable).read_text(encoding="utf-8"))

    failures: list[str] = []
    stability = kernel.get("refinement_stability", {})
    if kernel.get("transport_kind") != "conjugacy_class_family_kernel":
        failures.append("transport_kind is not conjugacy_class_family_kernel")
    if int(stability.get("stable_projector_count", 0)) != 3:
        failures.append("stable_projector_count != 3")
    if not bool(stability.get("projector_order_preserved", False)):
        failures.append("projector_order_preserved is false")
    if not bool(stability.get("uniform_three_cluster_gap", False)):
        failures.append("uniform_three_cluster_gap is false")

    gauge = kernel.get("gauge_invariance_checks", {})
    if not bool(gauge.get("conjugation_class_stable", False)):
        failures.append("conjugation_class_stable is false")
    projector_certificate = dict(observable.get("persistent_projector_certificate", {}))
    if int(projector_certificate.get("stable_projector_count", 0)) != 3:
        failures.append("observable persistent_projector_certificate lost the three-projector count")
    if not bool(projector_certificate.get("conjugation_class_stable", False)):
        failures.append("observable persistent_projector_certificate lost conjugation stability")

    spectral_gaps = [float(item) for item in observable.get("spectral_gaps", [])]
    if spectral_gaps and any(gap <= 0.0 for gap in spectral_gaps):
        failures.append("non-positive spectral gap detected")

    cycle_phases = {str(k): float(v) for k, v in dict(observable.get("cycle_phases", {})).items()}
    if not _cycle_consistent(cycle_phases):
        failures.append("cycle orientation inconsistency detected")

    if failures:
        for failure in failures:
            print(failure, file=sys.stderr)
        return 1

    print("transport-kernel persistence checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
