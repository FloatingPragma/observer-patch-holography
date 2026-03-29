#!/usr/bin/env python3
"""Validate the intrinsic generation-bundle branch-generator candidate surface."""

from __future__ import annotations

import argparse
import json
import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "flavor" / "generation_bundle_branch_generator.json"


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate generation-bundle branch-generator candidate fields.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT))
    args = parser.parse_args()

    payload = json.loads(pathlib.Path(args.input).read_text(encoding="utf-8"))
    if payload.get("artifact") != "oph_intrinsic_generation_bundle_branch_generator":
        print("generation-bundle artifact kind missing", file=sys.stderr)
        return 1
    if payload.get("carrier_dimension") != 3 or payload.get("realized_generation_count") != 3:
        print("generation-bundle artifact does not expose the realized three-generation carrier", file=sys.stderr)
        return 1
    if payload.get("remaining_missing_theorem") != "oph_generation_bundle_branch_generator_splitting":
        print("generation-bundle artifact does not reduce the blocker to the branch-generator splitting theorem", file=sys.stderr)
        return 1
    if payload.get("operator_theorem_candidate") != "oph_generation_bundle_branch_generator_splitting":
        print("generation-bundle theorem candidate id missing", file=sys.stderr)
        return 1
    transfer = dict(payload.get("actual_generator_transfer_candidate", {}))
    if transfer.get("smaller_exact_missing_clause") != "compression_descendant_commutator_vanishes_or_is_uniformly_quadratic_small_after_central_split":
        print("generation-bundle artifact is missing the commutator-transfer bridge reduction", file=sys.stderr)
        return 1
    if transfer.get("actual_proxy_centered_residual_kind") != "compression_descendant_commutator":
        print("generation-bundle artifact is missing the descended commutator residual kind", file=sys.stderr)
        return 1
    if transfer.get("first_order_residual_after_central_split") != "vanishes_if_commutator_zero":
        print("generation-bundle artifact is missing the first-order residual vanishing claim", file=sys.stderr)
        return 1
    if transfer.get("descended_commutator_control_mode") != "exact_zero_or_uniform_quadratic":
        print("generation-bundle artifact is missing the zero-or-quadratic commutator control mode", file=sys.stderr)
        return 1
    if transfer.get("quadratic_factorization_claim") != "all surviving centered P->P corrections factor through P->Q->P":
        print("generation-bundle artifact is missing the quadratic factorization claim", file=sys.stderr)
        return 1
    if transfer.get("transfer_if_closed_effect") != "proxy_defect_vs_gap_estimate_lifts_to_actual_generator":
        print("generation-bundle artifact is missing the transfer effect statement", file=sys.stderr)
        return 1
    promotion_gate = dict(payload.get("promotion_gate", {}))
    if promotion_gate.get("exact_missing_ingredient") is None:
        print("generation-bundle promotion gate missing the exact missing ingredient", file=sys.stderr)
        return 1
    if promotion_gate.get("smaller_exact_missing_clause") != "compression_descendant_commutator_vanishes_or_is_uniformly_quadratic_small_after_central_split":
        print("generation-bundle promotion gate is missing the reduced transfer clause", file=sys.stderr)
        return 1
    spectrum = dict(payload.get("simple_spectrum_certificate", {}))
    if not spectrum.get("simple_spectrum"):
        print("generation-bundle artifact does not carry a simple-spectrum certificate on the current family", file=sys.stderr)
        return 1
    persistent_gap = dict(payload.get("persistent_gap_certificate", {}))
    if not persistent_gap.get("riesz_bound_passes"):
        print("generation-bundle artifact does not record the current proxy Riesz safety margin", file=sys.stderr)
        return 1
    print("generation-bundle branch-generator candidate guard passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
