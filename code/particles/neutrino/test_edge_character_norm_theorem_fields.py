#!/usr/bin/env python3
"""Validate the sharpened edge-character / norm theorem fields."""

from __future__ import annotations

import argparse
import json
import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "neutrino" / "majorana_overlap_defect_scalar_evaluator.json"


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate sharpened Majorana norm-theorem fields.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT))
    args = parser.parse_args()

    payload = json.loads(pathlib.Path(args.input).read_text(encoding="utf-8"))
    if payload.get("theorem_candidate_id") != "oph_majorana_scalar_from_centered_edge_norm":
        print("Majorana scalar evaluator is missing the sharpened theorem candidate id", file=sys.stderr)
        return 1
    if payload.get("sublemma_candidate_id") != "selector_centered_quadraticity_polarization_law_on_edge_bundle":
        print("Majorana scalar evaluator is missing the sharpened quadraticity sublemma id", file=sys.stderr)
        return 1
    if payload.get("smallest_exact_missing_clause") != "selector_overlap_phase_coboundary_trivializes_same_label_edge_transport":
        print("Majorana scalar evaluator did not reduce the gap to the phase-cocycle triviality gate beneath the attachment route", file=sys.stderr)
        return 1
    if payload.get("remaining_theorem_object") != "oph_majorana_scalar_from_centered_edge_norm":
        print("Majorana scalar evaluator still points at the broader theorem object", file=sys.stderr)
        return 1
    if payload.get("exact_remaining_ingredient") != "selector-centered phase cocycle triviality for same-label edge transport":
        print("Majorana scalar evaluator does not expose the exact remaining phase-cocycle ingredient", file=sys.stderr)
        return 1
    if payload.get("phase_cocycle_triviality_status") != "candidate_only":
        print("Majorana scalar evaluator does not expose the phase-cocycle candidate status", file=sys.stderr)
        return 1
    print("edge-character norm theorem guard passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
