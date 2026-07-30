#!/usr/bin/env python3
"""Independent fail-closed checker for the certified zero-exclusion receipt.

This checker deliberately does not import the interval evaluator or
the certified producer.  It validates the v2 schema, recomputes the
artifact and module digests, pins the rational fixture and the exact
finite dimensional-prefactor corrections, replays the exact rational
interior cut-exclusion facts from the fixture alone, and recomputes
every certification verdict from the recorded row data instead of
trusting any stored flag.

A passing checker result means the committed receipt is internally
bound, accurately scoped, and that its principal-sheet zero-exclusion
verdict follows from its recorded gates.  It is not a pole, Laurent,
current-pole, physical-unit, or OPH-native certificate.
"""

from __future__ import annotations

import hashlib
import json
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any

import sympy as sp
from jsonschema import Draft202012Validator

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT_PATH = ROOT / "outputs" / "certified_wz_contours.json"
VECTOR_PATH = ROOT / "outputs" / "fj_direct_vector_blocks.json"
SCHEMA_PATH = ROOT / "schemas" / "certified_wz_contours_v2.schema.json"

PASS_STATUS = "CERTIFIED_CONTOUR_CHECK_PASS"

EXPECTED_CORRECTIONS = {
    "W": "(312000*p2 + 90163)/1944000",
    "Z": "(292020000*p2 + 79268959)/1944000000",
}

BOXES = {
    "W": ((Fraction(109, 1000), Fraction(114, 1000)),
          (Fraction(1, 5000), Fraction(1, 500))),
    "Z": ((Fraction(172, 1000), Fraction(176, 1000)),
          (Fraction(1, 5000), Fraction(3, 2000))),
}
TREE_MASSES = {"W": Fraction(1, 9), "Z": Fraction(25, 144)}


def fail(message: str) -> None:
    print(f"CERTIFIED_CONTOUR_CHECK_FAIL: {message}", file=sys.stderr)
    raise SystemExit(1)


def sha256_bytes(data: bytes) -> str:
    return "sha256:" + hashlib.sha256(data).hexdigest()


def canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n"


def check_digests(receipt: dict) -> None:
    recomputed = sha256_bytes(
        canonical_json(
            {k: v for k, v in receipt.items() if k != "receipt_sha256"}
        ).encode("utf-8")
    )
    if receipt["receipt_sha256"] != recomputed:
        fail("receipt digest mismatch")
    pins = receipt["pins"]
    for key, path in (
        ("vector_blocks_sha256", VECTOR_PATH),
        ("diagnostic_module_sha256", ROOT / "producers" / "wz_pole_receipts.py"),
        ("interval_module_sha256", ROOT / "producers" / "complex_interval.py"),
        ("verifier_module_sha256", ROOT / "producers" / "certified_wz_contours.py"),
    ):
        if pins[key] != sha256_bytes(path.read_bytes()):
            fail(f"pin {key} does not match the working tree")


def check_corrections(receipt: dict) -> None:
    p2 = sp.Symbol("p2")
    for name, expected in EXPECTED_CORRECTIONS.items():
        recorded = sp.sympify(
            receipt["dimensional_prefactor_finite_correction"][name],
            locals={"p2": p2},
        )
        if sp.simplify(recorded - sp.sympify(expected, locals={"p2": p2})) != 0:
            fail(f"{name} prefactor correction drifted")


def replay_interior(receipt: dict) -> None:
    """Exact rational replay of the cut-exclusion facts per particle."""

    for name, ((re_lo, re_hi), (im_lo, im_hi)) in BOXES.items():
        if not (im_lo > 0 and re_lo > 0):
            fail(f"{name} box does not satisfy the base facts")
        for row in receipt["results"][name].values():
            interior = row["interior_holomorphy"]
            if not interior["base_facts"]["box_in_open_upper_half_plane"]:
                fail(f"{name} recorded base fact contradicts the box")
            for chart in interior["loop_charts"]:
                if chart.get("kind") == "denominator":
                    fail(f"{name} records a denominator failure")
                m1 = Fraction(chart["m1"])
                m2 = Fraction(chart["m2"])
                if m1 > 0 and m2 > 0:
                    outside = (m1 + m2) < re_lo or (m1 + m2) > re_hi
                    if bool(chart["threshold_sum_outside_box_re_range"]) != outside:
                        fail(
                            f"{name} threshold-sum comparison for "
                            f"({m1},{m2}) does not replay"
                        )
                    if bool(chart["holomorphic"]) != outside:
                        fail(f"{name} two-massive chart verdict does not replay")
                elif not chart["holomorphic"]:
                    fail(f"{name} massless chart recorded non-holomorphic")


def recompute_verdicts(receipt: dict) -> None:
    all_certified = True
    for name in ("W", "Z"):
        for row in receipt["results"][name].values():
            winding = row["boundary_winding"]
            expected = bool(
                row["interior_holomorphy"]["holomorphic_on_box"]
                and winding["certified"]
                and winding.get("winding") == 0
            )
            if bool(row["zero_exclusion_certified"]) != expected:
                fail(f"{name} row verdict does not follow from its gates")
            all_certified = all_certified and expected
        nesting = receipt["precision_nesting"][name]
        all_certified = (
            all_certified
            and nesting["enclosures_nested_with_precision"]
            and all(nesting["per_quantity_probe_nesting"].values())
            and nesting["per_segment_enclosure_nesting"]["all_nested"]
        )
        rows = receipt["results"][name]
        partitions = [rows[p]["boundary_winding"]["partition"] for p in ("128", "192", "256")]
        if partitions != ["adaptive", "replayed_base_partition", "replayed_base_partition"]:
            fail(f"{name} partition ladder is not base-then-replay")
        segment_counts = {rows[p]["boundary_winding"]["segments"] for p in rows}
        if len(segment_counts) != 1:
            fail(f"{name} segment counts differ across the fixed partition")
    expected_status = (
        "PRINCIPAL_SHEET_ZERO_EXCLUSION_CERTIFIED"
        if all_certified
        else "CERTIFICATION_INCOMPLETE"
    )
    if receipt["status"] != expected_status:
        fail("status does not follow from the recorded gates")
    promotion = receipt["promotion"]
    for key in (
        "complex_ball_certified",
        "sheet_certified_on_declared_boxes",
        "principal_sheet_zero_exclusion_certified",
        "root_count_certified_on_declared_boxes",
    ):
        if bool(promotion[key]) != all_certified:
            fail(f"promotion flag {key} does not follow from the gates")


def main() -> int:
    receipt = json.loads(ARTIFACT_PATH.read_text(encoding="utf-8"))
    schema = json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))
    errors = sorted(
        Draft202012Validator(schema).iter_errors(receipt),
        key=lambda e: e.json_path,
    )
    if errors:
        fail(f"schema: {errors[0].json_path}: {errors[0].message}")
    check_digests(receipt)
    check_corrections(receipt)
    replay_interior(receipt)
    recompute_verdicts(receipt)
    for name, tree in TREE_MASSES.items():
        for row in receipt["results"][name].values():
            if Fraction(row["tree_mass_sq"]) != tree:
                fail(f"{name} tree mass drifted")
    print(PASS_STATUS)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
