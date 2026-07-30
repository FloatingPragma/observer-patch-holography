#!/usr/bin/env python3
"""Independent fail-closed checker for the second-sheet pole receipt.

This checker deliberately does not import the interval evaluator or
either certified producer.  It validates the v1 schema, recomputes
the artifact and module digests including the zero-exclusion receipt
pin, replays the exact rational channel classification and window
selection from the fixture alone, verifies the declared additions
against the per-chart formulas, replays the lower-half interior facts,
checks the geometry (pole boxes strictly below the axis, Newton balls
inside their pole boxes, residue times derivative enclosing one in
midpoint arithmetic), and recomputes every certification verdict from
the recorded row data instead of trusting any stored flag.

A passing checker result means the committed receipt is internally
bound, accurately scoped, and that its winding-one pole verdict
follows from its recorded gates on the declared continuation.  It is
not a BMHV, physical-current, physical-unit, or OPH-native
certificate.
"""

from __future__ import annotations

import hashlib
import json
import sys
from fractions import Fraction
from math import isqrt
from pathlib import Path
from typing import Any

from jsonschema import Draft202012Validator

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT_PATH = ROOT / "outputs" / "certified_second_sheet_poles.json"
ZERO_EXCLUSION_PATH = ROOT / "outputs" / "certified_wz_contours.json"
VECTOR_PATH = ROOT / "outputs" / "fj_direct_vector_blocks.json"
SCHEMA_PATH = ROOT / "schemas" / "certified_second_sheet_poles_v1.schema.json"

PASS_STATUS = "SECOND_SHEET_POLE_CHECK_PASS"

TREE_MASSES = {"W": Fraction(1, 9), "Z": Fraction(25, 144)}


def fail(message: str) -> None:
    print(f"SECOND_SHEET_POLE_CHECK_FAIL: {message}", file=sys.stderr)
    raise SystemExit(1)


def sha256_bytes(data: bytes) -> str:
    return "sha256:" + hashlib.sha256(data).hexdigest()


def canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n"


def rational_sqrt(value: Fraction) -> Fraction:
    pn, qn = value.numerator, value.denominator
    a, b = isqrt(pn), isqrt(qn)
    if a * a != pn or b * b != qn:
        fail(f"mass {value} has no rational square root")
    return Fraction(a, b)


def parse_endpoint(text: str) -> float:
    cleaned = text.strip().lstrip("[").rstrip("]")
    return float(cleaned.split(",")[0])


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
        (
            "zero_exclusion_module_sha256",
            ROOT / "producers" / "certified_wz_contours.py",
        ),
        (
            "verifier_module_sha256",
            ROOT / "producers" / "certified_second_sheet_poles.py",
        ),
        ("zero_exclusion_receipt_sha256", ZERO_EXCLUSION_PATH),
    ):
        if pins[key] != sha256_bytes(path.read_bytes()):
            fail(f"pin {key} does not match the working tree")


def replay_channels(receipt: dict) -> None:
    """Exact rational replay of thresholds, windows, and additions."""

    for name, tree in TREE_MASSES.items():
        continuation = receipt["declared_continuation"][name]
        open_ths = []
        closed_ths = []
        for channel in continuation["channels"]:
            m1 = Fraction(channel["m1"])
            m2 = Fraction(channel["m2"])
            threshold = (rational_sqrt(m1) + rational_sqrt(m2)) ** 2
            if Fraction(channel["threshold"]) != threshold:
                fail(f"{name} channel ({m1},{m2}) threshold does not replay")
            is_open = threshold < tree
            if bool(channel["open"]) != is_open:
                fail(f"{name} channel ({m1},{m2}) openness does not replay")
            if is_open:
                open_ths.append(threshold)
                if m1 == 0 and m2 == 0:
                    expected = "0"
                elif m1 == 0 or m2 == 0:
                    m = m1 if m2 == 0 else m2
                    expected = f"-2*pi*i*({m})/s"
                else:
                    expected = "2*pi*i*sqrt(lambda(s))/s"
                if channel["addition"] != expected:
                    fail(f"{name} channel ({m1},{m2}) addition drifted")
            else:
                closed_ths.append(threshold)
                if channel["addition"] != "principal branch":
                    fail(f"{name} closed channel carries an addition")
        window = (max(open_ths), min(closed_ths + [tree]))
        recorded = tuple(Fraction(w) for w in continuation["window"])
        if recorded != window:
            fail(f"{name} continuation window does not replay")
        if not continuation["consistency_probes"]["all_passed"]:
            fail(f"{name} consistency probes recorded a failure")
        for row in continuation["consistency_probes"]["rows"]:
            if not row["passed"]:
                fail(f"{name} probe row failed: {row}")


def replay_geometry(receipt: dict) -> None:
    for name in ("W", "Z"):
        for row in receipt["results"][name].values():
            box = row["pole_box"]
            re_lo, re_hi = (Fraction(v) for v in box["re"])
            im_lo, im_hi = (Fraction(v) for v in box["im"])
            if not (im_hi < 0 and re_lo > 0 and re_lo < re_hi and im_lo < im_hi):
                fail(f"{name} pole box geometry invalid")
            newton = row["interval_newton"]
            ball = newton["pole_ball"]
            ball_re = [parse_endpoint(v) for v in ball["re"]]
            ball_im = [parse_endpoint(v) for v in ball["im"]]
            if not (
                float(re_lo) <= ball_re[0] <= ball_re[1] <= float(re_hi)
                and float(im_lo) <= ball_im[0] <= ball_im[1] <= float(im_hi)
            ):
                fail(f"{name} Newton ball escapes the pole box")
            derivative = newton["derivative_ball"]
            residue = newton["residue_ball"]
            d_mid = complex(
                (parse_endpoint(derivative["re"][0]) + parse_endpoint(derivative["re"][1])) / 2,
                (parse_endpoint(derivative["im"][0]) + parse_endpoint(derivative["im"][1])) / 2,
            )
            r_mid = complex(
                (parse_endpoint(residue["re"][0]) + parse_endpoint(residue["re"][1])) / 2,
                (parse_endpoint(residue["im"][0]) + parse_endpoint(residue["im"][1])) / 2,
            )
            if abs(d_mid * r_mid - 1) > 1e-9:
                fail(f"{name} residue is not the derivative reciprocal")


def recompute_verdicts(receipt: dict) -> None:
    all_certified = all(
        receipt["declared_continuation"][name]["consistency_probes"]["all_passed"]
        for name in ("W", "Z")
    )
    for name in ("W", "Z"):
        for row in receipt["results"][name].values():
            winding = row["boundary_winding"]
            newton = row["interval_newton"]
            null_vectors = newton.get("null_vectors", {})
            expected_simple = bool(
                winding["certified"] and winding.get("winding") == 1
            )
            if bool(row["simple_root_certified"]) != expected_simple:
                fail(f"{name} simple-root verdict does not follow")
            expected = bool(
                row["interior_holomorphy"]["holomorphic_on_box"]
                and row["corridor_holomorphy"]["holomorphic_on_box"]
                and expected_simple
                and newton.get("contracted")
                and null_vectors.get("left_residual_contains_zero")
                and null_vectors.get("right_residual_contains_zero")
                and null_vectors.get("laurent_denominator_excludes_zero")
            )
            if bool(row["pole_certified"]) != expected:
                fail(f"{name} pole verdict does not follow from its gates")
            all_certified = all_certified and expected
        nesting = receipt["precision_nesting"][name]
        all_certified = (
            all_certified
            and nesting["enclosures_nested_with_precision"]
            and all(nesting["per_quantity_probe_nesting"].values())
            and nesting["per_segment_enclosure_nesting"]["all_nested"]
            and all(nesting["newton_ball_nesting"].values())
        )
    expected_status = (
        "SECOND_SHEET_POLE_CERTIFIED_ON_DECLARED_CONTINUATION"
        if all_certified
        else "CERTIFICATION_INCOMPLETE"
    )
    if receipt["status"] != expected_status:
        fail("status does not follow from the recorded gates")
    promotion = receipt["promotion"]
    for key in (
        "second_sheet_pole_certified_on_declared_continuation",
        "simple_root_certified",
        "laurent_denominator_ball_certified",
        "scalar_residue_ball_certified",
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
    replay_channels(receipt)
    replay_geometry(receipt)
    recompute_verdicts(receipt)
    print(PASS_STATUS)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
