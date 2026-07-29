#!/usr/bin/env python3
"""Independent fail-closed checker for the sampled W/Z boundary diagnostic.

This checker deliberately does not import the loop evaluator or the diagnostic
producer.  It validates the strict diagnostic schema, recomputes the artifact
and source digests, pins the complete rational fixture, replays the exact
``d = 4 - 2 eps`` finite corrections against fixed symbolic results, and
requires every certification and promotion field to remain false.

A passing checker result means that the committed diagnostic is internally
bound and accurately scoped.  It is not a root-count, Laurent, current-pole,
physical-unit, or OPH-native certificate.
"""

from __future__ import annotations

import copy
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Any

from jsonschema import Draft202012Validator

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT_PATH = ROOT / "outputs" / "wz_pole_receipts.json"
VECTOR_PATH = ROOT / "outputs" / "fj_direct_vector_blocks.json"
SCHEMA_PATH = ROOT / "schemas" / "wz_pole_boundary_diagnostic_v3.schema.json"

PASS_STATUS = "DIAGNOSTIC_CONTRACT_PASS__NO_ROOT_CERTIFICATE"

EXPECTED_FIXTURE = {
    "g1": "1/4",
    "g2": "1/3",
    "v": "2",
    "lam": "1/8",
    "xi": "1",
    "mu2": "1/2",
    "mfu1": "1/50",
    "mfu2": "1/20",
    "mfu3": "1/5",
    "mfd1": "1/60",
    "mfd2": "1/25",
    "mfd3": "1/10",
    "mfe1": "1/80",
    "mfe2": "1/30",
    "mfe3": "1/15",
    "mu_ren2": "1",
}
EXPECTED_CKM = {
    f"V{i}{j}": str(int(i == j))
    for i in (1, 2, 3)
    for j in (1, 2, 3)
}
EXPECTED_CORRECTIONS = {
    "W": "(312000*p2 + 90163)/1944000",
    "Z": "(292020000*p2 + 79268959)/1944000000",
    "AZ": "(8880000*p2 + 11104571)/648000000",
}
EXPECTED_RECEIPT_KEYS = {
    f"{boson}_{precision}"
    for boson in ("W", "Z")
    for precision in (128, 192, 256)
}


def canonical_json_bytes(value: Any) -> bytes:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
    ).encode("utf-8")


def expected_contour(boson: str) -> tuple[str, dict[str, list[str]]]:
    tree = Fraction(1, 9) if boson == "W" else Fraction(25, 144)
    radius = tree / 8
    return str(tree), {
        "lower_left": [str(tree - radius), str(-radius)],
        "upper_right": [str(tree + radius), str(radius)],
    }


def check(
    payload: dict[str, Any] | None = None,
    *,
    vector_path: Path = VECTOR_PATH,
) -> dict[str, Any]:
    if payload is None:
        payload = json.loads(ARTIFACT_PATH.read_text(encoding="utf-8"))
    schema = json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))
    problems: list[str] = []

    schema_errors = sorted(
        Draft202012Validator(schema).iter_errors(payload),
        key=lambda error: list(error.absolute_path),
    )
    problems.extend(
        "schema at "
        + (".".join(str(part) for part in error.absolute_path) or "<root>")
        + f": {error.message}"
        for error in schema_errors
    )

    artifact_copy = copy.deepcopy(payload)
    declared_artifact_digest = artifact_copy.pop("artifact_sha256", None)
    recomputed_artifact_digest = (
        "sha256:" + hashlib.sha256(canonical_json_bytes(artifact_copy)).hexdigest()
    )
    if declared_artifact_digest != recomputed_artifact_digest:
        problems.append("artifact self-digest mismatch")

    vector_raw = vector_path.read_bytes()
    vector = json.loads(vector_raw)
    declared_input = payload.get("input", {})
    expected_input = {
        "path": "outputs/fj_direct_vector_blocks.json",
        "bytes": len(vector_raw),
        "sha256": hashlib.sha256(vector_raw).hexdigest(),
        "schema": vector.get("schema"),
        "target": vector.get("target"),
    }
    if declared_input != expected_input:
        problems.append("source vector path, metadata, or byte digest mismatch")
    if vector.get("schema") != "fj_direct_vector_blocks.v1":
        problems.append("source vector schema is not fj_direct_vector_blocks.v1")
    if vector.get("target") != "FJ_DIRECT_1":
        problems.append("source vector target is not FJ_DIRECT_1")

    if payload.get("fixture") != EXPECTED_FIXTURE:
        problems.append("rational scalar fixture differs from the pinned diagnostic fixture")
    if payload.get("ckm_fixture") != EXPECTED_CKM:
        problems.append("CKM fixture differs from the pinned exact identity matrix")
    corrections = payload.get("conventions", {}).get(
        "dimensional_prefactor_finite_correction",
        {},
    )
    if corrections != EXPECTED_CORRECTIONS:
        problems.append("exact dimensional-prefactor correction map mismatch")

    receipts = payload.get("receipts", {})
    if set(receipts) != EXPECTED_RECEIPT_KEYS:
        problems.append("receipt key set mismatch")
    for key in sorted(EXPECTED_RECEIPT_KEYS & set(receipts)):
        receipt = receipts[key]
        boson, precision_text = key.split("_")
        expected_tree, contour = expected_contour(boson)
        if receipt.get("boson") != boson:
            problems.append(f"{key}: boson label mismatch")
        if receipt.get("precision_bits") != int(precision_text):
            problems.append(f"{key}: precision label mismatch")
        if receipt.get("tree_mass_sq") != expected_tree:
            problems.append(f"{key}: exact tree-mass coordinate mismatch")
        if receipt.get("contour_exact") != contour:
            problems.append(f"{key}: exact contour mismatch")
        expected_binary64 = {
            name: [float(Fraction(value)) for value in pair]
            for name, pair in contour.items()
        }
        if receipt.get("contour_binary64") != expected_binary64:
            problems.append(f"{key}: binary64 contour does not match the exact contour")
        if receipt.get("dimensional_prefactor_finite_correction") != (
            EXPECTED_CORRECTIONS[boson]
        ):
            problems.append(f"{key}: exact finite correction mismatch")
        if receipt.get("samples_per_edge") != 32:
            problems.append(f"{key}: samples-per-edge mismatch")
        if receipt.get("edge_sample_counts") != {
            "bottom": 32,
            "right": 32,
            "top": 32,
            "left": 32,
        }:
            problems.append(f"{key}: four-edge sample counts mismatch")
        if not receipt.get("diagnostic_checks_passed", False):
            problems.append(f"{key}: committed sampled diagnostic did not pass")
        if receipt.get("sampled_boundary_winding") != 1:
            problems.append(f"{key}: sampled winding differs from one")
        if receipt.get("sampled_points_exclude_zero") is not True:
            problems.append(f"{key}: sampled-point exclusion flag is false")
        minimum = receipt.get("sampled_boundary_min_modulus")
        radius = receipt.get("roundoff_heuristic_max_radius")
        if (
            not isinstance(minimum, (int, float))
            or not isinstance(radius, (int, float))
            or minimum <= 8 * radius
        ):
            problems.append(f"{key}: sampled modulus/heuristic-radius relation fails")

    receipt_diagnostics = [
        bool(receipt.get("diagnostic_checks_passed", False))
        for receipt in receipts.values()
        if isinstance(receipt, dict)
    ]
    expected_diagnostic_summary = (
        len(receipt_diagnostics) == len(EXPECTED_RECEIPT_KEYS)
        and all(receipt_diagnostics)
    )
    if payload.get("diagnostic_checks_passed") is not expected_diagnostic_summary:
        problems.append("top-level diagnostic summary does not match the receipts")

    return {
        "schema": "wz_pole_boundary_diagnostic_check.v1",
        "status": PASS_STATUS if not problems else "FAIL",
        "diagnostic_contract_valid": not problems,
        "root_count_certified": False,
        "acceptance_row_discharged": False,
        "promotion_allowed": False,
        "problems": problems,
    }


def main() -> int:
    verdict = check()
    print(json.dumps(verdict, sort_keys=True))
    return 0 if verdict["diagnostic_contract_valid"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
