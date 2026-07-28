#!/usr/bin/env python3
"""Fail-closed exact certificate for a finite collar-projection gap.

This utility verifies a declared finite *calibration* family.  It deliberately
does not turn such a family into a receipt for the physical compact-gauge OPH
tower: that promotion requires a complete physical source table and the
separate continuum, transfer, and OS/noncollapse receipts.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Any


SCHEMA = "oph.yang_mills.collar_gap_certificate.v1"


def rational(value: Any) -> Fraction:
    """Read an exact rational; binary floats are never proof inputs."""
    if isinstance(value, bool) or isinstance(value, float):
        raise ValueError("binary floating-point values are forbidden in certificates")
    if isinstance(value, int):
        return Fraction(value)
    if not isinstance(value, str):
        raise ValueError(f"expected integer or rational string, got {value!r}")
    try:
        return Fraction(value)
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError(f"invalid rational {value!r}") from exc


def canonical(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")


def expand_types(manifest: dict[str, Any]) -> list[dict[str, Any]]:
    """Expand the compact deterministic catalogue used by the calibration fixture."""
    table = manifest.get("type_table")
    if isinstance(table, list):
        return table
    family = manifest.get("calibration_type_family")
    if not isinstance(family, dict):
        raise ValueError("manifest needs an explicit type_table or calibration_type_family")
    count = family.get("count")
    prefix = family.get("id_prefix")
    template = family.get("template")
    if not isinstance(count, int) or count <= 0 or not isinstance(prefix, str) or not isinstance(template, dict):
        raise ValueError("malformed deterministic calibration type family")
    types: list[dict[str, Any]] = []
    for index in range(count):
        entry = dict(template)
        entry["id"] = f"{prefix}{index:03d}"
        entry["refinement_targets"] = [f"{prefix}{index:03d}"]
        for influence in entry.get("influences", []):
            influence = dict(influence)
            influence["target_type"] = entry["id"]
            entry["influences"] = [influence]
        types.append(entry)
    return types


def validate(manifest: dict[str, Any]) -> dict[str, Any]:
    if manifest.get("schema") != SCHEMA:
        raise ValueError("unsupported manifest schema")
    scope = manifest.get("scope")
    if scope not in {"calibration", "physical_source_receipt"}:
        raise ValueError("scope must be calibration or physical_source_receipt")
    types = expand_types(manifest)
    if not types:
        raise ValueError("the active type table is empty")
    ids = [entry.get("id") for entry in types]
    if any(not isinstance(type_id, str) for type_id in ids) or len(set(ids)) != len(ids):
        raise ValueError("type ids must be present and unique")
    known = set(ids)
    row_sums: dict[str, Fraction] = {}
    rates: list[Fraction] = []
    type_hashes: dict[str, str] = {}
    required = {
        "id", "rooted_neighborhood", "boundary_sector_flags", "local_alphabets",
        "potential_templates", "complete_repaired_readback", "conditional_kernel",
        "rate_lower", "influences", "refinement_targets",
    }
    for entry in types:
        missing = sorted(required - set(entry))
        if missing:
            raise ValueError(f"{entry.get('id', '<unknown>')}: incomplete source signature: {missing}")
        rate = rational(entry["rate_lower"])
        if rate <= 0:
            raise ValueError(f"{entry['id']}: rate lower bound must be positive")
        rates.append(rate)
        targets = entry["refinement_targets"]
        if not isinstance(targets, list) or not targets or any(target not in known for target in targets):
            raise ValueError(f"{entry['id']}: refinement leaves the finite active type table")
        total = Fraction(0)
        for influence in entry["influences"]:
            if influence.get("target_type") not in known:
                raise ValueError(f"{entry['id']}: influence references an unknown type")
            upper = rational(influence.get("upper"))
            multiplicity = influence.get("multiplicity", 1)
            if not isinstance(multiplicity, int) or multiplicity < 0 or upper < 0:
                raise ValueError(f"{entry['id']}: invalid influence bound")
            rows = influence.get("conditional_rows")
            if not isinstance(rows, list) or len(rows) != 2 or not all(isinstance(row, list) for row in rows):
                raise ValueError(f"{entry['id']}: influence needs two conditional rows")
            if len(rows[0]) != len(rows[1]) or not rows[0]:
                raise ValueError(f"{entry['id']}: conditional rows have incompatible support")
            left, right = ([rational(value) for value in row] for row in rows)
            if min(left + right) < 0 or sum(left) != 1 or sum(right) != 1:
                raise ValueError(f"{entry['id']}: conditional rows must be probability laws")
            tv = sum(abs(a - b) for a, b in zip(left, right)) / 2
            if tv > upper:
                raise ValueError(f"{entry['id']}: stated influence is below exact TV distance")
            total += multiplicity * upper
        row_sums[entry["id"]] = total
        type_hashes[entry["id"]] = hashlib.sha256(canonical(entry)).hexdigest()
    c_floor = min(rates)
    eta_upper = max(row_sums.values())
    if eta_upper >= 1:
        raise ValueError(f"Dobrushin upper bound must be < 1, got {eta_upper}")
    if scope == "physical_source_receipt":
        missing_physical = sorted({"physical_source_provenance", "continuum_receipts"} - set(manifest))
        if missing_physical:
            raise ValueError("physical receipt is incomplete: " + ", ".join(missing_physical))
    return {
        "schema": SCHEMA,
        "scope": scope,
        "physical_clay_receipt": scope == "physical_source_receipt",
        "active_type_count": len(types),
        "c_floor": str(c_floor),
        "eta_upper": str(eta_upper),
        "approximate_tensorization_upper": str(1 / (1 - eta_upper)),
        "gap_lower": str(c_floor * (1 - eta_upper)),
        "row_sums": {key: str(value) for key, value in sorted(row_sums.items())},
        "type_table_sha256": hashlib.sha256(canonical(types)).hexdigest(),
        "type_hashes": type_hashes,
    }


def certify(manifest_path: Path, output_path: Path) -> None:
    manifest = json.loads(manifest_path.read_text())
    receipt = validate(manifest)
    receipt["manifest_sha256"] = hashlib.sha256(canonical(manifest)).hexdigest()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")


def verify(manifest_path: Path, receipt_path: Path) -> None:
    manifest = json.loads(manifest_path.read_text())
    expected = json.loads(receipt_path.read_text())
    actual = validate(manifest)
    actual["manifest_sha256"] = hashlib.sha256(canonical(manifest)).hexdigest()
    if actual != expected:
        raise ValueError("receipt does not exactly recompute from its manifest")


def main() -> int:
    parser = argparse.ArgumentParser()
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("certify", "verify"):
        command = commands.add_parser(name)
        command.add_argument("--manifest", type=Path, required=True)
        command.add_argument("--output" if name == "certify" else "--receipt", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "certify":
        certify(args.manifest, args.output)
    else:
        verify(args.manifest, args.receipt)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
