#!/usr/bin/env python3
"""Verify the exact-rational finite collar-type gap contract from issue #306."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path
from typing import Any


def fraction(value: Any) -> Fraction:
    if isinstance(value, int):
        return Fraction(value)
    if not isinstance(value, str):
        raise ValueError(f"expected an integer or rational string, got {value!r}")
    return Fraction(value)


def verify(payload: dict[str, Any]) -> dict[str, Any]:
    if payload.get("schema") != "oph.yang_mills.collar_gap.v1":
        raise ValueError("unsupported certificate schema")
    if payload.get("scope") not in {"theorem_contract_witness", "physical_source_receipt"}:
        raise ValueError("scope must distinguish a theorem witness from a physical receipt")

    types = payload.get("types")
    if not isinstance(types, list) or not types:
        raise ValueError("types must be a nonempty list")
    by_id = {entry["id"]: entry for entry in types}
    if len(by_id) != len(types):
        raise ValueError("collar type ids must be unique")

    rates: list[Fraction] = []
    row_sums: dict[str, Fraction] = {}
    for type_id, entry in by_id.items():
        required = {
            "rooted_neighborhood",
            "boundary_sector_flags",
            "local_alphabets",
            "potential_templates",
            "complete_repaired_readback",
            "conditional_kernel",
            "rate_lower",
            "influences",
            "refinement_targets",
        }
        missing = sorted(required - set(entry))
        if missing:
            raise ValueError(f"{type_id}: incomplete source signature: {missing}")
        rate = fraction(entry["rate_lower"])
        if rate <= 0:
            raise ValueError(f"{type_id}: rate lower bound must be positive")
        rates.append(rate)

        influences = entry["influences"]
        if not isinstance(influences, list):
            raise ValueError(f"{type_id}: influences must be a list")
        row_sum = Fraction(0)
        for influence in influences:
            if influence["target_type"] not in by_id:
                raise ValueError(f"{type_id}: influence has unknown target type")
            bound = fraction(influence["upper"])
            multiplicity = int(influence.get("multiplicity", 1))
            if bound < 0 or multiplicity < 0:
                raise ValueError(f"{type_id}: influence bounds must be nonnegative")
            rows = influence.get("conditional_rows")
            if not isinstance(rows, list) or len(rows) != 2 or len(rows[0]) != len(rows[1]):
                raise ValueError(f"{type_id}: influence needs two equal-length conditional rows")
            left = [fraction(value) for value in rows[0]]
            right = [fraction(value) for value in rows[1]]
            if any(value < 0 for value in left + right):
                raise ValueError(f"{type_id}: conditional probabilities must be nonnegative")
            if sum(left) != 1 or sum(right) != 1:
                raise ValueError(f"{type_id}: conditional rows must each sum exactly to one")
            exact_tv = sum(abs(a - b) for a, b in zip(left, right)) / 2
            if exact_tv > bound:
                raise ValueError(
                    f"{type_id}: claimed influence {bound} is below exact TV {exact_tv}"
                )
            row_sum += multiplicity * bound
        row_sums[type_id] = row_sum

        targets = entry["refinement_targets"]
        if not targets or any(target not in by_id for target in targets):
            raise ValueError(f"{type_id}: refinement transition leaves the finite type table")

    c_floor = min(rates)
    eta_upper = max(row_sums.values())
    if eta_upper >= 1:
        raise ValueError(f"influence upper bound is not contractive: eta={eta_upper}")
    gap_lower = c_floor * (1 - eta_upper)

    expected = payload.get("expected", {})
    checks = {
        "c_floor": c_floor,
        "eta_upper": eta_upper,
        "gap_lower": gap_lower,
    }
    for key, actual in checks.items():
        if key in expected and fraction(expected[key]) != actual:
            raise ValueError(f"expected {key}={expected[key]}, computed {actual}")

    return {
        "valid": True,
        "scope": payload["scope"],
        "physical_clay_receipt": payload["scope"] == "physical_source_receipt",
        "type_count": len(types),
        "c_floor": str(c_floor),
        "eta_upper": str(eta_upper),
        "approximate_tensorization_upper": str(1 / (1 - eta_upper)),
        "gap_lower": str(gap_lower),
        "row_sums": {key: str(value) for key, value in sorted(row_sums.items())},
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("certificate", type=Path)
    args = parser.parse_args()
    payload = json.loads(args.certificate.read_text())
    print(json.dumps(verify(payload), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
