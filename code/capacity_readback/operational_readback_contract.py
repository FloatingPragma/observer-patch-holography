#!/usr/bin/env python3
"""Independent operational-resolution estimator.

This module does not define the capacity map F.  It evaluates a frozen public
scale-discrimination experiment and compares its geometric area-law estimate
with a separately constructed public stable-record capacity.
"""
from __future__ import annotations

from fractions import Fraction
import json
import math
from pathlib import Path
from typing import Iterable, Mapping, Sequence


HERE = Path(__file__).resolve().parent
FAIL_STATES = {
    "NO_STABLE_OBSERVER_SECTOR",
    "NO_SCALE_CALIBRATION",
    "NON_CANONICAL_READBACK_PROTOCOL",
    "NO_OPERATIONAL_RHO_ESTIMATE",
    "INCOMPLETE_TERMINAL_FIBER",
    "AMBIGUOUS_RHO_ESTIMATE",
    "RHO_DERIVED_FROM_CAPACITY",
    "CAPACITY_USED_TO_SELECT_RHO_PROTOCOL",
    "TARGET_TAINTED",
}


def _validate_probability(value: Fraction, label: str) -> None:
    if not 0 <= value <= 1:
        raise ValueError(f"{label} must lie in [0,1]")


def discrimination_error(p0: Sequence[Fraction], p1: Sequence[Fraction]) -> Fraction:
    """Equal-prior optimal classical discrimination error."""
    if len(p0) != len(p1) or not p0:
        raise ValueError("record distributions must have one common nonzero size")
    if sum(p0) != 1 or sum(p1) != 1 or min((*p0, *p1)) < 0:
        raise ValueError("record distributions must be normalized and nonnegative")
    l1 = sum(abs(a - b) for a, b in zip(p0, p1, strict=True))
    return Fraction(1, 2) * (1 - Fraction(1, 2) * l1)


def public_error(
    observers: Iterable[Mapping[str, Fraction | Sequence[Fraction] | str]],
) -> Fraction:
    """Worst-observer capped union bound for a frozen scale test.

    Distribution stage is mandatory so checkpoint or reread failures cannot be
    counted both inside the distributions and as separate terms.
    """
    errors: list[Fraction] = []
    for row in observers:
        checkpoint_error = row["checkpoint_error"]
        reread_error = row["reread_error"]
        if not isinstance(checkpoint_error, Fraction) or not isinstance(
            reread_error, Fraction
        ):
            raise TypeError("checkpoint and reread errors must be Fractions")
        _validate_probability(checkpoint_error, "checkpoint_error")
        _validate_probability(reread_error, "reread_error")
        stage = row.get("distribution_stage")
        if stage not in {"pre_checkpoint", "post_checkpoint"}:
            raise ValueError("distribution_stage must be pre_checkpoint or post_checkpoint")
        if stage == "post_checkpoint" and (checkpoint_error or reread_error):
            raise ValueError(
                "post-checkpoint distributions require zero separately added channel errors"
            )
        errors.append(
            min(
                Fraction(1),
                discrimination_error(row["p0"], row["p1"])
                + checkpoint_error
                + reread_error,
            )
        )
    if not errors:
        raise ValueError("NO_STABLE_OBSERVER_SECTOR")
    return max(errors)


def finest_stable_scale(
    scale_errors: Mapping[Fraction, Fraction], epsilon_read: Fraction
) -> Fraction | None:
    """Return the finest scale for which it and every coarser scale pass."""
    _validate_probability(epsilon_read, "epsilon_read")
    if any(scale <= 0 for scale in scale_errors):
        raise ValueError("all scale coordinates must be positive")
    for error in scale_errors.values():
        _validate_probability(error, "scale error")
    scales = sorted(scale_errors)
    for scale in scales:
        if all(
            scale_errors[coarser] <= epsilon_read
            for coarser in scales
            if coarser >= scale
        ):
            return scale
    return None


def scalarize_rho_fiber(
    rhos: Sequence[Fraction | None], *, fiber_manifest_complete: bool
) -> dict:
    """Scalarize an estimator only on a complete, nonempty agreeing fiber."""
    if not fiber_manifest_complete:
        return {"status": "INCOMPLETE_TERMINAL_FIBER"}
    if not rhos or any(rho is None for rho in rhos):
        return {"status": "NO_OPERATIONAL_RHO_ESTIMATE"}
    defined = [rho for rho in rhos if rho is not None]
    if any(rho <= 0 for rho in defined):
        return {"status": "NO_SCALE_CALIBRATION"}
    distinct = sorted(set(defined))
    if len(distinct) != 1:
        return {
            "status": "AMBIGUOUS_RHO_ESTIMATE",
            "rho_fiber_diameter": str(distinct[-1] - distinct[0]),
            "rho_values": [str(value) for value in distinct],
        }
    return {
        "status": "PASS",
        "rho_op": str(distinct[0]),
        "rho_fiber_diameter": "0",
        "defines_F": False,
    }


def compare_with_public_capacity(
    stable_public_code_size: int,
    rho_op: Fraction,
    *,
    rho_derived_from_capacity: bool,
    capacity_used_to_select_protocol: bool = False,
) -> dict:
    """Compare independent code and geometric estimates without defining F."""
    if not isinstance(stable_public_code_size, int) or stable_public_code_size < 1:
        raise ValueError("stable_public_code_size must be a positive integer")
    if rho_op <= 0:
        return {"status": "NO_SCALE_CALIBRATION"}
    if rho_derived_from_capacity:
        return {"status": "RHO_DERIVED_FROM_CAPACITY"}
    if capacity_used_to_select_protocol:
        return {"status": "CAPACITY_USED_TO_SELECT_RHO_PROTOCOL"}
    pi_coefficient = 1 / (rho_op * rho_op)
    code_value = math.log(stable_public_code_size)
    rho_value = math.pi * float(pi_coefficient)
    return {
        "status": "PASS",
        "defines_F": False,
        "code_capacity_expression": f"log({stable_public_code_size})",
        "rho_capacity_expression": f"{pi_coefficient}*pi",
        "area_record_residual_expression": (
            f"log({stable_public_code_size}) - {pi_coefficient}*pi"
        ),
        "area_record_residual_approximate": code_value - rho_value,
    }


def example_receipt() -> dict:
    """Build a deterministic, nonphysical independent-estimator control."""
    epsilon = Fraction(1, 10)
    scale_rows = {
        Fraction(1, 10): [
            {
                "p0": [Fraction(1, 2), Fraction(1, 2)],
                "p1": [Fraction(1, 2), Fraction(1, 2)],
                "distribution_stage": "pre_checkpoint",
                "checkpoint_error": Fraction(0),
                "reread_error": Fraction(0),
            }
        ],
        Fraction(1, 5): [
            {
                "p0": [Fraction(1), Fraction(0)],
                "p1": [Fraction(0), Fraction(1)],
                "distribution_stage": "pre_checkpoint",
                "checkpoint_error": Fraction(1, 100),
                "reread_error": Fraction(1, 100),
            }
        ],
        Fraction(2, 5): [
            {
                "p0": [Fraction(1), Fraction(0)],
                "p1": [Fraction(0), Fraction(1)],
                "distribution_stage": "pre_checkpoint",
                "checkpoint_error": Fraction(0),
                "reread_error": Fraction(0),
            }
        ],
    }
    errors = {scale: public_error(rows) for scale, rows in scale_rows.items()}
    rho = finest_stable_scale(errors, epsilon)
    assert rho == Fraction(1, 5)
    comparison = compare_with_public_capacity(
        4, rho, rho_derived_from_capacity=False
    )
    return {
        "artifact": "independent_operational_rho_estimator_example_v2",
        "status": "SCHEMA_ONLY",
        "physical_content": False,
        "defines_F": False,
        "moves_cl7": False,
        "cl7_status": "open",
        "capacity_input_accessed_by_estimator": False,
        "source_ancestry": ["declared synthetic record distributions"],
        "epsilon_read": str(epsilon),
        "scale_errors": {str(scale): str(error) for scale, error in errors.items()},
        "rho_fiber_result": scalarize_rho_fiber(
            [rho, rho], fiber_manifest_complete=True
        ),
        "ambiguous_rho_control": scalarize_rho_fiber(
            [rho, Fraction(2, 5)], fiber_manifest_complete=True
        ),
        "incomplete_fiber_control": scalarize_rho_fiber(
            [rho], fiber_manifest_complete=False
        ),
        "independent_capacity_comparison": comparison,
        "rho_from_capacity_negative_control": compare_with_public_capacity(
            4, rho, rho_derived_from_capacity=True
        ),
        "fail_states": sorted(FAIL_STATES),
        "open_estimator_gates": [
            "quotient-derived scale tests reachable from the committed trial family",
            "complete terminal-fiber and stable-observer registries",
            "canonical protocol uniqueness or equivalence theorem",
            "constructor independence from capacity and target labels",
            "physical area-law scale calibration",
            "scale-covariant refinement",
        ],
    }


if __name__ == "__main__":
    payload = example_receipt()
    out = HERE / "runtime" / "operational_capacity_readback_contract_example.json"
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
