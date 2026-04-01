#!/usr/bin/env python3
"""Emit the theorem-grade same-family D12 quark mass ray object.

Chain role: package the exact one-parameter same-family D12 mass-side object
already exposed by the current continuation diagnostics, without claiming an
intrinsic scale law on that ray.

Mathematics: the present D12 same-family branch emits the ray
``(Delta_ud_overlap, eta_Q_centered) = ray_modulus * (1/5, -((1 - x2^2) / 27))``
with the induced odd/even transport formulas inherited from the overlap and
quadratic shells.

OPH-derived inputs: the one-scalar D12 specialization and the current D12
mass-branch diagnostic.

Output: one exact emitted object ``D12_ud_mass_ray`` with unresolved coordinate
``ray_modulus`` and the smaller next object ``quark_d12_t1_value_law``.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
ONE_SCALAR_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_d12_one_scalar_specialization.json"
D12_MASS_BRANCH_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_d12_mass_branch_and_ckm_residual.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "quark_d12_mass_ray.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the emitted same-family D12 quark mass-ray artifact.")
    parser.add_argument("--one-scalar", default=str(ONE_SCALAR_JSON))
    parser.add_argument("--mass-branch", default=str(D12_MASS_BRANCH_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    one_scalar = _load_json(Path(args.one_scalar))
    mass_branch = _load_json(Path(args.mass_branch))

    sample_same_family_point = dict(one_scalar["sample_same_family_point"])
    sample_same_family_mass_point = dict(one_scalar["sample_same_family_mass_point"])
    specialization_formulas = dict(one_scalar["specialization_formulas"])
    comparison_only_best_same_family_point = dict(mass_branch["comparison_only_best_same_family_point"])
    ray_modulus = float(sample_same_family_point["ray_modulus"])
    x2 = float(sample_same_family_point["x2"])
    delta_per_ray_modulus = float(sample_same_family_point["Delta_ud_overlap"] / ray_modulus)
    eta_per_ray_modulus = float(sample_same_family_point["eta_Q_centered"] / ray_modulus)
    eta_per_delta = float(
        sample_same_family_point["eta_Q_centered"] / sample_same_family_point["Delta_ud_overlap"]
    )

    artifact = {
        "artifact": "oph_quark_d12_mass_ray",
        "generated_utc": _timestamp(),
        "scope": "D12_continuation_only",
        "proof_status": "same_family_mass_ray_emitted_modulus_open",
        "public_promotion_allowed": False,
        "theorem_statement": (
            "On the present same-family D12 continuation branch, the mass-side deformation "
            "already emits the one-parameter ray D12_ud_mass_ray with "
            "(Delta_ud_overlap, eta_Q_centered) = ray_modulus * "
            "(1/5, -((1 - x2^2) / 27)). "
            "The emitted ray already identifies its unresolved coordinate with one scalar, "
            "ray_modulus = t1, so the next exact object is the one-scalar value law on this "
            "emitted ray. The larger intrinsic_scale_law_D12 wrapper is derived from that scalar."
        ),
        "emitted_object": {
            "id": "D12_ud_mass_ray",
            "carrier": ["Delta_ud_overlap", "eta_Q_centered"],
            "ray_parameter": "ray_modulus",
            "unresolved_coordinate": "ray_modulus",
            "ray_modulus_equals_t1": True,
        },
        "same_family_ray": {
            "ray_name": "D12_ud_mass_ray",
            "ray_parameter": "ray_modulus",
            "ray_formulas": {
                "Delta_ud_overlap": specialization_formulas["Delta_ud_overlap"],
                "eta_Q_centered": specialization_formulas["eta_Q_centered"],
                "kappa_Q": specialization_formulas["kappa_Q"],
                "tau_u_log_per_side": specialization_formulas["tau_u_log_per_side"],
                "tau_d_log_per_side": specialization_formulas["tau_d_log_per_side"],
                "Lambda_ud_B_transport": specialization_formulas["Lambda_ud_B_transport"],
                "quadratic_even_log_formula_direct": specialization_formulas["quadratic_even_log_formula_direct"],
            },
            "ray_direction_values": {
                "delta_ud_overlap_per_ray_modulus": delta_per_ray_modulus,
                "eta_q_centered_per_ray_modulus": eta_per_ray_modulus,
                "eta_q_centered_per_delta_ud_overlap": eta_per_delta,
            },
            "ordered_family_coordinate_x2": x2,
        },
        "sample_same_family_point": sample_same_family_point,
        "sample_same_family_mass_point": sample_same_family_mass_point,
        "diagnostic_compare_only_best_same_family_point": comparison_only_best_same_family_point,
        "next_exact_missing_object": "quark_d12_t1_value_law",
        "t1_value_law_contract": {
            "id": "quark_d12_t1_value_law",
            "must_emit": "quark_d12_t1_value_law",
            "scalar_name": "t1",
            "identifies": "ray_modulus = t1",
            "unique_intersection_with": "D12_ud_mass_ray",
            "then_emits": ["t1", "ray_modulus", "Delta_ud_overlap", "eta_Q_centered"],
            "must_not_use_target_masses": True,
            "must_not_use_ckm_cp": True,
        },
        "intrinsic_scale_law_contract": {
            "id": "intrinsic_scale_law_D12",
            "derived_from": "quark_d12_t1_value_law",
            "unique_intersection_with": "D12_ud_mass_ray",
            "then_emits": ["ray_modulus", "Delta_ud_overlap", "eta_Q_centered"],
            "must_not_use_target_masses": True,
            "must_not_use_ckm_cp": True,
        },
        "notes": [
            "This artifact is the exact emitted same-family D12 mass object, not a compare-only point fit.",
            "The retained numerical same-family point is still sample-only; it witnesses one point on the emitted ray but does not fix ray_modulus intrinsically.",
            "The current D12 branch can therefore be discussed cleanly in three layers: the ray is emitted, the primitive open object is the one-scalar value law quark_d12_t1_value_law, and intrinsic_scale_law_D12 is the derived wrapper on top of that scalar.",
        ],
    }

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
