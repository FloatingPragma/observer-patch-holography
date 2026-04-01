#!/usr/bin/env python3
"""Emit the minimal one-scalar D12 quark value-law scaffold.

Chain role: record the smallest target-free theorem object left on the emitted
same-family D12 mass ray.

Mathematics: the emitted ray already identifies the unresolved mass-side
coordinate with one scalar, `ray_modulus = t1`, and the full mass-side wrapper
is then induced by
`Delta_ud_overlap = t1 / 5` and
`eta_Q_centered = -((1 - x2^2) / 27) * t1`.

Output: the primitive one-scalar target `quark_d12_t1_value_law`, with
`intrinsic_scale_law_D12` retained only as a derived wrapper.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
MASS_RAY_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_d12_mass_ray.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "quark_d12_t1_value_law.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def build_payload(mass_ray: dict[str, Any]) -> dict[str, Any]:
    sample = dict(mass_ray["sample_same_family_point"])
    formulas = dict(mass_ray["same_family_ray"]["ray_formulas"])
    return {
        "artifact": "oph_quark_d12_t1_value_law",
        "generated_utc": _timestamp(),
        "scope": "D12_continuation_only",
        "proof_status": "emitted_mass_ray_one_scalar_value_law_open",
        "public_promotion_allowed": False,
        "exact_missing_object": "quark_d12_t1_value_law",
        "parent_emitted_object": {
            "artifact": mass_ray["artifact"],
            "id": mass_ray["emitted_object"]["id"],
        },
        "theorem_statement": (
            "On the emitted same-family D12 mass ray, the unresolved coordinate is already identified "
            "with one scalar t1 through ray_modulus = t1. The exact next target-free theorem object is "
            "therefore the one-scalar value law quark_d12_t1_value_law. The larger wrapper "
            "intrinsic_scale_law_D12 is derived from that single scalar by the emitted ray formulas."
        ),
        "one_scalar_reduction": {
            "scalar_name": "t1",
            "ray_coordinate_identity": "ray_modulus = t1",
            "induced_formulas": {
                "Delta_ud_overlap": "t1 / 5",
                "eta_Q_centered": "-((1 - x2^2) / 27) * t1",
                "kappa_Q": "-t1 / 54",
            },
            "emitted_ray_formulas": formulas,
        },
        "sample_same_family_point": {
            "t1": float(sample["t1_sample"]),
            "ray_modulus": float(sample["ray_modulus"]),
            "x2": float(sample["x2"]),
            "Delta_ud_overlap": float(sample["Delta_ud_overlap"]),
            "eta_Q_centered": float(sample["eta_Q_centered"]),
            "kappa_Q": float(sample["kappa_Q"]),
        },
        "minimal_new_theorem": {
            "id": "quark_d12_t1_value_law",
            "must_emit": "quark_d12_t1_value_law",
            "scalar_name": "t1",
            "unique_intersection_with": "D12_ud_mass_ray",
            "identifies": "ray_modulus = t1",
            "then_emits": ["t1", "ray_modulus", "Delta_ud_overlap", "eta_Q_centered"],
            "must_not_use_target_masses": True,
            "must_not_use_ckm_cp": True,
        },
        "derived_wrapper": {
            "id": "intrinsic_scale_law_D12",
            "derived_from": "quark_d12_t1_value_law",
            "meaning": "the induced mass-side wrapper on D12_ud_mass_ray after t1 is fixed intrinsically",
        },
        "notes": [
            "This artifact resizes the primitive quark frontier to the single unresolved scalar already exposed by the emitted ray.",
            "It does not claim a theorem-grade value for t1 is already emitted on the current corpus.",
            "The larger wrapper intrinsic_scale_law_D12 remains valid language, but only as the derived mass-ray wrapper above quark_d12_t1_value_law.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the D12 quark t1 value-law scaffold.")
    parser.add_argument("--mass-ray", default=str(MASS_RAY_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    payload = build_payload(_load_json(Path(args.mass_ray)))
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
