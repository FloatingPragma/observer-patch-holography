#!/usr/bin/env python3
"""Export the local quadratic action-germ boundary for the OPH-only Majorana route."""

from __future__ import annotations

import argparse
import json
import pathlib
from datetime import datetime, timezone

import numpy as np


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_LIFT = ROOT / "particles" / "runs" / "neutrino" / "majorana_holonomy_lift.json"
DEFAULT_FAMILY = ROOT / "particles" / "runs" / "neutrino" / "family_response_tensor.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "neutrino" / "majorana_overlap_defect_action_germ.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the Majorana overlap-defect action-germ artifact.")
    parser.add_argument("--lift", default=str(DEFAULT_LIFT))
    parser.add_argument("--family", default=str(DEFAULT_FAMILY))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    lift = json.loads(pathlib.Path(args.lift).read_text(encoding="utf-8"))
    family = json.loads(pathlib.Path(args.family).read_text(encoding="utf-8"))

    diag_entries = np.asarray(family.get("majorana_normalization_diag", [1.0, 1.0, 1.0]), dtype=float)
    e_nu = np.asarray(family["E_nu"], dtype=float)
    edge_pairs = {"psi12": (0, 1), "psi23": (1, 2), "psi31": (2, 0)}
    edge_coefficients = {
        key: float((diag_entries[i] * diag_entries[j] * e_nu[i, j]) ** 2)
        for key, (i, j) in edge_pairs.items()
    }
    candidate_scale = float(np.mean(list(edge_coefficients.values())))
    selector_point = {
        "psi12": float(lift.get("Omega_012", 0.0)) / 3.0,
        "psi23": float(lift.get("Omega_012", 0.0)) / 3.0,
        "psi31": float(lift.get("Omega_012", 0.0)) / 3.0,
    }

    payload = {
        "artifact": "oph_majorana_overlap_defect_action_germ",
        "generated_utc": _timestamp(),
        "source_artifacts": {
            "lift": str(pathlib.Path(args.lift)),
            "family": str(pathlib.Path(args.family)),
        },
        "domain": "affine_majorana_lift",
        "selector_center": "principal_equal_split",
        "selector_point": selector_point,
        "proof_status": "local_quadratic_germ_closed",
        "residual_symmetry_group": "S3_on_A2_reflection_representation",
        "quadratic_invariant_residual": "I2(u,v)=u^2+u*v+v^2",
        "cubic_invariant_obstruction": {
            "status": "open",
            "formula": "I3(u,v)=u*v*(u+v)",
            "eliminated": False,
        },
        "invariant_ring_status": "quadratic_jet_closed_cubic_freedom_open",
        "action_formula_residual": "mu_nu*(u^2 + u*v + v^2)",
        "hessian_class_residual_2x2": [[2.0, 1.0], [1.0, 2.0]],
        "scale_status": "unresolved_without_scalar_evaluator",
        "candidate_scale_from_edge_coefficients": candidate_scale,
        "edge_coefficients": edge_coefficients,
        "upstream_missing_object": "oph_majorana_overlap_defect_scalar_evaluator",
        "notes": [
            "This artifact closes the local quadratic action germ on the affine Majorana lift, but not the exact scalar evaluator that fixes the scale and finite-angle continuation.",
            "The exact remaining no-go is the surviving S3-invariant cubic freedom beyond the closed quadratic jet.",
        ],
    }

    out_path = pathlib.Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
