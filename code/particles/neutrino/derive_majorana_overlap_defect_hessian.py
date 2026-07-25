#!/usr/bin/env python3
"""Export the upstream OPH overlap-defect Hessian candidate for the Majorana lift."""

from __future__ import annotations

import argparse
import json
import pathlib
import sys
from datetime import datetime, timezone

import numpy as np


ROOT = pathlib.Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from particles.artifact_paths import portable_json_dumps
DEFAULT_LIFT = ROOT / "particles" / "runs" / "neutrino" / "majorana_holonomy_lift.json"
DEFAULT_FAMILY = ROOT / "particles" / "runs" / "neutrino" / "family_response_tensor.json"
DEFAULT_ACTION_GERM = ROOT / "particles" / "runs" / "neutrino" / "majorana_overlap_defect_action_germ.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "neutrino" / "majorana_overlap_defect_hessian.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _residual_basis_matrix(lift: dict[str, object]) -> np.ndarray:
    basis_payload = list(lift.get("residual_basis") or [])
    if len(basis_payload) != 2:
        raise ValueError("majorana lift must provide exactly two residual basis vectors")
    basis = np.zeros((3, 2), dtype=float)
    keys = ("psi12", "psi23", "psi31")
    for column, item in enumerate(basis_payload):
        vector = dict(item)
        for row, key in enumerate(keys):
            basis[row, column] = float(vector[key])
    return basis


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the Majorana overlap-defect Hessian candidate artifact.")
    parser.add_argument("--lift", default=str(DEFAULT_LIFT))
    parser.add_argument("--family", default=str(DEFAULT_FAMILY))
    parser.add_argument("--action-germ", default=str(DEFAULT_ACTION_GERM))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    lift = json.loads(pathlib.Path(args.lift).read_text(encoding="utf-8"))
    family = json.loads(pathlib.Path(args.family).read_text(encoding="utf-8"))
    action_germ_path = pathlib.Path(args.action_germ)
    action_germ = json.loads(action_germ_path.read_text(encoding="utf-8")) if action_germ_path.exists() else None

    basis = _residual_basis_matrix(lift)
    diag_entries = np.asarray(family.get("majorana_normalization_diag", [1.0, 1.0, 1.0]), dtype=float)
    e_nu = np.asarray(family["E_nu"], dtype=float)
    edge_keys = ("psi12", "psi23", "psi31")
    edge_pairs = ((0, 1), (1, 2), (2, 0))
    edge_coefficients = {
        key: float((diag_entries[i] * diag_entries[j] * e_nu[i, j]) ** 2)
        for key, (i, j) in zip(edge_keys, edge_pairs, strict=True)
    }
    mean_weight = float(np.mean(list(edge_coefficients.values())))
    ambient_hessian = mean_weight * np.eye(3, dtype=float)
    residual_hessian = basis.T @ ambient_hessian @ basis
    selector_point = dict(lift.get("canonical_selector_point") or {})
    if not selector_point:
        selector_point = {
            "psi12": float(lift.get("Omega_012", 0.0)) / 3.0,
            "psi23": float(lift.get("Omega_012", 0.0)) / 3.0,
            "psi31": float(lift.get("Omega_012", 0.0)) / 3.0,
        }

    if action_germ is not None:
        residual_hessian = mean_weight * np.asarray(action_germ.get("hessian_class_residual_2x2", [[2.0, 1.0], [1.0, 2.0]]), dtype=float)
    payload = {
        "artifact": "oph_majorana_overlap_defect_hessian",
        "generated_utc": _timestamp(),
        "source_artifacts": {
            "lift": str(pathlib.Path(args.lift)),
            "family": str(pathlib.Path(args.family)),
            "action_germ": str(action_germ_path) if action_germ is not None else None,
        },
        "oph_scalar_action_kind": "local_quadratic_action_germ",
        "selector_point": selector_point,
        "ambient_phase_basis": ["psi12", "psi23", "psi31"],
        "ambient_hessian_3x3": ambient_hessian.tolist(),
        "residual_basis_order": ["b1=(1,-1,0)", "b2=(1,0,-1)"],
        "residual_hessian_2x2": residual_hessian.tolist(),
        "edge_coefficients": edge_coefficients,
        "candidate_scale_from_edge_coefficients": mean_weight,
        "longitudinal_null_certificate": {
            "status": "closed",
            "vector": [1.0, 1.0, 1.0],
            "drops_on_residual_plane": True,
        },
        "isotropic_branch_status": "local_quadratic_germ_closed",
        "nonisotropic_formula_status": "open",
        "proof_status": "local_quadratic_germ_closed",
        "proof_scope": "conditional_algebra_on_declared_family_response_and_majorana_lift",
        "public_surface_candidate_allowed": False,
        "source_only_physical_input_eligible": False,
        "source_closure_status": {
            "closed": False,
            "missing_objects": [
                "source_derived_neutrino_operator_and_family_response",
                "source_derived_overlap_defect_action",
            ],
        },
        "oph_origin_status": "not_established_from_source",
        "upstream_missing_object": "source_closed_neutrino_operator_basis_and_mass_label_contract",
        "upstream_missing_objects": [
            "source_derived_neutrino_operator_and_family_response",
            "source_derived_overlap_defect_action",
            "oph_neutrino_attachment_bridge_invariant",
        ],
        "primitive_metric_source": "oph_overlap_defect_candidate",
        "notes": [
            "This artifact records the local quadratic action-germ/Hessian class conditional on the declared family-response and Majorana-lift inputs.",
            "It does not establish a source-derived physical neutrino operator or overlap-defect action. A positive attachment bridge becomes meaningful only after those source, basis, and mass-label gates close.",
        ],
    }

    out_path = pathlib.Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(portable_json_dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
