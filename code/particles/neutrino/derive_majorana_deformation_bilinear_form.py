#!/usr/bin/env python3
"""Build the OPH-side deformation bilinear form for the Majorana lane.

Chain role: expose the OPH-derived residual quadratic form behind the neutrino
phase dynamics, distinct from the standard-math ambient metric closure.

Mathematics: residual-basis pullback of the deformation/Hessian class on the
three-edge Majorana phase space.

OPH-derived inputs: the local scale anchor, family tensor, Majorana lift, and
the overlap-defect Hessian when it is available.

Output: the OPH-only bilinear-form boundary consumed by the pullback-metric
comparison and promotion-gate logic.
"""

from __future__ import annotations

import argparse
import json
import pathlib
from datetime import datetime, timezone

import numpy as np


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_SCALE_ANCHOR = ROOT / "particles" / "runs" / "neutrino" / "neutrino_scale_anchor.json"
DEFAULT_LIFT = ROOT / "particles" / "runs" / "neutrino" / "majorana_holonomy_lift.json"
DEFAULT_FAMILY = ROOT / "particles" / "runs" / "neutrino" / "family_response_tensor.json"
DEFAULT_HESSIAN = ROOT / "particles" / "runs" / "neutrino" / "majorana_overlap_defect_hessian.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "neutrino" / "majorana_deformation_bilinear_form.json"


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
    parser = argparse.ArgumentParser(description="Build the Majorana deformation bilinear-form artifact.")
    parser.add_argument("--scale-anchor", default=str(DEFAULT_SCALE_ANCHOR))
    parser.add_argument("--lift", default=str(DEFAULT_LIFT))
    parser.add_argument("--family", default=str(DEFAULT_FAMILY))
    parser.add_argument("--hessian", default=str(DEFAULT_HESSIAN))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    scale_anchor = json.loads(pathlib.Path(args.scale_anchor).read_text(encoding="utf-8"))
    lift = json.loads(pathlib.Path(args.lift).read_text(encoding="utf-8"))
    family = json.loads(pathlib.Path(args.family).read_text(encoding="utf-8"))
    hessian_path = pathlib.Path(args.hessian)
    hessian = json.loads(hessian_path.read_text(encoding="utf-8")) if hessian_path.exists() else None

    basis = _residual_basis_matrix(lift)
    diag_entries = np.asarray(family.get("majorana_normalization_diag", [1.0, 1.0, 1.0]), dtype=float)
    e_nu = np.asarray(family["E_nu"], dtype=float)
    edge_keys = ("psi12", "psi23", "psi31")
    edge_pairs = ((0, 1), (1, 2), (2, 0))
    general_branch_weights = {}
    for key, (i, j) in zip(edge_keys, edge_pairs, strict=True):
        general_branch_weights[key] = float((diag_entries[i] * diag_entries[j] * e_nu[i, j]) ** 2)

    isotropic_closed = bool((lift.get("edge_weight_isotropy_certificate") or {}).get("closed", False))
    isotropic_weight = float(np.mean(list(general_branch_weights.values())))
    ambient_class_matrix = isotropic_weight * np.eye(3, dtype=float)
    residual_class_matrix = basis.T @ ambient_class_matrix @ basis
    if hessian is not None:
        ambient_class_matrix = np.asarray(hessian.get("ambient_hessian_3x3", ambient_class_matrix.tolist()), dtype=float)
        residual_class_matrix = np.asarray(hessian.get("residual_hessian_2x2", residual_class_matrix.tolist()), dtype=float)
    longitudinal = np.ones((3, 1), dtype=float)
    longitudinal_drop = basis.T @ (longitudinal @ longitudinal.T) @ basis

    payload = {
        "artifact": "oph_majorana_deformation_bilinear_form",
        "generated_utc": _timestamp(),
        "source_artifacts": {
            "scale_anchor": str(pathlib.Path(args.scale_anchor)),
            "lift": str(pathlib.Path(args.lift)),
            "family": str(pathlib.Path(args.family)),
            "hessian": str(hessian_path) if hessian is not None else None,
        },
        "ambient_phase_basis": ["psi12", "psi23", "psi31"],
        "proof_status": "local_quadratic_germ_closed",
        "oph_origin_status": "missing_scalar_evaluator",
        "law_closure_scope": "local_quadratic_germ_closed",
        "isotropic_branch_rigidity": isotropic_closed,
        "general_branch_weight_formula": "(n_i n_j E_ij)^2",
        "general_branch_weights": general_branch_weights,
        "ambient_metric_class_parameters": {
            "a": isotropic_weight,
            "b": 0.0,
            "class": "aI_3 + bJ_3",
        },
        "ambient_metric_class_3x3": ambient_class_matrix.tolist(),
        "residual_metric_class_2x2": residual_class_matrix.tolist(),
        "residual_basis_order": ["b1=(1,-1,0)", "b2=(1,0,-1)"],
        "longitudinal_piece_on_residual_plane_2x2": longitudinal_drop.tolist(),
        "selector_law_target": "w_e*sin(psi_e)=lambda on psi12+psi23+psi31=Omega_012",
        "scale_anchor_gev": float(scale_anchor["anchors"]["m_star_gev"]),
        "upstream_overlap_defect_hessian": None if hessian is None else {
            "artifact": hessian.get("artifact"),
            "proof_status": hessian.get("proof_status"),
            "oph_origin_status": hessian.get("oph_origin_status"),
            "oph_scalar_action_kind": hessian.get("oph_scalar_action_kind"),
            "upstream_missing_object": hessian.get("upstream_missing_object"),
        },
        "notes": [
            "On the current isotropic branch, any S3-equivariant positive bilinear form restricts to the unique residual class proportional to [[2,1],[1,2]].",
            "The remaining OPH-only burden is no longer the local quadratic class itself; it is the intrinsic scalar evaluator that upgrades the action germ to an exact scaled/finite-angle OPH action.",
        ],
    }

    out_path = pathlib.Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
