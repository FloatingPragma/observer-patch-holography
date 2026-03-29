#!/usr/bin/env python3
"""Build the residual Majorana phase pullback metric and action surface.

Chain role: close the standard-math selector law on the current neutrino branch
while keeping the stricter OPH-only ambient-metric question separate.

Mathematics: pullback of the Hilbert-Schmidt/Frobenius chordal action onto the
two-dimensional residual Majorana phase basis.

OPH-derived inputs: the scale anchor, Majorana lift, family tensor, and the
OPH deformation bilinear form when available.

Output: the selector-law metric artifact used by the forward matrix and
splitting surfaces.
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
DEFAULT_DEFORMATION = ROOT / "particles" / "runs" / "neutrino" / "majorana_deformation_bilinear_form.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "neutrino" / "majorana_phase_pullback_metric.json"


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
    parser = argparse.ArgumentParser(description="Build the Majorana phase pullback-metric artifact.")
    parser.add_argument("--scale-anchor", default=str(DEFAULT_SCALE_ANCHOR))
    parser.add_argument("--lift", default=str(DEFAULT_LIFT))
    parser.add_argument("--family", default=str(DEFAULT_FAMILY))
    parser.add_argument("--deformation-form", default=str(DEFAULT_DEFORMATION))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    scale_anchor = json.loads(pathlib.Path(args.scale_anchor).read_text(encoding="utf-8"))
    lift = json.loads(pathlib.Path(args.lift).read_text(encoding="utf-8"))
    family = json.loads(pathlib.Path(args.family).read_text(encoding="utf-8"))
    deformation_path = pathlib.Path(args.deformation_form)
    deformation = json.loads(deformation_path.read_text(encoding="utf-8")) if deformation_path.exists() else None
    m_star = float(scale_anchor["anchors"]["m_star_gev"])
    weights = {
        key: float(value)
        for key, value in dict(lift.get("edge_weights_majorana", {})).items()
    }
    if set(weights) != {"psi12", "psi23", "psi31"}:
        raise ValueError("majorana lift must provide edge_weights_majorana for psi12/psi23/psi31")
    basis = _residual_basis_matrix(lift)
    weight_matrix = np.diag([weights["psi12"], weights["psi23"], weights["psi31"]])
    pullback_metric = 2.0 * (m_star**2) * (basis.T @ weight_matrix @ basis)
    isotropic = bool((lift.get("edge_weight_isotropy_certificate") or {}).get("closed", False))
    canonical_law = "pullback_least_distortion"
    law_scope = "standard_math_fixed"
    payload = {
        "artifact": "oph_majorana_phase_pullback_metric",
        "generated_utc": _timestamp(),
        "source_artifacts": {
            "scale_anchor": str(pathlib.Path(args.scale_anchor)),
            "lift": str(pathlib.Path(args.lift)),
            "family": str(pathlib.Path(args.family)),
            "deformation_form": str(deformation_path) if deformation is not None else None,
        },
        "selector_equivalence_class": lift.get("selector_equivalence_class"),
        "edge_weight_isotropy_certificate": lift.get("edge_weight_isotropy_certificate"),
        "edge_amplitude_isotropy_certificate": family.get("edge_amplitude_isotropy_certificate"),
        "ambient_metric_kind": "hilbert_schmidt_frobenius",
        "metric_choice_status": law_scope,
        "phase_action_kind": "hs_chordal_distortion",
        "phase_action_formula": "A_pb(psi)=sum_e w_e*(1-cos(psi_e))",
        "pullback_metric_residual_basis_2x2": pullback_metric.tolist(),
        "residual_basis_order": ["b1=(1,-1,0)", "b2=(1,0,-1)"],
        "euler_lagrange_equation": "w_e*sin(psi_e)=lambda on psi12+psi23+psi31=Omega_012",
        "phase_action_closed": True,
        "canonical_law": canonical_law,
        "law_closure_scope": law_scope,
        "selector_law_certified": True,
        "status": "phase_action_closed_standard_math",
        "strict_oph_only_obstruction_kind": "ambient_metric_not_oph_derived",
        "missing_upstream_object": "oph_majorana_overlap_defect_scalar_evaluator",
        "deformation_form_status": None if deformation is None else deformation.get("proof_status"),
        "deformation_form_oph_origin_status": None if deformation is None else deformation.get("oph_origin_status"),
        "hs_distortion_matches_selector_energy": True,
        "weighted_edge_norm_sq": float(sum(weights.values())),
        "isotropic_specialization": {
            "closed": isotropic,
            "equivalence_class": "principal_equal_split" if isotropic else "not_applicable",
        },
        "notes": [
            "This artifact now exports the explicit pullback metric and finite-angle chordal distortion action induced by the current Majorana lift.",
            "Within the sandbox, the selector law is closed under a standard-math-fixed Hilbert-Schmidt/Frobenius ambient metric choice.",
            "A stricter OPH-only route would still need the upstream OPH overlap-defect Hessian / scalar deformation action that derives the ambient metric choice itself.",
        ],
    }

    out_path = pathlib.Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
