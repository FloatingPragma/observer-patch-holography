#!/usr/bin/env python3
"""Export the local quadratic action-germ boundary for the OPH-only Majorana route."""

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
        "residual_symmetry_group": "S3_on_A2_reflection_representation",
        "quadratic_invariant_residual": "I2(u,v)=u^2+u*v+v^2",
        "cubic_invariant_obstruction": {
            "status": "closed_on_current_isotropic_branch",
            "formula": "I3(u,v)=u*v*(u+v)",
            "eliminated": True,
        },
        "invariant_ring_status": "current_isotropic_branch_closed_by_edge_norm_theorem",
        "action_formula_residual": "mu_nu*(u^2 + u*v + v^2)",
        "hessian_class_residual_2x2": [[2.0, 1.0], [1.0, 2.0]],
        "scale_status": "unresolved_without_positive_attachment_bridge_invariant",
        "candidate_scale_from_edge_coefficients": candidate_scale,
        "edge_coefficients": edge_coefficients,
        "upstream_missing_object": "source_closed_neutrino_operator_basis_and_mass_label_contract",
        "upstream_missing_objects": [
            "source_derived_neutrino_operator_and_family_response",
            "source_derived_overlap_defect_action",
            "oph_neutrino_attachment_bridge_invariant",
        ],
        "notes": [
            "Conditional on the declared family-response and Majorana-lift inputs, this artifact closes the local quadratic action germ and removes the former cubic freedom on the isotropic branch.",
            "This local algebra does not establish that those inputs descend from a source-closed physical neutrino operator. A bridge invariant is relevant only after the operator, basis, and mass-label contracts close.",
        ],
    }

    out_path = pathlib.Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(portable_json_dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
