#!/usr/bin/env python3
"""Map source-side charged coefficients into the active excitation-gap carrier.

Chain role: translate the ordered charged package into the two-scalar gap data
used by the charged log-spectrum readout.

Mathematics: ordered-family linear/quadratic basis vectors, centered diagonal
shape data, and gap-pair reconstruction.

OPH-derived inputs: the charged mean certificate, family coordinate data, and
the source-side ordered-package coefficients.

Output: the charged excitation-gap map and the next missing charged scalar if
the support extension has not yet been emitted.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MEAN_CERT = ROOT / "particles" / "runs" / "flavor" / "charged_mean_eigenvalue_certificate.json"
DEFAULT_FAMILY = ROOT / "particles" / "runs" / "flavor" / "family_excitation_evaluator.json"
DEFAULT_ORDERED_PACKAGE_VALUE_LAW = ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_ordered_package_value_law.json"
DEFAULT_SUPPORT_EXTENSION_EMITTER = ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_minimal_source_support_extension_emitter.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "leptons" / "lepton_excitation_gap_map.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(
    mean_certificate: dict,
    family: dict,
    ordered_package_value_law: dict | None = None,
    support_extension_emitter: dict | None = None,
) -> dict:
    x2 = float(family["family_coordinate_x"][1])
    l_ord = [
        -(3.0 + x2) / 3.0,
        (2.0 * x2) / 3.0,
        (3.0 - x2) / 3.0,
    ]
    q_ord = [
        (1.0 - x2 * x2) / 3.0,
        -(2.0 * (1.0 - x2 * x2)) / 3.0,
        (1.0 - x2 * x2) / 3.0,
    ]
    ordered_package_value_law = ordered_package_value_law or {}
    support_extension_emitter = support_extension_emitter or {}
    source_coefficients_available = (
        ordered_package_value_law.get("sigma_source_total_log_per_side_readback") is not None
        and ordered_package_value_law.get("eta_source_split_log_per_side_readback") is not None
    )
    sigma_emitted = ordered_package_value_law.get("sigma_source_total_log_per_side_readback")
    eta_emitted = ordered_package_value_law.get("eta_source_split_log_per_side_readback")
    ordered_package = ordered_package_value_law.get("source_side_ordered_package_log_per_side_emitted", [None, None, None])
    mean_log = ordered_package_value_law.get("mu_source_log_per_side_readback")
    centered_logs = (
        [float(value) - float(mean_log) for value in ordered_package]
        if source_coefficients_available and mean_log is not None
        else [None, None, None]
    )
    shape_singular_values = (
        [float(math.exp(value)) for value in centered_logs]
        if source_coefficients_available
        else [None, None, None]
    )
    y_e_shape = (
        {
            "real": [
                [shape_singular_values[0], 0.0, 0.0],
                [0.0, shape_singular_values[1], 0.0],
                [0.0, 0.0, shape_singular_values[2]],
            ],
            "imag": [[0.0, 0.0, 0.0] for _ in range(3)],
        }
        if source_coefficients_available
        else None
    )
    return {
        "artifact": "oph_charged_lepton_excitation_gap_map",
        "generated_utc": _timestamp(),
        "parent_artifact": mean_certificate.get("artifact"),
        "input_kind": "ordered_branch_generator_spectral_package",
        "theorem_scope": "current_family_only",
        "proof_status": (
            ordered_package_value_law.get("proof_status")
            if source_coefficients_available
            else "carrier_family_closed_coefficients_missing"
        ),
        "carrier_family_status": "closed",
        "source_side_ordered_package_value_law_artifact": ordered_package_value_law.get("artifact"),
        "source_side_ordered_package_value_law_status": ordered_package_value_law.get("proof_status"),
        "ordered_family_coordinate": [-1.0, x2, 1.0],
        "rho_ord": float(family["rho_ord"]),
        "linear_basis_vector_centered": l_ord,
        "quadratic_basis_vector": q_ord,
        "a_e_log_coeff": (
            None if sigma_emitted is None else float(sigma_emitted) / 2.0
        ),
        "b_e_log_coeff": (
            None
            if eta_emitted is None
            else float(eta_emitted) / (2.0 * (1.0 - x2 * x2))
        ),
        "sigma_e_total_log_per_side_emitted": sigma_emitted,
        "eta_e_split_log_per_side_emitted": eta_emitted,
        "gamma21_log_per_side_emitted": (
            (((1.0 + x2) * float(sigma_emitted) - float(eta_emitted)) / 2.0)
        ) if source_coefficients_available else None,
        "gamma32_log_per_side_emitted": (
            (((1.0 - x2) * float(sigma_emitted) + float(eta_emitted)) / 2.0)
            if source_coefficients_available
            else None
        ),
        "E_e_log_centered_emitted": centered_logs,
        "shape_singular_values_emitted": shape_singular_values,
        "Y_e_shape_emitted": y_e_shape,
        "sigma_e_total_log_per_side_emitted_formula": "2 * a_e_log_coeff",
        "eta_e_split_log_per_side_emitted_formula": "2 * (1 - x2^2) * b_e_log_coeff",
        "centered_linear_quadratic_formula": "E_e_log_centered_emitted = a_e_log_coeff * L_ord + b_e_log_coeff * Q_ord",
        "gamma21_log_per_side_emitted_formula": "((1 + x2) * sigma_e_total_log_per_side_emitted - eta_e_split_log_per_side_emitted) / 2",
        "gamma32_log_per_side_emitted_formula": "((1 - x2) * sigma_e_total_log_per_side_emitted + eta_e_split_log_per_side_emitted) / 2",
        "E_e_log_centered_emitted_formula": [
            "-(2*gamma21_log_per_side_emitted + gamma32_log_per_side_emitted)/3",
            "(gamma21_log_per_side_emitted - gamma32_log_per_side_emitted)/3",
            "(gamma21_log_per_side_emitted + 2*gamma32_log_per_side_emitted)/3",
        ],
        "shape_singular_values_emitted_formula": "exp(E_e_log_centered_emitted)",
        "coefficient_recovery": {
            "a_e_log_coeff": "sigma_e_total_log_per_side_emitted / 2",
            "b_e_log_coeff": "eta_e_split_log_per_side_emitted / (2 * (1 - x2^2))",
        },
        "Y_e_shape_formula": "diag(exp(E_e_log_centered_emitted))",
        "downstream_targets": [
            "oph_lepton_log_spectrum_readout",
            "oph_forward_charged_leptons",
        ],
        "smallest_constructive_missing_object": (
            support_extension_emitter.get("smallest_constructive_missing_object")
            or ordered_package_value_law.get("smallest_constructive_missing_object")
            if source_coefficients_available
            else "oph_charged_lepton_excitation_gap_map"
        ),
        "notes": [
            "The ordered three-point support and the centered affine-quadratic parameterization are fixed.",
            "This artifact carries the current source-side value shell when the sector-local ordered package is available.",
            "No fitted charged masses are consumed here.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the charged-lepton excitation-gap map artifact.")
    parser.add_argument("--mean-certificate", default=str(DEFAULT_MEAN_CERT))
    parser.add_argument("--family", default=str(DEFAULT_FAMILY))
    parser.add_argument("--ordered-package-value-law", default=str(DEFAULT_ORDERED_PACKAGE_VALUE_LAW))
    parser.add_argument("--support-extension-emitter", default=str(DEFAULT_SUPPORT_EXTENSION_EMITTER))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    mean_certificate = json.loads(Path(args.mean_certificate).read_text(encoding="utf-8"))
    family = json.loads(Path(args.family).read_text(encoding="utf-8"))
    ordered_package_value_law_path = Path(args.ordered_package_value_law)
    ordered_package_value_law = (
        json.loads(ordered_package_value_law_path.read_text(encoding="utf-8"))
        if ordered_package_value_law_path.exists()
        else None
    )
    support_extension_path = Path(args.support_extension_emitter)
    support_extension_emitter = (
        json.loads(support_extension_path.read_text(encoding="utf-8"))
        if support_extension_path.exists()
        else None
    )
    artifact = build_artifact(mean_certificate, family, ordered_package_value_law, support_extension_emitter)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
