#!/usr/bin/env python3
"""Convert ordered charged coefficients into the diagonal log-spectrum branch.

Chain role: build the shape-only charged spectrum from the current ordered gap
carrier before the absolute scale and final hierarchy completion are attached.

Mathematics: ordered gap-pair reconstruction, centered log-spectrum assembly,
and diagonal shape-matrix emission.

OPH-derived inputs: family coordinate data plus the charged source-emission or
gap-map coefficients produced upstream in `/particles`.

Output: the charged log-spectrum artifact consumed by the forward charged
builder and the charged exactness audit.
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
DEFAULT_SOURCE_EMISSION = ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_ordered_package_source_emission.json"
DEFAULT_GAP_MAP = ROOT / "particles" / "runs" / "leptons" / "lepton_excitation_gap_map.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "leptons" / "lepton_log_spectrum_readout.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _diag_payload(diagonal_values: list[float]) -> dict[str, list[list[float]]]:
    real = [[0.0, 0.0, 0.0] for _ in range(3)]
    imag = [[0.0, 0.0, 0.0] for _ in range(3)]
    for idx, value in enumerate(diagonal_values):
        real[idx][idx] = float(value)
    return {"real": real, "imag": imag}


def build_artifact(
    mean_certificate: dict,
    family: dict,
    source_emission: dict | None = None,
    gap_map: dict | None = None,
) -> dict:
    source_emission = source_emission or {}
    if source_emission.get("artifact") == "oph_charged_sector_local_ordered_package_source_emission":
        ordered_eigenvalues = [
            float(value) for value in source_emission["source_side_ordered_package_log_per_side_emitted"]
        ]
        mean_value = float(source_emission["source_side_ordered_package_mean_log_per_side_emitted"])
        source_artifact = source_emission.get("artifact")
    else:
        ordered_eigenvalues = sorted(float(value) for value in mean_certificate["current_family_eigenvalues"])
        mean_value = sum(ordered_eigenvalues) / 3.0
        source_artifact = mean_certificate.get("artifact")
    x2 = float(family["family_coordinate_x"][1])
    rho_ord = float(family["rho_ord"])
    sigma_e = ordered_eigenvalues[2] - ordered_eigenvalues[0]
    gap_map = gap_map or {}
    sigma_emitted = gap_map.get("sigma_e_total_log_per_side_emitted")
    eta_emitted = gap_map.get("eta_e_split_log_per_side_emitted")
    gap_map_consumed = sigma_emitted is not None and eta_emitted is not None
    active_sigma = float(sigma_emitted) if gap_map_consumed else sigma_e
    eta_e_rigid_fallback = (((1.0 + x2) - (rho_ord * (1.0 - x2))) / (1.0 + rho_ord)) * active_sigma
    active_eta = float(eta_emitted) if gap_map_consumed else eta_e_rigid_fallback
    gamma21 = (((1.0 + x2) * active_sigma) - active_eta) / 2.0
    gamma32 = (((1.0 - x2) * active_sigma) + active_eta) / 2.0
    e_e_log_centered = [
        -((2.0 * gamma21) + gamma32) / 3.0,
        (gamma21 - gamma32) / 3.0,
        (gamma21 + (2.0 * gamma32)) / 3.0,
    ]
    shape_singular_values = [math.exp(value) for value in e_e_log_centered]
    q_ord = [
        (1.0 - x2 * x2) / 3.0,
        -(2.0 * (1.0 - x2 * x2)) / 3.0,
        (1.0 - x2 * x2) / 3.0,
    ]

    smallest_missing_object = (
        gap_map.get("smallest_constructive_missing_object")
        if gap_map_consumed
        else "oph_charged_lepton_excitation_gap_map"
    )

    return {
        "artifact": "oph_lepton_log_spectrum_readout",
        "generated_utc": _timestamp(),
        "parent_artifact": mean_certificate.get("artifact"),
        "proof_status": (
            "current_family_ordered_package_readback_candidate"
            if gap_map_consumed
            else "current_family_ordered_ratio_candidate"
        ),
        "predictive_promotion_allowed": False,
        "theorem_scope": "current_family_only",
        "consumer_builder": "build_forward_charged_leptons.py",
        "source_ordered_package_artifact": source_artifact,
        "excitation_gap_map_artifact": gap_map.get("artifact"),
        "excitation_gap_map_status": gap_map.get("proof_status"),
        "source_package_value_law_artifact": gap_map.get("source_side_ordered_package_value_law_artifact"),
        "source_package_value_law_status": gap_map.get("source_side_ordered_package_value_law_status"),
        "hierarchy_mode": "two_scalar_ordered_gap_pair_candidate",
        "ordered_family_coordinate": [-1.0, x2, 1.0],
        "rho_ord": rho_ord,
        "current_family_eigenvalues_sorted": ordered_eigenvalues,
        "current_mean_eigenvalue": mean_value,
        "sigma_e_total_log_per_side": active_sigma,
        "eta_e_split_log_per_side": (float(eta_emitted) if gap_map_consumed else None),
        "eta_e_rigid_fallback": eta_e_rigid_fallback,
        "gap_pair_e": {
            "gamma21_log_per_side": gamma21,
            "gamma32_log_per_side": gamma32,
        },
        "E_e_log_centered": e_e_log_centered,
        "shape_log_shift_e": 0.0 if gap_map_consumed else None,
        "gap_pair_formula_two_scalar": {
            "gamma21_log_per_side": "((1 + x2) * sigma_e_total_log_per_side - eta_e_split_log_per_side) / 2",
            "gamma32_log_per_side": "((1 - x2) * sigma_e_total_log_per_side + eta_e_split_log_per_side) / 2",
        },
        "gap_pair_formula_rigid_fallback": {
            "eta_e_rigid_fallback": "(((1 + x2) - rho_ord * (1 - x2)) / (1 + rho_ord)) * sigma_e_total_log_per_side",
            "gamma21_log_per_side": "((1 + x2) * sigma_e_total_log_per_side - eta_e_rigid_fallback) / 2",
            "gamma32_log_per_side": "((1 - x2) * sigma_e_total_log_per_side + eta_e_rigid_fallback) / 2",
        },
        "centered_shape_formula": "E_e_log_centered = [-(2*gamma21+gamma32)/3, (gamma21-gamma32)/3, (gamma21+2*gamma32)/3]",
        "quadratic_basis_vector": q_ord,
        "quadratic_coeff_recovery": {
            "a_e": "(gamma21_log_per_side + gamma32_log_per_side) / 2",
            "b_e": "((1 + x2) * gamma32_log_per_side - (1 - x2) * gamma21_log_per_side) / (2 * (1 - x2^2))",
        },
        "diagonal_emitter_rule": "Y_e_shape = diag(exp(shape_log_shift_e + E_e_log_centered))",
        "shape_singular_values": shape_singular_values,
        "Y_e_shape": _diag_payload(shape_singular_values),
        "shape_closed": gap_map_consumed,
        "expected_closure_state": "absolute_scale_pending" if gap_map_consumed else "hierarchy_split_missing",
        "smallest_constructive_missing_object": smallest_missing_object,
        "metadata": {
            "note": "Current-family charged-lepton hierarchy on the ordered three-point family. The upstream ordered package is explicit; the remaining charged gap sits in the source-side value law on top of that package.",
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build a current-family ordered-log charged-lepton spectrum candidate.")
    parser.add_argument("--mean-certificate", default=str(DEFAULT_MEAN_CERT))
    parser.add_argument("--family", default=str(DEFAULT_FAMILY))
    parser.add_argument("--source-emission", default=str(DEFAULT_SOURCE_EMISSION))
    parser.add_argument("--gap-map", default=str(DEFAULT_GAP_MAP))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    mean_certificate = json.loads(Path(args.mean_certificate).read_text(encoding="utf-8"))
    family = json.loads(Path(args.family).read_text(encoding="utf-8"))
    source_emission_path = Path(args.source_emission)
    source_emission = json.loads(source_emission_path.read_text(encoding="utf-8")) if source_emission_path.exists() else None
    gap_map_path = Path(args.gap_map)
    gap_map = json.loads(gap_map_path.read_text(encoding="utf-8")) if gap_map_path.exists() else None
    artifact = build_artifact(mean_certificate, family, source_emission, gap_map)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
