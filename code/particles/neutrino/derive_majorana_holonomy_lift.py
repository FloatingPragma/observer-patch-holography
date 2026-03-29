#!/usr/bin/env python3
"""Build the Majorana holonomy-lift boundary on the residual phase cycle.

Chain role: expose the phase-cycle constraint, residual basis, and canonical
selector candidates that sit between the real family tensor and the phase law.

Mathematics: weighted cycle-constrained least-distortion selectors on the
three-edge Majorana phase graph.

OPH-derived inputs: the shared sector transport pushforward that fixes the
cycle weights and residual basis data used in the neutrino lane.

Output: the Majorana lift artifact consumed by the pullback metric, phase
envelope, and forward matrix builders.
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib
from datetime import datetime, timezone

import numpy as np


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "flavor" / "sector_transport_pushforward.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "neutrino" / "majorana_holonomy_lift.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _selector_energy(psi: dict[str, float], weights: dict[str, float]) -> float:
    return float(
        sum(weights[key] * (1.0 - math.cos(float(psi[key]))) for key in ("psi12", "psi23", "psi31"))
    )


def _clip_ratio(value: float) -> float:
    return max(-1.0, min(1.0, value))


def _harmonic_selector(omega_012: float, weights: dict[str, float]) -> dict[str, float]:
    inverse_weights = {key: 1.0 / max(value, 1.0e-30) for key, value in weights.items()}
    normalizer = sum(inverse_weights.values())
    return {key: omega_012 * value / normalizer for key, value in inverse_weights.items()}


def _least_distortion_selector(omega_012: float, weights: dict[str, float]) -> dict[str, float]:
    edge_order = ("psi12", "psi23", "psi31")
    weight_values = [max(float(weights[key]), 1.0e-30) for key in edge_order]
    min_weight = min(weight_values)
    epsilon = min_weight * 1.0e-12

    def phase_sum(lambda_value: float) -> float:
        return sum(math.asin(_clip_ratio(lambda_value / weight)) for weight in weight_values)

    lo = -min_weight + epsilon
    hi = min_weight - epsilon
    flo = phase_sum(lo) - omega_012
    fhi = phase_sum(hi) - omega_012
    if flo * fhi > 0.0:
        return _harmonic_selector(omega_012, weights)

    for _ in range(90):
        mid = 0.5 * (lo + hi)
        fmid = phase_sum(mid) - omega_012
        if flo * fmid <= 0.0:
            hi = mid
            fhi = fmid
        else:
            lo = mid
            flo = fmid

    lambda_star = 0.5 * (lo + hi)
    return {
        key: math.asin(_clip_ratio(lambda_star / max(float(weights[key]), 1.0e-30)))
        for key in edge_order
    }


def _isotropic_certificate(values: dict[str, float], tolerance: float = 1.0e-30) -> dict[str, object]:
    ordered = [float(values[key]) for key in sorted(values)]
    spread = max(ordered) - min(ordered)
    return {
        "closed": bool(spread <= tolerance),
        "spread": spread,
        "tolerance": tolerance,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the neutrino Majorana holonomy-lift artifact.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT), help="Input sector-response JSON path.")
    parser.add_argument("--output", default=str(DEFAULT_OUT), help="Output JSON path.")
    args = parser.parse_args()

    payload = json.loads(pathlib.Path(args.input).read_text(encoding="utf-8"))
    nu = dict(payload.get("sector_response_object", {}).get("nu", {}))
    if not nu:
        raise ValueError("sector_response_object['nu'] is required")

    omega_012 = float(dict(nu.get("omega_cycle", {})).get("012", 0.0))
    k_sym = np.asarray(nu.get("K_core_majorana_sym", nu.get("K_core")), dtype=float)
    if k_sym.shape != (3, 3):
        raise ValueError("neutrino symmetric kernel must be 3x3")
    normalization_diag = np.sqrt(np.clip(np.diag(k_sym), 1.0e-15, None))
    inv_diag = np.diag([1.0 / value for value in normalization_diag])
    e_nu = inv_diag @ k_sym @ inv_diag

    edge_weights = {
        "psi12": float((normalization_diag[0] * normalization_diag[1] * e_nu[0, 1]) ** 2),
        "psi23": float((normalization_diag[1] * normalization_diag[2] * e_nu[1, 2]) ** 2),
        "psi31": float((normalization_diag[2] * normalization_diag[0] * e_nu[2, 0]) ** 2),
    }
    balanced = {"psi12": omega_012 / 3.0, "psi23": omega_012 / 3.0, "psi31": omega_012 / 3.0}
    harmonic = _harmonic_selector(omega_012, edge_weights)
    least_distortion = _least_distortion_selector(omega_012, edge_weights)
    selector_candidates = {
        "balanced": {
            **balanced,
            "selector": "balanced",
            "energy_majorana": _selector_energy(balanced, edge_weights),
        },
        "harmonic": {
            **harmonic,
            "selector": "harmonic",
            "energy_majorana": _selector_energy(harmonic, edge_weights),
        },
        "least_distortion": {
            **least_distortion,
            "selector": "least_distortion",
            "energy_majorana": _selector_energy(least_distortion, edge_weights),
        },
    }
    weight_isotropy = _isotropic_certificate(edge_weights)
    principal_equal_split = {
        "psi12": balanced["psi12"],
        "psi23": balanced["psi23"],
        "psi31": balanced["psi31"],
        "selector": "principal_equal_split",
        "status": "closed_equal_split" if weight_isotropy["closed"] else "candidate_only",
    }
    canonical_lift_closed = bool(weight_isotropy["closed"])
    canonical_selector_status = "closed_equal_split" if canonical_lift_closed else "candidate_only"

    artifact = {
        "artifact": "oph_majorana_holonomy_lift",
        "generated_utc": _timestamp(),
        "labels": payload.get("labels"),
        "symmetry_type": "majorana",
        "phase_space": ["psi12", "psi23", "psi31"],
        "phase_status": "residual_affine_plane",
        "canonical_lift_closed": canonical_lift_closed,
        "canonical_selector_status": canonical_selector_status,
        "selector_law_status": "candidate_only",
        "selector_equivalence_class": "principal_equal_split" if canonical_lift_closed else "unresolved",
        "selector_closure_reason": "s3_fixed_point" if canonical_lift_closed else "candidate_only",
        "principal_branch_certificate": {
            "closed": True,
            "branch": "principal",
            "range": "(-pi, pi]",
        },
        "selector_family": ["balanced", "harmonic", "least_distortion"],
        "affine_constraint_matrix": [[1.0, 1.0, 1.0]],
        "affine_constraint_rhs": [omega_012],
        "cycle_constraint": {
            "omega_012": omega_012,
            "equation": "psi12+psi23+psi31=Omega_012",
        },
        "edge_weights_majorana": edge_weights,
        "edge_weight_isotropy_certificate": weight_isotropy,
        "selector_candidates": selector_candidates,
        "canonical_selector_point": principal_equal_split,
        "selector_candidate_psi": {
            **selector_candidates["least_distortion"],
            "status": "candidate_only",
        },
        "balanced_lift_candidate": {
            "psi12": balanced["psi12"],
            "psi23": balanced["psi23"],
            "psi31": balanced["psi31"],
            "selector": "balanced_cycle_split_candidate_only",
        },
        "residual_basis": [
            {"psi12": 1.0, "psi23": -1.0, "psi31": 0.0},
            {"psi12": 1.0, "psi23": 0.0, "psi31": -1.0},
        ],
        "majorana_phase_obstruction_class": {
            "kind": "congruence_gauge_affine_plane",
            "residual_dimension": 2,
            "closed": False,
        },
        "notes": [
            "The neutrino Majorana phase lift is kept separate from the real symmetric amplitude branch.",
            "Balanced, harmonic, and least-distortion selectors are exported as selector-law candidates.",
            "On the current isotropic branch, the selector point closes as the principal equal split by S3 symmetry and the affine cycle constraint.",
            "The least-distortion selector law is not yet promoted beyond candidate status on this branch.",
        ],
    }

    out_path = pathlib.Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
