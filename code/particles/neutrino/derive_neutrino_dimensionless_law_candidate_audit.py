#!/usr/bin/env python3
"""Audit simple neutrino weighted-cycle law candidates around the live selector."""

from __future__ import annotations

import argparse
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from particles.artifact_paths import portable_json_dumps
DEFAULT_CERTIFICATE = ROOT / "particles" / "runs" / "neutrino" / "same_label_scalar_certificate.json"
DEFAULT_COCYCLE = ROOT / "particles" / "runs" / "flavor" / "overlap_edge_transport_cocycle.json"
DEFAULT_PHASE_SOURCE = ROOT / "particles" / "runs" / "neutrino" / "intrinsic_neutrino_mass_eigenstate_bundle_from_scalar_certificate.json"
DEFAULT_SELECTOR = ROOT / "particles" / "runs" / "neutrino" / "neutrino_transport_load_segment_selector.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "neutrino" / "neutrino_dimensionless_law_candidate_audit.json"
EDGE_ORDER = ("psi12", "psi23", "psi31")


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _pmns_parameters(unitary: np.ndarray) -> dict[str, float]:
    s13 = abs(unitary[0, 2])
    theta13 = math.asin(np.clip(s13, 0.0, 1.0))
    c13 = math.cos(theta13)

    s12 = abs(unitary[0, 1]) / max(c13, 1.0e-30)
    s12 = float(np.clip(s12, 0.0, 1.0))
    theta12 = math.asin(s12)

    s23 = abs(unitary[1, 2]) / max(c13, 1.0e-30)
    s23 = float(np.clip(s23, 0.0, 1.0))
    theta23 = math.asin(s23)

    jarlskog = float(np.imag(unitary[0, 0] * unitary[1, 1] * np.conjugate(unitary[0, 1]) * np.conjugate(unitary[1, 0])))

    c12 = math.cos(theta12)
    c23 = math.cos(theta23)
    denom = 2.0 * s12 * c12 * s23 * c23 * s13
    if abs(denom) <= 1.0e-30:
        delta = 0.0
    else:
        cos_delta = (
            (s12 * s23) ** 2 + (c12 * c23 * s13) ** 2 - abs(unitary[2, 0]) ** 2
        ) / denom
        cos_delta = float(np.clip(cos_delta, -1.0, 1.0))
        den_j = c12 * s12 * c23 * s23 * (c13**2) * s13
        sin_delta = 0.0 if abs(den_j) <= 1.0e-30 else float(np.clip(jarlskog / den_j, -1.0, 1.0))
        delta = math.atan2(sin_delta, cos_delta) % (2.0 * math.pi)

    return {
        "theta12_deg": math.degrees(theta12),
        "theta23_deg": math.degrees(theta23),
        "theta13_deg": math.degrees(theta13),
        "delta_deg": math.degrees(delta),
        "J": jarlskog,
    }


def _candidate_surface(q: dict[str, float], psi: dict[str, float], p: float, chi: float) -> dict[str, Any]:
    weights = {edge: float(q[edge] ** p) for edge in EDGE_ORDER}
    phases = {edge: -psi[edge] for edge in EDGE_ORDER}
    cycle_matrix = np.array(
        [
            [
                -chi * weights["psi31"],
                weights["psi31"] * np.exp(1j * phases["psi31"]),
                weights["psi23"] * np.exp(1j * phases["psi23"]),
            ],
            [
                weights["psi31"] * np.exp(1j * phases["psi31"]),
                -chi * weights["psi12"],
                weights["psi12"] * np.exp(1j * phases["psi12"]),
            ],
            [
                weights["psi23"] * np.exp(1j * phases["psi23"]),
                weights["psi12"] * np.exp(1j * phases["psi12"]),
                -chi * weights["psi23"],
            ],
        ],
        dtype=complex,
    )
    hermitian = cycle_matrix.conjugate().T @ cycle_matrix
    evals, unitary = np.linalg.eigh(hermitian)
    evals = np.asarray(np.real_if_close(evals), dtype=float)
    dm21 = float(evals[1] - evals[0])
    dm31 = float(evals[2] - evals[0])
    dm32 = float(evals[2] - evals[1])
    return {
        "p": float(p),
        "chi": float(chi),
        "ratio_dm21_over_dm32": float(dm21 / dm32),
        "dimensionless_dm2": {"21": dm21, "31": dm31, "32": dm32},
        "pmns": _pmns_parameters(unitary),
    }


def _segment_candidate(
    q: dict[str, float],
    psi: dict[str, float],
    *,
    gamma: float,
    eps: float,
    chi: float,
    denominator: float,
) -> dict[str, Any]:
    payload = _candidate_surface(q, psi, 1.0 + gamma + eps / denominator, chi)
    payload["segment_denominator"] = float(denominator)
    return payload


def main() -> int:
    parser = argparse.ArgumentParser(description="Audit compare-only neutrino weighted-cycle law candidates.")
    parser.add_argument("--certificate", default=str(DEFAULT_CERTIFICATE))
    parser.add_argument("--cocycle", default=str(DEFAULT_COCYCLE))
    parser.add_argument("--phase-source", default=str(DEFAULT_PHASE_SOURCE))
    parser.add_argument("--selector", default=str(DEFAULT_SELECTOR))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    certificate = _load_json(Path(args.certificate))
    cocycle = _load_json(Path(args.cocycle))
    phase_source = _load_json(Path(args.phase_source))
    selector = _load_json(Path(args.selector))

    q = {edge: float(certificate["q_e"][edge]) for edge in EDGE_ORDER}
    psi = {edge: float(phase_source["selector_point_absolute"][edge]) for edge in EDGE_ORDER}
    gamma = float(cocycle["theorem_gap_gamma"])
    eps = float(cocycle["defect_gap_ratio"])
    gamma_half = float(cocycle["hermitian_descendant_riesz_margin"]["gamma_half"])
    chi = 1.0 + eps
    representative_ratio = 7.49e-5 / 2.438e-3

    segment_arithmetic = (chi + 1.0 + gamma_half) / 2.0
    segment_geometric = math.sqrt(chi * (1.0 + gamma_half))
    segment_harmonic = 2.0 / ((1.0 / chi) + (1.0 / (1.0 + gamma_half)))

    candidates = {
        "current_law": _candidate_surface(q, psi, 1.0 + gamma + eps, chi),
        "normalized_eps_over_chi": _candidate_surface(q, psi, 1.0 + gamma + eps / chi, chi),
        "midpoint_normalized_gap_defect": _segment_candidate(
            q,
            psi,
            gamma=gamma,
            eps=eps,
            chi=chi,
            denominator=segment_arithmetic,
        ),
        "segment_geometric_mean": _segment_candidate(
            q,
            psi,
            gamma=gamma,
            eps=eps,
            chi=chi,
            denominator=segment_geometric,
        ),
        "segment_harmonic_mean": _segment_candidate(
            q,
            psi,
            gamma=gamma,
            eps=eps,
            chi=chi,
            denominator=segment_harmonic,
        ),
    }
    for payload in candidates.values():
        payload["representative_ratio_error"] = float(payload["ratio_dm21_over_dm32"] - representative_ratio)
        payload["absolute_ratio_error"] = float(abs(payload["representative_ratio_error"]))

    ranking = sorted(
        (
            {
                "candidate": name,
                "absolute_ratio_error": payload["absolute_ratio_error"],
                "theta12_deg": payload["pmns"]["theta12_deg"],
                "theta23_deg": payload["pmns"]["theta23_deg"],
                "theta13_deg": payload["pmns"]["theta13_deg"],
            }
            for name, payload in candidates.items()
        ),
        key=lambda item: item["absolute_ratio_error"],
    )

    artifact = {
        "artifact": "oph_neutrino_dimensionless_law_candidate_audit",
        "generated_utc": _timestamp(),
        "status": "law_space_audit_around_live_selector",
        "public_promotion_allowed": False,
        "source_artifacts": {
            "same_label_scalar_certificate": str(Path(args.certificate)),
            "overlap_edge_transport_cocycle": str(Path(args.cocycle)),
            "selector_phase_source": str(Path(args.phase_source)),
            "transport_load_selector": str(Path(args.selector)),
        },
        "live_parameters": {
            "gamma": gamma,
            "gamma_half": gamma_half,
            "defect_gap_ratio": eps,
            "chi": chi,
        },
        "midpoint_equivalent_law": {
            "denominator_formula": "D_mid = 1 + eps/2 + gamma/4 = 1 + eps/2 + gamma_half/2 = (chi + 1 + gamma_half)/2",
            "selector_family": "D_tau = (1-tau_nu) * chi + tau_nu * (1 + gamma_half)",
            "midpoint_choice": "tau_nu = 1/2",
            "theorem_status": "live_standard_math_fixed_selector",
            "exact_missing_object": None,
        },
        "live_promoted_selector": {
            "selector_id": selector.get("selected_selector"),
            "D_nu": selector.get("selected_D_nu"),
            "weight_exponent_value": selector.get("derived_quantities", {}).get("weight_exponent_value"),
            "closure_scope": selector.get("closure_scope"),
        },
        "segment_selector_candidates": {
            "left_endpoint_chi": chi,
            "right_endpoint_one_plus_gamma_half": 1.0 + gamma_half,
            "arithmetic_mean": segment_arithmetic,
            "geometric_mean": segment_geometric,
            "harmonic_mean": segment_harmonic,
            "note": "These are selector candidates on the positive load segment between chi and 1 + gamma_half. The arithmetic midpoint is the live standard-math-fixed selector because balanced and least-distortion coincide on a one-dimensional affine segment. The geometric and harmonic means remain compare-only alternatives, and all three mean selectors outperform the legacy raw-defect law by a wide margin.",
        },
        "representative_pdg_central_ratio": {
            "delta_m21_sq_eV2": 7.49e-5,
            "delta_m32_sq_eV2": 2.438e-3,
            "ratio": representative_ratio,
            "note": "Representative PDG 2025 central values used only for compare-only candidate-law auditing.",
        },
        "candidates": candidates,
        "ranking_by_absolute_ratio_error": ranking,
        "leading_observation": {
            "candidate": ranking[0]["candidate"] if ranking else None,
            "note": "The full positive-load segment family between chi and 1 + gamma_half is explicit in the audit. The arithmetic midpoint is the live standard-math-fixed selector on the repaired branch, while the geometric mean remains the numerical leader by raw central-ratio error on this compare-only audit surface.",
        },
    }

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(portable_json_dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
