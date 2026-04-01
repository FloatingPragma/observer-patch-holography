#!/usr/bin/env python3
"""Emit a full exact neutrino compare-only continuation adapter.

Chain role: expose the strongest exact compare-only neutrino fit currently
available on disk by moving along the already-emitted positive selector segment
and then rescaling the resulting dimensionless branch by one positive
normalization.

Mathematics: the live midpoint selector fixes one theorem-grade weighted-cycle
branch and leaves the exact missing theorem object `C_nu` open. This adapter
does not change that theorem status. It uses the broader already-explicit
segment family

    D_tau = (1 - tau_nu) * chi + tau_nu * (1 + gamma_half)

to solve the representative PDG central ratio exactly, then uses one positive
`lambda_nu` to hit the atmospheric splitting exactly as well. The result is an
exact two-observable compare-only continuation adapter, not a promoted OPH
mass theorem.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CERTIFICATE = ROOT / "particles" / "runs" / "neutrino" / "same_label_scalar_certificate.json"
DEFAULT_COCYCLE = ROOT / "particles" / "runs" / "flavor" / "overlap_edge_transport_cocycle.json"
DEFAULT_PHASE_SOURCE = ROOT / "particles" / "runs" / "neutrino" / "intrinsic_neutrino_mass_eigenstate_bundle_from_scalar_certificate.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "neutrino" / "neutrino_two_parameter_exact_adapter.json"
EDGE_ORDER = ("psi12", "psi23", "psi31")

# PDG 2025 neutrino review Table 14.7, Ref. [193], normal ordering,
# SK-ATM and IC24 representative central values used only for compare-only fits.
PDG_2025_NO_CENTRAL = {
    "source": "PDG 2025 neutrino review Table 14.7, Ref. [193] with SK-ATM and IC24, normal ordering, representative central values",
    "delta_m21_sq_eV2": 7.49e-5,
    "delta_m32_sq_eV2": 2.438e-3,
}


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


def _surface_at_tau(
    *,
    tau_nu: float,
    q: dict[str, float],
    psi: dict[str, float],
    gamma: float,
    eps: float,
    gamma_half: float,
    chi: float,
) -> dict[str, Any]:
    d_nu = (1.0 - tau_nu) * chi + tau_nu * (1.0 + gamma_half)
    p_nu = 1.0 + gamma + eps / d_nu
    weights = {edge: float(q[edge] ** p_nu) for edge in EDGE_ORDER}
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
    masses_dimless = np.sqrt(np.maximum(evals, 0.0))
    return {
        "tau_nu": float(tau_nu),
        "D_nu": float(d_nu),
        "p_nu": float(p_nu),
        "m_hat": [float(x) for x in masses_dimless],
        "delta_hat_m_sq_eV2": {
            "21": dm21,
            "31": dm31,
            "32": dm32,
        },
        "ratio_21_over_32": float(dm21 / dm32),
        "pmns_observables": _pmns_parameters(unitary),
    }


def _solve_tau_exact_ratio(
    *,
    q: dict[str, float],
    psi: dict[str, float],
    gamma: float,
    eps: float,
    gamma_half: float,
    chi: float,
    target_ratio: float,
) -> dict[str, Any]:
    left = _surface_at_tau(tau_nu=0.0, q=q, psi=psi, gamma=gamma, eps=eps, gamma_half=gamma_half, chi=chi)
    right = _surface_at_tau(tau_nu=1.0, q=q, psi=psi, gamma=gamma, eps=eps, gamma_half=gamma_half, chi=chi)
    f_left = left["ratio_21_over_32"] - target_ratio
    f_right = right["ratio_21_over_32"] - target_ratio
    if abs(f_left) <= 1.0e-18:
        return left
    if abs(f_right) <= 1.0e-18:
        return right
    if f_left * f_right > 0.0:
        raise SystemExit("positive selector segment does not bracket the representative PDG central ratio")

    lo = 0.0
    hi = 1.0
    f_lo = f_left
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        mid_payload = _surface_at_tau(tau_nu=mid, q=q, psi=psi, gamma=gamma, eps=eps, gamma_half=gamma_half, chi=chi)
        f_mid = mid_payload["ratio_21_over_32"] - target_ratio
        if abs(f_mid) <= 1.0e-18 or (hi - lo) <= 1.0e-18:
            return mid_payload
        if f_lo * f_mid <= 0.0:
            hi = mid
        else:
            lo = mid
            f_lo = f_mid
    return _surface_at_tau(tau_nu=0.5 * (lo + hi), q=q, psi=psi, gamma=gamma, eps=eps, gamma_half=gamma_half, chi=chi)


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the compare-only exact neutrino two-parameter continuation adapter.")
    parser.add_argument("--certificate", default=str(DEFAULT_CERTIFICATE))
    parser.add_argument("--cocycle", default=str(DEFAULT_COCYCLE))
    parser.add_argument("--phase-source", default=str(DEFAULT_PHASE_SOURCE))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    certificate = _load_json(Path(args.certificate))
    cocycle = _load_json(Path(args.cocycle))
    phase_source = _load_json(Path(args.phase_source))

    q = {edge: float(certificate["q_e"][edge]) for edge in EDGE_ORDER}
    psi = {edge: float(phase_source["selector_point_absolute"][edge]) for edge in EDGE_ORDER}
    gamma = float(cocycle["theorem_gap_gamma"])
    eps = float(cocycle["defect_gap_ratio"])
    gamma_half = float(cocycle["hermitian_descendant_riesz_margin"]["gamma_half"])
    chi = 1.0 + eps
    target_21 = float(PDG_2025_NO_CENTRAL["delta_m21_sq_eV2"])
    target_32 = float(PDG_2025_NO_CENTRAL["delta_m32_sq_eV2"])
    target_31 = target_21 + target_32
    target_ratio = target_21 / target_32

    midpoint_surface = _surface_at_tau(tau_nu=0.5, q=q, psi=psi, gamma=gamma, eps=eps, gamma_half=gamma_half, chi=chi)
    exact_surface = _solve_tau_exact_ratio(
        q=q,
        psi=psi,
        gamma=gamma,
        eps=eps,
        gamma_half=gamma_half,
        chi=chi,
        target_ratio=target_ratio,
    )

    dh = exact_surface["delta_hat_m_sq_eV2"]
    z = target_32 / float(dh["32"])
    lambda_nu = math.sqrt(z)
    exact_masses = [lambda_nu * value for value in exact_surface["m_hat"]]
    exact_dm2 = {
        "21": float(dh["21"] * z),
        "31": float(dh["31"] * z),
        "32": float(dh["32"] * z),
    }

    payload = {
        "artifact": "oph_neutrino_two_parameter_exact_adapter",
        "generated_utc": _timestamp(),
        "proof_status": "compare_only_exact_two_parameter_continuation_adapter",
        "scope": "compare_only_two_parameter_segment_adapter",
        "promotable": False,
        "source_artifacts": {
            "same_label_scalar_certificate": str(Path(args.certificate)),
            "overlap_edge_transport_cocycle": str(Path(args.cocycle)),
            "selector_phase_source": str(Path(args.phase_source)),
        },
        "theorem_boundary": {
            "status": "non_promotable_compare_only_segment_and_scale_inverse_adapter",
            "statement": (
                "This adapter solves exact representative PDG central splittings only by moving along the already-explicit "
                "positive selector segment and then rescaling the resulting scale-free branch. The live theorem lane still "
                "waits on the reduced bridge-correction invariant C_nu and therefore does not emit a theorem-grade absolute neutrino spectrum."
            ),
            "forbidden_feedback": "compare_only_segment_adapter_must_not_feed_back_into_theorem_state_or_C_nu_emission",
        },
        "reference_central_values": {
            **PDG_2025_NO_CENTRAL,
            "delta_m31_sq_eV2": target_31,
            "ratio_21_over_32": target_ratio,
        },
        "selector_family": {
            "formula": "D_tau = (1 - tau_nu) * chi + tau_nu * (1 + gamma_half)",
            "segment_endpoints": {
                "chi": chi,
                "one_plus_gamma_half": 1.0 + gamma_half,
            },
            "midpoint_live_selector": {
                "tau_nu": 0.5,
                "D_nu": midpoint_surface["D_nu"],
                "p_nu": midpoint_surface["p_nu"],
                "ratio_21_over_32": midpoint_surface["ratio_21_over_32"],
            },
        },
        "exact_solution": {
            "tau_nu": exact_surface["tau_nu"],
            "D_nu": exact_surface["D_nu"],
            "p_nu": exact_surface["p_nu"],
            "lambda_nu": lambda_nu,
            "relative_tau_shift_from_live_midpoint": exact_surface["tau_nu"] - 0.5,
        },
        "exact_outputs": {
            "masses_eV": exact_masses,
            "delta_m_sq_eV2": exact_dm2,
            "ratio_21_over_32": exact_dm2["21"] / exact_dm2["32"],
            "pmns_observables": exact_surface["pmns_observables"],
        },
        "exact_fit_residuals_eV2": {
            "21": exact_dm2["21"] - target_21,
            "31": exact_dm2["31"] - target_31,
            "32": exact_dm2["32"] - target_32,
        },
        "notes": [
            "The exact fit is achieved on the existing positive selector segment; no new theorem object is claimed.",
            "The exact match uses two compare-only degrees of freedom: tau_nu fixes the dimensionless ratio and lambda_nu fixes the overall positive scale.",
            "This exact adapter is stronger than the older one-observable atmospheric-only and solar-only slices, but it remains non-promotable until the live theorem lane emits C_nu.",
        ],
    }

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
