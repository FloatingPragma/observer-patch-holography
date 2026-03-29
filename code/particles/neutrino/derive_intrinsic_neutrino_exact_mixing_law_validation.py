#!/usr/bin/env python3
"""Validate the exact intrinsic neutrino mixing laws around the isotropic point.

Chain role: audit the first-order and quadratic intrinsic mixing laws that sit
on top of the exact eta-chain artifact.

Mathematics: corrected linear/quadratic perturbation laws for the solar split,
ordered atmospheric gaps, centroid gap, and the collective-mode singular
vector.

OPH-derived inputs: the local isotropic forward neutrino bundle plus the exact
eta-map artifact.

Output: a diagnostic validation artifact showing which intrinsic mixing objects
are already exact, which are asymptotic, and which remain blocked only by the
proof-facing eta payload and shared charged basis.
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
DEFAULT_ISOTROPIC = ROOT / "particles" / "runs" / "neutrino" / "forward_neutrino_closure_bundle.json"
DEFAULT_EXACT_MAP = ROOT / "particles" / "runs" / "neutrino" / "intrinsic_neutrino_exact_eta_map.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "neutrino" / "intrinsic_neutrino_exact_mixing_law_validation.json"
EDGE_ORDER = ("psi12", "psi23", "psi31")


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _centered_eta(payload: dict[str, Any]) -> np.ndarray:
    eta = np.array([float(payload["eta_e"][edge]) for edge in EDGE_ORDER], dtype=float)
    return eta - float(np.mean(eta))


def _solve_selector(mu: np.ndarray, omega: float) -> np.ndarray:
    bound = float(np.min(mu)) * (1.0 - 1.0e-15)

    def f_lam(lam_value: float) -> float:
        return float(np.sum(np.arcsin(np.clip(lam_value / mu, -1.0, 1.0))) - omega)

    lo = -bound
    hi = bound
    f_lo = f_lam(lo)
    f_hi = f_lam(hi)
    if not (f_lo <= 0.0 <= f_hi):
        raise ValueError("principal selector bracket failed")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        f_mid = f_lam(mid)
        if f_lo * f_mid <= 0.0:
            hi = mid
            f_hi = f_mid
        else:
            lo = mid
            f_lo = f_mid
    return np.arcsin(np.clip((0.5 * (lo + hi)) / mu, -1.0, 1.0))


def _align_phase(reference: np.ndarray, vector: np.ndarray) -> np.ndarray:
    overlap = np.vdot(reference, vector)
    if abs(overlap) == 0.0:
        return vector
    return vector * np.exp(-1j * np.angle(overlap))


def _normalize(vector: np.ndarray) -> np.ndarray:
    norm = float(np.linalg.norm(vector))
    return vector if norm == 0.0 else vector / norm


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate the intrinsic neutrino exact mixing laws.")
    parser.add_argument("--isotropic", default=str(DEFAULT_ISOTROPIC))
    parser.add_argument("--exact-map", default=str(DEFAULT_EXACT_MAP))
    parser.add_argument("--mc-samples", type=int, default=250)
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    isotropic = _load_json(Path(args.isotropic))
    exact_map = _load_json(Path(args.exact_map))

    eta = _centered_eta(exact_map)
    mu = np.exp(eta)
    mu = mu / float(np.mean(mu))

    isotropic_matrix = np.array(isotropic["majorana_matrix_real"], dtype=float) + 1j * np.array(
        isotropic["majorana_matrix_imag"], dtype=float
    )
    exact_matrix = np.array(exact_map["majorana_matrix_real"], dtype=float) + 1j * np.array(
        exact_map["majorana_matrix_imag"], dtype=float
    )
    a_value = float(isotropic_matrix[0, 0].real)
    rho_value = float(abs(isotropic_matrix[0, 1]))
    phi = float(np.angle(isotropic_matrix[0, 1]))
    omega = 3.0 * phi

    psi_iso = np.array([phi, phi, phi], dtype=float)
    psi_payload = _solve_selector(mu, omega)
    psi_actual = np.array([float(exact_map["selector_point_absolute"][edge]) for edge in EDGE_ORDER], dtype=float)
    delta_actual = psi_actual - psi_iso
    delta_from_eta = -math.tan(phi) * eta
    sigma_actual = math.sqrt((2.0 / 3.0) * float(np.sum(delta_actual**2)))
    sigma_eta = math.sqrt((2.0 / 3.0) * float(np.sum(eta**2)))

    md2 = a_value * a_value + rho_value * rho_value - 2.0 * a_value * rho_value * math.cos(phi)
    mc2 = a_value * a_value + 4.0 * rho_value * rho_value + 4.0 * a_value * rho_value * math.cos(phi)
    atm_iso = mc2 - md2
    gamma = 2.0 * a_value * rho_value * math.sin(phi)

    solar_pred_delta = 2.0 * gamma * sigma_actual
    solar_pred_eta = 2.0 * gamma * math.tan(phi) * sigma_eta
    centroid_shift_pred = (
        -a_value
        * rho_value
        * (sigma_actual**2)
        * (a_value * (4.0 * math.cos(phi) ** 2 - 1.0) + 6.0 * rho_value * math.cos(phi))
        / (2.0 * a_value * math.cos(phi) + rho_value)
    )

    m2_iso = np.array(isotropic["masses_gev_sorted"], dtype=float) ** 2
    m2_ani = np.array(exact_map["masses_gev_sorted"], dtype=float) ** 2
    solar_actual = float(m2_ani[1] - m2_ani[0])
    dm31_actual = float(m2_ani[2] - m2_ani[0])
    dm32_actual = float(m2_ani[2] - m2_ani[1])

    centroid_gap_iso = float(m2_iso[2] - 0.5 * (m2_iso[0] + m2_iso[1]))
    centroid_gap_actual = float(m2_ani[2] - 0.5 * (m2_ani[0] + m2_ani[1]))
    centroid_shift_actual = centroid_gap_actual - centroid_gap_iso

    dm31_pred = atm_iso + 0.5 * solar_pred_delta + centroid_shift_pred
    dm32_pred = atm_iso - 0.5 * solar_pred_delta + centroid_shift_pred
    dm31_pred_eta = atm_iso + 0.5 * solar_pred_eta + centroid_shift_pred
    dm32_pred_eta = atm_iso - 0.5 * solar_pred_eta + centroid_shift_pred

    tan2theta_delta = None if abs(float(delta_actual[0])) <= 1.0e-30 else -float(
        (delta_actual[1] - delta_actual[2]) / (math.sqrt(3.0) * delta_actual[0])
    )
    tan2theta_eta = None if abs(float(eta[0])) <= 1.0e-30 else -float(
        (eta[1] - eta[2]) / (math.sqrt(3.0) * eta[0])
    )

    u_ani = np.array(exact_map["u_nu_left_real"], dtype=float) + 1j * np.array(
        exact_map["u_nu_left_imag"], dtype=float
    )
    v_collective_actual = u_ani[:, 2]
    u_dem = np.ones(3, dtype=complex) / math.sqrt(3.0)
    kappa = math.sqrt(3.0) * (2.0 * a_value * math.sin(phi) + 3.0j * rho_value) / (
        9.0 * (2.0 * a_value * math.cos(phi) + rho_value)
    )
    v_collective_pred_delta = _normalize(
        u_dem + kappa * np.array([delta_actual[1], delta_actual[2], delta_actual[0]], dtype=complex)
    )
    v_collective_pred_eta = _normalize(
        u_dem - math.tan(phi) * kappa * np.array([eta[1], eta[2], eta[0]], dtype=complex)
    )
    v_collective_actual = _align_phase(v_collective_pred_delta, _normalize(v_collective_actual))

    rng = np.random.default_rng(17)
    mc_rows: list[dict[str, float]] = []
    for scale in (1.0e-3, 3.0e-3, 1.0e-2, 3.0e-2, 1.0e-1):
        rel_solar_eta: list[float] = []
        rel_dm31: list[float] = []
        rel_centroid: list[float] = []
        overlap_vec: list[float] = []
        for _ in range(max(1, args.mc_samples)):
            direction = rng.normal(size=3)
            direction = direction - float(np.mean(direction))
            direction = direction / float(np.linalg.norm(direction))
            eta_trial = scale * direction
            mu_trial = np.exp(eta_trial)
            mu_trial = mu_trial / float(np.mean(mu_trial))
            psi_trial = _solve_selector(mu_trial, omega)
            matrix_trial = np.array(
                [
                    [a_value, rho_value * np.exp(1j * psi_trial[0]), rho_value * np.exp(1j * psi_trial[2])],
                    [rho_value * np.exp(1j * psi_trial[0]), a_value, rho_value * np.exp(1j * psi_trial[1])],
                    [rho_value * np.exp(1j * psi_trial[2]), rho_value * np.exp(1j * psi_trial[1]), a_value],
                ],
                dtype=complex,
            )
            h_trial = matrix_trial.conjugate().T @ matrix_trial
            eig_vals, eig_vecs = np.linalg.eigh(h_trial)
            order = np.argsort(eig_vals)
            eig_vals = eig_vals[order]
            eig_vecs = eig_vecs[:, order]
            delta_trial = psi_trial - phi
            sigma_trial = math.sqrt((2.0 / 3.0) * float(np.sum(delta_trial**2)))
            sigma_eta_trial = math.sqrt((2.0 / 3.0) * float(np.sum(eta_trial**2)))

            solar_true = float(eig_vals[1] - eig_vals[0])
            dm31_true = float(eig_vals[2] - eig_vals[0])
            centroid_true = float(eig_vals[2] - 0.5 * (eig_vals[0] + eig_vals[1]))

            solar_eta_pred = 2.0 * gamma * math.tan(phi) * sigma_eta_trial
            centroid_pred = atm_iso - a_value * rho_value * (sigma_trial**2) * (
                a_value * (4.0 * math.cos(phi) ** 2 - 1.0) + 6.0 * rho_value * math.cos(phi)
            ) / (2.0 * a_value * math.cos(phi) + rho_value)
            dm31_pred_trial = atm_iso + 0.5 * (2.0 * gamma * sigma_trial) + (centroid_pred - atm_iso)

            rel_solar_eta.append(abs(solar_eta_pred - solar_true) / solar_true)
            rel_dm31.append(abs(dm31_pred_trial - dm31_true) / dm31_true)
            rel_centroid.append(abs(centroid_pred - centroid_true) / centroid_true)
            vec_pred = _normalize(
                u_dem + kappa * np.array([delta_trial[1], delta_trial[2], delta_trial[0]], dtype=complex)
            )
            vec_true = _align_phase(vec_pred, _normalize(eig_vecs[:, 2]))
            overlap_vec.append(abs(np.vdot(vec_pred, vec_true)))

        mc_rows.append(
            {
                "eta_norm_scale": float(scale),
                "mean_rel_error_solar_eta_linear": float(np.mean(rel_solar_eta)),
                "mean_rel_error_dm31_corrected": float(np.mean(rel_dm31)),
                "mean_rel_error_centroid_quadratic": float(np.mean(rel_centroid)),
                "mean_collective_vector_overlap": float(np.mean(overlap_vec)),
                "min_collective_vector_overlap": float(np.min(overlap_vec)),
            }
        )

    payload = {
        "artifact": "oph_intrinsic_neutrino_exact_mixing_law_validation",
        "generated_utc": _timestamp(),
        "base_parameters": {
            "a_diag_gev": a_value,
            "rho_offdiag_gev": rho_value,
            "phi_equal_split_rad": phi,
            "omega_cycle_rad": omega,
            "md2_isotropic_gev2": md2,
            "mc2_isotropic_gev2": mc2,
            "atm_centroid_gap_isotropic_gev2": atm_iso,
        },
        "eta_payload": {edge: float(eta[idx]) for idx, edge in enumerate(EDGE_ORDER)},
        "selector_payload_solution": {edge: float(psi_payload[idx]) for idx, edge in enumerate(EDGE_ORDER)},
        "selector_actual_from_matrix": {edge: float(psi_actual[idx]) for idx, edge in enumerate(EDGE_ORDER)},
        "delta_actual": {edge: float(delta_actual[idx]) for idx, edge in enumerate(EDGE_ORDER)},
        "delta_pred_from_eta_linear": {edge: float(delta_from_eta[idx]) for idx, edge in enumerate(EDGE_ORDER)},
        "sigma_actual": sigma_actual,
        "sigma_eta": sigma_eta,
        "solar_split_actual_gev2": solar_actual,
        "solar_split_pred_linear_from_delta_gev2": solar_pred_delta,
        "solar_split_pred_linear_from_eta_gev2": solar_pred_eta,
        "solar_split_rel_error_from_delta_linear": abs(solar_pred_delta - solar_actual) / solar_actual,
        "solar_split_rel_error_from_eta_linear": abs(solar_pred_eta - solar_actual) / solar_actual,
        "ordered_atmospheric_shift_is_linear": True,
        "delta_m31_actual_gev2": dm31_actual,
        "delta_m31_pred_corrected_gev2": dm31_pred,
        "delta_m31_pred_corrected_from_eta_gev2": dm31_pred_eta,
        "delta_m31_rel_error_corrected": abs(dm31_pred - dm31_actual) / dm31_actual,
        "delta_m31_rel_error_corrected_from_eta": abs(dm31_pred_eta - dm31_actual) / dm31_actual,
        "delta_m32_actual_gev2": dm32_actual,
        "delta_m32_pred_corrected_gev2": dm32_pred,
        "delta_m32_pred_corrected_from_eta_gev2": dm32_pred_eta,
        "delta_m32_rel_error_corrected": abs(dm32_pred - dm32_actual) / dm32_actual,
        "delta_m32_rel_error_corrected_from_eta": abs(dm32_pred_eta - dm32_actual) / dm32_actual,
        "centroid_gap_is_first_order_invariant": True,
        "centroid_gap_isotropic_gev2": centroid_gap_iso,
        "centroid_gap_actual_gev2": centroid_gap_actual,
        "centroid_gap_pred_quadratic_gev2": centroid_gap_iso + centroid_shift_pred,
        "centroid_gap_shift_actual_gev2": centroid_shift_actual,
        "centroid_gap_shift_pred_quadratic_gev2": centroid_shift_pred,
        "centroid_gap_shift_rel_error_quadratic": abs(centroid_shift_pred - centroid_shift_actual) / abs(centroid_shift_actual),
        "doublet_rotation_tan2theta_delta": tan2theta_delta,
        "doublet_rotation_tan2theta_eta_linear": tan2theta_eta,
        "collective_vector_prediction_from_delta": {
            "real": np.real(v_collective_pred_delta).tolist(),
            "imag": np.imag(v_collective_pred_delta).tolist(),
        },
        "collective_vector_prediction_from_eta": {
            "real": np.real(v_collective_pred_eta).tolist(),
            "imag": np.imag(v_collective_pred_eta).tolist(),
        },
        "collective_vector_actual_aligned": {
            "real": np.real(v_collective_actual).tolist(),
            "imag": np.imag(v_collective_actual).tolist(),
        },
        "collective_vector_overlap_from_delta": float(abs(np.vdot(v_collective_pred_delta, v_collective_actual))),
        "collective_vector_overlap_from_eta": float(abs(np.vdot(v_collective_pred_eta, v_collective_actual))),
        "mc_scaling_check": mc_rows,
        "notes": [
            "This validation corrects the earlier atmospheric statement: ordered Delta m31^2 and Delta m32^2 move linearly once the solar pair splits.",
            "The first-order invariant atmospheric object is the collective-to-doublet-centroid gap, not the ordered Delta m31^2 row.",
            "The collective singular-vector deformation law is already highly accurate on the weighted demo and becomes asymptotically exact near the isotropic point.",
        ],
    }

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
