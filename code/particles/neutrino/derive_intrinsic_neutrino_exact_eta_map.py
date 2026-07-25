#!/usr/bin/env python3
"""Build the exact intrinsic neutrino eta-chain artifact.

Chain role: turn a centered same-label eta-class into the exact intrinsic
Majorana matrix, ascending singular spectrum, squared gaps, and Takagi
subspaces. It does not select physical mass-eigenstate labels or ordering.

Mathematics: principal-branch selector from the scale-free eta-class, exact
complex symmetric Majorana matrix construction, depressed cubic for
`H = M^dagger M`, and exact singular-spectrum extraction.

Declared inputs: the local isotropic neutrino forward bundle plus a payload
carrying the centered same-label eta-class or any scale-equivalent positive
family. Source provenance must be checked separately before physical use.

Output: exact conditional spectral algebra once a centered same-label eta-class
is supplied.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from particles.artifact_paths import portable_json_dumps
DEFAULT_ISOTROPIC = ROOT / "particles" / "runs" / "neutrino" / "forward_majorana_matrix.json"
DEFAULT_PAYLOAD = ROOT / "particles" / "runs" / "neutrino" / "same_label_scalar_certificate.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "neutrino" / "intrinsic_neutrino_exact_eta_map.json"
EDGE_ORDER = ("psi12", "psi23", "psi31")


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _unwrap_near(value: float, reference: float) -> float:
    two_pi = 2.0 * math.pi
    k_shift = round((reference - value) / two_pi)
    return value + two_pi * k_shift


def _phase_vector_from_matrix(matrix: np.ndarray) -> np.ndarray:
    return np.array(
        [
            float(np.angle(matrix[0, 1])),
            float(np.angle(matrix[1, 2])),
            float(np.angle(matrix[0, 2])),
        ],
        dtype=float,
    )


def _centered_eta_from_payload(payload: dict[str, Any]) -> np.ndarray:
    if isinstance(payload.get("eta_e"), dict):
        eta = np.array([float(payload["eta_e"][edge]) for edge in EDGE_ORDER], dtype=float)
        return eta - float(np.mean(eta))
    if isinstance(payload.get("q_e"), dict):
        q_value = np.array([float(payload["q_e"][edge]) for edge in EDGE_ORDER], dtype=float)
        if np.any(q_value <= 0.0):
            raise ValueError("q_e values must be positive")
        eta = np.log(q_value)
        return eta - float(np.mean(eta))
    if isinstance(payload.get("mu_e"), dict):
        mu_value = np.array([float(payload["mu_e"][edge]) for edge in EDGE_ORDER], dtype=float)
        if np.any(mu_value <= 0.0):
            raise ValueError("mu_e values must be positive")
        eta = np.log(mu_value)
        return eta - float(np.mean(eta))
    raise ValueError("payload must provide eta_e, q_e, or mu_e")


def _normalized_mu_from_eta(eta: np.ndarray) -> np.ndarray:
    mu = np.exp(np.asarray(eta, dtype=float))
    return mu / float(np.mean(mu))


def _principal_selector_domain(mu: np.ndarray) -> float:
    mu = np.asarray(mu, dtype=float)
    if np.any(mu <= 0.0):
        raise ValueError("mu must be positive")
    mu_min = float(np.min(mu))
    return float(np.sum(np.arcsin(np.clip(mu_min / mu, -1.0, 1.0))))


def _solve_principal_selector(mu: np.ndarray, omega: float) -> tuple[float, np.ndarray]:
    mu = np.asarray(mu, dtype=float)
    if np.any(mu <= 0.0):
        raise ValueError("mu must be positive")
    domain_half_width = _principal_selector_domain(mu)
    if not abs(float(omega)) < domain_half_width:
        raise ValueError(
            "principal selector requires |Omega| < "
            f"sum_e arcsin(mu_min / mu_e) = {domain_half_width:.17g}"
        )
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
        if f_mid == 0.0:
            lo = hi = mid
            break
        if f_lo * f_mid <= 0.0:
            hi = mid
            f_hi = f_mid
        else:
            lo = mid
            f_lo = f_mid
    lam_value = 0.5 * (lo + hi)
    psi = np.arcsin(np.clip(lam_value / mu, -1.0, 1.0))
    return lam_value, psi


def _takagi_unitary(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return ascending singular values and U satisfying U.T @ M @ U > 0."""
    if np.max(np.abs(matrix - matrix.T)) > 1.0e-12:
        raise ValueError("Majorana matrix must be complex symmetric")
    eigenvalues, unitary = np.linalg.eigh(matrix.conjugate().T @ matrix)
    order = np.argsort(eigenvalues)
    eigenvalues = np.maximum(np.real(eigenvalues[order]), 0.0)
    unitary = unitary[:, order]
    congruence = unitary.T @ matrix @ unitary
    offdiag = congruence - np.diag(np.diag(congruence))
    tolerance = 1.0e-10 * max(1.0e-30, float(np.max(np.sqrt(eigenvalues))))
    if np.max(np.abs(offdiag)) > tolerance:
        raise ValueError("Takagi eigenspaces require a degenerate-block congruence resolution")
    unitary = unitary @ np.diag(np.exp(-0.5j * np.angle(np.diag(congruence))))
    diagonalized = unitary.T @ matrix @ unitary
    if np.max(np.abs(diagonalized - np.diag(np.diag(diagonalized)))) > tolerance:
        raise ValueError("Takagi congruence is not diagonal")
    if np.max(np.abs(np.imag(np.diag(diagonalized)))) > tolerance:
        raise ValueError("Takagi diagonal is not real")
    if np.any(np.real(np.diag(diagonalized)) < -tolerance):
        raise ValueError("Takagi diagonal is not positive")
    return np.sqrt(eigenvalues), unitary


@dataclass
class ExactEtaMap:
    a: float
    rho: float
    omega: float
    eta: np.ndarray
    mu: np.ndarray
    lam: float
    psi: np.ndarray
    majorana: np.ndarray
    h_shift: np.ndarray
    d_trace_shift: float
    p_invariant: float
    q_invariant: float
    cubic_roots: np.ndarray
    masses_squared: np.ndarray
    masses: np.ndarray
    u_takagi: np.ndarray


def _build_exact_eta_map(a_value: float, rho_value: float, omega: float, eta: np.ndarray) -> ExactEtaMap:
    mu = _normalized_mu_from_eta(eta)
    lam_value, psi = _solve_principal_selector(mu, omega)

    matrix = np.array(
        [
            [a_value, rho_value * np.exp(1j * psi[0]), rho_value * np.exp(1j * psi[2])],
            [rho_value * np.exp(1j * psi[0]), a_value, rho_value * np.exp(1j * psi[1])],
            [rho_value * np.exp(1j * psi[2]), rho_value * np.exp(1j * psi[1]), a_value],
        ],
        dtype=complex,
    )
    h_mat = matrix.conjugate().T @ matrix
    d_value = float(a_value * a_value + 2.0 * rho_value * rho_value)
    t_mat = h_mat - np.eye(3, dtype=complex) * d_value

    x12 = t_mat[0, 1]
    x23 = t_mat[1, 2]
    x13 = t_mat[0, 2]
    p_invariant = float(abs(x12) ** 2 + abs(x23) ** 2 + abs(x13) ** 2)
    q_invariant = float(np.real(x12 * x23 * np.conjugate(x13)))

    if p_invariant <= 0.0:
        cubic_roots = np.array([0.0, 0.0, 0.0], dtype=float)
    else:
        argument = 3.0 * math.sqrt(3.0) * q_invariant / (p_invariant ** 1.5)
        argument = max(-1.0, min(1.0, float(argument)))
        if 1.0 - abs(argument) < 1.0e-12:
            cubic_roots = np.linalg.eigvalsh(t_mat).astype(float)
        else:
            theta = math.acos(argument)
            scale = 2.0 * math.sqrt(p_invariant / 3.0)
            cubic_roots = np.array(
                [
                    scale * math.cos(theta / 3.0),
                    scale * math.cos(theta / 3.0 - 2.0 * math.pi / 3.0),
                    scale * math.cos(theta / 3.0 - 4.0 * math.pi / 3.0),
                ],
                dtype=float,
            )

    masses_squared = np.sort(d_value + cubic_roots)
    masses = np.sqrt(np.clip(masses_squared, 0.0, None))
    singular_values, u_mat = _takagi_unitary(matrix)

    return ExactEtaMap(
        a=a_value,
        rho=rho_value,
        omega=omega,
        eta=np.asarray(eta, dtype=float),
        mu=mu,
        lam=lam_value,
        psi=psi,
        majorana=matrix,
        h_shift=t_mat,
        d_trace_shift=d_value,
        p_invariant=p_invariant,
        q_invariant=q_invariant,
        cubic_roots=np.sort(cubic_roots),
        masses_squared=masses_squared,
        masses=masses,
        u_takagi=u_mat,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the exact intrinsic neutrino eta-chain artifact.")
    parser.add_argument("--isotropic", default=str(DEFAULT_ISOTROPIC))
    parser.add_argument("--payload", default=str(DEFAULT_PAYLOAD))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    isotropic = _load_json(Path(args.isotropic))
    payload = _load_json(Path(args.payload))

    isotropic_matrix = np.array(isotropic["majorana_matrix_real"], dtype=float) + 1j * np.array(
        isotropic["majorana_matrix_imag"], dtype=float
    )
    eta = _centered_eta_from_payload(payload)
    base_phase = float(np.angle(isotropic_matrix[0, 1]))
    omega = 3.0 * base_phase
    exact_map = _build_exact_eta_map(
        a_value=float(isotropic_matrix[0, 0].real),
        rho_value=float(abs(isotropic_matrix[0, 1])),
        omega=omega,
        eta=eta,
    )
    payload_is_live_certificate = bool(payload.get("sufficient_for_intrinsic_mass_eigenstates"))
    payload_source_closure = dict(payload.get("source_closure_status") or {})
    source_only_physical_input_eligible = (
        payload.get("source_only_physical_input_eligible") is True
        and payload_source_closure.get("closed") is True
    )

    actual_phase = _phase_vector_from_matrix(exact_map.majorana)
    for idx, reference in enumerate(exact_map.psi):
        actual_phase[idx] = _unwrap_near(float(actual_phase[idx]), float(reference))
    det_formula = (
        exact_map.a**3
        - exact_map.a * exact_map.rho**2 * np.sum(np.exp(2j * exact_map.psi))
        + 2.0 * exact_map.rho**3 * np.exp(1j * np.sum(exact_map.psi))
    )

    result = {
        "artifact": "oph_intrinsic_neutrino_exact_eta_map",
        "generated_utc": _timestamp(),
        "theorem_surface_status": (
            "intrinsic_builder_complete_exact_from_source_closed_inputs"
            if source_only_physical_input_eligible
            else "intrinsic_builder_exact_conditional_on_source_open_inputs"
        ),
        "source_only_physical_input_eligible": source_only_physical_input_eligible,
        "source_closure_status": payload_source_closure or {"closed": False},
        "builder_facing_exact_object": "centered_log_pullback_class_[log_q_e]",
        "proof_facing_residual_object": None if payload_is_live_certificate else "realized_arrow_pullback_from_flavor_gap_and_defect_certificates",
        "public_flavor_rows_gate": "pmns_and_flavor_rows_formed_downstream_from_shared_charged_basis" if payload_is_live_certificate else "blocked_pending_proof_facing_eta_emission_and_shared_charged_lepton_left_basis",
        "pmns_status": "not_formed_here",
        "payload_source": str(Path(args.payload)),
        "payload_kind": "eta_e" if payload.get("eta_e") else "q_or_mu_equivalent_payload",
        "isotropic_reference_bundle": str(Path(args.isotropic)),
        "base_parameters": {
            "a_diag_gev": exact_map.a,
            "rho_offdiag_gev": exact_map.rho,
            "phi_equal_split_rad": base_phase,
            "omega_cycle_rad": exact_map.omega,
            "trace_h_fixed_gev2": float(np.real(np.trace(exact_map.majorana.conjugate().T @ exact_map.majorana))),
        },
        "eta_e": {edge: float(exact_map.eta[idx]) for idx, edge in enumerate(EDGE_ORDER)},
        "mu_e_normalized": {edge: float(exact_map.mu[idx]) for idx, edge in enumerate(EDGE_ORDER)},
        "selector_equation": "sum_e arcsin(lam / mu_e) = Omega on the principal branch",
        "selector_domain": {
            "condition": "|Omega| < F_max, F_max = sum_e arcsin(mu_min / mu_e)",
            "F_max_rad": _principal_selector_domain(exact_map.mu),
            "abs_omega_rad": abs(float(exact_map.omega)),
            "margin_rad": _principal_selector_domain(exact_map.mu) - abs(float(exact_map.omega)),
            "satisfied": abs(float(exact_map.omega)) < _principal_selector_domain(exact_map.mu),
        },
        "selector_scale_free_multiplier": float(exact_map.lam),
        "selector_common_scale_invariant": True,
        "selector_point_absolute": {edge: float(actual_phase[idx]) for idx, edge in enumerate(EDGE_ORDER)},
        "selector_sum_check": float(np.sum(actual_phase)),
        "majorana_matrix_real": np.real(exact_map.majorana).tolist(),
        "majorana_matrix_imag": np.imag(exact_map.majorana).tolist(),
        "exact_cubic_invariants": {
            "depressed_cubic": "lambda^3 - P lambda - 2Q = 0 for H - d I",
            "d_trace_shift_gev2": exact_map.d_trace_shift,
            "p_invariant_gev4": exact_map.p_invariant,
            "q_invariant_gev6": exact_map.q_invariant,
            "cubic_argument_clamped": (
                0.0
                if exact_map.p_invariant <= 0.0
                else max(-1.0, min(1.0, 3.0 * math.sqrt(3.0) * exact_map.q_invariant / (exact_map.p_invariant ** 1.5)))
            ),
        },
        "cubic_roots_closed_form_gev2_sorted": [float(value) for value in exact_map.cubic_roots.tolist()],
        "h_shift_eigenvalues_gev2_sorted": [float(value) for value in np.linalg.eigvalsh(exact_map.h_shift).tolist()],
        "masses_gev_sorted": [float(value) for value in exact_map.masses.tolist()],
        "masses_from_h_eigvalsh_gev_sorted": [
            float(math.sqrt(max(value, 0.0)))
            for value in np.sort(np.linalg.eigvalsh(exact_map.majorana.conjugate().T @ exact_map.majorana)).tolist()
        ],
        "intrinsic_mass_eigenstates": {
            f"s{idx}": float(value) for idx, value in enumerate(exact_map.masses.tolist())
        },
        "mass_eigenstate_label_status": "ascending_singular_states_only",
        "physical_ordering_assignments": {
            "normal_ordering_hypothesis": {"nu1": "s0", "nu2": "s1", "nu3": "s2"},
            "inverted_ordering_hypothesis": {"nu3": "s0", "nu1": "s1", "nu2": "s2"},
            "selected": None,
            "missing_source_object": "solar_pair_and_atmospheric_sign_eigenstate_label_rule",
        },
        "mass_audit_max_abs_gev": float(
            np.max(
                np.abs(
                    exact_map.masses
                    - np.sqrt(
                        np.clip(
                            np.sort(np.linalg.eigvalsh(exact_map.majorana.conjugate().T @ exact_map.majorana)),
                            0.0,
                            None,
                        )
                    )
                )
            )
        ),
        "cubic_root_audit_max_abs_gev2": float(
            np.max(np.abs(np.sort(np.linalg.eigvalsh(exact_map.h_shift)) - exact_map.cubic_roots))
        ),
        "delta_m21_sq_gev2": float(exact_map.masses_squared[1] - exact_map.masses_squared[0]),
        "delta_m31_sq_gev2": float(exact_map.masses_squared[2] - exact_map.masses_squared[0]),
        "delta_m32_sq_gev2": float(exact_map.masses_squared[2] - exact_map.masses_squared[1]),
        "delta_mij_field_status": "ascending_gap_coordinates_under_declared_normal_ordering_hypothesis",
        "splitting_ratio_r": float(
            (exact_map.masses_squared[1] - exact_map.masses_squared[0])
            / (exact_map.masses_squared[2] - exact_map.masses_squared[0])
        ),
        "ordering_phase_certified": None,
        "ordering_status": "unresolved_without_mass_eigenstate_label_rule",
        "u_nu_takagi_real": np.real(exact_map.u_takagi).tolist(),
        "u_nu_takagi_imag": np.imag(exact_map.u_takagi).tolist(),
        "u_nu_left_real": np.real(exact_map.u_takagi).tolist(),
        "u_nu_left_imag": np.imag(exact_map.u_takagi).tolist(),
        "u_nu_field_semantics": "Takagi U satisfying U^T M U = diag(m_i) > 0; legacy u_nu_left keys are aliases",
        "takagi_congruence_max_offdiag_gev": float(
            np.max(
                np.abs(
                    exact_map.u_takagi.T
                    @ exact_map.majorana
                    @ exact_map.u_takagi
                    - np.diag(np.diag(exact_map.u_takagi.T @ exact_map.majorana @ exact_map.u_takagi))
                )
            )
        ),
        "determinant_formula_complex_gev3": {"real": float(np.real(det_formula)), "imag": float(np.imag(det_formula))},
        "determinant_direct_complex_gev3": {
            "real": float(np.real(np.linalg.det(exact_map.majorana))),
            "imag": float(np.imag(np.linalg.det(exact_map.majorana))),
        },
        "determinant_audit_abs_gev3": float(abs(np.linalg.det(exact_map.majorana) - det_formula)),
        "notes": [
            "This artifact computes the intrinsic builder-facing neutrino chain exactly from the centered eta-class.",
            "Once eta_e is emitted at the flavor boundary, no further selector ambiguity remains on the principal branch.",
            (
                "The same-label scalar certificate is numerically complete, but its upstream family kernel and overlap-line lift remain source-open."
                if payload_is_live_certificate and not source_only_physical_input_eligible
                else "The proof-facing eta provenance is supplied by a source-closed same-label scalar certificate; PMNS is formed downstream from the shared charged-lepton left basis."
                if source_only_physical_input_eligible
                else "The remaining supported blockers are proof-facing eta provenance and the shared charged-lepton left basis required for PMNS."
            ),
        ],
    }

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(portable_json_dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
