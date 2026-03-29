#!/usr/bin/env python3
"""Export the explicit odd-response-law boundary for the quark lane."""

from __future__ import annotations

import argparse
import json
import pathlib
from datetime import datetime, timezone
from typing import Any

import numpy as np


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_OBSERVABLE = ROOT / "particles" / "runs" / "flavor" / "flavor_observable_artifact.json"
DEFAULT_CHARGED_BUDGET = ROOT / "particles" / "runs" / "flavor" / "charged_budget_transport.json"
DEFAULT_TENSORS = ROOT / "particles" / "runs" / "flavor" / "suppression_phase_tensors.json"
DEFAULT_ODD_FORM = ROOT / "particles" / "runs" / "flavor" / "charged_dirac_odd_deformation_form.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "quark_odd_response_law.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _decode_complex_matrix(payload: Any) -> np.ndarray:
    if isinstance(payload, dict) and "real" in payload and "imag" in payload:
        return np.asarray(payload["real"], dtype=float) + 1j * np.asarray(payload["imag"], dtype=float)
    return np.asarray(payload, dtype=complex)


def _kappa_max_admissible(s_ch: np.ndarray, phi_ch: np.ndarray, delta_s: np.ndarray, delta_phi: np.ndarray) -> float | None:
    bounds: list[float] = []
    for i in range(3):
        for j in range(i + 1, 3):
            s_val = float(abs(delta_s[i, j]))
            phi_val = float(abs(delta_phi[i, j]))
            if s_val > 0.0:
                bounds.append(float(s_ch[i, j] / s_val))
            if phi_val > 0.0:
                bounds.append(float((np.pi - abs(phi_ch[i, j])) / phi_val))
    finite_bounds = [value for value in bounds if np.isfinite(value) and value > 0.0]
    if not finite_bounds:
        return None
    return float(min(finite_bounds))


def build_artifact(
    observable: dict[str, Any],
    charged_budget: dict[str, Any],
    tensors: dict[str, Any],
    odd_form: dict[str, Any] | None = None,
) -> dict[str, Any]:
    projectors = [_decode_complex_matrix(item) for item in observable.get("family_projectors", [])]
    if len(projectors) != 3:
        raise ValueError("flavor observable must provide three family projectors")
    eigenvalues = np.asarray(observable.get("family_eigenvalues", []), dtype=float)
    if eigenvalues.shape != (3,):
        raise ValueError("flavor observable must provide three family eigenvalues")

    centered = eigenvalues - float(np.mean(eigenvalues))
    degenerate_fallback = False
    norm = float(np.linalg.norm(centered))
    if norm <= 0.0:
        centered = np.asarray([1.0, 0.0, -1.0], dtype=float)
        norm = float(np.linalg.norm(centered))
        degenerate_fallback = True
    coefficients = centered / norm
    odd_splitter = sum(coefficients[idx] * projectors[idx] for idx in range(3))
    odd_splitter = 0.5 * (odd_splitter + odd_splitter.conj().T)

    s_u = np.asarray(tensors["S_u"], dtype=float)
    s_d = np.asarray(tensors["S_d"], dtype=float)
    phi_u = np.asarray(tensors["Phi_u"], dtype=float)
    phi_d = np.asarray(tensors["Phi_d"], dtype=float)
    s_ch = 0.5 * (s_u + s_d)
    phi_ch = 0.5 * (phi_u + phi_d)

    xi_real = np.real(odd_splitter)
    xi_imag = np.imag(odd_splitter)
    delta_s = xi_real - np.diag(np.diag(xi_real))
    delta_phi = xi_imag
    diag_part = np.real(np.diag(odd_splitter))
    delta_logd_left = diag_part
    delta_logd_right = -diag_part
    delta_logg_q = 0.0

    budget_total = list(charged_budget.get("B_ch_by_refinement", []))
    g_ch = float(budget_total[-1]["value"] / 3.0) if budget_total else float(tensors.get("u_raw_channel_norm_candidate", 1.0))
    charged_certificate = dict(charged_budget.get("charged_dirac_scalarization_certificate", {}))
    shared_scalarization_law_status = str(charged_certificate.get("law_status", "candidate_only"))
    corollary_status = str((odd_form or {}).get("quark_zero_odd_scalar_corollary", "open"))
    delta_logg_q_status = "closed_zero_corollary" if corollary_status == "closed" else "quotient_corollary_pending_odd_form"
    y_ch = g_ch * np.exp(-s_ch) * np.exp(1j * phi_ch)
    edge_weights = {
        "12": float(2.0 * abs(y_ch[0, 1]) ** 2),
        "13": float(2.0 * abs(y_ch[0, 2]) ** 2),
        "23": float(2.0 * abs(y_ch[1, 2]) ** 2),
    }
    diag_laplacian = [
        [edge_weights["12"] + edge_weights["13"], -edge_weights["12"], -edge_weights["13"]],
        [-edge_weights["12"], edge_weights["12"] + edge_weights["23"], -edge_weights["23"]],
        [-edge_weights["13"], -edge_weights["23"], edge_weights["13"] + edge_weights["23"]],
    ]
    kappa_max = _kappa_max_admissible(s_ch, phi_ch, delta_s, delta_phi)

    return {
        "artifact": "oph_quark_odd_response_law",
        "generated_utc": _timestamp(),
        "response_map_id": "projector_witness_pullback_riesz_candidate",
        "proof_status": "candidate_only",
        "coefficient_free": True,
        "orientation_fixed": True,
        "projector_covariant": True,
        "budget_neutral": True,
        "degenerate_placeholder_fallback_used": degenerate_fallback,
        "shared_scalarization_law_status": shared_scalarization_law_status,
        "delta_logg_q_status": delta_logg_q_status,
        "upstream_missing_object": "oph_charged_dirac_odd_deformation_form",
        "projector_split_coefficients": [float(x) for x in coefficients.tolist()],
        "Xi_q_real": xi_real.tolist(),
        "Xi_q_imag": xi_imag.tolist(),
        "shared_charged_response_seed": {
            "seed_kind": "averaged_u_d_candidate",
            "S_ch": s_ch.tolist(),
            "Phi_ch": phi_ch.tolist(),
            "logD_ch_left": [0.0, 0.0, 0.0],
            "logD_ch_right": [0.0, 0.0, 0.0],
            "g_ch_candidate": g_ch,
        },
        "lift_parameterization_kind": "single_kappa_signed_projector_ray",
        "charged_dirac_odd_deformation_form": None if odd_form is None else {
            "artifact": odd_form.get("artifact"),
            "proof_status": odd_form.get("proof_status"),
            "hidden_normalization_free": odd_form.get("hidden_normalization_free"),
            "quark_zero_odd_scalar_corollary": odd_form.get("quark_zero_odd_scalar_corollary"),
            "odd_scalar_slot_present": odd_form.get("odd_scalar_slot_present"),
        },
        "lift_constants": None,
        "cone_choice": "seed_admissibility_tangent_cone",
        "kappa": 1.0,
        "kappa_max_admissible": kappa_max,
        "edge_weights": edge_weights,
        "diag_laplacian": diag_laplacian,
        "Delta_S_q": delta_s.tolist(),
        "Delta_Phi_q": delta_phi.tolist(),
        "Delta_logD_left_q": [float(x) for x in delta_logd_left.tolist()],
        "Delta_logD_right_q": [float(x) for x in delta_logd_right.tolist()],
        "delta_logg_q": float(delta_logg_q),
        "metadata": {
            "note": "Explicit odd-response-law boundary for the quark lane. The current constructive candidate is the single-kappa signed projector ray on the shape-only odd codomain, with the dual covector taken from the shared-seed pullback Hilbert-Schmidt form. The common charged norm remains a separate shared-law question.",
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the quark odd-response-law artifact.")
    parser.add_argument("--observable", default=str(DEFAULT_OBSERVABLE))
    parser.add_argument("--charged-budget", default=str(DEFAULT_CHARGED_BUDGET))
    parser.add_argument("--tensors", default=str(DEFAULT_TENSORS))
    parser.add_argument("--odd-form", default=str(DEFAULT_ODD_FORM))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    observable = json.loads(pathlib.Path(args.observable).read_text(encoding="utf-8"))
    charged_budget = json.loads(pathlib.Path(args.charged_budget).read_text(encoding="utf-8"))
    tensors = json.loads(pathlib.Path(args.tensors).read_text(encoding="utf-8"))
    odd_form_path = pathlib.Path(args.odd_form)
    odd_form = json.loads(odd_form_path.read_text(encoding="utf-8")) if odd_form_path.exists() else None
    artifact = build_artifact(observable, charged_budget, tensors, odd_form=odd_form)

    out_path = pathlib.Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
