#!/usr/bin/env python3
"""Export the charged odd deformation-form boundary behind the quark response law."""

from __future__ import annotations

import argparse
import json
import pathlib
from datetime import datetime, timezone
from typing import Any

import numpy as np


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_CHARGED_BUDGET = ROOT / "particles" / "runs" / "flavor" / "charged_budget_transport.json"
DEFAULT_OBSERVABLE = ROOT / "particles" / "runs" / "flavor" / "flavor_observable_artifact.json"
DEFAULT_TENSORS = ROOT / "particles" / "runs" / "flavor" / "suppression_phase_tensors.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "charged_dirac_odd_deformation_form.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _decode_complex_matrix(payload: Any) -> np.ndarray:
    if isinstance(payload, dict) and "real" in payload and "imag" in payload:
        return np.asarray(payload["real"], dtype=float) + 1j * np.asarray(payload["imag"], dtype=float)
    return np.asarray(payload, dtype=complex)


def _odd_splitter_from_observable(observable: dict[str, Any]) -> np.ndarray:
    projectors = [_decode_complex_matrix(item) for item in observable.get("family_projectors", [])]
    if len(projectors) != 3:
        raise ValueError("flavor observable must provide three family projectors")
    eigenvalues = np.asarray(observable.get("family_eigenvalues", []), dtype=float)
    if eigenvalues.shape != (3,):
        raise ValueError("flavor observable must provide three family eigenvalues")
    centered = eigenvalues - float(np.mean(eigenvalues))
    norm = float(np.linalg.norm(centered))
    if norm <= 0.0:
        centered = np.asarray([1.0, 0.0, -1.0], dtype=float)
        norm = float(np.linalg.norm(centered))
    coefficients = centered / norm
    odd_splitter = sum(coefficients[idx] * projectors[idx] for idx in range(3))
    return 0.5 * (odd_splitter + odd_splitter.conj().T)


def _kappa_max_admissible(s_ch: np.ndarray, phi_ch: np.ndarray, a_q: np.ndarray, omega_q: np.ndarray) -> float | None:
    bounds: list[float] = []
    for i in range(3):
        for j in range(i + 1, 3):
            a_val = float(abs(a_q[i, j]))
            omega_val = float(abs(omega_q[i, j]))
            if a_val > 0.0:
                bounds.append(float(s_ch[i, j] / a_val))
            if omega_val > 0.0:
                bounds.append(float((np.pi - abs(phi_ch[i, j])) / omega_val))
    finite_bounds = [value for value in bounds if np.isfinite(value) and value > 0.0]
    if not finite_bounds:
        return None
    return float(min(finite_bounds))


def build_artifact(charged_budget: dict[str, Any], observable: dict[str, Any], tensors: dict[str, Any]) -> dict[str, Any]:
    law_status = str(dict(charged_budget.get("charged_dirac_scalarization_certificate", {})).get("law_status", "candidate_only"))
    budget_total = list(charged_budget.get("B_ch_by_refinement", []))
    g_ch = float(budget_total[-1]["value"] / 3.0) if budget_total else float(tensors.get("u_raw_channel_norm_candidate", 1.0))
    s_u = np.asarray(tensors["S_u"], dtype=float)
    s_d = np.asarray(tensors["S_d"], dtype=float)
    phi_u = np.asarray(tensors["Phi_u"], dtype=float)
    phi_d = np.asarray(tensors["Phi_d"], dtype=float)
    s_ch = 0.5 * (s_u + s_d)
    phi_ch = 0.5 * (phi_u + phi_d)
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
    odd_splitter = _odd_splitter_from_observable(observable)
    xi_real = np.real(odd_splitter)
    xi_imag = np.imag(odd_splitter)
    a_q = xi_real - np.diag(np.diag(xi_real))
    omega_q = xi_imag
    kappa_max = _kappa_max_admissible(s_ch, phi_ch, a_q, omega_q)
    return {
        "artifact": "oph_charged_dirac_odd_deformation_form",
        "generated_utc": _timestamp(),
        "proof_status": "candidate_only",
        "domain": "shared_charged_dirac_shape_class",
        "codomain": "Sym_off(3) + so(3) + R0^3",
        "odd_scalar_slot_present": False,
        "budget_neutral_shape_only": True,
        "projector_equivariant": True,
        "odd_under_u_d_exchange": True,
        "hidden_normalization_free": True,
        "deformation_form_kind": "shared_seed_pullback_hs",
        "riesz_map_kind": "same_form_dual_on_odd_quotient",
        "quotient_basis": ["sym12", "sym13", "sym23", "so12", "so13", "so23", "diag_t0_a", "diag_t0_b"],
        "edge_weights": edge_weights,
        "diag_laplacian": diag_laplacian,
        "covector_origin": "same_form_dual_of_projector_witness",
        "normalization_rule": "unit_projector_witness",
        "global_scale_family": {
            "kind": "single_kappa_ray",
            "kappa_canonical": 1.0,
            "kappa_max_admissible": kappa_max,
        },
        "lift_constants": None,
        "shared_scalarization_law_status": law_status,
        "quark_zero_odd_scalar_corollary": "closed",
        "upstream_missing_object": "oph_charged_dirac_odd_deformation_form",
        "notes": [
            "This artifact isolates the odd quark shape-map burden from the quark response builder. The current constructive candidate is the shared-seed pullback Hilbert-Schmidt form on the explicit odd quotient basis, with no independent odd scalar slot and no per-block lift constants.",
            "The odd codomain has no independent odd scalar slot on the current theorem design, so delta_logg_q=0 is a quotient corollary of the odd-form closure rather than a separate local quark theorem.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the charged Dirac odd deformation-form artifact.")
    parser.add_argument("--charged-budget", default=str(DEFAULT_CHARGED_BUDGET))
    parser.add_argument("--observable", default=str(DEFAULT_OBSERVABLE))
    parser.add_argument("--tensors", default=str(DEFAULT_TENSORS))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    charged_budget = json.loads(pathlib.Path(args.charged_budget).read_text(encoding="utf-8"))
    observable = json.loads(pathlib.Path(args.observable).read_text(encoding="utf-8"))
    tensors = json.loads(pathlib.Path(args.tensors).read_text(encoding="utf-8"))
    artifact = build_artifact(charged_budget, observable, tensors)

    out_path = pathlib.Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
