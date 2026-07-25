#!/usr/bin/env python3
"""Assemble the blind forward Majorana matrix from the local neutrino chain.

Chain role: turn the scale anchor, family tensor, and phase selector data into
an explicit complex Majorana mass matrix.

Mathematics: matrix assembly across `real_seed`, `canonical_selector`, and
`residual_envelope` modes plus Takagi spectral extraction.

Declared pipeline inputs: the local scale anchor, family-response tensor,
Majorana lift, and current phase-law/envelope artifacts. This matrix builder
does not itself establish their physical source closure.

Output: the forward Majorana matrix used for ascending singular-value gaps and
the exported neutrino closure bundle. Physical mass ordering remains a separate
source-label problem.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys
from typing import Any

import numpy as np


ROOT = pathlib.Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from particles.artifact_paths import portable_json_dumps


def load_json(path: pathlib.Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _phase_matrix_from_selector(selector: dict[str, Any] | None) -> np.ndarray:
    phase_matrix = np.ones((3, 3), dtype=complex)
    if not selector:
        return phase_matrix
    phase_matrix[0, 1] = np.exp(1j * float(selector["psi12"]))
    phase_matrix[1, 0] = phase_matrix[0, 1]
    phase_matrix[1, 2] = np.exp(1j * float(selector["psi23"]))
    phase_matrix[2, 1] = phase_matrix[1, 2]
    phase_matrix[2, 0] = np.exp(1j * float(selector["psi31"]))
    phase_matrix[0, 2] = phase_matrix[2, 0]
    return phase_matrix


def _sorted_takagi(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
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
    return np.sqrt(eigenvalues), unitary


def main() -> int:
    ap = argparse.ArgumentParser(description="Build a blind forward Majorana matrix artifact.")
    ap.add_argument("--scale-anchor", default="particles/runs/neutrino/neutrino_scale_anchor.json")
    ap.add_argument("--family", default="particles/runs/neutrino/family_response_tensor.json")
    ap.add_argument("--lift", default="particles/runs/neutrino/majorana_holonomy_lift.json")
    ap.add_argument("--pullback-metric", default="particles/runs/neutrino/majorana_phase_pullback_metric.json")
    ap.add_argument("--envelope", default="particles/runs/neutrino/majorana_phase_envelope.json")
    ap.add_argument(
        "--mode",
        default="canonical_selector",
        choices=["real_seed", "canonical_selector", "residual_envelope"],
        help="Which neutrino phase mode to export.",
    )
    ap.add_argument("--out", default="particles/runs/neutrino/forward_majorana_matrix.json")
    args = ap.parse_args()

    scale_anchor_path = pathlib.Path(args.scale_anchor)
    if not scale_anchor_path.is_absolute():
        scale_anchor_path = ROOT / args.scale_anchor
    family_path = pathlib.Path(args.family)
    if not family_path.is_absolute():
        family_path = ROOT / args.family
    lift_path = pathlib.Path(args.lift)
    if not lift_path.is_absolute():
        lift_path = ROOT / args.lift
    pullback_metric_path = pathlib.Path(args.pullback_metric)
    if not pullback_metric_path.is_absolute():
        pullback_metric_path = ROOT / args.pullback_metric
    envelope_path = pathlib.Path(args.envelope)
    if not envelope_path.is_absolute():
        envelope_path = ROOT / args.envelope
    out_path = pathlib.Path(args.out)
    if not out_path.is_absolute():
        out_path = ROOT / args.out
    out_path.parent.mkdir(parents=True, exist_ok=True)

    scale_anchor = load_json(scale_anchor_path)
    family = load_json(family_path)
    lift = load_json(lift_path) if lift_path.exists() else None
    pullback_metric = load_json(pullback_metric_path) if pullback_metric_path.exists() else None
    envelope = load_json(envelope_path) if envelope_path.exists() else None
    family_source_closure = dict(family.get("source_closure_status") or {"closed": False})
    lift_source_closure = dict((lift or {}).get("source_closure_status") or {"closed": False})
    source_closed = (
        family_source_closure.get("closed") is True
        and lift_source_closure.get("closed") is True
    )
    m_star = float(scale_anchor["anchors"]["m_star_gev"])
    c_nu_hat_real = np.array(family["C_nu_hat_real"], dtype=float)
    delta_nu = np.array(family["Delta_nu"], dtype=float)
    e_nu = np.array(family["E_nu"], dtype=float)
    diag_entries = np.array(family["majorana_normalization_diag"], dtype=float)
    n_diag = np.diag(diag_entries)
    selector = None
    selector_status = None if lift is None else lift.get("canonical_selector_status")
    phase_mode = str(args.mode)
    certification_status = "real_seed_surrogate"
    selector_used = None
    envelope_used = None
    selector_point_certified = False
    selector_law_certified = False

    if phase_mode == "real_seed":
        majorana_matrix = m_star * c_nu_hat_real.astype(complex)
    elif phase_mode == "canonical_selector":
        selector = None if lift is None else dict(lift.get("canonical_selector_point") or lift.get("selector_candidate_psi", {}))
        if not selector:
            raise ValueError("canonical_selector mode requires a canonical selector point in the lift artifact")
        majorana_matrix = m_star * (n_diag @ (e_nu * _phase_matrix_from_selector(selector)) @ n_diag)
        selector_used = selector.get("selector")
        selector_point_certified = selector_status in {"closed_equal_split", "closed_least_distortion"}
        selector_law_certified = bool((pullback_metric or {}).get("phase_action_closed", False)) and selector_point_certified
        if selector_status == "closed_equal_split":
            certification_status = "selector_closed_equal_split"
        elif selector_status == "closed_least_distortion":
            certification_status = "selector_closed_least_distortion"
        else:
            certification_status = "selector_candidate_only"
        if selector_law_certified:
            certification_status = "selector_law_closed_standard_math"
    else:
        majorana_matrix = m_star * c_nu_hat_real.astype(complex)
        envelope_used = str(envelope_path) if envelope is not None else None
        certification_status = "envelope_only_no_single_matrix"

    if phase_mode == "real_seed" and lift is not None and float(lift["cycle_constraint"]["omega_012"]) != 0.0:
        certification_status = "real_seed_phase_unresolved"

    if np.allclose(np.imag(majorana_matrix), 0.0, atol=1.0e-18):
        eigenvalues, eigenvectors = np.linalg.eigh(np.real(majorana_matrix))
        order = np.argsort(np.abs(eigenvalues))
        eigenvalues = eigenvalues[order]
        eigenvectors = eigenvectors[:, order]
        masses = np.abs(eigenvalues)
        takagi_vectors = eigenvectors.astype(complex)
        takagi_vectors[:, eigenvalues < 0.0] *= 1j
        raw_eigenvalues = [float(value) for value in eigenvalues.tolist()]
    else:
        singular_values, takagi_vectors = _sorted_takagi(majorana_matrix)
        masses = singular_values
        raw_eigenvalues = None

    u_vector = np.asarray(scale_anchor["collective_mode"]["u_vector"], dtype=float)
    collective_overlaps = [float(abs(np.vdot(u_vector, takagi_vectors[:, idx])) ** 2) for idx in range(3)]
    principal_minors = [
        float(np.real(majorana_matrix[0, 0])),
        float(abs(np.linalg.det(majorana_matrix[:2, :2]))),
        float(abs(np.linalg.det(majorana_matrix))),
    ]
    payload = {
        "artifact": "oph_neutrino_forward_majorana_matrix",
        "status": "blind_forward_matrix",
        "proof_scope": "exact_matrix_and_takagi_algebra_conditional_on_declared_inputs",
        "source_only_physical_input_eligible": source_closed,
        "public_surface_candidate_allowed": False,
        "source_closure_status": {
            "closed": source_closed,
            "family_response": family_source_closure,
            "majorana_lift": lift_source_closure,
        },
        "majorana_mode": phase_mode,
        "inputs": {
            "scale_anchor_artifact": str(scale_anchor_path),
            "family_artifact": str(family_path),
            "lift_artifact": str(lift_path) if lift is not None else None,
            "pullback_metric_artifact": str(pullback_metric_path) if pullback_metric is not None else None,
            "envelope_artifact": str(envelope_path) if envelope is not None else None,
        },
        "m_star_gev": m_star,
        "phase_mode": phase_mode,
        "selector_used": selector_used,
        "envelope_used": envelope_used,
        "selector_law_status": None if lift is None else lift.get("selector_law_status"),
        "selector_closure_reason": None if lift is None else lift.get("selector_closure_reason"),
        "selector_point_certified": selector_point_certified,
        "selector_law_certified": selector_law_certified,
        "law_closure_scope": None if pullback_metric is None else pullback_metric.get("law_closure_scope"),
        "ambient_metric_kind": None if pullback_metric is None else pullback_metric.get("ambient_metric_kind"),
        "certification_status": certification_status,
        "majorana_matrix_real": np.real(majorana_matrix).tolist(),
        "majorana_matrix_imag": np.imag(majorana_matrix).tolist(),
        "C_nu_hat_real": c_nu_hat_real.tolist(),
        "Delta_nu": delta_nu.tolist(),
        "eigenvalues_raw_gev": raw_eigenvalues,
        "singular_values_raw_gev": [float(value) for value in masses.tolist()],
        "masses_sorted_gev": [float(value) for value in masses.tolist()],
        "U_nu_real": np.real(takagi_vectors).tolist(),
        "U_nu_imag": np.imag(takagi_vectors).tolist(),
        "u_nu_left_real": np.real(takagi_vectors).tolist(),
        "u_nu_left_imag": np.imag(takagi_vectors).tolist(),
        "U_nu_semantics": "Takagi U satisfying U^T M U = diag(m_i) > 0; legacy left-vector keys are aliases",
        "takagi_congruence_max_offdiag_gev": float(
            np.max(
                np.abs(
                    takagi_vectors.T @ majorana_matrix @ takagi_vectors
                    - np.diag(np.diag(takagi_vectors.T @ majorana_matrix @ takagi_vectors))
                )
            )
        ),
        "rank_proxy": int(np.linalg.matrix_rank(majorana_matrix)),
        "collective_mode_overlap_by_eigenvector": collective_overlaps,
        "phase_status": None if lift is None else lift.get("phase_status"),
        "invariants": {
            "trace_real": float(np.real(np.trace(majorana_matrix))),
            "trace_h": float(np.real(np.trace(majorana_matrix.conj().T @ majorana_matrix))),
            "det_abs": float(abs(np.linalg.det(majorana_matrix))),
            "principal_minors": principal_minors,
        },
        "notes": [
            "The real-seed branch is a surrogate when the Majorana selector is not theorem-closed.",
            "Canonical-selector output may close as a selector point or, under the pullback-distortion route, as a selector law with explicit scope.",
            "Selector/Takagi closure is conditional algebra and does not promote a source-open family response or Majorana lift.",
            "Residual-envelope mode carries certification metadata rather than claiming a unique complex Majorana matrix.",
        ],
    }
    out_path.write_text(portable_json_dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(out_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
