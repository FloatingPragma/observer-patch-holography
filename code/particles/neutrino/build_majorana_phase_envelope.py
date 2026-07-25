#!/usr/bin/env python3
"""Sweep the residual Majorana phase plane around the active selector branch.

Chain role: provide the ascending-spectrum stability envelope around the
current Majorana selector so that forward splittings can distinguish point
closure from law closure without inferring physical mass ordering.

Mathematics: residual-basis phase-plane sampling, complex Majorana assembly,
and ascending-gap/collective-mode-location stability certificates.

Declared pipeline inputs: the local scale anchor, family tensor, and Majorana
lift. The envelope certifies the stated matrix family, not its physical source
closure.

Output: the phase envelope used by the forward splittings and closure-bundle
export.
"""

from __future__ import annotations

import argparse
import json
import math
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


def _phase_vector_from_lift(lift: dict[str, Any]) -> np.ndarray:
    selector = dict(lift.get("selector_candidate_psi", {}))
    if not selector:
        selector = {
            "psi12": float(lift["cycle_constraint"]["omega_012"]) / 3.0,
            "psi23": float(lift["cycle_constraint"]["omega_012"]) / 3.0,
            "psi31": float(lift["cycle_constraint"]["omega_012"]) / 3.0,
        }
    return np.asarray([selector["psi12"], selector["psi23"], selector["psi31"]], dtype=float)


def _basis_matrix(lift: dict[str, Any]) -> np.ndarray:
    basis = lift.get("residual_basis", [])
    if len(basis) != 2:
        raise ValueError("majorana_holonomy_lift must provide two residual basis vectors")
    return np.asarray(
        [
            [basis[0]["psi12"], basis[0]["psi23"], basis[0]["psi31"]],
            [basis[1]["psi12"], basis[1]["psi23"], basis[1]["psi31"]],
        ],
        dtype=float,
    )


def _phase_matrix(psi: np.ndarray) -> np.ndarray:
    matrix = np.ones((3, 3), dtype=complex)
    matrix[0, 1] = np.exp(1j * psi[0])
    matrix[1, 0] = matrix[0, 1]
    matrix[1, 2] = np.exp(1j * psi[1])
    matrix[2, 1] = matrix[1, 2]
    matrix[2, 0] = np.exp(1j * psi[2])
    matrix[0, 2] = matrix[2, 0]
    return matrix


def _build_complex_majorana(m_star: float, e_nu: np.ndarray, diag_entries: np.ndarray, psi: np.ndarray) -> np.ndarray:
    n_diag = np.diag(diag_entries)
    return m_star * (n_diag @ (e_nu * _phase_matrix(psi)) @ n_diag)


def _sorted_takagi(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
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


def _collective_mode_location(overlaps: list[float]) -> str:
    dominant_index = int(max(range(3), key=lambda idx: overlaps[idx]))
    return f"s{dominant_index}"


def main() -> int:
    ap = argparse.ArgumentParser(description="Build the residual Majorana phase envelope artifact.")
    ap.add_argument("--scale-anchor", default="particles/runs/neutrino/neutrino_scale_anchor.json")
    ap.add_argument("--family", default="particles/runs/neutrino/family_response_tensor.json")
    ap.add_argument("--lift", default="particles/runs/neutrino/majorana_holonomy_lift.json")
    ap.add_argument("--majorana", default="particles/runs/neutrino/forward_majorana_matrix.json")
    ap.add_argument("--out", default="particles/runs/neutrino/majorana_phase_envelope.json")
    ap.add_argument("--grid", type=int, default=33, help="Residual-basis sweep resolution per axis.")
    args = ap.parse_args()

    scale_anchor_path = ROOT / args.scale_anchor if not pathlib.Path(args.scale_anchor).is_absolute() else pathlib.Path(args.scale_anchor)
    family_path = ROOT / args.family if not pathlib.Path(args.family).is_absolute() else pathlib.Path(args.family)
    lift_path = ROOT / args.lift if not pathlib.Path(args.lift).is_absolute() else pathlib.Path(args.lift)
    majorana_path = ROOT / args.majorana if not pathlib.Path(args.majorana).is_absolute() else pathlib.Path(args.majorana)
    out_path = ROOT / args.out if not pathlib.Path(args.out).is_absolute() else pathlib.Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    scale_anchor = load_json(scale_anchor_path)
    family = load_json(family_path)
    lift = load_json(lift_path)
    majorana = load_json(majorana_path)
    family_source_closure = dict(family.get("source_closure_status") or {"closed": False})
    lift_source_closure = dict(lift.get("source_closure_status") or {"closed": False})
    matrix_source_closure = dict(majorana.get("source_closure_status") or {"closed": False})
    source_closed = all(
        closure.get("closed") is True
        for closure in (family_source_closure, lift_source_closure, matrix_source_closure)
    )

    m_star = float(scale_anchor["anchors"]["m_star_gev"])
    e_nu = np.asarray(family["E_nu"], dtype=float)
    diag_entries = np.asarray(family["majorana_normalization_diag"], dtype=float)
    u_vector = np.asarray(scale_anchor["collective_mode"]["u_vector"], dtype=float)
    base_psi = _phase_vector_from_lift(lift)
    basis = _basis_matrix(lift)
    samples = np.linspace(-math.pi, math.pi, int(args.grid), dtype=float)

    mass_samples: list[list[float]] = []
    gap_s10_samples: list[float] = []
    gap_s20_samples: list[float] = []
    collective_mode_location_samples: list[str] = []

    for a_value in samples:
        for b_value in samples:
            psi = base_psi + a_value * basis[0] + b_value * basis[1]
            matrix = _build_complex_majorana(m_star, e_nu, diag_entries, psi)
            singular_values, left_vectors = _sorted_takagi(matrix)
            masses = [float(value) for value in singular_values.tolist()]
            overlaps = [float(abs(np.vdot(u_vector, left_vectors[:, idx])) ** 2) for idx in range(3)]
            mass_samples.append(masses)
            gap_s10_samples.append((masses[1] ** 2) - (masses[0] ** 2))
            gap_s20_samples.append((masses[2] ** 2) - (masses[0] ** 2))
            collective_mode_location_samples.append(_collective_mode_location(overlaps))

    masses_array = np.asarray(mass_samples, dtype=float)
    real_seed_masses = np.asarray(majorana["masses_sorted_gev"], dtype=float)
    real_seed_gaps = [float(real_seed_masses[1] - real_seed_masses[0]), float(real_seed_masses[2] - real_seed_masses[1])]
    weighted_edge_norm_sq = float(family["weighted_edge_norm_sq"])
    delta_sigma_radius = float(2.0 * math.sqrt(2.0) * m_star * math.sqrt(weighted_edge_norm_sq))
    m_max = float(real_seed_masses[-1])
    splitting_radius = float(4.0 * m_max * delta_sigma_radius + 2.0 * (delta_sigma_radius ** 2))
    projector_gap_seed = float(family["projector_gap_seed"])
    collective_projector_phase_stable = bool(projector_gap_seed > 2.0 * splitting_radius)
    collective_mode_location_phase_stable = bool(
        all(gap > 2.0 * delta_sigma_radius for gap in real_seed_gaps)
        and collective_projector_phase_stable
        and len(set(collective_mode_location_samples)) == 1
    )

    payload = {
        "artifact": "oph_neutrino_majorana_phase_envelope",
        "status": "residual_phase_envelope",
        "proof_scope": "spectral_envelope_conditional_on_declared_matrix_family",
        "source_only_physical_input_eligible": source_closed,
        "public_surface_candidate_allowed": False,
        "source_closure_status": {
            "closed": source_closed,
            "family_response": family_source_closure,
            "majorana_lift": lift_source_closure,
            "forward_matrix": matrix_source_closure,
        },
        "phase_mode": "residual_envelope",
        "inputs": {
            "scale_anchor_artifact": str(scale_anchor_path),
            "family_artifact": str(family_path),
            "lift_artifact": str(lift_path),
            "majorana_artifact": str(majorana_path),
        },
        "sample_grid_per_axis": int(args.grid),
        "sample_count": int(len(mass_samples)),
        "selector_reference": lift.get("selector_candidate_psi"),
        "mass_bounds_gev": [
            {"state": "s0", "min": float(np.min(masses_array[:, 0])), "max": float(np.max(masses_array[:, 0]))},
            {"state": "s1", "min": float(np.min(masses_array[:, 1])), "max": float(np.max(masses_array[:, 1]))},
            {"state": "s2", "min": float(np.min(masses_array[:, 2])), "max": float(np.max(masses_array[:, 2]))},
        ],
        "ascending_gap_s10_sq_bounds_gev2": {
            "min": float(min(gap_s10_samples)),
            "max": float(max(gap_s10_samples)),
        },
        "ascending_gap_s20_sq_bounds_gev2": {
            "min": float(min(gap_s20_samples)),
            "max": float(max(gap_s20_samples)),
        },
        "collective_mode_location_samples": sorted(set(collective_mode_location_samples)),
        "collective_mode_location_phase_stable": collective_mode_location_phase_stable,
        "physical_mass_ordering_status": "not_established_mass_label_rule_absent",
        "collective_projector_phase_stable": collective_projector_phase_stable,
        "delta_sigma_radius_gev": delta_sigma_radius,
        "gap_vs_radius_certificate": {
            "real_seed_gaps_gev": real_seed_gaps,
            "projector_gap_seed_gev2": projector_gap_seed,
            "delta_sigma_radius_gev": delta_sigma_radius,
            "delta_h_radius_gev2": splitting_radius,
            "collective_location_gap_test_passes": bool(all(gap > 2.0 * delta_sigma_radius for gap in real_seed_gaps)),
            "projector_gap_test_passes": collective_projector_phase_stable,
        },
        "notes": [
            "This is the compulsory support-check surface for phase dependence in the neutrino lane.",
            "The envelope bounds ascending singular states and gaps, not physically labeled neutrino splittings.",
            "Phase stability of a collective spectral mode cannot supply the missing mass-eigenstate label rule.",
        ],
    }
    out_path.write_text(portable_json_dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(out_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
