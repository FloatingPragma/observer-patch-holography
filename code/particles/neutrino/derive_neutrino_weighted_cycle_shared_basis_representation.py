#!/usr/bin/env python3
"""Transport the repaired weighted-cycle branch into the shared same-label basis.

Chain role: expose the repaired weighted-cycle Majorana branch on the explicit
`[f1, f2, f3]` basis closed by `U_e_left`, so the physical PMNS path is written
as `U_PMNS = U_e_left^dagger * U_nu_shared` rather than left inside row-order
metadata.

Mathematics:
1. The repaired weighted-cycle artifact already emits one canonical Takagi
   charged-basis surface `U_wc` with positive masses `m_i`.
2. Once the shared charged-lepton left basis closes, the same Majorana matrix is
   transported to the shared same-label basis by the exact Majorana congruence
   `M_shared = U_e_left^* M_wc U_e_left^dagger`.
3. The corresponding same-label neutrino Takagi unitary is
   `U_nu_shared = U_e_left U_wc`.
4. Therefore `U_e_left^dagger U_nu_shared = U_wc`, so the physical PMNS surface
   is rebuilt explicitly through the closed shared basis without introducing any
   new external convention.

OPH-derived inputs: the repaired weighted-cycle branch and the shared
charged-lepton left basis.

Output: an exact shared-basis representation of the repaired weighted-cycle
physical PMNS/Majorana branch.
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
DEFAULT_WEIGHTED_CYCLE = ROOT / "particles" / "runs" / "neutrino" / "neutrino_weighted_cycle_repair.json"
DEFAULT_SHARED_CHARGED_LEFT = ROOT / "particles" / "runs" / "neutrino" / "shared_charged_lepton_left_basis.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "neutrino" / "neutrino_weighted_cycle_shared_basis_representation.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _complex_matrix(payload: dict[str, Any], real_key: str, imag_key: str) -> np.ndarray:
    return np.array(payload[real_key], dtype=float) + 1j * np.array(payload[imag_key], dtype=float)


def _wrap_signed_phase(angle_rad: float) -> float:
    wrapped = float((angle_rad + math.pi) % (2.0 * math.pi) - math.pi)
    if wrapped <= -math.pi:
        return wrapped + 2.0 * math.pi
    return wrapped


def _standard_pmns_parameters(unitary: np.ndarray) -> dict[str, float]:
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
        sin_delta = 0.0
        den_j = c12 * s12 * c23 * s23 * (c13**2) * s13
        if abs(den_j) > 1.0e-30:
            sin_delta = float(np.clip(jarlskog / den_j, -1.0, 1.0))
        delta = math.atan2(sin_delta, cos_delta) % (2.0 * math.pi)

    return {
        "theta12_rad": float(theta12),
        "theta23_rad": float(theta23),
        "theta13_rad": float(theta13),
        "delta_rad": float(delta),
        "theta12_deg": math.degrees(theta12),
        "theta23_deg": math.degrees(theta23),
        "theta13_deg": math.degrees(theta13),
        "delta_deg": math.degrees(delta),
        "J": jarlskog,
    }


def _majorana_pair_from_pmns(unitary: np.ndarray, delta_rad: float) -> dict[str, float]:
    electron_row_phase = float(np.angle(unitary[0, 0]))
    row_gauged = np.diag([np.exp(-1j * electron_row_phase), 1.0, 1.0]) @ unitary
    alpha21_signed = _wrap_signed_phase(2.0 * float(np.angle(row_gauged[0, 1])))
    alpha31_signed = _wrap_signed_phase(2.0 * (float(np.angle(row_gauged[0, 2])) + delta_rad))
    alpha21_mod = float(alpha21_signed % (2.0 * math.pi))
    alpha31_mod = float(alpha31_signed % (2.0 * math.pi))
    return {
        "electron_row_gauge_phase_rad": electron_row_phase,
        "electron_row_gauge_phase_deg": math.degrees(electron_row_phase),
        "alpha21_rad": alpha21_signed,
        "alpha21_deg": math.degrees(alpha21_signed),
        "alpha21_rad_0_to_2pi": alpha21_mod,
        "alpha21_deg_0_to_360": math.degrees(alpha21_mod),
        "alpha31_rad": alpha31_signed,
        "alpha31_deg": math.degrees(alpha31_signed),
        "alpha31_rad_0_to_2pi": alpha31_mod,
        "alpha31_deg_0_to_360": math.degrees(alpha31_mod),
    }


def _canonical_takagi_unitary(matrix: np.ndarray) -> dict[str, Any]:
    if np.max(np.abs(matrix - matrix.T)) > 1.0e-12:
        raise ValueError("weighted-cycle matrix must remain complex symmetric")

    left_vectors, singular_values, vh = np.linalg.svd(matrix)
    order = np.argsort(singular_values)
    left_vectors = left_vectors[:, order]
    singular_values = singular_values[order]
    right_vectors = vh.conj().T[:, order]

    phase_relation = left_vectors.T @ right_vectors
    offdiag_relation = phase_relation - np.diag(np.diag(phase_relation))
    relation_offdiag_max = float(np.max(np.abs(offdiag_relation)))
    if relation_offdiag_max > 1.0e-8:
        raise ValueError("SVD phase relation is not diagonal enough to define the canonical Takagi gauge")

    congruence_half_angles = np.angle(np.diag(phase_relation))
    takagi_vectors = left_vectors @ np.diag(np.exp(-0.5j * congruence_half_angles))
    unitary = np.conjugate(takagi_vectors)
    diagonalized = unitary.T @ matrix @ unitary
    diagonalized_offdiag = diagonalized - np.diag(np.diag(diagonalized))
    diagonal_masses = np.real(np.diag(diagonalized))
    if np.any(diagonal_masses <= 0.0):
        raise ValueError("canonical Takagi diagonal must stay positive")

    return {
        "unitary": unitary,
        "singular_values": [float(x) for x in singular_values.tolist()],
        "congruence_half_angles_rad": [float(0.5 * angle) for angle in congruence_half_angles.tolist()],
        "congruence_half_angles_deg": [math.degrees(float(0.5 * angle)) for angle in congruence_half_angles.tolist()],
        "phase_relation_offdiag_max_abs": relation_offdiag_max,
        "diagonalized_real_masses": [float(x) for x in diagonal_masses.tolist()],
        "diagonalized_imag_max_abs": float(np.max(np.abs(np.imag(np.diag(diagonalized))))),
        "diagonalized_offdiag_max_abs": float(np.max(np.abs(diagonalized_offdiag))),
    }


def build_payload(
    weighted_cycle: dict[str, Any],
    shared_charged_left: dict[str, Any],
    *,
    weighted_cycle_path: Path,
    shared_charged_left_path: Path,
) -> dict[str, Any]:
    basis_contract = dict(shared_charged_left.get("basis_contract", {}))
    if shared_charged_left.get("status") != "closed":
        raise ValueError("shared charged-lepton left basis must be closed")
    if not basis_contract.get("orientation_preserved", False):
        raise ValueError("shared charged-lepton basis must preserve orientation")
    if weighted_cycle.get("physical_window_status") != "pmns_and_hierarchy_repaired":
        raise ValueError("weighted-cycle repair must close the physical PMNS window first")

    weighted_cycle_matrix = _complex_matrix(weighted_cycle, "repaired_cycle_matrix_real", "repaired_cycle_matrix_imag")
    takagi = _canonical_takagi_unitary(weighted_cycle_matrix)
    pmns = takagi["unitary"]
    weighted_cycle_observables = dict(weighted_cycle["pmns_observables"])

    u_e_left = _complex_matrix(shared_charged_left["U_e_left"], "real", "imag")
    shared_basis_matrix = u_e_left.conj() @ weighted_cycle_matrix @ u_e_left.conj().T
    u_nu_shared = u_e_left @ pmns
    recovered_pmns = np.conjugate(u_e_left).T @ u_nu_shared

    shared_diagonalized = u_nu_shared.T @ shared_basis_matrix @ u_nu_shared
    shared_offdiag = shared_diagonalized - np.diag(np.diag(shared_diagonalized))
    shared_diag_real = np.real(np.diag(shared_diagonalized))
    if np.any(shared_diag_real <= 0.0):
        raise ValueError("shared-basis transport must preserve a positive Takagi diagonal")

    observables = _standard_pmns_parameters(recovered_pmns)
    observable_match = {
        "theta12_deg_abs_delta": abs(observables["theta12_deg"] - float(weighted_cycle_observables["theta12_deg"])),
        "theta23_deg_abs_delta": abs(observables["theta23_deg"] - float(weighted_cycle_observables["theta23_deg"])),
        "theta13_deg_abs_delta": abs(observables["theta13_deg"] - float(weighted_cycle_observables["theta13_deg"])),
        "delta_deg_abs_delta": abs(((observables["delta_deg"] - float(weighted_cycle_observables["delta_deg"]) + 180.0) % 360.0) - 180.0),
        "J_abs_delta": abs(observables["J"] - float(weighted_cycle_observables["J"])),
    }
    if max(observable_match.values()) > 1.0e-8:
        raise ValueError("shared-basis transported branch must recover the weighted-cycle PMNS observables exactly")

    majorana_pair = _majorana_pair_from_pmns(recovered_pmns, observables["delta_rad"])
    transport_checks = {
        "shared_basis_symmetry_max_abs": float(np.max(np.abs(shared_basis_matrix - shared_basis_matrix.T))),
        "shared_basis_diagonalized_offdiag_max_abs": float(np.max(np.abs(shared_offdiag))),
        "shared_basis_diagonalized_imag_max_abs": float(np.max(np.abs(np.imag(np.diag(shared_diagonalized))))),
        "shared_basis_diagonalized_real_masses": [float(x) for x in shared_diag_real.tolist()],
        "pmns_recovery_max_abs": float(np.max(np.abs(recovered_pmns - pmns))),
    }

    return {
        "artifact": "oph_neutrino_weighted_cycle_shared_basis_representation",
        "generated_utc": _timestamp(),
        "status": "theorem_grade_emitted",
        "proof_chain_role": "active_theorem_lane",
        "theorem_object": "shared_basis_weighted_cycle_physical_pmns",
        "theorem_surface": "weighted_cycle_shared_basis_transport",
        "public_surface_candidate_allowed": True,
        "physical_branch_closed": True,
        "statement": (
            "The repaired weighted-cycle Majorana matrix is transported exactly into the closed shared same-label basis, "
            "and the physical PMNS surface is recovered there as U_e_left^dagger * U_nu_shared. This closes the "
            "weighted-cycle theorem lane on the shared basis without identifying it with the separate intrinsic/shared-"
            "basis PMNS diagnostic surface."
        ),
        "construction": {
            "shared_basis_matrix": "M_shared = U_e_left^* M_wc U_e_left^dagger",
            "shared_basis_neutrino_unitary": "U_nu_shared = U_e_left U_wc",
            "physical_pmns_recovery": "U_PMNS = U_e_left^dagger U_nu_shared",
        },
        "basis_contract": basis_contract,
        "basis_labels": list(shared_charged_left.get("labels") or []),
        "source_artifacts": {
            "weighted_cycle_branch": str(weighted_cycle_path),
            "shared_charged_left_basis": str(shared_charged_left_path),
        },
        "depends_on": [
            "oph_neutrino_weighted_cycle_repair",
            "oph_shared_charged_lepton_left_basis",
        ],
        "weighted_cycle_takagi_congruence": {
            "singular_values": takagi["singular_values"],
            "congruence_half_angles_rad": takagi["congruence_half_angles_rad"],
            "congruence_half_angles_deg": takagi["congruence_half_angles_deg"],
            "phase_relation_offdiag_max_abs": takagi["phase_relation_offdiag_max_abs"],
            "diagonalized_real_masses": takagi["diagonalized_real_masses"],
            "diagonalized_imag_max_abs": takagi["diagonalized_imag_max_abs"],
            "diagonalized_offdiag_max_abs": takagi["diagonalized_offdiag_max_abs"],
        },
        "charged_basis_matrix_real": np.real(weighted_cycle_matrix).tolist(),
        "charged_basis_matrix_imag": np.imag(weighted_cycle_matrix).tolist(),
        "shared_basis_matrix_real": np.real(shared_basis_matrix).tolist(),
        "shared_basis_matrix_imag": np.imag(shared_basis_matrix).tolist(),
        "u_nu_shared_real": np.real(u_nu_shared).tolist(),
        "u_nu_shared_imag": np.imag(u_nu_shared).tolist(),
        "pmns_matrix_real": np.real(recovered_pmns).tolist(),
        "pmns_matrix_imag": np.imag(recovered_pmns).tolist(),
        "pmns_observables": observables,
        "weighted_cycle_observables_match": observable_match,
        "transport_checks": transport_checks,
        "emitted_parameters": majorana_pair,
        "notes": [
            "This artifact does not use the old intrinsic/shared-basis PMNS surface; it transports the repaired weighted-cycle theorem lane itself into the closed same-label basis.",
            "Because the transport is a Majorana congruence by the already-closed U_e_left, the shared-basis matrix remains complex symmetric and is diagonalized by the transported U_nu_shared with the same positive masses.",
            "The physical PMNS matrix recovered from U_e_left^dagger * U_nu_shared matches the repaired weighted-cycle PMNS surface exactly, so Majorana readout can be promoted on this anchored branch.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the shared-basis representation of the repaired weighted-cycle branch.")
    parser.add_argument("--weighted-cycle", default=str(DEFAULT_WEIGHTED_CYCLE))
    parser.add_argument("--shared-charged-left", default=str(DEFAULT_SHARED_CHARGED_LEFT))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    weighted_cycle_path = Path(args.weighted_cycle)
    shared_charged_left_path = Path(args.shared_charged_left)
    payload = build_payload(
        _load_json(weighted_cycle_path),
        _load_json(shared_charged_left_path),
        weighted_cycle_path=weighted_cycle_path,
        shared_charged_left_path=shared_charged_left_path,
    )

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
