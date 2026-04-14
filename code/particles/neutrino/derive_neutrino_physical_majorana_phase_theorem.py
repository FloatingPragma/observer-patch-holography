#!/usr/bin/env python3
"""Build the physical Majorana-phase theorem on the repaired weighted-cycle lane.

Chain role: promote the repaired weighted-cycle Majorana pair when an explicit
shared-basis physical PMNS representation is available, and otherwise keep the
stable Takagi pair on its honest candidate-only fallback surface.

Mathematics:
1. The repaired weighted-cycle artifact already emits one complex symmetric
   cycle matrix and one canonical Takagi PMNS surface.
2. A readout from PMNS columns alone is not honest because arbitrary
   computational column phases move `alpha_21` / `alpha_31` while leaving the
   oscillation observables unchanged.
3. The canonical congruence gauge is fixed by the Takagi condition
   `U_PMNS^T M_nu U_PMNS = diag(m_i)` with `diag(m_i)` real and non-negative.
4. In the electron-row gauge `U_e1 in R_{>0}`, the canonical pair is then read
   out as `alpha_21 = 2 arg(U_e2)` and
   `alpha_31 = 2 (arg(U_e3) + delta_PMNS)`.
5. If the repaired weighted-cycle branch is transported explicitly into the
   closed shared same-label basis, the physical PMNS path becomes
   `U_PMNS = U_e_left^dagger * U_nu_shared` on an anchored branch and the same
   canonical Majorana pair is promoted publicly there.
6. If that anchored representation is absent, the weighted-cycle Takagi readout
   remains only a candidate.

OPH-derived inputs: the repaired weighted-cycle neutrino branch, the closed
shared charged-lepton basis contract, and optionally an explicit shared-basis
representation of the repaired branch.

Output: one public physical Majorana phase pair when the shared-basis weighted-
cycle representation closes, or otherwise one candidate-only pair plus the
exact blocker for public promotion.
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
DEFAULT_SHARED_BASIS_REPRESENTATION = (
    ROOT / "particles" / "runs" / "neutrino" / "neutrino_weighted_cycle_shared_basis_representation.json"
)
DEFAULT_OUT = ROOT / "particles" / "runs" / "neutrino" / "neutrino_physical_majorana_phase_theorem.json"


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


def _canonical_pmns_from_weighted_cycle_matrix(matrix: np.ndarray) -> dict[str, Any]:
    if np.max(np.abs(matrix - matrix.T)) > 1.0e-12:
        raise ValueError("weighted-cycle matrix must remain complex symmetric for Majorana readout")

    left_vectors, singular_values, vh = np.linalg.svd(matrix)
    order = np.argsort(singular_values)
    left_vectors = left_vectors[:, order]
    singular_values = singular_values[order]
    right_vectors = vh.conj().T[:, order]

    phase_relation = left_vectors.T @ right_vectors
    offdiag_relation = phase_relation - np.diag(np.diag(phase_relation))
    relation_offdiag_max = float(np.max(np.abs(offdiag_relation)))
    if relation_offdiag_max > 1.0e-8:
        raise ValueError("SVD phase relation is not diagonal enough to define a canonical Takagi gauge")

    congruence_half_angles = np.angle(np.diag(phase_relation))
    takagi_vectors = left_vectors @ np.diag(np.exp(-0.5j * congruence_half_angles))
    pmns = np.conjugate(takagi_vectors)
    diagonalized = pmns.T @ matrix @ pmns

    diagonalized_offdiag = diagonalized - np.diag(np.diag(diagonalized))
    diag_real = np.real(np.diag(diagonalized))
    if np.any(diag_real <= 0.0):
        raise ValueError("canonical Takagi readout must produce a positive diagonal mass row")

    return {
        "pmns": pmns,
        "singular_values": [float(x) for x in singular_values.tolist()],
        "congruence_half_angles_rad": [float(0.5 * angle) for angle in congruence_half_angles.tolist()],
        "congruence_half_angles_deg": [math.degrees(float(0.5 * angle)) for angle in congruence_half_angles.tolist()],
        "phase_relation_offdiag_max_abs": relation_offdiag_max,
        "diagonalized_real_masses": [float(x) for x in diag_real.tolist()],
        "diagonalized_imag_max_abs": float(np.max(np.abs(np.imag(np.diag(diagonalized))))),
        "diagonalized_offdiag_max_abs": float(np.max(np.abs(diagonalized_offdiag))),
    }


def _observable_match(observables: dict[str, float], weighted_cycle_observables: dict[str, Any]) -> dict[str, float]:
    return {
        "theta12_deg_abs_delta": abs(observables["theta12_deg"] - float(weighted_cycle_observables["theta12_deg"])),
        "theta23_deg_abs_delta": abs(observables["theta23_deg"] - float(weighted_cycle_observables["theta23_deg"])),
        "theta13_deg_abs_delta": abs(observables["theta13_deg"] - float(weighted_cycle_observables["theta13_deg"])),
        "delta_deg_abs_delta": abs(((observables["delta_deg"] - float(weighted_cycle_observables["delta_deg"]) + 180.0) % 360.0) - 180.0),
        "J_abs_delta": abs(observables["J"] - float(weighted_cycle_observables["J"])),
    }


def build_payload(
    weighted_cycle: dict[str, Any],
    shared_charged_left: dict[str, Any],
    shared_basis_representation: dict[str, Any] | None = None,
    *,
    weighted_cycle_path: Path,
    shared_charged_left_path: Path,
    shared_basis_representation_path: Path | None = None,
) -> dict[str, Any]:
    basis_contract = dict(shared_charged_left.get("basis_contract", {}))
    if shared_charged_left.get("status") != "closed":
        raise ValueError("shared charged-lepton left basis must be closed before Majorana phases can be emitted")
    if not basis_contract.get("orientation_preserved", False):
        raise ValueError("shared charged-lepton basis must preserve orientation")
    if weighted_cycle.get("pmns_row_order_for_pdg") != ["e", "mu", "tau"]:
        raise ValueError("weighted-cycle branch must already be in PDG row order")

    matrix = _complex_matrix(weighted_cycle, "repaired_cycle_matrix_real", "repaired_cycle_matrix_imag")
    canonical = _canonical_pmns_from_weighted_cycle_matrix(matrix)
    pmns = canonical["pmns"]
    observables = _standard_pmns_parameters(pmns)
    majorana_pair = _majorana_pair_from_pmns(pmns, observables["delta_rad"])
    weighted_cycle_observables = dict(weighted_cycle["pmns_observables"])
    observable_match = _observable_match(observables, weighted_cycle_observables)
    if max(observable_match.values()) > 1.0e-8:
        raise ValueError("canonical Takagi PMNS readout must agree with the emitted weighted-cycle oscillation observables")

    public_promotion_status = "blocked_missing_shared_basis_representation"
    public_promotion_blocker = (
        "The canonical Takagi readout is stable on the repaired weighted-cycle matrix, but public promotion remains "
        "open until that repaired branch is represented explicitly on the closed shared same-label basis through "
        "U_e_left."
    )
    shared_basis_use_status = "weighted_cycle_takagi_readout_not_yet_promoted_without_shared_basis_transport"
    shared_basis_representation_summary = None
    status = "candidate_only"
    proof_chain_role = "candidate_surface"
    public_surface_candidate_allowed = False
    theorem_object = "canonical_majorana_phase_pair_candidate"
    theorem_surface = "weighted_cycle_takagi_readout_candidate"
    statement = (
        "On the repaired weighted-cycle branch, canonical Takagi congruence fixes one stable Majorana-phase "
        "pair candidate, but the public physical promotion remains open until the shared-basis weighted-cycle "
        "representation is emitted."
    )
    emitted_parameters: dict[str, Any] | None = None
    takagi_congruence_payload = {
        "singular_values": canonical["singular_values"],
        "congruence_half_angles_rad": canonical["congruence_half_angles_rad"],
        "congruence_half_angles_deg": canonical["congruence_half_angles_deg"],
        "phase_relation_offdiag_max_abs": canonical["phase_relation_offdiag_max_abs"],
        "diagonalized_real_masses": canonical["diagonalized_real_masses"],
        "diagonalized_imag_max_abs": canonical["diagonalized_imag_max_abs"],
        "diagonalized_offdiag_max_abs": canonical["diagonalized_offdiag_max_abs"],
    }

    if shared_basis_representation is not None:
        if shared_basis_representation.get("artifact") != "oph_neutrino_weighted_cycle_shared_basis_representation":
            raise ValueError("unexpected shared-basis weighted-cycle artifact id")
        if not bool(shared_basis_representation.get("physical_branch_closed", False)):
            raise ValueError("shared-basis weighted-cycle representation must explicitly close the physical branch")
        transport_checks = dict(shared_basis_representation.get("transport_checks") or {})
        if transport_checks:
            if float(transport_checks.get("shared_basis_symmetry_max_abs", 0.0)) > 1.0e-8:
                raise ValueError("shared-basis weighted-cycle representation must remain complex symmetric")
            if float(transport_checks.get("shared_basis_diagonalized_offdiag_max_abs", 0.0)) > 1.0e-8:
                raise ValueError("shared-basis weighted-cycle representation must diagonalize exactly")
            if float(transport_checks.get("shared_basis_diagonalized_imag_max_abs", 0.0)) > 1.0e-8:
                raise ValueError("shared-basis weighted-cycle representation must keep a real positive Takagi diagonal")
            if float(transport_checks.get("pmns_recovery_max_abs", 0.0)) > 1.0e-8:
                raise ValueError("shared-basis weighted-cycle representation must recover the physical PMNS exactly")
        pmns = _complex_matrix(shared_basis_representation, "pmns_matrix_real", "pmns_matrix_imag")
        observables = _standard_pmns_parameters(pmns)
        majorana_pair = _majorana_pair_from_pmns(pmns, observables["delta_rad"])
        observable_match = _observable_match(observables, weighted_cycle_observables)
        if max(observable_match.values()) > 1.0e-8:
            raise ValueError("shared-basis weighted-cycle representation must recover the weighted-cycle oscillation observables")
        takagi_congruence_payload = dict(shared_basis_representation.get("weighted_cycle_takagi_congruence") or takagi_congruence_payload)
        public_promotion_status = "closed_on_weighted_cycle_shared_basis_representation"
        public_promotion_blocker = None
        shared_basis_use_status = "public_emission_closed_on_weighted_cycle_shared_basis_representation"
        shared_basis_representation_summary = {
            "artifact": shared_basis_representation.get("artifact"),
            "status": shared_basis_representation.get("status"),
            "physical_branch_closed": bool(shared_basis_representation.get("physical_branch_closed", False)),
            "statement": shared_basis_representation.get("statement"),
            "transport_checks": transport_checks,
        }
        if shared_basis_representation_path is not None:
            shared_basis_representation_summary["source_path"] = str(shared_basis_representation_path)
        status = "theorem_grade_emitted"
        proof_chain_role = "active_theorem_lane"
        public_surface_candidate_allowed = True
        theorem_object = "physical_majorana_phase_pair"
        theorem_surface = "weighted_cycle_shared_basis_transport"
        statement = (
            "On the repaired weighted-cycle branch, the explicit shared-basis weighted-cycle representation closes the "
            "physical PMNS path through U_e_left^dagger * U_nu_shared, so canonical Takagi congruence emits one public "
            "Majorana phase pair on that weighted-cycle transport branch."
        )
        emitted_parameters = dict(majorana_pair)

    return {
        "artifact": "oph_neutrino_physical_majorana_phase_theorem",
        "generated_utc": _timestamp(),
        "status": status,
        "proof_chain_role": proof_chain_role,
        "public_surface_candidate_allowed": public_surface_candidate_allowed,
        "theorem_object": theorem_object,
        "theorem_surface": theorem_surface,
        "statement": statement,
        "readout_convention": {
            "takagi_condition": "U_PMNS^T M_nu U_PMNS = diag(m_i) with diag(m_i) in R_{>0}",
            "row_gauge": "U_e1 in R_{>0}",
            "phase_parameterization": "U_PMNS = V(theta12, theta23, theta13, delta_PMNS) * diag(1, exp(i alpha21 / 2), exp(i alpha31 / 2))",
        },
        "public_promotion_status": public_promotion_status,
        "public_promotion_blocker": public_promotion_blocker,
        "shared_basis_use_status": shared_basis_use_status,
        "source_artifacts": {
            "weighted_cycle_branch": str(weighted_cycle_path),
            "shared_charged_left_basis": str(shared_charged_left_path),
            "shared_basis_representation": str(shared_basis_representation_path) if shared_basis_representation_path is not None else None,
        },
        "depends_on": [
            dependency
            for dependency in [
                "oph_neutrino_weighted_cycle_repair",
                "oph_shared_charged_lepton_left_basis",
                "oph_neutrino_weighted_cycle_shared_basis_representation" if shared_basis_representation is not None else None,
            ]
            if dependency is not None
        ],
        "basis_contract": basis_contract,
        "takagi_congruence": takagi_congruence_payload,
        "pmns_matrix_real": np.real(pmns).tolist(),
        "pmns_matrix_imag": np.imag(pmns).tolist(),
        "pmns_observables": observables,
        "weighted_cycle_observables_match": observable_match,
        "shared_basis_representation": shared_basis_representation_summary,
        "candidate_parameters": majorana_pair,
        "emitted_parameters": emitted_parameters,
        "notes": [
            "A strict readout from PMNS columns alone would remain computational-gauge dependent, because arbitrary intermediate column phases leave the oscillation observables unchanged while moving alpha21 and alpha31.",
            "The readout fixes that ambiguity by reading the phases only after the canonical Takagi congruence of the emitted symmetric weighted-cycle matrix.",
            (
                "The resulting pair is promoted publicly because the repaired weighted-cycle branch is now represented explicitly on the closed shared basis and the physical PMNS path is recovered there exactly. This promotion does not identify that branch with the separate intrinsic/shared-basis PMNS diagnostic surface."
                if shared_basis_representation is not None
                else "The resulting pair is stable under computational column rephasings of the intermediate unitary, but it is not promoted here as a public physical emission because the repaired weighted-cycle branch has not yet been represented explicitly on the closed shared basis."
            ),
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the Majorana-phase theorem artifact.")
    parser.add_argument("--weighted-cycle", default=str(DEFAULT_WEIGHTED_CYCLE))
    parser.add_argument("--shared-charged-left", default=str(DEFAULT_SHARED_CHARGED_LEFT))
    parser.add_argument("--shared-basis-representation", default=str(DEFAULT_SHARED_BASIS_REPRESENTATION))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    weighted_cycle_path = Path(args.weighted_cycle)
    shared_charged_left_path = Path(args.shared_charged_left)
    shared_basis_representation_path = Path(args.shared_basis_representation)
    payload = build_payload(
        _load_json(weighted_cycle_path),
        _load_json(shared_charged_left_path),
        _load_json(shared_basis_representation_path) if shared_basis_representation_path.exists() else None,
        weighted_cycle_path=weighted_cycle_path,
        shared_charged_left_path=shared_charged_left_path,
        shared_basis_representation_path=shared_basis_representation_path if shared_basis_representation_path.exists() else None,
    )

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
