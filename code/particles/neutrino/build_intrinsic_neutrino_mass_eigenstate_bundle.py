#!/usr/bin/env python3
"""Build intrinsic neutrino mass eigenstates from a scalar certificate.

Chain role: consume the proof-facing same-label scalar certificate and export
the exact intrinsic neutrino mass-eigenstate bundle.

Mathematics: reuse the exact principal-branch selector and depressed-cubic
spectral law from the centered eta-class builder.

Declared inputs: the isotropic neutrino forward bundle plus a numerically
complete same-label scalar certificate. Physical promotion additionally
requires that certificate's upstream source graph to close.

Output: the intrinsic neutrino mass-eigenstate bundle, kept separate from PMNS
and flavor-labeled rows until the shared charged basis exists.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import numpy as np
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from particles.artifact_paths import portable_json_dumps

EXACT_ETA_MAP_SCRIPT = ROOT / "particles" / "neutrino" / "derive_intrinsic_neutrino_exact_eta_map.py"
DEFAULT_ISOTROPIC = ROOT / "particles" / "runs" / "neutrino" / "forward_majorana_matrix.json"
DEFAULT_CERTIFICATE = ROOT / "particles" / "runs" / "neutrino" / "same_label_scalar_certificate.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "neutrino" / "intrinsic_neutrino_mass_eigenstate_bundle_from_scalar_certificate.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _load_exact_eta_module():
    spec = importlib.util.spec_from_file_location("oph_intrinsic_eta_map", EXACT_ETA_MAP_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError("could not load exact eta-map module")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def main() -> int:
    parser = argparse.ArgumentParser(description="Build intrinsic neutrino mass eigenstates from a scalar certificate.")
    parser.add_argument("--isotropic", default=str(DEFAULT_ISOTROPIC))
    parser.add_argument("--scalar-certificate", default=str(DEFAULT_CERTIFICATE))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    isotropic = _load_json(Path(args.isotropic))
    certificate = _load_json(Path(args.scalar_certificate))
    if certificate.get("sufficient_for_intrinsic_mass_eigenstates") is not True:
        raise ValueError("scalar certificate is incomplete; intrinsic mass eigenstates require a complete certificate")

    exact = _load_exact_eta_module()
    matrix_iso = np.array(isotropic["majorana_matrix_real"], dtype=float) + 1j * np.array(isotropic["majorana_matrix_imag"], dtype=float)
    a_value = float(matrix_iso[0, 0].real)
    rho_value = float(abs(matrix_iso[0, 1]))
    phase = exact._phase_vector_from_matrix(matrix_iso)
    omega = float(3.0 * phase[0])
    eta = exact._centered_eta_from_payload(certificate)
    mu_from_eta = np.exp(eta)
    mu_normalized = mu_from_eta / float(np.mean(mu_from_eta))
    exact_map = exact._build_exact_eta_map(a_value, rho_value, omega, eta)

    masses = exact_map.masses
    masses_sq = exact_map.masses_squared
    ordering = "unresolved_without_mass_eigenstate_label_rule"
    source_closure_status = dict(certificate.get("source_closure_status") or {})
    source_only_physical_input_eligible = (
        certificate.get("source_only_physical_input_eligible") is True
        and source_closure_status.get("closed") is True
    )
    rows = [
        {"state": "s0", "ascending_index": 0, "mass_gev": float(masses[0]), "mass_sq_gev2": float(masses_sq[0])},
        {"state": "s1", "ascending_index": 1, "mass_gev": float(masses[1]), "mass_sq_gev2": float(masses_sq[1])},
        {"state": "s2", "ascending_index": 2, "mass_gev": float(masses[2]), "mass_sq_gev2": float(masses_sq[2])},
    ]

    payload = {
        "artifact": "oph_intrinsic_neutrino_mass_eigenstate_bundle",
        "generated_utc": _timestamp(),
        "source_isotropic": str(Path(args.isotropic)),
        "source_scalar_certificate": str(Path(args.scalar_certificate)),
        "completion_scope": "intrinsic_mass_eigenstates_only",
        "source_only_physical_input_eligible": source_only_physical_input_eligible,
        "source_closure_status": source_closure_status or {"closed": False},
        "not_completed_items": [
            "source_closed_same_label_scalar_certificate"
            if not source_only_physical_input_eligible
            else None,
            "shared_charged_lepton_left_basis_for_pmns",
            "source_derived_mass_eigenstate_label_and_ordering_rule",
        ],
        "a_gev": float(a_value),
        "rho_gev": float(rho_value),
        "omega": omega,
        "eta_e": dict(certificate["eta_e"]),
        "mu_e_certificate": dict(certificate["mu_e"]),
        "mu_e_normalized": {
            edge: float(mu_normalized[idx]) for idx, edge in enumerate(exact.EDGE_ORDER)
        },
        "selector_lambda": float(exact_map.lam),
        "selector_point_absolute": {edge: float(exact_map.psi[idx]) for idx, edge in enumerate(exact.EDGE_ORDER)},
        "majorana_matrix_real": [[float(v.real) for v in row] for row in exact_map.majorana],
        "majorana_matrix_imag": [[float(v.imag) for v in row] for row in exact_map.majorana],
        "cubic_invariants": {
            "d": float(exact_map.d_trace_shift),
            "P": float(exact_map.p_invariant),
            "Q": float(exact_map.q_invariant),
            "cubic_roots": [float(x) for x in exact_map.cubic_roots],
            "equation": "lambda^3 - P lambda - 2Q = 0 for eigenvalues of H - d I",
            "spectral_crosscheck_max_abs_residual_gev2": float(np.max(np.abs((masses * masses) - masses_sq))),
        },
        "mass_eigenstates": rows,
        "mass_eigenstate_label_status": "ascending_singular_states_only",
        "physical_ordering_assignments": {
            "normal_ordering_hypothesis": {"nu1": "s0", "nu2": "s1", "nu3": "s2"},
            "inverted_ordering_hypothesis": {"nu3": "s0", "nu1": "s1", "nu2": "s2"},
            "selected": None,
            "missing_source_object": "solar_pair_and_atmospheric_sign_eigenstate_label_rule",
        },
        "masses_gev_sorted": [float(x) for x in masses],
        "masses_sq_gev2_sorted": [float(x) for x in masses_sq],
        "delta_m21_sq_gev2": float(masses_sq[1] - masses_sq[0]),
        "delta_m31_sq_gev2": float(masses_sq[2] - masses_sq[0]),
        "delta_m32_sq_gev2": float(masses_sq[2] - masses_sq[1]),
        "delta_mij_field_status": "ascending_gap_coordinates_under_declared_normal_ordering_hypothesis",
        "ascending_mass_sq_gaps_gev2": {
            "s1_minus_s0": float(masses_sq[1] - masses_sq[0]),
            "s2_minus_s0": float(masses_sq[2] - masses_sq[0]),
            "s2_minus_s1": float(masses_sq[2] - masses_sq[1]),
        },
        "ordering": ordering,
        "u_nu_real": [[float(v.real) for v in row] for row in exact_map.u_takagi],
        "u_nu_imag": [[float(v.imag) for v in row] for row in exact_map.u_takagi],
        "u_nu_semantics": "Takagi U satisfying U^T M U = diag(m_i) > 0",
        "row_basis_labels": ["f1", "f2", "f3"],
        "same_label_basis_contract": {
            "labels": ["f1", "f2", "f3"],
            "orientation_preserved": True,
        },
        "paper_export_policy": {
            "recommended_particle_rows": [],
            "diagnostic_singular_rows": ["s0", "s1", "s2"],
            "flavor_rows_status": "keep_gated_until_shared_charged_lepton_left_basis_closes",
            "pmns_status": "not_formed_here",
        },
        "notes": [
            "This bundle computes an exact ascending singular spectrum conditional on the same-label scalar certificate.",
            (
                "The certificate is source-closed."
                if source_only_physical_input_eligible
                else "The current certificate is numerically complete but inherits template/candidate source inputs."
            ),
            "Sorting does not select normal versus inverted physical mass labels; a source-side eigenstate label rule is still required.",
            "It does not by itself close PMNS or flavor-labeled neutrino rows.",
        ],
    }
    payload["not_completed_items"] = [item for item in payload["not_completed_items"] if item is not None]

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(portable_json_dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
