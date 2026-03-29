#!/usr/bin/env python3
"""Write the rho phase-shift / resonance-readout scaffold artifact."""

from __future__ import annotations

import argparse
import json
import pathlib
from typing import Any


def load_json(path: pathlib.Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> int:
    ap = argparse.ArgumentParser(description="Write a rho phase-shift-fit placeholder artifact.")
    ap.add_argument("--basis", default="particles/runs/hadron/rho_operator_basis.json")
    ap.add_argument("--levels", default="particles/runs/hadron/rho_levels.json")
    ap.add_argument("--stable", default="particles/runs/hadron/stable_channel_groundstate_readout.json")
    ap.add_argument("--out", default="particles/runs/hadron/rho_phase_shift_fit.json")
    args = ap.parse_args()

    basis_path = pathlib.Path(args.basis)
    if not basis_path.is_absolute():
        basis_path = pathlib.Path(__file__).resolve().parents[3] / basis_path
    levels_path = pathlib.Path(args.levels)
    if not levels_path.is_absolute():
        levels_path = pathlib.Path(__file__).resolve().parents[3] / levels_path
    stable_path = pathlib.Path(args.stable)
    if not stable_path.is_absolute():
        stable_path = pathlib.Path(__file__).resolve().parents[3] / stable_path
    out_path = pathlib.Path(args.out)
    if not out_path.is_absolute():
        out_path = pathlib.Path(__file__).resolve().parents[3] / out_path
    out_path.parent.mkdir(parents=True, exist_ok=True)

    basis = load_json(basis_path)
    levels = load_json(levels_path)
    stable = load_json(stable_path)
    payload = {
        "artifact": "oph_hadron_rho_scattering_readout",
        "status": "candidate_only",
        "basis_source": str(basis_path),
        "levels_source": str(levels_path),
        "stable_channel_source": str(stable_path),
        "basis_status": basis.get("status"),
        "levels_status": levels.get("status"),
        "stable_channel_status": stable.get("proof_status"),
        "channel": {
            "name": "rho",
            "elastic_channel": "pi_pi",
            "partial_wave": 1,
            "quantum_numbers": {"I": 1, "JPC": "1--"},
        },
        "moving_frames": basis.get("moving_frames", []),
        "irreps": basis.get("irreps", []),
        "operator_ids": basis.get("operator_ids", []),
        "pion_input": {
            "artifact": stable.get("artifact"),
            "channel": "pi_iso",
            "mass_source_rule": "stable_channel_groundstate_readout.pi_iso",
        },
        "phase_shift_model": {
            "kind": "elastic_p_wave_breit_wigner",
            "declared_chain": [
                "E_lab",
                "E_cm",
                "k",
                "q",
                "Luescher_phi",
                "delta1",
                "resonance_readout",
            ],
            "hidden_fit_choices_forbidden": True,
            "lab_to_cm_formula": "E_cm = sqrt(E_lab^2 - |P|^2)",
            "momentum_formula": "k(E_cm) = sqrt(E_cm^2 / 4 - m_pi^2)",
            "dimensionless_momentum_formula": "q = k L / (2*pi)",
            "quantization_condition": "delta1(E_cm) + phi_Luescher(q) = n*pi",
            "breit_wigner": {
                "tan_delta_formula": "E * Gamma(E) / (m_rho^2 - E^2)",
                "width_formula": "(g_rho_pipi^2 / (6*pi)) * k(E)^3 / E^2",
                "linearized_formula": "(k(E)^3 / E) * cot(delta1(E)) = (6*pi / g_rho_pipi^2) * (m_rho^2 - E^2)",
            },
        },
        "level_points": levels.get("level_points", []),
        "resonance_readout": {
            "readout_kind": "real_axis_pi_over_2_or_pole",
            "m_rho_gev": None,
            "Gamma_rho_gev": None,
            "pole_s_gev2": None,
            "real_axis_rule": "delta1(E_rho) = pi / 2",
        },
        "promotion_gate": "RhoElasticPhaseShiftReadoutConsistency",
        "anti_cheat": {
            "local_mass_ansatz_used": False,
            "measured_rho_mass_used_as_input": False,
            "local_rho_effective_mass_promoted": False,
        },
        "notes": [
            "Use this surface for Luescher / pole-fit metadata once finite-volume levels exist.",
            "The current local-rho effective mass is not the closure observable.",
            "This artifact consumes the operator basis, finite-volume levels, and the pion stable-channel input without promoting any local rho effective mass.",
        ],
    }
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(out_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
