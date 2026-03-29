#!/usr/bin/env python3
"""Build the finite-volume level scaffold for the rho scattering program."""

from __future__ import annotations

import argparse
import json
import pathlib
from typing import Any


def load_json(path: pathlib.Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> int:
    ap = argparse.ArgumentParser(description="Write a rho level-extraction placeholder artifact.")
    ap.add_argument("--basis", default="particles/runs/hadron/rho_operator_basis.json")
    ap.add_argument("--out", default="particles/runs/hadron/rho_levels.json")
    args = ap.parse_args()

    basis_path = pathlib.Path(args.basis)
    if not basis_path.is_absolute():
        basis_path = pathlib.Path(__file__).resolve().parents[3] / basis_path
    out_path = pathlib.Path(args.out)
    if not out_path.is_absolute():
        out_path = pathlib.Path(__file__).resolve().parents[3] / out_path
    out_path.parent.mkdir(parents=True, exist_ok=True)

    basis = load_json(basis_path)
    payload = {
        "artifact": "oph_hadron_rho_scattering_levels",
        "status": "candidate_only",
        "basis_source": str(basis_path),
        "basis_status": basis.get("status"),
        "moving_frames": basis.get("moving_frames", []),
        "irreps": basis.get("irreps", []),
        "channel": {
            "name": "rho",
            "elastic_channel": "pi_pi",
            "partial_wave": 1,
            "quantum_numbers": {"I": 1, "JPC": "1--"},
        },
        "level_points": [],
        "level_point_schema": {
            "required_fields": ["frame_d", "irrep", "E_lab_gev"],
            "derived_fields": ["E_cm_gev", "k_gev", "q", "phi_luescher", "delta1_rad"],
        },
        "declared_chain": [
            "E_n(L,P,irrep)",
            "E_cm",
            "k2",
            "q2",
            "Luescher_phi",
            "delta1",
        ],
        "notes": [
            "No finite-volume spectrum has been extracted yet.",
            "Populate this artifact after generalized-eigenvalue / spectrum extraction work exists.",
            "The rho closure observable starts from pi-pi finite-volume levels, not a local rho effective mass.",
        ],
        "anti_cheat": {
            "local_rho_effective_mass_promoted": False,
            "measured_rho_mass_used_as_input": False,
        },
    }
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(out_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
