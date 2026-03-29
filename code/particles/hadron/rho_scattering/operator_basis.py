#!/usr/bin/env python3
"""Emit the minimal rho-scattering operator-basis artifact."""

from __future__ import annotations

import argparse
import json
import pathlib


def main() -> int:
    ap = argparse.ArgumentParser(description="Write a rho operator-basis artifact.")
    ap.add_argument("--out", default="particles/runs/hadron/rho_operator_basis.json")
    args = ap.parse_args()

    out_path = pathlib.Path(args.out)
    if not out_path.is_absolute():
        out_path = pathlib.Path(__file__).resolve().parents[3] / out_path
    out_path.parent.mkdir(parents=True, exist_ok=True)

    payload = {
        "artifact": "oph_hadron_rho_operator_basis",
        "status": "candidate_only",
        "target": "rho resonance via pi-pi finite-volume spectra",
        "channel": {
            "name": "rho",
            "elastic_channel": "pi_pi",
            "partial_wave": 1,
            "quantum_numbers": {"I": 1, "JPC": "1--"},
        },
        "operator_ids": [
            "V_local",
            "pipi_d000_T1u",
            "pipi_d001_A1",
            "pipi_d110_A1",
        ],
        "operators": [
            "local vector interpolator",
            "two-pion operators in moving frames",
            "irrep-resolved operator sets",
        ],
        "moving_frames": [
            {"d": [0, 0, 0]},
            {"d": [0, 0, 1]},
            {"d": [1, 1, 0]},
        ],
        "irreps": [
            {"frame_d": [0, 0, 0], "irrep": "T1u"},
            {"frame_d": [0, 0, 1], "irrep": "A1"},
            {"frame_d": [1, 1, 0], "irrep": "A1"},
        ],
        "next_outputs": [
            "basis description",
            "finite-volume level extraction inputs",
            "Luescher-fit-ready metadata",
        ],
        "anti_cheat": {
            "local_mass_ansatz_used": False,
            "measured_rho_mass_used_as_input": False,
        },
    }
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(out_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
