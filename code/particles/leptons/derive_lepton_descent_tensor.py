#!/usr/bin/env python3
"""Derive a sandbox charged-lepton channel tensor from the sector response object."""

from __future__ import annotations

import argparse
import json
import pathlib
from datetime import datetime, timezone
from typing import Any

import numpy as np


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "flavor" / "sector_transport_pushforward.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "leptons" / "lepton_descent_tensor.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _encode_complex_matrix(matrix: np.ndarray) -> dict[str, Any]:
    return {
        "real": np.real(matrix).tolist(),
        "imag": np.imag(matrix).tolist(),
    }


def _dirac_phase_lift(omega: float) -> np.ndarray:
    scale = float(omega) / 3.0
    return np.asarray(
        [
            [0.0, scale, -scale],
            [-scale, 0.0, scale],
            [scale, -scale, 0.0],
        ],
        dtype=float,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Build a charged-lepton channel tensor artifact.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT), help="Input sector-pushforward JSON path.")
    parser.add_argument("--output", default=str(DEFAULT_OUT), help="Output JSON path.")
    args = parser.parse_args()

    payload = json.loads(pathlib.Path(args.input).read_text(encoding="utf-8"))
    sector = dict(payload.get("sector_response_object", {}).get("e", {}))
    if not sector:
        raise ValueError("sector_response_object['e'] is required")

    s_e = np.asarray(sector.get("S_core"), dtype=float)
    if s_e.shape != (3, 3):
        raise ValueError("sector_response_object['e'].S_core must be a 3x3 matrix")
    omega_012 = float(dict(sector.get("omega_cycle", {})).get("012", 0.0))
    phi_e = _dirac_phase_lift(omega_012)
    k_e = np.exp(-s_e)
    trace_norm = float(np.trace(k_e))
    if trace_norm <= 0.0:
        trace_norm = float(np.linalg.norm(k_e)) or 1.0
    c_e_hat_matrix = (k_e * np.exp(1j * phi_e)) / trace_norm

    artifact = {
        "artifact": "oph_lepton_channel_tensor",
        "generated_utc": _timestamp(),
        "labels": payload.get("labels"),
        "channel": sector.get("channel", "L-H-E"),
        "family_projectors": payload.get("family_projectors"),
        "spectral_gaps": payload.get("spectral_gaps"),
        "projector_basis_provenance": {
            "observable_artifact": payload.get("metadata", {}).get("observable_artifact"),
            "sector_response_artifact": payload.get("artifact"),
            "orientation_preserved": bool(sector.get("orientation_preserved", False)),
        },
        "sector_response_summary": {
            "symmetry_type": sector.get("symmetry_type"),
            "positive_kernel_kind": sector.get("positive_kernel_kind"),
            "normalization_class": sector.get("normalization_class"),
            "raw_channel_norm_candidate": sector.get("raw_channel_norm_candidate"),
            "raw_channel_norm_by_refinement": sector.get("raw_channel_norm_by_refinement", []),
            "scale_scope_candidate": sector.get("scale_scope_candidate"),
            "shared_budget_key": sector.get("shared_budget_key"),
            "functoriality_certificate": sector.get("functoriality_certificate", {}),
        },
        "C_e_hat": {
            "K_e_plus": k_e.tolist(),
            "S_e": s_e.tolist(),
            "Omega_e": {"012": omega_012, "021": -omega_012},
            "Phi_e": phi_e.tolist(),
            "C_e_hat_matrix": _encode_complex_matrix(c_e_hat_matrix),
            "trace_normalization": trace_norm,
        },
        "phase_normal_form_certificate": {
            "normal_form_available": False,
            "phase_can_be_demoted_for_singular_spectrum": False,
            "status": "open",
        },
        "shape_closed": False,
        "channel_norm_required": True,
        "absolute_scale_closed": False,
        "residual_burden": "lepton_channel_norm_g_e",
        "metadata": {
            "sector_response_artifact": payload.get("artifact", "unknown"),
            "note": "Sandbox charged-lepton channel tensor. Shape and absolute channel norm stay split until the descent and norm theorems close.",
            **dict(payload.get("metadata", {})),
        },
    }

    out_path = pathlib.Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
