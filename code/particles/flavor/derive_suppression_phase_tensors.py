#!/usr/bin/env python3
"""Lift a sector response object into suppression and phase tensors."""

from __future__ import annotations

import argparse
import json
import pathlib
from datetime import datetime, timezone
from typing import Any

import numpy as np


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "flavor" / "sector_transport_pushforward.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "suppression_phase_tensors.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


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


def build_artifact(payload: dict[str, Any]) -> dict[str, Any]:
    labels = payload.get("labels", ["f1", "f2", "f3"])
    if not isinstance(labels, list) or len(labels) != 3:
        raise ValueError("labels must be a length-3 list")

    response_object = dict(payload.get("sector_response_object", {}))
    if not response_object:
        raise ValueError("sector_response_object is required")

    artifact = {
        "artifact": "oph_suppression_phase_tensors",
        "generated_utc": _timestamp(),
        "status": str(payload.get("status", "lifted_from_sector_response")),
        "labels": [str(item) for item in labels],
        "metadata": {
            "observable_artifact": payload.get("artifact", "unknown"),
            "note": "Sector-response lift only. Dirac sectors get the canonical K3 lift; Majorana neutrino phase lift remains explicit if unresolved.",
            **dict(payload.get("metadata", {})),
        },
    }

    for sector, item in response_object.items():
        entry = dict(item)
        s_core = np.asarray(entry.get("S_core", np.zeros((3, 3))), dtype=float)
        k_core = np.asarray(entry.get("K_core", np.exp(-s_core)), dtype=float)
        omega = dict(entry.get("omega_cycle", {}))
        omega_012 = float(omega.get("012", 0.0))
        artifact[f"S_{sector}"] = s_core.tolist()
        artifact[f"Omega_{sector}"] = {"012": omega_012, "021": -omega_012}
        artifact[f"{sector}_symmetry_type"] = entry.get("symmetry_type")
        artifact[f"{sector}_normalization_class"] = entry.get("normalization_class")
        artifact[f"{sector}_residual_norms"] = entry.get("residual_norms")
        artifact[f"{sector}_raw_channel_norm_candidate"] = entry.get("raw_channel_norm_candidate")

        if entry.get("symmetry_type") == "majorana":
            artifact["S_nu_majorana_real"] = s_core.tolist()
            artifact["K_nu_majorana_real"] = k_core.tolist()
            artifact["phi_nu_unclosed"] = True
            artifact["majorana_lift_rule"] = "congruence_gauge_pending"
            artifact["majorana_phase_obstruction_class"] = {
                "kind": "normalized_edge_phase_affine_plane",
                "residual_dimension": 2,
                "constraint": "psi12+psi23+psi31=Omega_012",
            }
        else:
            artifact[f"Phi_{sector}"] = _dirac_phase_lift(omega_012).tolist()

    return artifact


def main() -> int:
    parser = argparse.ArgumentParser(description="Lift a sector response object into suppression and phase tensors.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT), help="Input sector-response JSON path.")
    parser.add_argument("--output", default=str(DEFAULT_OUT), help="Output tensor JSON path.")
    args = parser.parse_args()

    payload = json.loads(pathlib.Path(args.input).read_text(encoding="utf-8"))
    artifact = build_artifact(payload)

    out_path = pathlib.Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
