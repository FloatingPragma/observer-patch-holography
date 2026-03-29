#!/usr/bin/env python3
"""Build the real symmetric neutrino family-response tensor.

Chain role: combine the local neutrino scale anchor with the shared sector
response object to produce the family tensor behind the Majorana lane.

Mathematics: symmetric normalization, projector overlaps, and gap/isotropy
certificates on the neutrino response matrix.

OPH-derived inputs: the local D10-based scale anchor and the shared flavor
sector transport pushforward.

Output: the family-response tensor consumed by the Majorana lift, pullback
metric, and forward matrix builders.
"""

from __future__ import annotations

import argparse
import json
import pathlib
from datetime import datetime, timezone
from typing import Any

import numpy as np


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_SECTOR = ROOT / "particles" / "runs" / "flavor" / "sector_transport_pushforward.json"
DEFAULT_SCALE_ANCHOR = ROOT / "particles" / "runs" / "neutrino" / "neutrino_scale_anchor.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "neutrino" / "family_response_tensor.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: pathlib.Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _gap_certificate(values: np.ndarray) -> dict[str, Any]:
    ordered = np.sort(np.asarray(values, dtype=float))
    if ordered.size < 2:
        return {"min_gap": None, "stable": False}
    diffs = np.diff(ordered)
    min_gap = float(np.min(np.abs(diffs)))
    return {"min_gap": min_gap, "stable": bool(min_gap > 1.0e-9)}


def _projector_gap_seed(matrix: np.ndarray) -> float:
    h_seed = np.asarray(matrix, dtype=float).T @ np.asarray(matrix, dtype=float)
    eigvals = np.sort(np.linalg.eigvalsh(h_seed))
    if eigvals.size < 2:
        return 0.0
    return float(np.min(np.abs(np.diff(eigvals))))


def _value_isotropy_certificate(values: list[float], tolerance: float = 1.0e-30) -> dict[str, Any]:
    spread = float(max(values) - min(values))
    return {"closed": bool(spread <= tolerance), "spread": spread, "tolerance": tolerance}


def main() -> int:
    ap = argparse.ArgumentParser(description="Derive the neutrino family-response tensor.")
    ap.add_argument("--sector", default=str(DEFAULT_SECTOR), help="Input sector-response JSON path")
    ap.add_argument("--scale-anchor", default=str(DEFAULT_SCALE_ANCHOR), help="Input neutrino scale-anchor JSON path")
    ap.add_argument("--out", default=str(DEFAULT_OUT), help="Output JSON path")
    args = ap.parse_args()

    sector_path = pathlib.Path(args.sector)
    if not sector_path.is_absolute():
        sector_path = ROOT / args.sector
    scale_anchor_path = pathlib.Path(args.scale_anchor)
    if not scale_anchor_path.is_absolute():
        scale_anchor_path = ROOT / args.scale_anchor
    out_path = pathlib.Path(args.out)
    if not out_path.is_absolute():
        out_path = ROOT / args.out
    out_path.parent.mkdir(parents=True, exist_ok=True)

    sector_payload = _load_json(sector_path)
    scale_anchor_payload = _load_json(scale_anchor_path)
    nu = dict(sector_payload.get("sector_response_object", {}).get("nu", {}))
    if not nu:
        raise ValueError("sector_response_object['nu'] is required")

    directed = np.asarray(nu.get("K_core_directed", nu.get("K_core")), dtype=float)
    k_sym = np.asarray(nu.get("K_core_majorana_sym", nu.get("K_core")), dtype=float)
    if directed.shape != (3, 3) or k_sym.shape != (3, 3):
        raise ValueError("neutrino sector kernels must be 3x3")

    diag_entries = np.sqrt(np.clip(np.diag(k_sym), 1.0e-15, None))
    inv_diag = np.diag([1.0 / value for value in diag_entries])
    e_nu = inv_diag @ k_sym @ inv_diag
    s_nu = -np.log(np.clip(e_nu, 1.0e-15, None))

    u_vector = np.asarray(scale_anchor_payload["collective_mode"]["u_vector"], dtype=float)
    u_projector = np.asarray(scale_anchor_payload["collective_mode"]["u_uT"], dtype=float)
    collective_overlap = float(u_vector.T @ e_nu @ u_vector)
    if abs(collective_overlap) <= 1.0e-15:
        collective_overlap = 1.0
    c_nu_hat_real = e_nu / collective_overlap
    delta_nu = c_nu_hat_real - u_projector
    eigvals = np.linalg.eigvalsh(c_nu_hat_real)
    edge_amplitudes = {
        "a12": float(e_nu[0, 1]),
        "a23": float(e_nu[1, 2]),
        "a31": float(e_nu[2, 0]),
    }
    edge_amplitude_values = [edge_amplitudes["a12"], edge_amplitudes["a23"], edge_amplitudes["a31"]]
    weighted_edge_norm_sq = float(
        (diag_entries[0] * diag_entries[1] * e_nu[0, 1]) ** 2
        + (diag_entries[1] * diag_entries[2] * e_nu[1, 2]) ** 2
        + (diag_entries[2] * diag_entries[0] * e_nu[2, 0]) ** 2
    )

    payload = {
        "artifact": "oph_neutrino_family_response",
        "generated_utc": _timestamp(),
        "status": "real_majorana_response",
        "labels": sector_payload.get("labels"),
        "source_artifacts": {
            "sector_response": str(sector_path),
            "scale_anchor": str(scale_anchor_path),
        },
        "u_vector": u_vector.tolist(),
        "u_projector": u_projector.tolist(),
        "u_normalization": collective_overlap,
        "collective_mode_overlap": collective_overlap,
        "K_core_directed": directed.tolist(),
        "K_core_majorana_sym": k_sym.tolist(),
        "majorana_normalization_diag": diag_entries.tolist(),
        "E_nu": e_nu.tolist(),
        "S_nu": s_nu.tolist(),
        "C_nu_hat_real": c_nu_hat_real.tolist(),
        "Delta_nu": delta_nu.tolist(),
        "edge_amplitudes": edge_amplitudes,
        "edge_amplitude_isotropy_certificate": _value_isotropy_certificate(edge_amplitude_values),
        "weighted_edge_norm_sq": weighted_edge_norm_sq,
        "family_symmetry_certificate": {
            "s3_equivariant_data": bool(_value_isotropy_certificate(edge_amplitude_values)["closed"]),
            "reason": "equal_edge_amplitudes" if _value_isotropy_certificate(edge_amplitude_values)["closed"] else "nonisotropic_edges",
        },
        "rank_lower_bound": int(np.sum(np.abs(eigvals) > 1.0e-12)),
        "gap_certificate": _gap_certificate(eigvals),
        "projector_gap_seed": _projector_gap_seed(c_nu_hat_real),
        "phase_status": "unresolved",
        "notes": [
            "This artifact builds the real symmetric Majorana response from the sector-response object.",
            "The complex Majorana phase lift stays separate until the congruence-gauge theorem closes.",
            "Unique-edge amplitudes and the weighted edge norm are exported for the residual-envelope theorem burden.",
        ],
    }
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(out_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
