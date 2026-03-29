#!/usr/bin/env python3
"""Verify the exported pullback action matches Hilbert-Schmidt distortion samples."""

from __future__ import annotations

import argparse
import json
import math
import pathlib
import sys

import numpy as np


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_SCALE_ANCHOR = ROOT / "particles" / "runs" / "neutrino" / "neutrino_scale_anchor.json"
DEFAULT_FAMILY = ROOT / "particles" / "runs" / "neutrino" / "family_response_tensor.json"
DEFAULT_LIFT = ROOT / "particles" / "runs" / "neutrino" / "majorana_holonomy_lift.json"
DEFAULT_PULLBACK = ROOT / "particles" / "runs" / "neutrino" / "majorana_phase_pullback_metric.json"


def _phase_matrix(psi12: float, psi23: float, psi31: float) -> np.ndarray:
    out = np.ones((3, 3), dtype=complex)
    out[0, 1] = np.exp(1j * psi12)
    out[1, 0] = out[0, 1]
    out[1, 2] = np.exp(1j * psi23)
    out[2, 1] = out[1, 2]
    out[2, 0] = np.exp(1j * psi31)
    out[0, 2] = out[2, 0]
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description="Check the pullback action against HS distortion.")
    parser.add_argument("--scale-anchor", default=str(DEFAULT_SCALE_ANCHOR))
    parser.add_argument("--family", default=str(DEFAULT_FAMILY))
    parser.add_argument("--lift", default=str(DEFAULT_LIFT))
    parser.add_argument("--pullback", default=str(DEFAULT_PULLBACK))
    args = parser.parse_args()

    scale_anchor = json.loads(pathlib.Path(args.scale_anchor).read_text(encoding="utf-8"))
    family = json.loads(pathlib.Path(args.family).read_text(encoding="utf-8"))
    lift = json.loads(pathlib.Path(args.lift).read_text(encoding="utf-8"))
    pullback = json.loads(pathlib.Path(args.pullback).read_text(encoding="utf-8"))
    if not bool(pullback.get("phase_action_closed", False)):
        print("pullback action not closed; skip HS distortion identity test")
        return 0

    m_star = float(scale_anchor["anchors"]["m_star_gev"])
    e_nu = np.asarray(family["E_nu"], dtype=float)
    n_diag = np.diag(np.asarray(family["majorana_normalization_diag"], dtype=float))
    weights = {key: float(value) for key, value in dict(lift["edge_weights_majorana"]).items()}
    omega = float(lift["cycle_constraint"]["omega_012"])

    samples = [
        (omega / 3.0, omega / 3.0, omega / 3.0),
        (omega / 2.0, omega / 4.0, omega / 4.0),
        (omega - 0.1, 0.05, 0.05),
    ]
    anchor = m_star * (n_diag @ (e_nu * _phase_matrix(0.0, 0.0, 0.0)) @ n_diag)
    for psi12, psi23, psi31 in samples:
        lifted = m_star * (n_diag @ (e_nu * _phase_matrix(psi12, psi23, psi31)) @ n_diag)
        lhs = float(np.linalg.norm(lifted - anchor, ord="fro") ** 2 / (4.0 * (m_star**2)))
        rhs = (
            weights["psi12"] * (1.0 - math.cos(psi12))
            + weights["psi23"] * (1.0 - math.cos(psi23))
            + weights["psi31"] * (1.0 - math.cos(psi31))
        )
        if abs(lhs - rhs) > 1.0e-15:
            print(f"HS distortion mismatch: lhs={lhs} rhs={rhs}", file=sys.stderr)
            return 1

    print("pullback action matches HS distortion samples")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
