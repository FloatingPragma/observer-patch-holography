#!/usr/bin/env python3
"""Build PMNS only if a charged-lepton left basis is supplied explicitly."""

from __future__ import annotations

import argparse
import json
import pathlib
from typing import Any

import numpy as np


ROOT = pathlib.Path(__file__).resolve().parents[2]


def load_json(path: pathlib.Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_complex_matrix(real_rows: list[list[float]], imag_rows: list[list[float]]) -> np.ndarray:
    return np.array(real_rows, dtype=float) + 1j * np.array(imag_rows, dtype=float)


def main() -> int:
    ap = argparse.ArgumentParser(description="Build PMNS from a shared charged-lepton basis if available.")
    ap.add_argument("--majorana", default="particles/runs/neutrino/forward_majorana_matrix.json")
    ap.add_argument("--charged-left", default="", help="JSON artifact containing charged-lepton left diagonalizer")
    ap.add_argument("--out", default="particles/runs/neutrino/pmns_from_shared_basis.json")
    args = ap.parse_args()

    majorana_path = pathlib.Path(args.majorana)
    if not majorana_path.is_absolute():
        majorana_path = ROOT / args.majorana
    out_path = pathlib.Path(args.out)
    if not out_path.is_absolute():
        out_path = ROOT / args.out
    out_path.parent.mkdir(parents=True, exist_ok=True)

    majorana = load_json(majorana_path)
    u_nu = load_complex_matrix(majorana["U_nu_real"], majorana["U_nu_imag"])
    payload: dict[str, Any] = {
        "status": "conditional_blocked",
        "majorana_artifact": str(majorana_path),
        "notes": [
            "Physical PMNS requires a forward-derived charged-lepton left basis.",
        ],
    }
    if args.charged_left:
        charged_path = pathlib.Path(args.charged_left)
        if not charged_path.is_absolute():
            charged_path = ROOT / args.charged_left
        charged = load_json(charged_path)
        basis_contract = dict(charged.get("basis_contract", {}))
        if charged.get("labels") != majorana.get("labels") or not basis_contract.get("orientation_preserved", False):
            payload = {
                "status": "blocked_basis_mismatch",
                "majorana_artifact": str(majorana_path),
                "charged_left_artifact": str(charged_path),
                "notes": [
                    "Refusing PMNS because the charged-lepton basis contract is missing or mismatched.",
                ],
            }
        else:
            u_e = load_complex_matrix(charged["U_e_left"]["real"], charged["U_e_left"]["imag"])
            pmns = np.conjugate(u_e).T @ u_nu
            payload = {
                "status": "conditional_pmns",
                "majorana_artifact": str(majorana_path),
                "charged_left_artifact": str(charged_path),
                "pmns_real": np.real(pmns).tolist(),
                "pmns_imag": np.imag(pmns).tolist(),
                "notes": [
                    "This artifact is conditional on the charged-lepton left basis being forward-derived.",
                ],
            }
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(out_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
