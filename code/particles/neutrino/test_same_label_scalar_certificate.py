#!/usr/bin/env python3
"""Validate the same-label scalar certificate builder."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys
import tempfile


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "neutrino" / "derive_same_label_scalar_certificate.py"


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="oph_same_label_cert_") as tmpdir:
        tmp = pathlib.Path(tmpdir)
        source = tmp / "pullback.json"
        out = tmp / "certificate.json"
        source.write_text(
            json.dumps(
                {
                    "artifact": "oph_realized_same_label_gap_defect_readback",
                    "same_label": {"psi12": "f1", "psi23": "f2", "psi31": "f3"},
                    "rank_e": {"psi12": 1.0, "psi23": 1.0, "psi31": 1.0},
                    "same_label_overlap_sq": {"psi12": 0.81, "psi23": 0.64, "psi31": 0.49},
                    "defect_e": {"psi12": 0.19, "psi23": 0.36, "psi31": 0.51},
                    "gap_e": {"psi12": 2.0, "psi23": 3.0, "psi31": 5.0},
                    "base_mu_nu": 2.5e-12,
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        subprocess.run([sys.executable, str(SCRIPT), "--input", str(source), "--output", str(out)], check=True, cwd=ROOT)
        payload = json.loads(out.read_text(encoding="utf-8"))
        if payload.get("proof_status") != "fixed_cutoff_scalar_sufficient_downstream_certificate":
            print("certificate should close on complete same-label scalars", file=sys.stderr)
            return 1
        if payload.get("sufficient_for_intrinsic_mass_eigenstates") is not True:
            print("certificate should be sufficient for intrinsic mass eigenstates", file=sys.stderr)
            return 1
        if any(payload["q_e"][edge] is None for edge in ("psi12", "psi23", "psi31")):
            print("certificate should emit q_e on every edge", file=sys.stderr)
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
