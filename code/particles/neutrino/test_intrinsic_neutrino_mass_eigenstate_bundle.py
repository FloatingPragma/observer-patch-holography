#!/usr/bin/env python3
"""Validate the intrinsic neutrino mass-eigenstate export from a scalar certificate."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys
import tempfile


ROOT = pathlib.Path(__file__).resolve().parents[2]
CERT_SCRIPT = ROOT / "particles" / "neutrino" / "derive_same_label_scalar_certificate.py"
BUNDLE_SCRIPT = ROOT / "particles" / "neutrino" / "build_intrinsic_neutrino_mass_eigenstate_bundle.py"
ISOTROPIC = ROOT / "particles" / "runs" / "neutrino" / "forward_neutrino_closure_bundle.json"


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="oph_intrinsic_mass_bundle_") as tmpdir:
        tmp = pathlib.Path(tmpdir)
        source = tmp / "pullback.json"
        cert = tmp / "certificate.json"
        out = tmp / "bundle.json"
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
        subprocess.run([sys.executable, str(CERT_SCRIPT), "--input", str(source), "--output", str(cert)], check=True, cwd=ROOT)
        subprocess.run(
            [sys.executable, str(BUNDLE_SCRIPT), "--isotropic", str(ISOTROPIC), "--scalar-certificate", str(cert), "--output", str(out)],
            check=True,
            cwd=ROOT,
        )
        payload = json.loads(out.read_text(encoding="utf-8"))
        if payload.get("artifact") != "oph_intrinsic_neutrino_mass_eigenstate_bundle":
            print("unexpected intrinsic mass-eigenstate bundle artifact", file=sys.stderr)
            return 1
        if len(payload.get("mass_eigenstates", [])) != 3:
            print("bundle should emit three intrinsic neutrino mass eigenstates", file=sys.stderr)
            return 1
        if payload.get("paper_export_policy", {}).get("pmns_status") != "not_formed_here":
            print("bundle should keep PMNS out of scope", file=sys.stderr)
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
