#!/usr/bin/env python3
"""Ensure PMNS construction refuses mismatched charged-lepton basis artifacts."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys
import tempfile


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "neutrino" / "build_pmns_from_shared_flavor_basis.py"


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="oph_pmns_basis_guard_") as tmpdir:
        tmp = pathlib.Path(tmpdir)
        majorana_path = tmp / "majorana.json"
        charged_path = tmp / "charged.json"
        out_path = tmp / "pmns.json"

        majorana_path.write_text(
            json.dumps(
                {
                    "labels": ["f1", "f2", "f3"],
                    "U_nu_real": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                    "U_nu_imag": [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        charged_path.write_text(
            json.dumps(
                {
                    "labels": ["x1", "x2", "x3"],
                    "basis_contract": {"orientation_preserved": False},
                    "U_e_left": {
                        "real": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                        "imag": [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
                    },
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )

        result = subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--majorana",
                str(majorana_path),
                "--charged-left",
                str(charged_path),
                "--out",
                str(out_path),
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        if result.returncode != 0:
            print(result.stderr.strip() or "pmns builder failed", file=sys.stderr)
            return 1
        payload = json.loads(out_path.read_text(encoding="utf-8"))
        if payload.get("status") != "blocked_basis_mismatch":
            print("basis mismatch did not block PMNS construction", file=sys.stderr)
            return 1

    print("shared flavor-basis contract blocks mismatched PMNS inputs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
