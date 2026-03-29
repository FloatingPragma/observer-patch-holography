#!/usr/bin/env python3
"""Emit the typed charged absolute-scale coordinate shell.

Chain role: keep the charged absolute-scale lane explicit about the difference
between a linear mass scale `g_e` and its logarithmic coordinate `mu_e_abs`.

Mathematics: coordinate-shell declaration only; no value law is emitted here.

OPH-derived inputs: none beyond the declared charged-lane closure contract.

Output: a machine-readable shell stating that representation consistency is
closed while the actual absolute-scale value remains open.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = ROOT / "particles" / "runs" / "leptons" / "charged_lepton_absolute_scale_coordinate_shell.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact() -> dict[str, object]:
    return {
        "artifact": "oph_charged_lepton_absolute_scale_coordinate_shell",
        "generated_utc": _timestamp(),
        "g_e_linear_candidate": None,
        "mu_e_absolute_log_candidate": None,
        "conversion_formula": "g_e_linear_candidate = exp(mu_e_absolute_log_candidate)",
        "representation_consistency_closed": True,
        "proof_status": "coordinate_shell_only_values_open",
        "next_single_residual_object": "mu_e_absolute_log_candidate",
        "notes": [
            "The charged absolute-scale lane must not subtract log gaps directly from a linear scale.",
            "Any future charged absolute-scale emitter should produce either mu_e_absolute_log_candidate or g_e_linear_candidate and convert exactly once.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the charged absolute-scale coordinate shell.")
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    payload = build_artifact()
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
