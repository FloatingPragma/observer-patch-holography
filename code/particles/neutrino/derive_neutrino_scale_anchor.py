#!/usr/bin/env python3
"""Derive the active neutrino scale anchor from the local D10 core.

Chain role: provide the neutrino mass scale used by the local Majorana branch
without importing any legacy neutrino estimate.

Mathematics: the scale law `m_star = v^2 / mu_u` plus the democratic rank-one
collective projector that seeds the neutrino family lane.

OPH-derived inputs: `v`, `mu_u`, `P`, and `alpha_u` from the active D10 core.

Output: the local neutrino scale anchor consumed by the family-response,
Majorana, and splitting builders.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_D10 = ROOT / "particles" / "runs" / "calibration" / "d10_ew_observable_family.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "neutrino" / "neutrino_scale_anchor.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact(d10: dict[str, object]) -> dict[str, object]:
    core = dict(d10["core_source"])
    v = float(core["v"])
    mu_u = float(core["mu_u"])
    m_star = (v * v) / mu_u
    collective_mode = [
        [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
        [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
        [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
    ]
    return {
        "artifact": "oph_neutrino_scale_anchor",
        "generated_utc": _timestamp(),
        "status": "local_d10_scale_anchor",
        "proof_status": "derived_from_local_d10_core",
        "source_artifacts": {
            "d10_observable_family": d10.get("artifact"),
        },
        "inputs": {
            "p": float(core["p"]),
            "mu_u_used_gev": mu_u,
            "v_used_gev": v,
            "alpha_u": float(core["alpha_u"]),
        },
        "anchors": {
            "m_star_gev": m_star,
            "m_star_formula": "v^2 / mu_u",
        },
        "collective_mode": {
            "basis": "democratic_rank_one",
            "u_vector": [0.5773502691896258, 0.5773502691896258, 0.5773502691896258],
            "u_uT": collective_mode,
        },
        "notes": [
            "This artifact is the active local neutrino scale anchor.",
            "It uses the current D10 core only: m_star = v^2 / mu_u.",
            "No legacy predictor surface or imported neutrino mass estimate is used here.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the local neutrino scale-anchor artifact.")
    parser.add_argument("--d10", default=str(DEFAULT_D10))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    d10 = json.loads(Path(args.d10).read_text(encoding="utf-8"))
    artifact = build_artifact(d10)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
