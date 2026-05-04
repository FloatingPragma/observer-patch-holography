#!/usr/bin/env python3
"""Emit the constructive QCD spectral-measure contract for the Thomson endpoint."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = ROOT / "particles" / "runs" / "hadron" / "ward_projected_spectral_measure_contract.json"
SCHEMA = ROOT / "particles" / "hadron" / "ward_projected_spectral_measure.schema.json"


def _now_utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact() -> dict[str, Any]:
    return {
        "artifact": "oph_ward_projected_spectral_measure_contract",
        "generated_utc": _now_utc(),
        "status": "constructive_contract_emitted_not_production_data",
        "schema": str(SCHEMA.relative_to(ROOT)),
        "promotion_allowed": False,
        "constructive_next_artifact": "oph_qcd_ward_projected_hadronic_spectral_measure",
        "why_this_is_forward_progress": (
            "The stable-channel hadron backend is not enough for the Thomson endpoint, so the next "
            "implementation target is now a concrete production export schema for the required "
            "Ward-projected electromagnetic spectral measure."
        ),
        "minimum_payload": {
            "projection": {
                "lane": "U(1)_Q",
                "ward_projected": True,
                "zero_momentum_endpoint_compatible": True,
            },
            "spectral_requirements": [
                "finite_volume_levels",
                "ward_projected_residues",
                "current_normalization",
                "rho_had_or_primitive_measure",
                "pushforward_rule_to_rho_had(s;P)",
            ],
            "required_budgets": [
                "statistical_budget",
                "continuum_budget",
                "finite_volume_budget",
                "chiral_budget",
                "current_matching_budget",
                "quadrature_budget",
                "endpoint_remainder_budget",
            ],
        },
        "forbidden_promotions": [
            "stable_channel_only_backend_export",
            "surrogate_hadron_artifact",
            "free_quark_screened_ansatz",
            "compare_only_external_Thomson_endpoint",
        ],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Emit the Ward-projected spectral-measure contract.")
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    parser.add_argument("--print-json", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    payload = build_artifact()
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text, encoding="utf-8")
    if args.print_json:
        print(text, end="")
    else:
        print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
