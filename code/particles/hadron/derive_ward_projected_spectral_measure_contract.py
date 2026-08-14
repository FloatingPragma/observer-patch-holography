#!/usr/bin/env python3
"""Emit the scientific QCD spectral-measure contract for the Thomson endpoint."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = ROOT / "particles" / "runs" / "hadron" / "ward_projected_spectral_measure_contract.json"
SCHEMA = ROOT / "particles" / "hadron" / "ward_projected_spectral_measure.schema.json"


def build_artifact() -> dict[str, Any]:
    return {
        "artifact": "oph_ward_projected_spectral_measure_contract",
        "classification": "scientific_contract_without_production_data",
        "schema": SCHEMA.relative_to(ROOT).as_posix(),
        "production_data_supplied": False,
        "promotion_allowed": False,
        "production_boundary": {
            "required_artifact": "oph_qcd_ward_projected_hadronic_spectral_measure",
            "requires_working_oph_hadron_backend": True,
            "backend_class": "OPH hardware backend such as GLORB/Echosahedron",
            "local_surrogate_promotable": False,
            "promotion_rule": (
                "Promotion requires a schema-valid production payload emitted by a working OPH "
                "hadron backend with every declared systematic budget supplied."
            ),
            "no_go_without_production_payload": True,
        },
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
        "local_real_engine": {
            "package": "particles/hadron/lattice_backend/",
            "runner": "particles/hadron/run_lattice_diagnostic_backend.py",
            "evidence_artifact": "particles/runs/hadron/lattice_engine_lane_status.json",
            "classification": "real_lattice_diagnostic_toy_scale",
            "satisfies_production_contract": False,
            "role": (
                "Physics-true lattice engine executed locally at toy scale: real gauge "
                "ensembles, real Dirac solves, real contractions, and no target anchoring. "
                "Its toy-scale evidence does not supply the seeded 2+1 production family or "
                "the full systematic budgets required by this contract."
            ),
        },
        "empirical_companion": {
            "artifact": "oph_empirical_ward_projected_hadronic_spectral_measure",
            "schema": "particles/hadron/empirical_ward_projected_spectral_measure.schema.json",
            "emitted_payload": "particles/runs/hadron/empirical_ward_projected_spectral_measure.json",
            "builder": "particles/hadron/derive_empirical_ward_projected_spectral_measure.py",
            "classification": "oph_plus_empirical_hadron_closure",
            "satisfies_production_contract": False,
            "role": (
                "Declared-empirical hadronic spectral input for the empirical closure surface. "
                "It carries the e+e- R(s) compilation as an explicit positive spectral measure "
                "with both endpoint kernels and a requadrature consistency gate, so the empirical "
                "Thomson endpoint consumes a spectral object rather than a bare integral. It is "
                "not an OPH-source production payload."
            ),
        },
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
