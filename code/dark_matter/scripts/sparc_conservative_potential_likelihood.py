#!/usr/bin/env python3
"""Gate and contract for the conservative-potential SPARC likelihood."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "outputs" / "sparc_conservative_potential_likelihood.json"


def payload(args: argparse.Namespace) -> dict[str, Any]:
    return {
        "status": {
            "category": "SPARC conservative-potential likelihood gate",
            "ready": False,
            "paper_grade": False,
            "solver": "conservative_potential",
            "branch_selector_fixed": True,
            "reason": "axisymmetric disk-potential solver and hierarchical likelihood are not emitted by this code path",
        },
        "field_equation": {
            "primary": "nabla^2 Phi = nabla dot [nu_lambda(|nabla Phi_b|/a0_OPH) nabla Phi_b]",
            "acceleration": "g_obs = -nabla Phi",
            "algebraic_scope": "exact only in spherical, planar, or one-dimensional symmetry",
        },
        "required_nuisance_parameters": [
            "distance",
            "inclination",
            "disk mass-to-light ratio",
            "bulge mass-to-light ratio",
            "gas scale",
            "intrinsic scatter",
            "stellar-population hyperparameters",
            "velocity covariance",
            "grid convergence covariance",
        ],
        "required_outputs": [
            "global likelihood rows for OPH unit, OPH Z6/Poisson, empirical MOND, QUMOND, NFW, and Burkert",
            "per-galaxy posterior or profiled nuisance summaries",
            "grid convergence report",
            "quality-cut sensitivity report",
            "posterior predictive residuals",
        ],
        "compatibility_scripts": [
            "cosmology/scripts/sparc_rar_compare.py",
            "cosmology/scripts/sparc_rotation_curve_compare.py",
            "cosmology/scripts/sparc_profiled_ml_likelihood.py",
            "cosmology/scripts/sparc_systematic_likelihood.py",
        ],
        "acceptance_tests": [
            "spherical test mass reproduces the algebraic response",
            "disk mode refuses non-conservative vector response",
            "algebraic compatibility mode reproduces existing scaffold rows",
            "conservative mode emits separate publication-grade rows",
        ],
        "inputs": {
            "max_quality": args.max_quality,
            "min_inclination_deg": args.min_inclination_deg,
        },
    }


def print_markdown(data: dict[str, Any]) -> None:
    status = data["status"]
    print("# SPARC Conservative Potential Likelihood Gate")
    print()
    print(f"Ready: `{status['ready']}`")
    print(f"Paper grade: `{status['paper_grade']}`")
    print(f"Reason: `{status['reason']}`")
    print()
    print("## Primary Equation")
    print()
    print(f"`{data['field_equation']['primary']}`")
    print()
    print("## Required Nuisance Parameters")
    print()
    for item in data["required_nuisance_parameters"]:
        print(f"- {item}")
    print()
    print("## Acceptance Tests")
    print()
    for item in data["acceptance_tests"]:
        print(f"- {item}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-quality", type=int, default=2)
    parser.add_argument("--min-inclination-deg", type=float, default=30.0)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--json", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    data = payload(args)
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.json:
        print(json.dumps(data, indent=2, sort_keys=True))
    else:
        print_markdown(data)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
