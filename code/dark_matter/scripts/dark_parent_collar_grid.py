#!/usr/bin/env python3
"""Build a parent-collar grid for the OPH dark-sector likelihood interface."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import dark_repair_transition_matrix


DEFAULT_MU_EQ = 5.363470440729118
DEFAULT_A_GRID = "1.0,0.5,0.1,0.000908265"
DEFAULT_K_GRID = "0.01,0.05,0.1,0.2,0.5,1.0"


def parse_csv_floats(value: str) -> list[float]:
    return [float(item.strip()) for item in value.split(",") if item.strip()]


def load_samples(path: str | None) -> dict[str, Any] | None:
    if path is None:
        return None
    with Path(path).open("r", encoding="utf-8") as handle:
        return json.load(handle)


def from_samples(samples: dict[str, Any]) -> dict[str, Any]:
    collar_samples = samples.get("collar_samples", [])
    if not collar_samples:
        raise ValueError("sample file contains no collar_samples")

    weighted_r = 0.0
    weighted_dr = 0.0
    weighted_k = 0.0
    has_direct_k = False
    total_weight = 0.0
    for sample in collar_samples:
        weight = float(sample.get("measure", 1.0))
        weighted_r += weight * float(sample["R_C"])
        weighted_dr += weight * float(sample.get("dR_dlogrho_b", 0.0))
        if "K_A_density_response" in sample:
            weighted_k += weight * float(sample["K_A_density_response"])
            has_direct_k = True
        total_weight += weight
    if total_weight <= 0.0:
        raise ValueError("sample measure must be positive")

    rho_ratio = weighted_r / total_weight
    mean_dr = weighted_dr / total_weight
    if has_direct_k:
        density_response = weighted_k / total_weight
    else:
        # If R_C stores rho_A/rho_b, then rho_A = rho_b R_C and
        # d rho_A / d rho_b = R_C + dR_C/dlogrho_b.
        density_response = rho_ratio + mean_dr
    contrast_response = density_response / rho_ratio if rho_ratio else 0.0
    a_grid = samples.get("a_grid", [1.0])
    k_grid = samples.get("k_grid_hMpc", [0.1])
    return {
        "status": {
            "category": "finite-collar parent functional evaluation",
            "source": "collar sample file",
            "paper_grade": bool(samples.get("paper_grade", False)),
            "notes": samples.get("notes", []),
        },
        "rho_A_over_rho_b": rho_ratio,
        "K_A_density_response_grid": [
            {"a": float(a), "k_hMpc": float(k), "K_A": density_response}
            for a in a_grid
            for k in k_grid
        ],
        "B_A_grid": [
            {"a": float(a), "k_hMpc": float(k), "B_A": contrast_response}
            for a in a_grid
            for k in k_grid
        ],
        "a_grid": [float(a) for a in a_grid],
        "k_grid_hMpc": [float(k) for k in k_grid],
        "sample_summary": {
            "n_samples": len(collar_samples),
            "total_measure": total_weight,
            "mean_R_C": weighted_r / total_weight,
            "mean_dR_dlogrho_b": weighted_dr / total_weight,
            "mean_K_A_density_response": density_response,
            "mean_B_A_contrast_response": contrast_response,
        },
    }


def diagnostic_scalar_load(args: argparse.Namespace) -> dict[str, Any]:
    a_grid = parse_csv_floats(args.a_grid)
    k_grid = parse_csv_floats(args.k_grid)
    repair_args = argparse.Namespace(mu_eq=args.mu_eq, n_max=args.n_max, hold=args.hold)
    repair_payload = dark_repair_transition_matrix.compute(repair_args)
    return {
        "status": {
            "category": "diagnostic scalar-load parent grid",
            "source": "finite scalar-load schema",
            "paper_grade": False,
            "notes": [
                "The default mu_eq is the flat capacity-saturated Omega_A/Omega_b diagnostic ratio.",
                "This diagnostic tests the likelihood plumbing. A finite-collar derivation of the abundance is not emitted by this code path.",
                "B_A is set to a scale-independent contrast response for interface testing.",
            ],
        },
        "rho_A_over_rho_b": args.mu_eq,
        "K_A_density_response_grid": [
            {"a": a, "k_hMpc": k, "K_A": args.mu_eq * args.B_A}
            for a in a_grid
            for k in k_grid
        ],
        "B_A_grid": [
            {"a": a, "k_hMpc": k, "B_A": args.B_A}
            for a in a_grid
            for k in k_grid
        ],
        "a_grid": a_grid,
        "k_grid_hMpc": k_grid,
        "repair_matrix": repair_payload,
    }


def compute(args: argparse.Namespace) -> dict[str, Any]:
    samples = load_samples(args.samples_json)
    if samples is not None:
        return from_samples(samples)
    if args.diagnostic_scalar_load:
        return diagnostic_scalar_load(args)
    return {
        "status": {
            "category": "parent-collar grid absent",
            "source": "no finite-collar sample file supplied",
            "paper_grade": False,
            "notes": [
                "Run with --samples-json for finite-collar samples, or with --diagnostic-scalar-load for a plumbing diagnostic."
            ],
        },
        "rho_A_over_rho_b": None,
        "B_A_grid": [],
        "a_grid": [],
        "k_grid_hMpc": [],
    }


def print_markdown(payload: dict[str, Any]) -> None:
    status = payload["status"]
    print("# Dark Parent Collar Grid")
    print()
    print(f"Category: `{status['category']}`")
    print(f"Paper grade: `{status['paper_grade']}`")
    rho_ratio = payload.get("rho_A_over_rho_b")
    if rho_ratio is None:
        print("rho_A/rho_b: `absent`")
    else:
        print(f"rho_A/rho_b: `{rho_ratio:.9f}`")
    print(f"B_A grid points: `{len(payload.get('B_A_grid', []))}`")
    print(f"K_A grid points: `{len(payload.get('K_A_density_response_grid', []))}`")
    if "repair_matrix" in payload:
        gap = payload["repair_matrix"]["gap"]
        print(f"gamma_rec: `{gap['gamma_rec']:.9f}`")
    print()
    print("## Notes")
    print()
    for note in status.get("notes", []):
        print(f"- {note}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples-json", default=None)
    parser.add_argument("--diagnostic-scalar-load", action="store_true")
    parser.add_argument("--mu-eq", type=float, default=DEFAULT_MU_EQ)
    parser.add_argument("--B-A", type=float, default=1.0)
    parser.add_argument("--n-max", type=int, default=40)
    parser.add_argument("--hold", type=float, default=0.25)
    parser.add_argument("--a-grid", default=DEFAULT_A_GRID)
    parser.add_argument("--k-grid", default=DEFAULT_K_GRID)
    parser.add_argument("--out", default=None)
    parser.add_argument("--json", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    payload = compute(args)
    if args.out:
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print_markdown(payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
