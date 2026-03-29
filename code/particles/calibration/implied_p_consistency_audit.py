#!/usr/bin/env python3
"""Audit which D10 observables imply which pixel constant under the current code."""

from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Dict, Optional


REPO_ROOT = Path(__file__).resolve().parents[2]


def _find_ancillary_particle_code_dir() -> Path:
    for base in Path(__file__).resolve().parents:
        candidate = base / "arXiv" / "RC1" / "ancillary" / "code" / "particles"
        if candidate.exists():
            return candidate
    raise FileNotFoundError("Could not locate arXiv/RC1/ancillary/code/particles")


PARTICLE_CODE_DIR = _find_ancillary_particle_code_dir()
for candidate in [PARTICLE_CODE_DIR / "core", PARTICLE_CODE_DIR]:
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

from particle_masses_paper_d10_d11 import (  # type: ignore
    PAPER_D10_TARGETS,
    P_DEFAULT,
    D10Closure,
    build_paper_d10,
)


ObservableFn = Callable[[D10Closure], float]


def _alpha_u_inv(d10: D10Closure) -> float:
    return 1.0 / d10.alpha_u


def _alpha_em_inv_mz(d10: D10Closure) -> float:
    return 1.0 / d10.alpha_em_mz


OBSERVABLES: Dict[str, ObservableFn] = {
    "alpha_u": lambda d10: d10.alpha_u,
    "alpha_u_inv": _alpha_u_inv,
    "mz_run": lambda d10: d10.mz_run,
    "v": lambda d10: d10.v,
    "alpha1_mz": lambda d10: d10.alpha1_mz,
    "alpha2_mz": lambda d10: d10.alpha2_mz,
    "alpha3_mz": lambda d10: d10.alpha3_mz,
    "alpha_em_inv_mz": _alpha_em_inv_mz,
    "sin2w_mz": lambda d10: d10.sin2w_mz,
    "m_z_pole_stage3": lambda d10: d10.m_z_pole_stage3,
    "m_w_run": lambda d10: d10.m_w_run,
}

DEFAULT_KEYS = (
    "alpha_em_inv_mz",
    "sin2w_mz",
    "mz_run",
    "m_z_pole_stage3",
    "m_w_run",
    "v",
)


def evaluate_observable(pix_area: float, key: str) -> float:
    return OBSERVABLES[key](build_paper_d10(pix_area=pix_area))


def estimate_implied_p_from_local_slope(
    *,
    p_center: float,
    current_value: float,
    target_value: float,
    derivative: float,
    p_min: float,
    p_max: float,
) -> Dict[str, object]:
    delta = current_value - target_value
    if derivative == 0.0:
        return {
            "status": "zero_slope",
            "implied_p": None,
            "delta_p_from_default": None,
        }
    implied_p = p_center - (delta / derivative)
    status = "estimated"
    if not (p_min <= implied_p <= p_max):
        status = "estimated_out_of_scan_range"
    return {
        "status": status,
        "implied_p": implied_p,
        "delta_p_from_default": implied_p - p_center,
    }


def bracket_root(
    key: str,
    target_value: float,
    p_min: float,
    p_max: float,
    grid_points: int,
) -> Optional[tuple[float, float]]:
    previous_p: Optional[float] = None
    previous_residual: Optional[float] = None
    for i in range(grid_points):
        current_p = p_min + (p_max - p_min) * float(i) / float(grid_points - 1)
        residual = evaluate_observable(current_p, key) - target_value
        if previous_residual is not None and residual == 0.0:
            return current_p, current_p
        if previous_residual is not None and residual * previous_residual < 0.0:
            return previous_p, current_p
        previous_p = current_p
        previous_residual = residual
    return None


def solve_implied_p(
    key: str,
    target_value: float,
    p_min: float,
    p_max: float,
    grid_points: int,
    iterations: int,
) -> Dict[str, object]:
    bracket = bracket_root(key, target_value, p_min, p_max, grid_points)
    if bracket is None:
        return {
            "status": "unbracketed",
            "implied_p": None,
            "iterations": 0,
        }

    lo, hi = bracket
    if lo == hi:
        return {
            "status": "exact_grid_hit",
            "implied_p": lo,
            "iterations": 0,
        }

    flo = evaluate_observable(lo, key) - target_value
    fhi = evaluate_observable(hi, key) - target_value
    for _ in range(iterations):
        mid = 0.5 * (lo + hi)
        fmid = evaluate_observable(mid, key) - target_value
        if flo * fmid <= 0.0:
            hi, fhi = mid, fmid
        else:
            lo, flo = mid, fmid

    implied_p = 0.5 * (lo + hi)
    final_residual = evaluate_observable(implied_p, key) - target_value
    return {
        "status": "solved",
        "implied_p": implied_p,
        "iterations": iterations,
        "final_residual": final_residual,
        "bracket": [bracket[0], bracket[1]],
    }


def derivative_at(pix_area: float, key: str, step: float) -> float:
    left = evaluate_observable(pix_area - step, key)
    right = evaluate_observable(pix_area + step, key)
    return (right - left) / (2.0 * step)


def build_audit(
    *,
    p_center: float,
    p_span: float,
    grid_points: int,
    iterations: int,
    derivative_step: float,
    keys: tuple[str, ...],
    refine: bool,
) -> Dict[str, object]:
    current_d10 = build_paper_d10(pix_area=p_center)
    d10_left = build_paper_d10(pix_area=p_center - derivative_step)
    d10_right = build_paper_d10(pix_area=p_center + derivative_step)
    p_min = p_center - p_span
    p_max = p_center + p_span
    observables: Dict[str, object] = {}
    implied_values = []

    for key in keys:
        target_value = PAPER_D10_TARGETS[key]
        current_value = OBSERVABLES[key](current_d10)
        derivative = (
            OBSERVABLES[key](d10_right) - OBSERVABLES[key](d10_left)
        ) / (2.0 * derivative_step)
        estimate = estimate_implied_p_from_local_slope(
            p_center=p_center,
            current_value=current_value,
            target_value=target_value,
            derivative=derivative,
            p_min=p_min,
            p_max=p_max,
        )
        record = {
            "target_value": target_value,
            "current_value_at_p_default": current_value,
            "current_delta": current_value - target_value,
            "dp_sensitivity_near_default": derivative,
            **estimate,
        }
        if refine and record["status"] in {"estimated", "estimated_out_of_scan_range"}:
            record["refined_root"] = solve_implied_p(key, target_value, p_min, p_max, grid_points, iterations)
        implied_p = record.get("implied_p")
        if isinstance(implied_p, float):
            implied_values.append(implied_p)
        observables[key] = record

    implied_spread = None
    implied_mid = None
    if implied_values:
        implied_spread = max(implied_values) - min(implied_values)
        implied_mid = 0.5 * (max(implied_values) + min(implied_values))

    focus_pair = None
    if (
        isinstance(observables["m_w_run"]["implied_p"], float)
        and isinstance(observables["m_z_pole_stage3"]["implied_p"], float)
    ):
        focus_pair = {
            "m_w_run_vs_m_z_pole_stage3": observables["m_z_pole_stage3"]["implied_p"]
            - observables["m_w_run"]["implied_p"]
        }

    return {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "source_branch": "arXiv/RC1/ancillary/code/particles/particle_masses_paper_d10_d11.py",
        "p_default": p_center,
        "p_scan_min": p_min,
        "p_scan_max": p_max,
        "grid_points": grid_points,
        "iterations": iterations,
        "d10_snapshot_at_p_default": asdict(current_d10),
        "root_refinement_enabled": refine,
        "observables": observables,
        "summary": {
            "observable_count": len(keys),
            "estimated_count": sum(1 for record in observables.values() if record["status"].startswith("estimated")),
            "implied_p_midpoint": implied_mid,
            "implied_p_spread": implied_spread,
            "exact_single_p_candidate": bool(implied_spread is not None and implied_spread < 1e-12),
            "focus_pair_spreads": focus_pair,
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=REPO_ROOT / "particles" / "runs" / "calibration" / "implied_p_consistency.json",
    )
    parser.add_argument("--p-center", type=float, default=P_DEFAULT)
    parser.add_argument("--p-span", type=float, default=5.0e-4)
    parser.add_argument("--grid-points", type=int, default=401)
    parser.add_argument("--iterations", type=int, default=80)
    parser.add_argument("--derivative-step", type=float, default=1.0e-5)
    parser.add_argument("--keys", nargs="*", default=list(DEFAULT_KEYS))
    parser.add_argument("--refine", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    keys = tuple(args.keys)
    unknown = [key for key in keys if key not in OBSERVABLES]
    if unknown:
        raise SystemExit(f"Unknown observable keys: {unknown}")

    audit = build_audit(
        p_center=args.p_center,
        p_span=args.p_span,
        grid_points=args.grid_points,
        iterations=args.iterations,
        derivative_step=args.derivative_step,
        keys=keys,
        refine=args.refine,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    summary = audit["summary"]
    print(f"wrote {args.output}")
    print(
        "implied_p_spread="
        f"{summary['implied_p_spread']:.12g}" if summary["implied_p_spread"] is not None else "implied_p_spread=unavailable"
    )
    focus_pair = summary["focus_pair_spreads"]
    if focus_pair:
        print(
            "m_w_run_vs_m_z_pole_stage3="
            f"{focus_pair['m_w_run_vs_m_z_pole_stage3']:.12g}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
