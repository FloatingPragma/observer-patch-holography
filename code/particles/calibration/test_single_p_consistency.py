#!/usr/bin/env python3
"""Sanity checks for the local implied-P calibration audit."""

from __future__ import annotations

import math

from implied_p_consistency_audit import build_audit, estimate_implied_p_from_local_slope
from particle_masses_paper_d10_d11 import P_DEFAULT, build_paper_d10


def assert_close(label: str, left: float, right: float, tol: float = 1.0e-10) -> None:
    if not math.isclose(left, right, rel_tol=0.0, abs_tol=tol):
        raise AssertionError(f"{label}: expected {right}, got {left}")


def main() -> int:
    d10 = build_paper_d10(pix_area=P_DEFAULT)
    d10_left = build_paper_d10(pix_area=P_DEFAULT - 1.0e-5)
    d10_right = build_paper_d10(pix_area=P_DEFAULT + 1.0e-5)

    w_slope = (d10_right.m_w_run - d10_left.m_w_run) / (2.0e-5)
    w_estimate = estimate_implied_p_from_local_slope(
        p_center=P_DEFAULT,
        current_value=d10.m_w_run,
        target_value=d10.m_w_run,
        derivative=w_slope,
        p_min=P_DEFAULT - 5.0e-4,
        p_max=P_DEFAULT + 5.0e-4,
    )
    if w_estimate["status"] not in {"estimated", "estimated_out_of_scan_range"}:
        raise AssertionError(f"m_w_run estimate status was {w_estimate['status']}")
    assert_close("m_w_run roundtrip P", float(w_estimate["implied_p"]), P_DEFAULT)

    z_slope = (d10_right.m_z_pole_stage3 - d10_left.m_z_pole_stage3) / (2.0e-5)
    z_estimate = estimate_implied_p_from_local_slope(
        p_center=P_DEFAULT,
        current_value=d10.m_z_pole_stage3,
        target_value=d10.m_z_pole_stage3,
        derivative=z_slope,
        p_min=P_DEFAULT - 5.0e-4,
        p_max=P_DEFAULT + 5.0e-4,
    )
    if z_estimate["status"] not in {"estimated", "estimated_out_of_scan_range"}:
        raise AssertionError(f"m_z_pole_stage3 estimate status was {z_estimate['status']}")
    assert_close("m_z_pole_stage3 roundtrip P", float(z_estimate["implied_p"]), P_DEFAULT)

    audit = build_audit(
        p_center=P_DEFAULT,
        p_span=5.0e-4,
        grid_points=401,
        iterations=80,
        derivative_step=1.0e-5,
        keys=("alpha_em_inv_mz", "sin2w_mz", "m_z_pole_stage3", "m_w_run", "v"),
        refine=False,
    )
    if audit["summary"]["estimated_count"] < 5:
        raise AssertionError("expected all focus observables to estimate implied P")
    if "m_w_run" not in audit["observables"] or "m_z_pole_stage3" not in audit["observables"]:
        raise AssertionError("missing core electroweak observables from audit")

    print("calibration implied-P audit sanity checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
