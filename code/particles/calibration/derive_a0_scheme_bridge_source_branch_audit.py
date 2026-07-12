#!/usr/bin/env python3
"""Executed source-branch attempt at the a0 scheme bridge (#545).

The certified anchor gap [0.649, 0.855] names the shift the source side must
emit at a0(P) for the empirical-closure Thomson endpoint to meet the measured
value.  The only non-circular source route is to continue the declared
running convention to the next order and test whether the induced anchor
shift lands inside the certified interval.

This builder executes that test.  It reproduces the one-loop anchor from the
declared MSSM coefficients b = (33/5, 1, -3) and the certificate boundary
data, integrates the standard two-loop MSSM gauge system (with a top-Yukawa
convention scan, since no tan-beta or superpartner spectrum is declared), and
compares the induced shift against the certified gap.

Result encoded fail-closed: the two-loop continuation overshoots the gap by a
factor of two to three for every Yukawa convention, and the terms that could
restore the balance (threshold placements, scheme conversion) are exactly the
open scheme-lock and threshold-map interfaces of the RG matching contract,
which forbids using them as hidden fit parameters.  The source branch of
issue #545 therefore does not close under the currently declared conventions;
Route A (empirical class) remains the carrying lane.  The audit records the
quantitative balance requirement so the eventual scheme-lock object knows
what it must supply.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
REPO_ROOT = ROOT.parent
RUNS = ROOT / "particles" / "runs"
RUNTIME = ROOT / "P_derivation" / "runtime"

CERTIFICATE_JSON = RUNS / "calibration" / "d10_ew_forward_transmutation_certificate.json"
ENDPOINT_JSON = RUNTIME / "empirical_thomson_endpoint_current.json"
RG_CONTRACT_JSON = RUNTIME / "rg_matching_threshold_contract_current.json"
DEFAULT_OUT = RUNS / "calibration" / "a0_scheme_bridge_source_branch_audit.json"

B_ONE_LOOP = (33.0 / 5.0, 1.0, -3.0)
B_TWO_LOOP = (
    (199.0 / 25.0, 27.0 / 5.0, 88.0 / 5.0),
    (9.0 / 5.0, 25.0, 24.0),
    (11.0 / 5.0, 9.0, 14.0),
)
C_TOP = (26.0 / 5.0, 6.0, 4.0)
YUKAWA_SCAN = (0.0, 0.3, 0.5, 0.7, 1.0, 1.2)
MZ_GEV = 91.1876

# Known-physics scheme spread between the MS-bar-like and on-shell effective
# electromagnetic couplings at MZ (PDG: 127.93 vs 128.94).  Declared reference
# only; this audit does not derive it.
SCHEME_CONVERSION_MSBAR_TO_ON_SHELL = 1.01


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _artifact_ref(path: Path) -> str:
    try:
        return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()
    except ValueError:
        return str(path)


def integrate_two_loop(
    x_unified: float, efolds: float, yt_unified: float, two_loop: bool = True, steps: int = 4000
) -> tuple[list[float], float]:
    """RK4 descent of the MSSM gauge system from the unification boundary.

    State: x_i = alpha_i^-1 in GUT normalization, plus the top Yukawa.
    dx_i/dt = -(1/2pi)[ b_i + (1/4pi)(sum_j b_ij alpha_j - c_i yt^2/(4pi)*4pi) ]
    with the standard one-loop MSSM Yukawa running for yt.
    """

    x = [x_unified, x_unified, x_unified]
    yt = yt_unified
    dt = -efolds / steps

    def derivs(x: list[float], yt: float) -> tuple[list[float], float]:
        a = [1.0 / xi for xi in x]
        dx = []
        for i in range(3):
            two = sum(B_TWO_LOOP[i][j] * a[j] for j in range(3)) - C_TOP[i] * (
                yt * yt / (4.0 * math.pi)
            )
            dx.append(
                -(1.0 / (2.0 * math.pi))
                * (B_ONE_LOOP[i] + (two / (4.0 * math.pi)) * (1.0 if two_loop else 0.0))
            )
        g_sq = [4.0 * math.pi * ai for ai in a]
        dyt = (yt / (16.0 * math.pi**2)) * (
            6.0 * yt * yt - (16.0 / 3.0) * g_sq[2] - 3.0 * g_sq[1] - (13.0 / 15.0) * g_sq[0]
        )
        return dx, dyt

    for _ in range(steps):
        k1x, k1y = derivs(x, yt)
        s2 = [x[i] + 0.5 * dt * k1x[i] for i in range(3)]
        k2x, k2y = derivs(s2, yt + 0.5 * dt * k1y)
        s3 = [x[i] + 0.5 * dt * k2x[i] for i in range(3)]
        k3x, k3y = derivs(s3, yt + 0.5 * dt * k2y)
        s4 = [x[i] + dt * k3x[i] for i in range(3)]
        k4x, k4y = derivs(s4, yt + dt * k3y)
        x = [x[i] + dt / 6.0 * (k1x[i] + 2 * k2x[i] + 2 * k3x[i] + k4x[i]) for i in range(3)]
        yt = yt + dt / 6.0 * (k1y + 2 * k2y + 2 * k3y + k4y)
    return x, yt


def alpha_em_inv(x: list[float]) -> float:
    """alpha_em^-1 = alpha2^-1 + alphaY^-1 with xY = (5/3) x1."""

    return x[1] + (5.0 / 3.0) * x[0]


def _find_numeric(tree: Any, key: str) -> float:
    if isinstance(tree, dict):
        for k, v in tree.items():
            if k == key and isinstance(v, (int, float)):
                return float(v)
            found = _find_numeric(v, key)
            if found is not None:
                return found
    elif isinstance(tree, list):
        for item in tree:
            found = _find_numeric(item, key)
            if found is not None:
                return found
    return None


def build(out_path: Path = DEFAULT_OUT) -> dict[str, Any]:
    certificate = _load_json(CERTIFICATE_JSON)
    endpoint = _load_json(ENDPOINT_JSON)

    alpha_u = _find_numeric(certificate, "alpha_u")
    alpha2_mz = _find_numeric(certificate, "alpha2_mz")
    alpha_y_mz = _find_numeric(certificate, "alphaY_mz")
    gap_interval = [
        float(v)
        for v in endpoint["compare_only"]["same_scheme_anchor_gap_interval_inv_alpha"]
    ]

    x_unified = 1.0 / alpha_u
    x1_mz = (3.0 / 5.0) / alpha_y_mz
    efolds = (x1_mz - x_unified) * 2.0 * math.pi / B_ONE_LOOP[0]

    # One-loop reproduction: the SU(2) coupling is an independent check because
    # the e-fold count was inferred from the hypercharge line alone.
    x_one_loop, _ = integrate_two_loop(x_unified, efolds, 0.0, two_loop=False)
    a0_one_loop = alpha_em_inv(x_one_loop)
    a0_certificate = 1.0 / alpha2_mz + 1.0 / alpha_y_mz
    reproduction_residual_alpha2 = x_one_loop[1] - 1.0 / alpha2_mz

    scan_rows = []
    for yt_u in YUKAWA_SCAN:
        x_two_loop, yt_mz = integrate_two_loop(x_unified, efolds, yt_u, two_loop=True)
        shift = alpha_em_inv(x_two_loop) - a0_one_loop
        scan_rows.append(
            {
                "yt_at_unification": yt_u,
                "yt_at_mz": yt_mz,
                "alpha3_mz": 1.0 / x_two_loop[2],
                "anchor_shift_inv_alpha": shift,
                "inside_certified_gap": gap_interval[0] <= shift <= gap_interval[1],
            }
        )
    shifts = [row["anchor_shift_inv_alpha"] for row in scan_rows]
    shift_interval = [min(shifts), max(shifts)]

    # Balance requirement: what the undeclared threshold map would have to
    # supply for the two-loop continuation plus the declared scheme spread to
    # land at the certified gap center.  Naive delta-b accounting between the
    # MSSM and SM one-loop electromagnetic coefficients gives the equivalent
    # single effective threshold scale.  Diagnostic only; not a claim.
    gap_center = 0.5 * (gap_interval[0] + gap_interval[1])
    shift_center = 0.5 * (shift_interval[0] + shift_interval[1])
    required_threshold_shift = gap_center - shift_center - SCHEME_CONVERSION_MSBAR_TO_ON_SHELL
    delta_b_em = 12.0 - 11.0 / 3.0  # MSSM minus SM electromagnetic one-loop coefficient
    effective_threshold_gev = MZ_GEV * math.exp(
        -required_threshold_shift * 2.0 * math.pi / delta_b_em
    )

    result = {
        "artifact": "oph_a0_scheme_bridge_source_branch_audit",
        "issue": 545,
        "generated_utc": _timestamp(),
        "row_class": "source_branch_audit_fail_closed",
        "guards": {
            "source_only_bridge_emitted": False,
            "public_promotion_allowed": False,
            "measured_alpha_in_two_loop_computation": False,
            "certified_gap_used_only_as_membership_test": True,
            "new_axiom_introduced": False,
        },
        "inputs": {
            "alpha_u": alpha_u,
            "alpha2_mz_one_loop": alpha2_mz,
            "alphaY_mz_one_loop": alpha_y_mz,
            "boundary_ref": _artifact_ref(CERTIFICATE_JSON),
            "certified_gap_interval": gap_interval,
            "gap_ref": _artifact_ref(ENDPOINT_JSON),
            "rg_contract_ref": _artifact_ref(RG_CONTRACT_JSON),
            "efolds_unification_to_mz": efolds,
        },
        "one_loop_reproduction": {
            "a0_integrated": a0_one_loop,
            "a0_certificate": a0_certificate,
            "residual": a0_one_loop - a0_certificate,
            "independent_alpha2_check_residual": reproduction_residual_alpha2,
        },
        "two_loop_scan": scan_rows,
        "anchor_shift_interval_over_yukawa_conventions": shift_interval,
        "membership_test": {
            "certified_gap_interval": gap_interval,
            "any_convention_lands_inside": any(r["inside_certified_gap"] for r in scan_rows),
            "overshoot_factor_range": [
                shift_interval[0] / gap_interval[1],
                shift_interval[1] / gap_interval[0],
            ],
        },
        "balance_requirement": {
            "scheme_conversion_msbar_to_on_shell_declared_reference": (
                SCHEME_CONVERSION_MSBAR_TO_ON_SHELL
            ),
            "required_threshold_shift_inv_alpha": required_threshold_shift,
            "delta_b_em_mssm_minus_sm": delta_b_em,
            "equivalent_single_effective_threshold_gev": effective_threshold_gev,
            "reading": (
                "for the two-loop MSSM continuation plus the declared scheme spread "
                "to land at the certified gap center, the threshold map must remove "
                "about this much inverse-alpha, equivalent under naive delta-b "
                "accounting to a single effective superpartner threshold at the "
                "quoted scale; recorded as the quantitative content the scheme-lock "
                "object must supply, not as a physics claim"
            ),
        },
        "verdict": {
            "status": "source_branch_not_closable_under_declared_conventions",
            "statement": (
                "The two-loop MSSM continuation of the declared one-loop convention "
                "shifts the anchor by two to three times the certified gap for every "
                "top-Yukawa convention, and restoring the balance requires threshold "
                "and scheme choices that the RG matching contract classifies as "
                "hidden fit parameters when undeclared. Issue #545 therefore stays "
                "open on the source branch; the empirical class (Route A) remains "
                "the carrying lane, and the refined missing object is the declared "
                "scheme-lock plus threshold-map pair of the RG contract."
            ),
            "forbidden_shortcut_avoided": "using_threshold_choices_as_hidden_fit_parameters",
            "refined_missing_objects": ["scheme_lock", "threshold_map"],
            "route_a_status": "unchanged_empirical_class_carrying_lane",
        },
    }

    out_path.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    result = build(args.out)
    print(f"one-loop reproduction residual: {result['one_loop_reproduction']['residual']:+.2e}")
    print(f"two-loop anchor shift interval: {result['anchor_shift_interval_over_yukawa_conventions']}")
    print(f"certified gap: {result['membership_test']['certified_gap_interval']}")
    print(f"lands inside: {result['membership_test']['any_convention_lands_inside']}")
    print(f"verdict: {result['verdict']['status']}")


if __name__ == "__main__":
    main()
