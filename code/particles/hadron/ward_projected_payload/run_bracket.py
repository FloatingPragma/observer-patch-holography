#!/usr/bin/env python3
"""Strict first-principles bracket for the Ward-projected payload (G1).

LABEL (mandatory): development bracket, non-blind environment; the protocol
pass requires an isolated re-run.

Runs the declared variant grid through the payload harness and reports the
resulting interval for S_eff, c_Q, Delta_had, Delta_source, and the residual
R_Q against the implemented 1-x screen. The grid is declared, not tuned:

- parton_free: 1 run.
- pqcd: Lambda3 in {lane_lo, lane_central, lane_hi} x k in {2, 4, 8}
  x below-cutoff in {free, zero} x truncation order in {1, 2, 3}: 54 runs.
- constituent: kappa = 1, Lambda3 in {lane_lo, lane_central, lane_hi}: 3 runs.

The chain's implemented screen S = 1 - x is reported as a reference row and
does not enter the bracket. No CODATA/NIST value, measured hadronic cross
section, PDG hadron datum, or empirical endpoint interval enters any number
emitted here. The comparison block against the canon S_required scale is a
non-blind development comparison and is labeled as such.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path
from typing import Any

import payload_harness as ph
import spectral_modules as sm

HERE = Path(__file__).resolve().parent
DEFAULT_OUT = HERE / "runtime" / "ward_projected_payload_bracket_current.json"

PASS_TOLERANCE_ALPHA_INV = 2.1e-8  # ledger-scale pass tolerance in Delta_source units

# Canon values seen in non-blind canon documents (THOMSON_TRANSPORT_THEOREMS.md).
# They are compare-only development targets, never inputs.
CANON_NON_BLIND = {
    "S_required_at_P_C": 0.895400132647658797805800283181670641,
    "c_Q_at_P_C": 0.658025759927155435638230170232360050,
    "Delta_required_at_root": 8.727731131834786107447994009818221065,
    "Delta_impl_at_root": 8.686567119456435565397988595605414327,
    "Delta_missing_at_root": 0.041163999587062704710570535563112120,
    "note": (
        "development bracket, non-blind environment; the protocol pass "
        "requires an isolated re-run"
    ),
}


def declared_grid(fast: bool = False) -> list[Any]:
    modules: list[Any] = [sm.make_parton_free()]
    lambda_keys = ["lane_central"] if fast else ["lane_lo", "lane_central", "lane_hi"]
    k_cuts = [4.0] if fast else [2.0, 4.0, 8.0]
    orders = [3] if fast else [1, 2, 3]
    for key in lambda_keys:
        for k in k_cuts:
            for below in ("free", "zero"):
                for order in orders:
                    modules.append(sm.make_pqcd(key, k, below, order))
    for key in lambda_keys:
        modules.append(sm.make_constituent(key, kappa=1.0))
    return modules


def build_bracket(
    ep: ph.EvaluationPoint,
    *,
    fast: bool = False,
    gauss_n: int = 48,
    splits_per_decade: int = 4,
) -> dict[str, Any]:
    start = time.perf_counter()
    rows: list[dict[str, Any]] = []
    for module in declared_grid(fast=fast):
        payload = ph.emit_delta_source(
            module, ep, gauss_n=gauss_n, splits_per_decade=splits_per_decade
        )
        diag = payload["diagnostics"]
        rows.append(
            {
                "module_id": module.module_id,
                "declared_branch": module.declared_branch,
                "delta_had_alpha_inv": payload["components_alpha_inv"]["delta_had"],
                "delta_source_alpha_inv": payload["delta_source_alpha_inv"],
                "s_effective": diag["s_effective"],
                "c_q_implied": diag["c_q_implied"],
                "r_q_residual_vs_implemented_screen": diag[
                    "r_q_residual_vs_implemented_screen"
                ],
                "positivity_ok": diag["positivity_ok"],
                "content_sha256": payload["content_sha256"],
            }
        )

    def interval(field: str) -> dict[str, float]:
        values = [row[field] for row in rows]
        return {"lo": min(values), "hi": max(values), "width": max(values) - min(values)}

    delta_source_interval = interval("delta_source_alpha_inv")
    s_interval = interval("s_effective")

    x = ep.x_screen
    naive = ph.quark_naive_transport(ep)
    lepton = ph.lepton_transport(ep)
    screened_impl = ph.implemented_screen(ep) * naive

    canon_s = CANON_NON_BLIND["S_required_at_P_C"]
    wall = {
        "pass_tolerance_alpha_inv": PASS_TOLERANCE_ALPHA_INV,
        "bracket_width_alpha_inv": delta_source_interval["width"],
        "width_over_tolerance": delta_source_interval["width"] / PASS_TOLERANCE_ALPHA_INV,
        "required_relative_precision_on_delta_had": PASS_TOLERANCE_ALPHA_INV
        / max(row["delta_had_alpha_inv"] for row in rows),
        "statement": (
            "The first-principles bracket is wider than the pass tolerance by "
            "the factor width_over_tolerance. See PAYLOAD_STATUS.md for the "
            "precision-wall statement with literature scales."
        ),
    }

    payload = {
        "artifact": "oph_ward_projected_payload_first_principles_bracket",
        "label": ph.NON_BLIND_LABEL,
        "source_family_id": "d10_running_tree",
        "current": "U1_Q",
        "scheme": {
            "same_subtraction_as_a0": True,
            "scheme_id": "d10_ward_projected_once_subtracted_at_mZ2",
            "normalization_convention": "R_ratio_massless_parton_NcQ2",
            "kernel": "mZ^2/(3*pi*s*(s+mZ^2))",
        },
        "evaluation_point": ep.to_json(),
        "grid": {
            "fast_mode": fast,
            "declared": (
                "parton_free; pqcd Lambda3{lane_lo,lane_central,lane_hi} x "
                "k{2,4,8} x below{free,zero} x order{1,2,3}; constituent "
                "kappa=1 x Lambda3{lane_lo,lane_central,lane_hi}"
            ),
            "runs": len(rows),
            "gauss_n": gauss_n,
            "splits_per_decade": splits_per_decade,
        },
        "chain_reference": {
            "lepton_delta_alpha_inv": lepton,
            "quark_delta_alpha_inv_naive": naive,
            "implemented_screen_1_minus_x": ph.implemented_screen(ep),
            "quark_delta_alpha_inv_screened_impl": screened_impl,
            "delta_impl_total_alpha_inv": lepton + screened_impl,
            "x_screen": x,
        },
        "rows": rows,
        "bracket": {
            "delta_had_alpha_inv": interval("delta_had_alpha_inv"),
            "delta_source_alpha_inv": delta_source_interval,
            "s_effective": s_interval,
            "c_q_implied": interval("c_q_implied"),
            "r_q_residual_vs_implemented_screen": interval(
                "r_q_residual_vs_implemented_screen"
            ),
        },
        "precision_wall": wall,
        "non_blind_development_comparison": {
            **CANON_NON_BLIND,
            "canon_S_required_inside_s_effective_bracket": bool(
                s_interval["lo"] <= canon_s <= s_interval["hi"]
            ),
        },
        "delta_EW_branch": "declared_zero_branch_unproven (Theorem 4 open)",
        "promotion_allowed": False,
        "external_inputs_used_in_computation": False,
        "wall_time_seconds": round(time.perf_counter() - start, 3),
    }
    digest_source = {k: v for k, v in payload.items() if k != "wall_time_seconds"}
    payload["content_sha256"] = hashlib.sha256(
        json.dumps(digest_source, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    return payload


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    parser.add_argument("--fast", action="store_true")
    parser.add_argument("--precision", type=int, default=ph.DEFAULT_PRECISION)
    args = parser.parse_args()

    ep = ph.build_evaluation_point(precision=args.precision)
    payload = build_bracket(ep, fast=args.fast)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "runs": payload["grid"]["runs"],
                "delta_source_bracket": payload["bracket"]["delta_source_alpha_inv"],
                "s_effective_bracket": payload["bracket"]["s_effective"],
                "width_over_tolerance": payload["precision_wall"]["width_over_tolerance"],
                "canon_S_inside": payload["non_blind_development_comparison"][
                    "canon_S_required_inside_s_effective_bracket"
                ],
                "wall_time_seconds": payload["wall_time_seconds"],
                "output": str(out),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
