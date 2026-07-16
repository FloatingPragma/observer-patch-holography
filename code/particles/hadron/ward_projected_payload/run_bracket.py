#!/usr/bin/env python3
"""Emit a source-only Ward-projected diagnostic bracket (G1).

Runs the declared variant grid through the payload harness and reports the
resulting sampled-module envelopes for S_had, conditional S_QEW, c_Q,
Delta_had, Delta_source_total, and the residual against the implemented 1-x
screen. The grid is declared, not tuned:

- parton_free: 1 run.
- pqcd: Lambda3 in {lane_lo, lane_central, lane_hi} x k in {2, 4, 8}
  x below-cutoff in {free, zero} x truncation order in {1, 2, 3}: 54 runs.
- constituent: kappa = 1, Lambda3 in {lane_lo, lane_central, lane_hi}: 3 runs.

The chain's implemented screen S = 1 - x is reported as a reference row and
does not enter the bracket. This process contains and reads no comparison
target, comparison tolerance, measurement-located P, or scoring rule. It uses
only the source-derived internal Stage-5 root. Its output is a diagnostic
sampled-grid envelope, not a certified interval and not the P-domain function
required by the corrective target contract. Seal the output before passing it
to ``score_bracket.py`` in a separate process.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from decimal import Decimal
from pathlib import Path
from typing import Any

import payload_harness as ph
import spectral_modules as sm

HERE = Path(__file__).resolve().parent
DEFAULT_OUT = HERE / "runtime" / "ward_projected_payload_source_bracket_current.json"
ARTIFACT_SCHEMA_VERSION = 3


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
    if Decimal(ep.p) != Decimal(ph.P_STAR_INTERNAL):
        raise ValueError(
            "run_bracket is restricted to the source-derived internal P; "
            "measurement-located or target-supplied P is prohibited"
        )

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
                "delta_source_total_alpha_inv": payload["delta_source_total_alpha_inv"],
                "s_qew_effective": diag["s_qew_effective"],
                "s_hadronic": diag["s_hadronic"],
                "c_q_implied": diag["c_q_implied"],
                "delta_source_residual_vs_implemented_alpha_inv": diag[
                    "delta_source_residual_vs_implemented_alpha_inv"
                ],
                "positivity_ok": diag["positivity_ok"],
                "content_sha256": payload["content_sha256"],
            }
        )

    def interval(field: str) -> dict[str, float]:
        values = [row[field] for row in rows]
        return {
            "lo": min(values),
            "hi": max(values),
            "width": max(values) - min(values),
        }

    x = ep.x_screen
    naive = ph.quark_naive_transport(ep)
    lepton = ph.lepton_transport(ep)
    screened_impl = ph.implemented_screen(ep) * naive

    payload = {
        "artifact": "oph_ward_projected_payload_source_bracket",
        "schema_version": ARTIFACT_SCHEMA_VERSION,
        "label": ph.SOURCE_ONLY_LABEL,
        "source_family_id": "d10_running_tree",
        "current": "U1_Q",
        "current_definition_id": "ward_projected_U1_Q_once_subtracted",
        "scheme": {
            "same_subtraction_as_a0": True,
            "scheme_id": "d10_ward_projected_once_subtracted_at_mZ2",
            "normalization_convention": "R_ratio_massless_parton_NcQ2",
            "kernel": "mZ^2/(3*pi*s*(s+mZ^2))",
        },
        "evaluation_point": ep.to_json(),
        "p_domain": {
            "kind": "singleton_source_diagnostic",
            "lo": ep.p,
            "hi": ep.p,
            "units": "dimensionless",
            "eligible_as_registered_domain": False,
        },
        "payload_object": "sampled_module_grid_envelope_at_singleton_source_P",
        "coordinate_schema": {
            "delta_source_total_alpha_inv": {
                "kind": "total",
                "units": "inverse_alpha",
                "definition": "delta_lep + delta_had + delta_EW",
                "artifact_path": "bracket.delta_source_total_alpha_inv",
                "scoring_role": "map_input_only",
            },
            "delta_source_residual_vs_implemented_alpha_inv": {
                "kind": "residual",
                "units": "inverse_alpha",
                "definition": (
                    "delta_source_total - "
                    "(delta_lep + implemented_screen * delta_quark_naive)"
                ),
                "artifact_path": (
                    "bracket.delta_source_residual_vs_implemented_alpha_inv"
                ),
                "scoring_role": "diagnostic_only",
            },
            "s_qew_effective": {
                "kind": "screening_ratio_qew",
                "units": "dimensionless",
                "definition": ("(delta_had + delta_EW) / delta_quark_naive_one_loop"),
                "artifact_path": "bracket.s_qew_effective",
                "scoring_role": "diagnostic_only",
                "status": "conditional_on_unproven_delta_EW_zero_branch",
            },
            "s_hadronic": {
                "kind": "screening_ratio_hadronic",
                "units": "dimensionless",
                "definition": "delta_had / delta_quark_naive_one_loop",
                "artifact_path": "bracket.s_hadronic",
                "scoring_role": "diagnostic_only",
            },
        },
        "grid": {
            "fast_mode": fast,
            "declared": (
                "TEST-ONLY REDUCED GRID: parton_free; pqcd "
                "Lambda3{lane_central} x k{4} x below{free,zero} x "
                "order{3}; constituent kappa=1 x Lambda3{lane_central}"
                if fast
                else "parton_free; pqcd Lambda3{lane_lo,lane_central,lane_hi} x "
                "k{2,4,8} x below{free,zero} x order{1,2,3}; constituent "
                "kappa=1 x Lambda3{lane_lo,lane_central,lane_hi}"
            ),
            "eligible_for_scientific_interpretation": False,
            "runs": len(rows),
            "gauss_n": gauss_n,
            "splits_per_decade": splits_per_decade,
            "envelope_semantics": (
                "sampled module-variant extrema only; not numerical or theory "
                "error bounds and not an interval certificate"
            ),
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
            "delta_source_total_alpha_inv": interval("delta_source_total_alpha_inv"),
            "delta_source_residual_vs_implemented_alpha_inv": interval(
                "delta_source_residual_vs_implemented_alpha_inv"
            ),
            "s_qew_effective": interval("s_qew_effective"),
            "s_hadronic": interval("s_hadronic"),
            "c_q_implied": interval("c_q_implied"),
        },
        "certification": {
            "status": "uncertified_sampled_grid_envelope",
            "numerical_error_interval": None,
            "theory_error_interval": None,
            "derivative_or_lipschitz_bound_over_P_domain": None,
            "delta_EW_gate": "open_declared_zero_branch_unproven",
            "sampled_grid_extrema_interval_certificate": False,
        },
        "scoring_status": "NOT_EVALUABLE_SOURCE_DIAGNOSTIC",
        "scoring_process": (
            "score_bracket.py must run after artifact sealing; the canonical "
            "corrective target currently fails closed before scoring"
        ),
        "promotion_allowed": False,
        "promotion_reason": (
            "singleton sampled-grid envelope is not a certified P-domain payload"
        ),
        "target_or_measurement_inputs_used_in_computation": False,
        "wall_time_seconds": round(time.perf_counter() - start, 3),
    }
    digest_source = {k: v for k, v in payload.items() if k != "wall_time_seconds"}
    payload["content_sha256"] = hashlib.sha256(
        json.dumps(
            digest_source,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        ).encode("utf-8")
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
    out.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "runs": payload["grid"]["runs"],
                "delta_source_total_bracket": payload["bracket"][
                    "delta_source_total_alpha_inv"
                ],
                "s_qew_effective_bracket": payload["bracket"]["s_qew_effective"],
                "scoring_status": payload["scoring_status"],
                "content_sha256": payload["content_sha256"],
                "wall_time_seconds": payload["wall_time_seconds"],
                "output": str(out),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
