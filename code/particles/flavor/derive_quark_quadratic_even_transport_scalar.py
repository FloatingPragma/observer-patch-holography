#!/usr/bin/env python3
"""Record the scalarized quadratic-even D12 quark transport shell.

Chain role: expose the smallest even-transport scalar on the D12 continuation
mass branch without promoting the compare-derived value as an OPH law.

Mathematics: ordered-family quadratic scalarization on the fixed three-point
carrier, with `Q_ord` collapsing the even residual to one centered scalar
`eta_Q_centered`.

OPH-derived inputs: the ordered-family coordinate `x2` together with the closed
spread totals from the D12 continuation quark lane.

Output: a diagnostic-only scalar shell that reduces the D12 even transport to
one value law `eta_Q_centered`.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "quark_quadratic_even_transport_scalar.json"

X2 = -0.5175863354681689
ONE_MINUS_X2_SQ = 1.0 - X2 * X2
Q_ORD = [
    (ONE_MINUS_X2_SQ / 3.0),
    -(2.0 * ONE_MINUS_X2_SQ / 3.0),
    (ONE_MINUS_X2_SQ / 3.0),
]


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_artifact() -> dict[str, object]:
    return {
        "artifact": "oph_quark_quadratic_even_transport_scalar",
        "generated_utc": _timestamp(),
        "scope": "D12_continuation_only",
        "proof_status": "scalarized_transport_shape_closed_value_open",
        "x2": X2,
        "one_minus_x2_squared": ONE_MINUS_X2_SQ,
        "Q_ord": Q_ORD,
        "Q_ord_formula": "((1 - x2^2) / 3) * (1, -2, 1)",
        "eta_Q_centered_name": "eta_Q_centered",
        "eta_Q_projector_formula": "v1 - 2*v2 + v3",
        "kappa_Q_formula": "eta_Q_centered / (2 * (1 - x2^2))",
        "quadratic_even_log_formula_via_Q_ord": "(eta_Q_centered / (2 * (1 - x2^2))) * Q_ord",
        "quadratic_even_log_formula_direct": "(eta_Q_centered / 6) * (1, -2, 1)",
        "current_carrier_audit_lift_seed": {
            "Delta_ud_overlap": 0.14042767914755594,
            "kappa_Q": -0.012343623680372975,
            "eta_Q_centered": -0.018073642054692307,
            "source": "current-carrier audit lift",
        },
        "same_t1_candidate": {
            "t1": 0.6695617711471163,
            "kappa_Q_formula": "-t1 / 54",
            "eta_Q_centered_formula": "2 * (1 - x2^2) * (-t1 / 54) = -((1 - x2^2) / 27) * t1",
            "eta_Q_centered": -0.018155152181872827,
            "kappa_Q": -0.012399292058279932,
        },
        "next_single_residual_object": "eta_Q_centered_value_law",
        "notes": [
            "The D12 even quark residual on the ordered family collapses to one centered scalar eta_Q_centered.",
            "The values recorded here are continuation-level diagnostics only and are not promoted as OPH-forward outputs.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the scalarized quadratic-even D12 quark transport shell.")
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()
    payload = build_artifact()
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
