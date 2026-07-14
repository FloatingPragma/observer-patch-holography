#!/usr/bin/env python3
"""CAP-K candidate family for the capacity readback map F (construction run 2026-07-14).

Cell-count readback (derivation: F_CONSTRUCTION_2026-07-14.md, Step 7): the sector
reconstructs the equal-area chart count K = 4N/P (compact paper,
def:self-closure-density) and reads back N_hat = (P/4) * K_readable, with the
readable cell fraction set by the reserve reading:

    F(N) = s * N,  s in {e^(-P/24), 1 - P/24, 5/6, 1/2}.

Every row is linear through the origin with s < 1: a contraction on (0, inf) whose
unique fixed point is N = 0, outside the admissible interval. The family is recorded
as executed and excluded; no positive fixed point exists on any row.
"""

from __future__ import annotations

import json
from pathlib import Path

from mpmath import iv, mp

from toy_readback import _interval_json

from F_candidate_capL import P_LO, P_HI

ARTIFACT_NAME = "oph_capacity_readback_candidate_capK_2026-07-14"
PRECISION = 40


def build() -> dict:
    iv.dps = PRECISION
    mp.dps = PRECISION
    P = iv.mpf([P_LO, P_HI])
    one = iv.mpf(1)

    branches = [
        ("capK.s_poisson", "readable cell fraction e^(-P/24)", iv.exp(-P / 24)),
        ("capK.s_presence", "readable cell fraction 1 - P/24", one - P / 24),
        ("capK.s_nat_share", "readable nat share (P/4 - P/24)/(P/4) = 5/6", iv.mpf(5) / 6),
        ("capK.s_edge_share", "readable nat share per-edge reserve: 1/2", iv.mpf(1) / 2),
    ]

    rows = []
    for branch_id, desc, s in branches:
        rows.append({
            "branch": branch_id,
            "reading": desc,
            "map": "F(N) = s*N",
            "s": _interval_json(s),
            "status": "no_positive_fixed_point",
            "existence_note": (
                "linear map through the origin with certified s < 1; unique fixed point "
                "N = 0 outside (0, inf); F(a) < a for every a > 0, so the self-map "
                "bracketing pair of spec property P3 does not exist"
            ),
            "count_density_p4": {"verdict": "no readback fixed point on this row"},
        })

    return {
        "artifact": ARTIFACT_NAME,
        "specification": "F_READBACK_SPEC.md",
        "derivation": "F_CONSTRUCTION_2026-07-14.md",
        "family": "CAP-K cell-count readback: F(N) = (P/4) * K_readable(N) = s*N",
        "interval_backend": {
            "library": "mpmath.iv",
            "precision_decimal_digits": PRECISION,
            "rounding": "mpmath_interval_outward",
        },
        "inputs": {
            "P_certified_enclosure": {"lo": P_LO, "hi": P_HI},
            "P_source": "code/P_derivation/runtime/p_interval_contraction_certificate_2026-07-14.json",
            "structure_integers": [24, 6, 4],
        },
        "blindness": {
            "inputs": [
                "certified forward-closure P enclosure (source-side computed fixed point)",
                "declared structure integers",
            ],
            "reads_measured_lambda": False,
            "reads_sl4_estimate": False,
            "reads_cl3_bridge_value": False,
            "reads_alpha_U": False,
            "dependency_cone": ["mpmath", "toy_readback interval machinery", "P_derivation certified enclosure"],
        },
        "summary": {
            "rows": len(rows),
            "by_status": {"no_positive_fixed_point": len(rows)},
            "certified_fixed_point_min_nats": None,
            "certified_fixed_point_max_nats": None,
        },
        "rows": rows,
        "moves_cl7": False,
        "cl7_status": "open",
    }


def main() -> int:
    certificate = build()
    out = Path(__file__).resolve().parent / "runtime" / "F_candidate_capK_certificates.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(certificate["summary"], indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
