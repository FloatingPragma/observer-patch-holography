#!/usr/bin/env python3
"""CAP-P candidate family for the capacity readback map F (construction run 2026-07-14).

Port-inversion readback (derivation: F_CONSTRUCTION_2026-07-14.md, Steps 4 and 7):
the selected sector reconstructs the boundary capacity by inverting the twelve-port
read of the invariant screen load X = log(N/pi) (icosahedral screen-sieve theorem,
paper/screen_microphysics_and_observer_synchronization.tex), with the surviving
per-port load degraded by the Z6 reserve.

Multiplicative sub-branches (survival factor s on the load):
    F(N) = pi * (N/pi)^s,  s in {e^(-P/24), 1-P/24, e^(-P/12), (1-P/24)^2}
Additive sub-branches (reserve cost subtracted from the per-port load):
    F(N) = N * e^(-P/2)   (P/24 per oriented slot, 12 ports x mean over 2 slots)
    F(N) = N * e^(-P)     (P/12 per port)

Blindness: same cone as CAP-L; P enters only as the certified forward enclosure.
"""

from __future__ import annotations

import json
from pathlib import Path

from mpmath import iv, mp, mpf

from toy_readback import _endpoints, _interval_json

from F_candidate_capL import (
    ENCLOSURE_HALF_WIDTH,
    P_HI,
    P_LO,
    contraction_certificate_centered,
    fixed_point_enclosure_centered,
)

ARTIFACT_NAME = "oph_capacity_readback_candidate_capP_2026-07-14"
PRECISION = 40


def build() -> dict:
    iv.dps = PRECISION
    mp.dps = PRECISION
    P = iv.mpf([P_LO, P_HI])
    one = iv.mpf(1)

    mult_branches = [
        ("capP.s_poisson_port", "s = e^(-P/24), poisson survival per port", iv.exp(-P / 24)),
        ("capP.s_presence_port", "s = 1 - P/24, presence survival per port", one - P / 24),
        ("capP.s_poisson_pair", "s = e^(-P/12), poisson survival per oriented slot pair", iv.exp(-P / 12)),
        ("capP.s_presence_pair", "s = (1 - P/24)^2, presence survival per oriented slot pair", (one - P / 24) ** 2),
    ]
    add_branches = [
        ("capP.add_slot", "x_read = X/12 - P/24: F(N) = N e^(-P/2)", iv.exp(-P / 2)),
        ("capP.add_port", "x_read = X/12 - P/12: F(N) = N e^(-P)", iv.exp(-P)),
    ]

    rows = []
    for branch_id, desc, s in mult_branches:
        def F(x, s=s):
            return iv.pi * iv.exp(s * iv.log(x / iv.pi))

        def Fp(x, s=s):
            return s * iv.exp((s - 1) * iv.log(x / iv.pi))

        interval = iv.mpf(["2.8", "3.6"])
        cert = contraction_certificate_centered(F, Fp, interval)
        enclosure = fixed_point_enclosure_centered(F, Fp, interval, ENCLOSURE_HALF_WIDTH)
        lo = mpf(enclosure["enclosure"]["lo"])
        hi = mpf(enclosure["enclosure"]["hi"])
        rows.append({
            "branch": branch_id,
            "reading": desc,
            "map": "F(N) = pi * (N/pi)^s",
            "s": _interval_json(s),
            "status": "fixed_point_certified" if (cert["banach_pass"] and enclosure["located"]) else "certificate_failed",
            "contraction_certificate": cert,
            "fixed_point": enclosure,
            "fixed_point_nats": mp.nstr((lo + hi) / 2, 25),
            "fixed_point_exact": "N = pi (the unique solution of (N/pi)^s = N/pi with s != 1)",
            "count_density_p4": {
                "interior_stationary_point": "inherited from the paired count branch",
                "n_star_range_over_admissible_count_branches": {"lo": "1.93", "hi": "18"},
                "verdict": "registered_discrepancy (N_CRC = pi differs from every interior N_star of the lattice)",
            },
        })

    for branch_id, desc, s in add_branches:
        s_lo, s_hi = _endpoints(s)
        rows.append({
            "branch": branch_id,
            "reading": desc,
            "map": "F(N) = s*N (linear through the origin)",
            "s": _interval_json(s),
            "status": "no_positive_fixed_point",
            "existence_note": (
                "F' = s < 1 everywhere, so F is a contraction on (0, inf), and the unique "
                "fixed point is N = 0, outside the admissible interval; F(a) < a for every "
                "a > 0, so no certificate interval exists"
            ),
            "count_density_p4": {"verdict": "no readback fixed point on this row"},
        })

    fps = [mpf(r["fixed_point_nats"]) for r in rows if r["status"] == "fixed_point_certified"]
    summary = {
        "rows": len(rows),
        "by_status": {
            s: sum(1 for r in rows if r["status"] == s)
            for s in sorted({r["status"] for r in rows})
        },
        "certified_fixed_point_min_nats": mp.nstr(min(fps), 20) if fps else None,
        "certified_fixed_point_max_nats": mp.nstr(max(fps), 20) if fps else None,
    }

    return {
        "artifact": ARTIFACT_NAME,
        "specification": "F_READBACK_SPEC.md",
        "derivation": "F_CONSTRUCTION_2026-07-14.md",
        "family": "CAP-P port-inversion readback: F(N) = pi * exp(12 * x_read(N))",
        "interval_backend": {
            "library": "mpmath.iv",
            "precision_decimal_digits": PRECISION,
            "rounding": "mpmath_interval_outward",
            "promotion_backend_required": "arb_or_mpfi_directed_outward",
        },
        "inputs": {
            "P_certified_enclosure": {"lo": P_LO, "hi": P_HI},
            "P_source": "code/P_derivation/runtime/p_interval_contraction_certificate_2026-07-14.json",
            "structure_integers": [12, 24, 6],
        },
        "blindness": {
            "inputs": [
                "certified forward-closure P enclosure (source-side computed fixed point)",
                "declared structure integers 12/24/6 and pi",
            ],
            "reads_measured_lambda": False,
            "reads_sl4_estimate": False,
            "reads_cl3_bridge_value": False,
            "reads_alpha_U": False,
            "dependency_cone": ["mpmath", "toy_readback interval machinery", "P_derivation certified enclosure"],
        },
        "summary": summary,
        "rows": rows,
        "moves_cl7": False,
        "cl7_status": "open",
    }


def main() -> int:
    certificate = build()
    out = Path(__file__).resolve().parent / "runtime" / "F_candidate_capP_certificates.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(certificate["summary"], indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
