"""Evaluate the frozen flavor-orbit selector candidate menu.

Computes, for each of the twelve preregistered candidates in CANDIDATES.md, the fixed
spread moduli (sigma_u, sigma_d) as pure functions of the computed fixed point P_fwd
and declared integers, then the implied dimensionless Yukawa spread ratios, and writes
runtime/candidates_evaluated.json.

Blindness (SELECTOR_SPEC.md clause S3): no measured mass, no measured ratio, and no
reference artifact appears anywhere in this module. The only numeric inputs are the
computed fixed point, its certified source-audit companion, declared integer counts,
and the ray shape datum used exclusively for the downstream gap display.
"""

from __future__ import annotations

import json
import math
import pathlib

# Computed fixed point (code/P_derivation/FULL_DERIVATION.md, CL-6 forward value).
P_FWD = 1.630972095858897

# Certified source-audit companion alpha_U(P_fwd) (same derivation, same branch).
ALPHA_U_P_FWD = 0.041124247441816685

# Ray shape datum of the ordered profiles (adjacent-gap ratio of the up ray; the down
# ray carries its reciprocal). Domain data of the rays; used only for the implied
# per-generation gap display, never inside a candidate closed form.
RHO_ORD = 1.2942849363777058

# Declared integers (CANDIDATES.md, declared structure block).
VERTICES = 12
FACES = 20
EDGES = 30
FLAGS = 60
WRITE_CHECK_PORTS = 24
Z6_ORDER = 6
GENERATIONS = 3

FROZEN_CANDIDATE_COUNT = 12

ARTIFACT_PATH = pathlib.Path(__file__).resolve().parent / "runtime" / "candidates_evaluated.json"


def candidate_menu() -> list[dict]:
    """The frozen menu. Closed forms mirror CANDIDATES.md exactly."""
    p = P_FWD
    reserve = math.exp(-p / WRITE_CHECK_PORTS)
    depth = math.log(1.0 / ALPHA_U_P_FWD)
    quad = math.sqrt(8.0 / 5.0)
    return [
        {
            "id": "C-01",
            "name": "Z6-traced face/vertex budget",
            "closed_form": {"sigma_u": "P*F/6", "sigma_d": "P*V/6"},
            "sigma_u": p * FACES / Z6_ORDER,
            "sigma_d": p * VERTICES / Z6_ORDER,
        },
        {
            "id": "C-02",
            "name": "Z6-traced edge/face budget",
            "closed_form": {"sigma_u": "P*E/6", "sigma_d": "P*F/6"},
            "sigma_u": p * EDGES / Z6_ORDER,
            "sigma_d": p * FACES / Z6_ORDER,
        },
        {
            "id": "C-03",
            "name": "Z6-traced write-check/port budget",
            "closed_form": {"sigma_u": "P*N_wc/6", "sigma_d": "P*12/6"},
            "sigma_u": p * WRITE_CHECK_PORTS / Z6_ORDER,
            "sigma_d": p * VERTICES / Z6_ORDER,
        },
        {
            "id": "C-04",
            "name": "Pixel-budget orbit product",
            "closed_form": {"sigma_u": "(P/4)*F", "sigma_d": "(P/4)*V"},
            "sigma_u": (p / 4.0) * FACES,
            "sigma_d": (p / 4.0) * VERTICES,
        },
        {
            "id": "C-05",
            "name": "Port-total pixel budget, F:V split",
            "closed_form": {"sigma_u": "(P/4)*24*(5/8)", "sigma_d": "(P/4)*24*(3/8)"},
            "sigma_u": (p / 4.0) * WRITE_CHECK_PORTS * FACES / (FACES + VERTICES),
            "sigma_d": (p / 4.0) * WRITE_CHECK_PORTS * VERTICES / (FACES + VERTICES),
        },
        {
            "id": "C-06",
            "name": "Reserve-withheld traced budget",
            "closed_form": {"sigma_u": "(P*F/6)*e^(-P/24)", "sigma_d": "(P*V/6)*e^(-P/24)"},
            "sigma_u": p * FACES / Z6_ORDER * reserve,
            "sigma_d": p * VERTICES / Z6_ORDER * reserve,
        },
        {
            "id": "C-07",
            "name": "Reserve-restored traced budget",
            "closed_form": {"sigma_u": "(P*F/6)*e^(+P/24)", "sigma_d": "(P*V/6)*e^(+P/24)"},
            "sigma_u": p * FACES / Z6_ORDER / reserve,
            "sigma_d": p * VERTICES / Z6_ORDER / reserve,
        },
        {
            "id": "C-08",
            "name": "Orbit log-cardinality",
            "closed_form": {"sigma_u": "ln(60)", "sigma_d": "ln(20)"},
            "sigma_u": math.log(FLAGS),
            "sigma_d": math.log(FACES),
        },
        {
            "id": "C-09",
            "name": "P-scaled log-cardinality",
            "closed_form": {"sigma_u": "P*ln(60)", "sigma_d": "P*ln(20)"},
            "sigma_u": p * math.log(FLAGS),
            "sigma_d": p * math.log(FACES),
        },
        {
            "id": "C-10",
            "name": "Transmutation depth, F/V split",
            "closed_form": {"sigma_u": "(5/3)*ln(1/alpha_U)", "sigma_d": "ln(1/alpha_U)"},
            "sigma_u": (FACES / VERTICES) * depth,
            "sigma_d": depth,
        },
        {
            "id": "C-11",
            "name": "Quadrupole-read traced budget",
            "closed_form": {"sigma_u": "sqrt(8/5)*(P*F/6)", "sigma_d": "sqrt(8/5)*(P*V/6)"},
            "sigma_u": quad * p * FACES / Z6_ORDER,
            "sigma_d": quad * p * VERTICES / Z6_ORDER,
        },
        {
            "id": "C-12",
            "name": "Flag-total traced budget, F:V split",
            "closed_form": {"sigma_u": "(P*G/6)*(5/8)", "sigma_d": "(P*G/6)*(3/8)"},
            "sigma_u": p * FLAGS / Z6_ORDER * FACES / (FACES + VERTICES),
            "sigma_d": p * FLAGS / Z6_ORDER * VERTICES / (FACES + VERTICES),
        },
    ]


def implied_ratios(sigma_u: float, sigma_d: float) -> dict:
    """Implied dimensionless Yukawa spread ratios from the fixed moduli.

    Convention: sigma_q = (1/2) ln(y_q3 / y_q1). The unit-span profiles carry
    adjacent-gap ratios rho_ord (up) and 1/rho_ord (down), so
    ln(y_c/y_u) = 2 sigma_u rho/(1+rho), ln(y_t/y_c) = 2 sigma_u/(1+rho),
    ln(y_s/y_d) = 2 sigma_d/(1+rho),     ln(y_b/y_s) = 2 sigma_d rho/(1+rho).
    """
    rho = RHO_ORD
    return {
        "ln_yt_over_yu": 2.0 * sigma_u,
        "ln_yb_over_yd": 2.0 * sigma_d,
        "ln_yc_over_yu": 2.0 * sigma_u * rho / (1.0 + rho),
        "ln_yt_over_yc": 2.0 * sigma_u / (1.0 + rho),
        "ln_ys_over_yd": 2.0 * sigma_d / (1.0 + rho),
        "ln_yb_over_ys": 2.0 * sigma_d * rho / (1.0 + rho),
        "spread_ratio_sigma_u_over_sigma_d": sigma_u / sigma_d,
        "spread_sum_sigma_u_plus_sigma_d": sigma_u + sigma_d,
    }


def build_payload() -> dict:
    menu = candidate_menu()
    if len(menu) != FROZEN_CANDIDATE_COUNT:
        raise RuntimeError("menu count does not match the frozen count")
    candidates = []
    for entry in menu:
        su = float(entry["sigma_u"])
        sd = float(entry["sigma_d"])
        if not (su > 0.0 and sd > 0.0):
            raise RuntimeError(f"{entry['id']} does not fix both moduli in (R_>0)^2")
        candidates.append(
            {
                "id": entry["id"],
                "name": entry["name"],
                "closed_form": entry["closed_form"],
                "moduli": {"sigma_u": su, "sigma_d": sd},
                "implied": implied_ratios(su, sd),
            }
        )
    return {
        "artifact": "oph_flavor_orbit_selector_candidates_evaluated",
        "spec": "SELECTOR_SPEC.md",
        "menu": "CANDIDATES.md",
        "frozen_candidate_count": FROZEN_CANDIDATE_COUNT,
        "inputs": {
            "P_fwd": P_FWD,
            "alpha_U_P_fwd": ALPHA_U_P_FWD,
            "rho_ord_display_only": RHO_ORD,
            "declared_integers": {
                "vertices": VERTICES,
                "faces": FACES,
                "edges": EDGES,
                "flags": FLAGS,
                "write_check_ports": WRITE_CHECK_PORTS,
                "z6_order": Z6_ORDER,
                "generations": GENERATIONS,
            },
        },
        "blindness": {
            "reads_measured_masses": False,
            "reads_measured_ratios": False,
            "reads_reference_artifacts": False,
            "comparison_deferred_to": "COMPARISON.md",
        },
        "candidates": candidates,
        "sigma_convention": "sigma_q = (1/2) * ln(y_q3 / y_q1)",
        "status": "evaluated_no_selection",
    }


def serialize(payload: dict) -> str:
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    payload = build_payload()
    ARTIFACT_PATH.parent.mkdir(parents=True, exist_ok=True)
    ARTIFACT_PATH.write_text(serialize(payload), encoding="utf-8")
    print(f"wrote {ARTIFACT_PATH}")


if __name__ == "__main__":
    main()
