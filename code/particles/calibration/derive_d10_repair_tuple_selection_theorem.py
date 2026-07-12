#!/usr/bin/env python3
"""D10 repair-tuple selection theorem (conditional) and candidate discrimination.

The unconditional theorem on this surface leaves a free quadratic family
beneath the target-free D10 repair:

    tau2_exact    = -c * eta_source^2
    delta_n_exact =  d * (1 - beta_EW) * eta_source^2,   (c, d) free.

Two principled selections are on disk:

  A. running-tree law (EWTargetFreeRepairValueLaw_D10):
     tau2 = -lambda (1 + (2/3) eta + (1 - beta/6) eta^2),
     dn   =  lambda (1 + (4/3) eta + (2 - beta/6) eta^2),
     lambda = eta^2 / (4 beta); effective (c, d) follow by division.

  B. color-balanced quadratic repair descent
     (ColorBalancedQuadraticRepairDescent_D10):
     c = sqrt(N_c)/2, d = N_c/2 on the realized N_c = 3 branch.

This module emits the conditional selection theorem for B: under the three
descent axioms (quadratic leading order; charged contraction through the
color-singlet projection of the realized triplet carrier, weight
sqrt(N_c)/2; coherent neutral color sum, weight N_c/2), the pair (c, d) is
uniquely (sqrt(3)/2, 3/2). The axioms are named conditions, so the artifact
is a conditional selection theorem, and the source derivation of the axioms
is the remaining content of the Higgs-pipeline issue (#521). The
unconditional underdetermination stays in force; promotion stays blocked.

The artifact also records the discrimination table: the forward W/Z/H/top
observables under both candidates, their spread, and the compare-only
measured rows. The spread (about 8 MeV in m_W, 0.03 MeV in m_Z, 16 MeV in
m_H, 30 MeV in m_t) sits below current experimental resolution, so the
selection is a theorem task, not an empirical one.

Run:
    python3 code/particles/calibration/derive_d10_repair_tuple_selection_theorem.py
writes code/particles/runs/calibration/d10_repair_tuple_selection_theorem.json.
"""

from __future__ import annotations

import json
import math
import pathlib
from datetime import datetime, timezone

HERE = pathlib.Path(__file__).resolve().parent
RUNS = HERE.parent / "runs" / "calibration"
OUT_PATH = RUNS / "d10_repair_tuple_selection_theorem.json"

SQRT_PI = math.sqrt(math.pi)
N_C = 3


def load_basis() -> dict:
    with open(RUNS / "d10_ew_target_free_repair_value_law.json", encoding="utf-8") as f:
        law = json.load(f)
    basis = dict(law["basis"])
    basis["lambda_EW"] = basis["eta_source"] ** 2 / (4.0 * basis["beta_EW"])
    return basis


def candidate_a_tuple(basis: dict) -> dict:
    eta, beta = basis["eta_source"], basis["beta_EW"]
    lam = basis["lambda_EW"]
    tau2 = -lam * (1.0 + (2.0 / 3.0) * eta + (1.0 - beta / 6.0) * eta * eta)
    dn = lam * (1.0 + (4.0 / 3.0) * eta + (2.0 - beta / 6.0) * eta * eta)
    return {
        "id": "running_tree_value_law",
        "tau2_exact": tau2,
        "delta_n_exact": dn,
        "c_effective": -tau2 / (eta * eta),
        "d_effective": dn / ((1.0 - beta) * eta * eta),
    }


def candidate_b_tuple(basis: dict) -> dict:
    eta, beta = basis["eta_source"], basis["beta_EW"]
    c = math.sqrt(N_C) / 2.0
    d = N_C / 2.0
    return {
        "id": "color_balanced_quadratic_descent",
        "tau2_exact": -c * eta * eta,
        "delta_n_exact": d * (1.0 - beta) * eta * eta,
        "c_effective": c,
        "d_effective": d,
    }


def forward_wz(basis: dict, tau2: float, dn: float) -> dict:
    """Closed forward map of the value-law artifact."""
    a2, a_y = basis["alpha2_mz"], basis["alphaY_mz"]
    eta, v = basis["eta_source"], basis["v_report_gev"]
    a2p = a2 * (1.0 + tau2)
    a_y_star = a_y * (1.0 - 2.0 * eta)
    d_par = a_y * (8.0 * eta * tau2 * tau2 - tau2) / (1.0 + 4.0 * tau2 * tau2)
    d_perp = (a2 + a_y) * dn
    a_yp = a_y_star + d_par + d_perp
    return {
        "MW_pole_gev": v * math.sqrt(math.pi * a2p),
        "MZ_pole_gev": v * math.sqrt(math.pi * (a2p + a_yp)),
        "sin2w_eff": a_yp / (a2p + a_yp),
        "alpha_em_eff_inv": (a2p + a_yp) / (a2p * a_yp),
    }


def forward_ht(basis: dict, tau2: float, dn: float, surface: dict) -> dict:
    """Closed D11 split map with the declared Jacobian cores."""
    eta, beta, lam = basis["eta_source"], basis["beta_EW"], basis["lambda_EW"]
    rho = math.log(1.0 + tau2)
    top_res = (-tau2 * eta ** 2 + (1.0 + beta / 28.0) * eta ** 6
               + eta ** 8 / 14.0 + eta ** 9 / 27.0)
    higgs_res = (eta ** 5 - (3.0 / 25.0) * eta ** 6
                 + lam * eta ** 6 / 18.0 + eta ** 8 / (2.0 * beta))
    pi_y = (eta + (1.5 + beta / 4.0) * rho + top_res) / SQRT_PI
    pi_lam = (eta - (4.0 / 3.0 - beta / 54.0) * rho + higgs_res) / SQRT_PI
    core, jac = surface["core"], surface["jacobian"]
    d_y = pi_y * core["y_t_core_mt"]
    d_lam = -(16.0 / 9.0) * pi_lam * core["lambda_core_mt"]
    return {
        "pi_y": pi_y,
        "pi_lambda": pi_lam,
        "mt_pole_gev": core["mt_pole_core_gev"] + jac["d_mt_pole_d_y_t"] * d_y,
        "mH_gev": core["mH_core_gev"] + jac["d_mH_d_lambda"] * d_lam,
    }


def build() -> dict:
    basis = load_basis()
    with open(RUNS / "d11_declared_calibration_surface.json", encoding="utf-8") as f:
        surface = json.load(f)

    cands = [candidate_a_tuple(basis), candidate_b_tuple(basis)]
    rows = {}
    for cand in cands:
        wz = forward_wz(basis, cand["tau2_exact"], cand["delta_n_exact"])
        ht = forward_ht(basis, cand["tau2_exact"], cand["delta_n_exact"], surface)
        rows[cand["id"]] = {**cand, **wz, **ht}

    a, b = rows["running_tree_value_law"], rows["color_balanced_quadratic_descent"]
    spread = {key: abs(a[key] - b[key])
              for key in ("MW_pole_gev", "MZ_pole_gev", "sin2w_eff",
                          "mt_pole_gev", "mH_gev")}

    return {
        "artifact": "oph_d10_repair_tuple_selection_theorem",
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "object_id": "D10RepairTupleSelectionTheorem",
        "status": "conditional_selection_theorem_axioms_named",
        "promotion_allowed": False,
        "free_family": {
            "tau2_exact": "-c * eta_source^2",
            "delta_n_exact": "d * (1 - beta_EW) * eta_source^2",
            "statement": "the current corpus underdetermines (c, d); this is "
                         "the unconditional boundary beneath every candidate",
        },
        "selection_axioms": {
            "A1_quadratic_leading_order": "the first nonzero beyond-carrier "
                "repair is quadratic in eta_source (unconditional on this "
                "surface)",
            "A2_charged_contraction": "the charged leg contracts through the "
                "color-singlet projection of the realized triplet carrier, "
                "weight sqrt(N_c)/2, so c = sqrt(N_c)/2",
            "A3_neutral_color_sum": "the neutral leg sums the N_c colors "
                "coherently against the (1 - beta_EW) hypercharge complement, "
                "so d = N_c/2",
        },
        "conditional_theorem": {
            "statement": "Under A1-A3 on the realized N_c = 3 branch, the "
                         "repair pair is uniquely (c, d) = (sqrt(3)/2, 3/2); "
                         "the running-tree law is the lambda-generated "
                         "companion selection and differs from it at relative "
                         "order one in (c, d) while moving every forward "
                         "observable by less than current experimental "
                         "resolution.",
            "proof": "A1 fixes the monomial degree; A2 fixes c by evaluating "
                     "the singlet projection norm on the triplet, "
                     "|P_singlet T|/|T| = 1/sqrt(N_c) against the half-weight "
                     "charged normalization, giving c = sqrt(N_c)/2; A3 fixes "
                     "d as the coherent color multiplicity N_c against the "
                     "half-weight neutral normalization, d = N_c/2. Uniqueness "
                     "is immediate since (c, d) are determined scalars.",
            "remaining_content": "the source derivation of A2 and A3 from the "
                                 "realized carrier package is the open step "
                                 "of the Higgs-pipeline issue; until then the "
                                 "selection is conditional and promotion "
                                 "stays blocked",
        },
        "basis": basis,
        "candidates": rows,
        "selection_spread": spread,
        "discrimination": {
            "statement": "the candidate spread is below current experimental "
                         "resolution in every observable (m_W spread ~8 MeV "
                         "vs sigma 13 MeV; m_Z spread ~0.03 MeV; m_H spread "
                         "~16 MeV vs sigma 110 MeV; m_t spread ~30 MeV vs "
                         "sigma 600 MeV), so the selection cannot be decided "
                         "empirically today and must be decided by the source "
                         "descent theorem",
            "compare_only_references": {
                "MW_pole_gev": [80.3692, 0.0133],
                "MZ_pole_gev": [91.1876, 0.0021],
                "mH_gev": [125.13, 0.11],
                "mt_pole_gev": [172.1, 0.6],
            },
        },
    }


def main() -> int:
    report = build()
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_PATH, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=1)
        f.write("\n")
    for cid, row in report["candidates"].items():
        print(f"{cid}:")
        print(f"  MW={row['MW_pole_gev']:.5f}  MZ={row['MZ_pole_gev']:.5f}  "
              f"sin2w={row['sin2w_eff']:.6f}")
        print(f"  mH={row['mH_gev']:.5f}  mt={row['mt_pole_gev']:.5f}")
    print("spread:", {k: round(v, 6) for k, v in report["selection_spread"].items()})
    print(f"wrote {OUT_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
