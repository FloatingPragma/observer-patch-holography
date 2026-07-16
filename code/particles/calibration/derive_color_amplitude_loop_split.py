#!/usr/bin/env python3
"""A2/A3 derivation attempt: the color amplitude/loop split (#521).

Tests whether the Doplicher-Roberts / large-N_c structure of the realized
carrier package fixes the charged and neutral repair coefficients to
c = sqrt(N_c)/2 (A2) and d = N_c/2 (A3), and confronts the prediction with
the cleanest tree-level data.

Framework. The compact-gauge reconstruction (compact paper Theorem 6.1,
DR/Tannaka) returns the gauge group and the color triplet sector rho_3 with
statistical dimension d(rho_3) = N_c = 3. The conjugate intertwiners
R in Hom(iota, rho_3-bar rho_3) satisfy R*R = d(rho_3) = N_c, so a single
intertwiner isometry carries norm sqrt(N_c), while a closed color loop
(R* (rho x 1) R) carries the full dimension N_c.

The two repair channels carry different normalizations:

  charged leg (tau_2, SU(2)_L, the BROKEN mass-generating channel):
      the W mass is driven by the color-singlet condensate amplitude
      <q-bar q>, whose single-channel amplitude scales as the large-N_c
      decay-constant f ~ sqrt(N_c) (one intertwiner isometry R). Hence
      tau_2 = -(sqrt(N_c)/2) eta^2.

  neutral leg (delta_n, U(1)_Y, the UNBROKEN screening channel):
      hypercharge screening is a vacuum-polarization loop closing the color
      line, R* (rho x 1) R = d = N_c. Hence
      delta_n = (N_c/2) (1 - beta_EW) eta^2.

The factor 1/2 in both is the half-weight current normalization shared by
the two legs; it cancels in every relative statement and is the one common
convention.

Historical chart diagnostic. The forward map with a free charged coefficient
was profiled against stale W/Z Breit--Wigner coordinates. The emitted values
are running/tree chart coordinates with no complete map to those references.
The profile is retained for arithmetic provenance only: it cannot select or
exclude a physical coefficient.

Status. This module records a retrospective comparison of the color-balanced
quadratic candidate with the stored W/Z rows and with other constant-c
quadratic proxies.  The proxy c=1/(4 beta) is only the leading constant term of
the complete archived value law; the complete law has eta-dependent return and
central-subtraction terms and is not excluded by this proxy comparison.  The
module supplies a conditional DR/large-N_c mechanism for the sqrt(N_c) versus
raw N_c split, but it does not emit the explicit collar-transport operators or
prove channel identification.  Promotion stays blocked.

Run:
    python3 code/particles/calibration/derive_color_amplitude_loop_split.py
writes code/particles/runs/calibration/color_amplitude_loop_split.json.
"""

from __future__ import annotations

import json
import math
import pathlib
import sys
from datetime import datetime, timezone

HERE = pathlib.Path(__file__).resolve().parent
RUNS = HERE.parent / "runs" / "calibration"
OUT_PATH = RUNS / "color_amplitude_loop_split.json"
sys.path.insert(0, str(HERE))

from derive_d10_repair_tuple_selection_theorem import (  # noqa: E402
    forward_ht,
    forward_wz,
    load_basis,
)

N_C = 3

# cleanest compare-only reference masses (compare-only; never in a solve path)
COMPARE_ONLY_REFERENCES = {
    "MW_chart_gev": (80.3692, 0.0133),
    "MZ_chart_gev": (91.1876, 0.0021),
    "mH_gev": (125.13, 0.11),
    "mt_pole_gev": (172.1, 0.6),
}


def obs(basis: dict, surface: dict, c: float, d: float) -> dict:
    eta, beta = basis["eta_source"], basis["beta_EW"]
    tau2 = -c * eta * eta
    dn = d * (1.0 - beta) * eta * eta
    wz = forward_wz(basis, tau2, dn)
    ht = forward_ht(basis, tau2, dn, surface)
    return {
        "MW_chart_gev": wz["MW_chart_gev"],
        "MZ_chart_gev": wz["MZ_chart_gev"],
        "mH_gev": ht["mH_gev"],
        "mt_pole_gev": ht["mt_pole_gev"],
    }


def profile_min(basis, surface, key, d_fixed, c_grid):
    ref, sig = COMPARE_ONLY_REFERENCES[key]
    best_c, best_chi = None, None
    for c in c_grid:
        val = obs(basis, surface, c, d_fixed)[key]
        chi = ((val - ref) / sig) ** 2
        if best_chi is None or chi < best_chi:
            best_c, best_chi = c, chi
    return best_c, best_chi


def build() -> dict:
    basis = load_basis()
    surface = json.load(open(RUNS / "d11_declared_calibration_surface.json"))
    beta = basis["beta_EW"]

    c_pred = math.sqrt(N_C) / 2.0
    d_pred = N_C / 2.0
    c_grid = [i * 0.0005 for i in range(6001)]     # 0 .. 3
    d_grid = [i * 0.0005 for i in range(6001)]

    mw_best_c, mw_best_chi = profile_min(basis, surface, "MW_chart_gev", d_pred, c_grid)
    # d-profile on MZ (neutral-sensitive) at c = c_pred
    ref, sig = COMPARE_ONLY_REFERENCES["MZ_chart_gev"]
    mz_best_d, mz_best_chi = None, None
    for d in d_grid:
        val = obs(basis, surface, c_pred, d)["MZ_chart_gev"]
        chi = ((val - ref) / sig) ** 2
        if mz_best_chi is None or chi < mz_best_chi:
            mz_best_d, mz_best_chi = d, chi

    # MW+MZ joint 1-sigma band on c (at d = d_pred)
    def chi2_mwmz(c):
        o = obs(basis, surface, c, d_pred)
        return sum(
            ((o[k] - COMPARE_ONLY_REFERENCES[k][0]) / COMPARE_ONLY_REFERENCES[k][1]) ** 2
            for k in ("MW_chart_gev", "MZ_chart_gev")
        )
    joint = [(c, chi2_mwmz(c)) for c in c_grid]
    c_star, chi_star = min(joint, key=lambda t: t[1])
    band = [c for c, x in joint if x <= chi_star + 1.0]
    band_lo, band_hi = min(band), max(band)

    competitors = {
        "color_amplitude_loop  c=sqrt(Nc)/2": c_pred,
        "loop_symmetric        c=Nc/2": N_C / 2.0,
        "loop_symmetric        c=Nc": float(N_C),
        "leading_constant_proxy c=1/(4 beta)": 1.0 / (4.0 * beta),
    }
    competitor_rows = {
        name: {
            "c": c,
            "inside_legacy_reference_error_profile": band_lo <= c <= band_hi,
        }
        for name, c in competitors.items()
    }

    pred = obs(basis, surface, c_pred, d_pred)
    comparison = {}
    for k, (ref_v, sig_v) in COMPARE_ONLY_REFERENCES.items():
        row = {
            "predicted": pred[k],
            "reference": ref_v,
            "reference_experimental_sigma": sig_v,
        }
        if k in {"MW_chart_gev", "MZ_chart_gev"}:
            row.update({
                "physical_comparison_status": "NOT_EVALUABLE",
                "physical_delta": None,
                "physical_pull": None,
                "legacy_experimental_error_only_standardized_difference": (
                    pred[k] - ref_v
                ) / sig_v,
            })
        else:
            row.update({
                "physical_comparison_status": "COMPARE_ONLY",
                "delta_over_sigma": (pred[k] - ref_v) / sig_v,
            })
        comparison[k] = row

    return {
        "artifact": "oph_color_amplitude_loop_split",
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "object_id": "ColorAmplitudeLoopSplit_A2A3",
        "row_class": "retrospective_data_comparison_conditional_mechanism_candidate",
        "guards": {
            "measured_values_in_any_oph_solve_path": False,
            "source_only_theorem_emitted": False,
            "public_promotion_allowed": False,
            "complete_archived_value_law_tested_as_constant_c_competitor": False,
            "wz_physical_comparison_status": "NOT_EVALUABLE",
        },
        "N_c": N_C,
        "dr_mechanism": {
            "statistical_dimension_d_rho3": N_C,
            "intertwiner_norm": "R*R = d(rho_3) = N_c, single isometry ||R|| = sqrt(N_c)",
            "charged_leg": "SU(2)_L broken mass channel = color-singlet condensate "
                           "amplitude, large-N_c decay-constant scaling sqrt(N_c) "
                           "(one intertwiner isometry): c = sqrt(N_c)/2",
            "neutral_leg": "U(1)_Y unbroken screening channel = vacuum-polarization "
                           "loop closing the color line, dimension N_c: d = N_c/2",
            "shared_half_weight": "the 1/2 is the common current normalization; it "
                                  "cancels in every relative statement",
        },
        "prediction": {"c": c_pred, "d": d_pred},
        "data_selection": {
            "MW_best_c": mw_best_c,
            "MW_best_c_chi2": mw_best_chi,
            "MW_prediction_c": c_pred,
            "MW_agreement": f"data optimum {mw_best_c:.4f} vs prediction "
                            f"{c_pred:.4f}",
            "MZ_best_d": mz_best_d,
            "MZ_best_d_chi2": mz_best_chi,
            "MZ_prediction_d": d_pred,
            "legacy_MWMZ_reference_error_profile_band": [band_lo, band_hi],
            "competitors": competitor_rows,
        },
        "comparison_compare_only": comparison,
        "verdict": {
            "loop_symmetric_law": "not physically evaluated; outside only the legacy reference-error profile",
            "leading_constant_proxy": "outside the legacy reference-error profile; no physical exclusion follows",
            "complete_archived_value_law": "not a constant-c competitor and not excluded by this comparison",
            "color_amplitude_loop_law": "retrospective chart candidate; W/Z physical comparison not evaluable",
            "remaining_rigorous_step": "prove from the explicit collar-transport "
                                       "operators that the charged repair channel "
                                       "factors through a single intertwiner "
                                       "(amplitude, sqrt(N_c)) and the neutral "
                                       "through a closed loop (N_c); the current "
                                       "reconstruction fixes N_c and the DR "
                                       "intertwiner norms but does not emit those "
                                       "channel operators",
            "status": "A2A3_retrospective_conditional_mechanism_"
                      "channel_identification_open",
        },
    }


def main() -> int:
    report = build()
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_PATH, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=1)
        f.write("\n")
    ds = report["data_selection"]
    print(f"prediction: c = sqrt(3)/2 = {report['prediction']['c']:.4f}, "
          f"d = 3/2 = {report['prediction']['d']:.4f}")
    print(f"m_W selects c = {ds['MW_best_c']:.4f}  (prediction {report['prediction']['c']:.4f})")
    print(f"m_Z selects d = {ds['MZ_best_d']:.4f}  (prediction {report['prediction']['d']:.4f})")
    profile_band = ds["legacy_MWMZ_reference_error_profile_band"]
    print(
        "legacy MW+MZ reference-error profile on c = "
        f"[{profile_band[0]:.3f}, {profile_band[1]:.3f}]"
    )
    for name, row in ds["competitors"].items():
        print(
            f"  {name:36s} c={row['c']:.4f}  "
            "inside_legacy_profile="
            f"{row['inside_legacy_reference_error_profile']}"
        )
    print(f"status: {report['verdict']['status']}")
    print(f"wrote {OUT_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
