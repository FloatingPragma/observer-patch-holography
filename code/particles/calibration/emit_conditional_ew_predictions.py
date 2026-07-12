#!/usr/bin/env python3
"""Conditional electroweak predictions: H, top, W, Z under the named freedoms.

Emits the four-mass prediction surface with every conditionality explicit.
Row class `conditional_on_P_and_repair_selection`: the rows carry (i) the
repair-selection freedom between the two on-disk candidate laws (running
tree and color-balanced descent, see the selection-theorem artifact), and
(ii) the pixel freedom between the calibration P and the empirical-closure
P interval from the Thomson endpoint artifact. Nothing here is promoted;
the promotion gates of the D10/D11 lane are unchanged.

P sensitivity is the eta-channel partial derivative: eta_source = alpha_U
times beta_EW with d alpha_U / dP taken from the hierarchy certificate's
derivative interval; the running couplings, v, and the Jacobian cores are
held at their repo-P values. The full-surface dP row requires re-deriving
the running package at shifted P and is recorded as the open instruction.

Run:
    python3 code/particles/calibration/emit_conditional_ew_predictions.py
writes code/particles/runs/calibration/conditional_ew_predictions_current.json.
"""

from __future__ import annotations

import json
import pathlib
from datetime import datetime, timezone

HERE = pathlib.Path(__file__).resolve().parent
RUNS = HERE.parent / "runs" / "calibration"
OUT_PATH = RUNS / "conditional_ew_predictions_current.json"

import sys
sys.path.insert(0, str(HERE))
from derive_d10_repair_tuple_selection_theorem import (  # noqa: E402
    candidate_a_tuple,
    candidate_b_tuple,
    forward_ht,
    forward_wz,
    load_basis,
)

# pixel inputs
P_CALIBRATION = 1.630968209403959          # public endpoint (hierarchy lane)
P_EMPIRICAL_INTERVAL = (1.631031463127, 1.631051577591)  # Thomson endpoint artifact
D_ALPHA_U_D_P = -10.9905                   # hierarchy certificate derivative (midpoint)

MEASURED = {
    "mH_gev": (125.13, 0.11, "PDG 2025"),
    "mt_pole_gev": (172.1, 0.6, "PDG 2025 direct-average context row"),
    "MW_pole_gev": (80.3692, 0.0133, "PDG 2025 world average"),
    "MZ_pole_gev": (91.1876, 0.0021, "PDG 2025"),
}


def eta_shift(basis: dict, delta_p: float) -> float:
    """eta-channel first-order shift of eta_source for a pixel shift."""
    return basis["beta_EW"] * D_ALPHA_U_D_P * delta_p


def rows_for_basis(basis: dict, surface: dict) -> dict:
    out = {}
    for cand in (candidate_a_tuple(basis), candidate_b_tuple(basis)):
        wz = forward_wz(basis, cand["tau2_exact"], cand["delta_n_exact"])
        ht = forward_ht(basis, cand["tau2_exact"], cand["delta_n_exact"], surface)
        out[cand["id"]] = {
            "MW_pole_gev": wz["MW_pole_gev"],
            "MZ_pole_gev": wz["MZ_pole_gev"],
            "sin2w_eff": wz["sin2w_eff"],
            "mH_gev": ht["mH_gev"],
            "mt_pole_gev": ht["mt_pole_gev"],
        }
    return out


def build() -> dict:
    basis = load_basis()
    with open(RUNS / "d11_declared_calibration_surface.json", encoding="utf-8") as f:
        surface = json.load(f)

    base_rows = rows_for_basis(basis, surface)

    # eta-channel P sensitivity to the empirical-closure interval endpoints
    p_rows = {}
    for label, p_val in (
        ("empirical_closure_P_lo", P_EMPIRICAL_INTERVAL[0]),
        ("empirical_closure_P_hi", P_EMPIRICAL_INTERVAL[1]),
    ):
        shifted = dict(basis)
        shifted["eta_source"] = basis["eta_source"] + eta_shift(
            basis, p_val - P_CALIBRATION)
        shifted["lambda_EW"] = shifted["eta_source"] ** 2 / (4.0 * shifted["beta_EW"])
        p_rows[label] = rows_for_basis(shifted, surface)

    # per-observable envelope over selection and P
    observables = ("mH_gev", "mt_pole_gev", "MW_pole_gev", "MZ_pole_gev")
    envelope, comparison = {}, {}
    for obs in observables:
        values = [base_rows[c][obs] for c in base_rows]
        for pr in p_rows.values():
            values += [pr[c][obs] for c in pr]
        lo, hi = min(values), max(values)
        envelope[obs] = [lo, hi]
        central = 0.5 * (lo + hi)
        ref, sigma, src = MEASURED[obs]
        band_lo, band_hi = ref - sigma, ref + sigma
        comparison[obs] = {
            "conditional_central": central,
            "conditional_envelope": [lo, hi],
            "measured": ref,
            "measured_sigma": sigma,
            "measured_source": src,
            "delta": central - ref,
            "delta_over_sigma": (central - ref) / sigma,
            "envelope_overlaps_one_sigma_band": not (hi < band_lo or lo > band_hi),
            "envelope_inside_one_sigma_band": band_lo <= lo and hi <= band_hi,
        }

    return {
        "artifact": "oph_conditional_ew_predictions",
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "row_class": "conditional_on_P_and_repair_selection",
        "guards": {
            "public_promotion_allowed": False,
            "conditional_display_allowed": True,
            "promotion_blockers": [
                "current_corpus_underdetermination_of_forward_d10_repair_law",
                "d10.same-scheme-anchor-bridge (issue #545)",
                "source derivation of the selection axioms (issue #521)",
            ],
        },
        "inputs": {
            "basis_source": "d10_ew_target_free_repair_value_law.json",
            "surface_source": "d11_declared_calibration_surface.json",
            "P_calibration": P_CALIBRATION,
            "P_empirical_interval": list(P_EMPIRICAL_INTERVAL),
            "d_alpha_U_d_P": D_ALPHA_U_D_P,
            "sensitivity_scope": "eta-channel partial only; running couplings, "
                                 "v, and Jacobian cores held at repo P; the "
                                 "full-surface dP row is the open instruction",
        },
        "rows_at_calibration_P": base_rows,
        "rows_at_empirical_P": p_rows,
        "conditional_envelope": envelope,
        "comparison_compare_only": comparison,
    }


def main() -> int:
    report = build()
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_PATH, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=1)
        f.write("\n")
    print("conditional envelope (selection x P), compare-only deltas:")
    for obs, row in report["comparison_compare_only"].items():
        lo, hi = row["conditional_envelope"]
        print(f"  {obs:14s} [{lo:.4f}, {hi:.4f}]  measured {row['measured']} "
              f"+- {row['measured_sigma']}  delta/sigma {row['delta_over_sigma']:+.2f}  "
              f"overlap_1s={row['envelope_overlaps_one_sigma_band']} "
              f"inside_1s={row['envelope_inside_one_sigma_band']}")
    print(f"wrote {OUT_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
