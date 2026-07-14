#!/usr/bin/env python3
"""Preregistered CL-5 electroweak repair sweep (generator G3, 2026-07-14).

Executes every entry of falsification/preregistered/ew_repair_menu_2026-07-14.json
with the repository's own W/Z lane objects and emits
falsification/preregistered/ew_repair_results_2026-07-14.json.

The forward chain per entry is the lane's own:

1. ``build_paper_d10(pix_area=P, n_c=n_c)`` from the pinned legacy D10 module,
   with the gauge one-loop triple taken from the menu axis (``B_MSSM`` or the
   declared ``B_SM_ALPHA_GUT``) by patching the module constant the lane
   already reads.
2. The source-pair reduction ``eta_source = alpha_U cos(2 theta_W0)`` exactly
   as in ``derive_d10_ew_source_transport_pair.build_artifact``.
3. ``strict_branch_two_law_evaluation`` from
   ``derive_wzh_residual_elimination_boundary`` emits both branches of the
   recorded two-way value-law selection on that basis.
4. ``v/E_star = exp(-2 pi / (beta_EW alpha_U)) / sqrt(P)`` (the lane's
   transmutation law in E_star units) and the fixed clock display adapter
   ``E_star = h f_Cs / epsilon_Cs`` convert to GeV; the Stage-3 ``delta_rho``
   Z-only pole surrogate is the single declared scheme switch.

No physics is forked; every formula above is imported from or identical to the
lane code. All axes are discrete; nothing continuous is tuned.
"""

from __future__ import annotations

import argparse
import json
import math
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REPO = ROOT.parent

try:
    from calibration._legacy_d10 import require_legacy_d10_path
except ModuleNotFoundError:
    from _legacy_d10 import require_legacy_d10_path

require_legacy_d10_path()

import particle_masses_paper_d10_d11 as pm  # type: ignore

try:
    from calibration.derive_wzh_residual_elimination_boundary import (
        strict_branch_two_law_evaluation,
    )
except ModuleNotFoundError:
    from derive_wzh_residual_elimination_boundary import (
        strict_branch_two_law_evaluation,
    )

MENU = REPO / "falsification" / "preregistered" / "ew_repair_menu_2026-07-14.json"
DEFAULT_OUT = (
    REPO / "falsification" / "preregistered" / "ew_repair_results_2026-07-14.json"
)
SOURCE_CERT = (
    ROOT / "particles" / "hierarchy" / "certificates"
    / "R_P_source_audit_pixel_certificate.json"
)
CLOCK_CERT = (
    ROOT / "particles" / "hierarchy" / "certificates"
    / "R_gamma_noG_DAG_certificate.json"
)
BOUNDARY_ARTIFACT = (
    ROOT / "particles" / "runs" / "calibration"
    / "wzh_residual_elimination_boundary.json"
)

MW_TARGET, MW_SIGMA = 80.3692, 0.0133
MZ_TARGET, MZ_SIGMA = 91.1880, 0.0020

B_TRIPLES = {
    "mssm": pm.B_MSSM,
    "sm": pm.B_SM_ALPHA_GUT,
}
BETA_AXIS = {"beta3": 2, "beta4": 3, "beta5": 4}  # id -> n_c (beta_EW = n_c + 1)
LAW_AXIS = ("zero_selector", "nonzero_carrier")
Z_AXIS = ("z_tree", "z_stage3")
DELTA_RHO_STAGE3 = 3.0 / (32.0 * math.pi**2)


def e_star_display_gev(epsilon_cs: float) -> float:
    """Fixed clock display adapter, identical to the boundary script."""

    return 6.62607015e-34 * 9192631770.0 / epsilon_cs / 1.602176634e-10


def run_chain(b_id: str, beta_id: str, p_value: float) -> dict[str, object]:
    """One D10 chain solve with the menu's discrete axis values."""

    original = pm.B_MSSM
    pm.B_MSSM = B_TRIPLES[b_id]
    try:
        d10 = pm.build_paper_d10(pix_area=p_value, n_c=BETA_AXIS[beta_id])
    finally:
        pm.B_MSSM = original
    alpha_y = 0.6 * d10.alpha1_mz
    alpha_2 = d10.alpha2_mz
    sin2 = alpha_y / (alpha_y + alpha_2)
    eta_source = d10.alpha_u * (1.0 - 2.0 * sin2)  # alpha_U cos(2 theta_W0)
    beta_ew = pm.beta_ew(BETA_AXIS[beta_id])
    v_over_e = math.exp(-2.0 * math.pi / (beta_ew * d10.alpha_u)) / math.sqrt(p_value)
    two_law = strict_branch_two_law_evaluation(alpha_2, alpha_y, eta_source, v_over_e)
    return {
        "alpha_u": d10.alpha_u,
        "alpha2_mz": alpha_2,
        "alphaY_mz": alpha_y,
        "eta_source": eta_source,
        "beta_EW_exponent": beta_ew,
        "v_over_E_star": v_over_e,
        "mu_u_gev": d10.mu_u,
        "two_law": two_law,
    }


def _wide_mz_fixed_point(
    alpha_u: float, pix_area: float, n_c: int, mu_u: float
) -> tuple[float, float, float, float, float]:
    """The lane's mz fixed-point equations on a wider diagnostic mu window.

    Identical mathematics to ``pm.solve_mz_fixed_point_tree`` (same
    ``v_from_transmutation``, ``run_alphas_from_unification``, and
    ``mz_tree_from_v_and_couplings``); only the bracketing window is widened
    from logspace(0, 5) to logspace(-10, 10). Diagnostic use only.
    """

    import numpy as np

    v_ev = pm.v_from_transmutation(alpha_u, pix_area, n_c)

    def f(mu: float) -> float:
        a1, a2, a3 = pm.run_alphas_from_unification(alpha_u, mu, mu_u)
        return pm.mz_tree_from_v_and_couplings(v_ev, a1, a2) - mu

    grid = np.logspace(-10, 10, 600)
    prev_mu = prev_f = None
    for mu in grid:
        mu = float(mu)
        try:
            val = f(mu)
        except (ValueError, OverflowError, ZeroDivisionError):
            # coupling unphysical at this mu (past a one-loop pole); skip
            prev_mu = prev_f = None
            continue
        if prev_f is not None and val * prev_f < 0:
            lo, hi = float(prev_mu), mu
            flo = float(prev_f)
            for _ in range(90):
                mid = math.sqrt(lo * hi)
                fm = f(mid)
                if flo * fm > 0:
                    lo, flo = mid, fm
                else:
                    hi = mid
            mz_run = 0.5 * (lo + hi)
            a1, a2, a3 = pm.run_alphas_from_unification(alpha_u, mz_run, mu_u)
            return mz_run, v_ev, a1, a2, a3
        prev_mu, prev_f = mu, float(val)
    raise RuntimeError("wide window: could not bracket the mZ fixed point")


def diagnose_chain(
    b_id: str, beta_id: str, p_value: float, e_star: float
) -> dict[str, object]:
    """Quantify why a menu chain failed. Diagnostic only; entries stay failed.

    Reports (a) the alpha_3 Landau ceiling on alpha_U for the entry's gauge
    triple, (b) a wide-window alpha_U scan of the lane's pixel residual, and
    (c) if a wide-window root exists, the implied forward (MW, MZ) displays.
    """

    import numpy as np

    n_c = BETA_AXIS[beta_id]
    beta_exp = pm.beta_ew(n_c)
    mu_u = pm.unification_scale_gev(p_value)
    triple = B_TRIPLES[b_id]
    diag: dict[str, object] = {
        "gauge_triple": list(triple),
        "beta_EW_exponent": beta_exp,
        "mu_u_gev": mu_u,
    }

    # alpha_3 (and alpha_2) positivity ceiling at mu = 91.19 GeV reference.
    log_ratio = math.log(mu_u / 91.19)
    ceilings = {}
    for name, b in zip(("alpha1", "alpha2", "alpha3"), triple):
        if b < 0.0:
            inv_u_min = (-b / (2.0 * math.pi)) * log_ratio
            ceilings[name] = {
                "inverse_alpha_U_must_exceed": inv_u_min,
                "alpha_U_ceiling": 1.0 / inv_u_min if inv_u_min > 0 else None,
            }
    diag["positivity_ceilings_at_91GeV"] = ceilings

    original = pm.B_MSSM
    pm.B_MSSM = triple
    scan = []
    root_entry: dict[str, object] | None = None
    try:
        prev = None
        for alpha_u in np.linspace(0.005, 0.15, 146):
            alpha_u = float(alpha_u)
            row: dict[str, object] = {"alpha_U": alpha_u}
            try:
                mz_run, v_ev, a1, a2, a3 = _wide_mz_fixed_point(
                    alpha_u, p_value, n_c, mu_u
                )
                residual = pm.pixel_residual(a2, a3, p_value)
                if not math.isfinite(residual):
                    raise OverflowError("pixel residual not finite")
                row.update(
                    {
                        "mz_fixed_point_gev": mz_run,
                        "pixel_residual": residual,
                        "alpha3_at_fixed_point": a3,
                    }
                )
                if prev is not None and prev[1] * residual < 0:
                    lo, hi = prev[0], alpha_u
                    for _ in range(80):
                        mid = 0.5 * (lo + hi)
                        _, _, _, a2m, a3m = _wide_mz_fixed_point(
                            mid, p_value, n_c, mu_u
                        )
                        if pm.pixel_residual(a2m, a3m, p_value) * prev[1] > 0:
                            lo = mid
                        else:
                            hi = mid
                    root = 0.5 * (lo + hi)
                    mz_r, v_r, a1_r, a2_r, a3_r = _wide_mz_fixed_point(
                        root, p_value, n_c, mu_u
                    )
                    ay = 0.6 * a1_r
                    eta = root * (1.0 - 2.0 * ay / (ay + a2_r))
                    v_over_e = math.exp(
                        -2.0 * math.pi / (beta_exp * root)
                    ) / math.sqrt(p_value)
                    two = strict_branch_two_law_evaluation(a2_r, ay, eta, v_over_e)
                    root_entry = {
                        "alpha_U_root_wide_window": root,
                        "mz_fixed_point_gev": mz_r,
                        "zero_selector_MW_GeV": two["zero_selector_law"][
                            "MW_over_E_star"
                        ]
                        * e_star,
                        "zero_selector_MZ_GeV": two["zero_selector_law"][
                            "MZ_over_E_star"
                        ]
                        * e_star,
                        "nonzero_carrier_MW_GeV": two["nonzero_carrier_law"][
                            "MW_over_E_star"
                        ]
                        * e_star,
                        "nonzero_carrier_MZ_GeV": two["nonzero_carrier_law"][
                            "MZ_over_E_star"
                        ]
                        * e_star,
                    }
                prev = (alpha_u, residual)
            except OverflowError:
                row["status"] = "overflow_coupling_past_landau_pole"
                prev = None
            except (RuntimeError, ValueError) as exc:
                row["status"] = f"{type(exc).__name__}: {exc}"
                prev = None
            scan.append(row)
    finally:
        pm.B_MSSM = original

    statuses = [r.get("status") for r in scan]
    n_valid = sum(1 for s in statuses if s is None)
    diag["alpha_U_scan_window"] = [0.005, 0.15]
    diag["alpha_U_scan_points"] = len(scan)
    diag["alpha_U_scan_valid_points"] = n_valid
    diag["alpha_U_scan_overflow_points"] = sum(
        1 for s in statuses if s == "overflow_coupling_past_landau_pole"
    )
    valid = [r for r in scan if r.get("status") is None]
    if valid:
        diag["pixel_residual_range_on_valid_points"] = [
            min(r["pixel_residual"] for r in valid),
            max(r["pixel_residual"] for r in valid),
        ]
        diag["valid_alpha_U_interval"] = [
            min(r["alpha_U"] for r in valid),
            max(r["alpha_U"] for r in valid),
        ]
    diag["wide_window_root"] = root_entry
    if root_entry is None:
        diag["verdict"] = (
            "no pixel-residual root exists for this triple/beta pair even on "
            "the widened diagnostic windows; the failure is structural, not "
            "a solver-window artifact"
        )
    else:
        diag["verdict"] = (
            "a wide-window root exists outside the lane's declared solver "
            "window; the implied masses are recorded above and are far from "
            "the targets; the entry remains failed per the preregistered rule"
        )
    return diag


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--menu", type=Path, default=MENU)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    started = time.time()
    menu = json.loads(args.menu.read_text(encoding="utf-8"))
    cert = json.loads(SOURCE_CERT.read_text(encoding="utf-8"))
    clock = json.loads(CLOCK_CERT.read_text(encoding="utf-8"))
    boundary = json.loads(BOUNDARY_ARTIFACT.read_text(encoding="utf-8"))
    epsilon_cs = float(clock["candidate_values_from_sources"]["epsilon_Cs"])
    e_star = e_star_display_gev(epsilon_cs)

    p_axis = {
        item["id"]: float(item["P"])
        for item in menu["axes"]["R4_evaluation_pixel"]["values"]
    }

    chains: dict[tuple[str, str, str], dict[str, object]] = {}
    chain_errors: dict[tuple[str, str, str], str] = {}
    for b_id in B_TRIPLES:
        for beta_id in BETA_AXIS:
            for p_id, p_value in p_axis.items():
                key = (b_id, beta_id, p_id)
                try:
                    chains[key] = run_chain(b_id, beta_id, p_value)
                except Exception as exc:  # bracket loss is a recorded outcome
                    chain_errors[key] = f"{type(exc).__name__}: {exc}"

    failure_diagnoses: dict[str, object] = {}
    for key, error in chain_errors.items():
        b_id, beta_id, p_id = key
        diag_key = "|".join(key)
        diag = diagnose_chain(b_id, beta_id, p_axis[p_id], e_star)
        diag["lane_solver_exception"] = error
        failure_diagnoses[diag_key] = diag

    # Baseline consistency block, evaluated before any verdict is read.
    baseline_key = ("mssm", "beta4", "P_sl3")
    baseline_report: dict[str, object] = {"chain_key": "|".join(baseline_key)}
    if baseline_key in chains:
        chain = chains[baseline_key]
        recorded_basis = boundary["strict_branch_two_law_boundary"]["basis"]
        law0 = chain["two_law"]["zero_selector_law"]
        law1 = chain["two_law"]["nonzero_carrier_law"]
        recorded0 = boundary["two_law_display"]["zero_selector_law"]
        recorded1 = boundary["two_law_display"]["nonzero_carrier_law"]
        baseline_report.update(
            {
                "alpha_u_recomputed": chain["alpha_u"],
                "alpha_u_certificate": float(cert["alpha_U_P_cand"]),
                "alpha2_mz_recomputed_minus_recorded": chain["alpha2_mz"]
                - float(recorded_basis["alpha2_mz"]),
                "alphaY_mz_recomputed_minus_recorded": chain["alphaY_mz"]
                - float(recorded_basis["alphaY_mz"]),
                "eta_source_recomputed_minus_recorded": chain["eta_source"]
                - float(recorded_basis["eta_source"]),
                "v_over_E_star_transmutation": chain["v_over_E_star"],
                "v_over_E_star_certificate": float(cert["v_over_E_star_P_cand"]),
                "MW_zero_selector_recomputed_GeV": law0["MW_over_E_star"] * e_star,
                "MW_zero_selector_recorded_GeV": recorded0["MW_GeV"],
                "MZ_zero_selector_recomputed_GeV": law0["MZ_over_E_star"] * e_star,
                "MZ_zero_selector_recorded_GeV": recorded0["MZ_GeV"],
                "MW_nonzero_recomputed_GeV": law1["MW_over_E_star"] * e_star,
                "MW_nonzero_recorded_GeV": recorded1["MW_GeV"],
                "MZ_nonzero_recomputed_GeV": law1["MZ_over_E_star"] * e_star,
                "MZ_nonzero_recorded_GeV": recorded1["MZ_GeV"],
            }
        )

    entries = []
    landings = []
    for b_id in B_TRIPLES:
        for beta_id in BETA_AXIS:
            for p_id in p_axis:
                key = (b_id, beta_id, p_id)
                for law_id in LAW_AXIS:
                    for z_id in Z_AXIS:
                        entry_id = (
                            f"B={b_id}|beta={beta_id}|P={p_id}"
                            f"|law={law_id}|z={z_id}"
                        )
                        is_baseline = key == baseline_key and (
                            law_id == "zero_selector" and z_id == "z_tree"
                        )
                        if key in chain_errors:
                            entries.append(
                                {
                                    "id": entry_id,
                                    "baseline": is_baseline,
                                    "status": "chain_failed",
                                    "error": chain_errors[key],
                                    "diagnosis_key": "|".join(key),
                                    "landing": False,
                                }
                            )
                            continue
                        chain = chains[key]
                        law = chain["two_law"][f"{law_id}_law"] if law_id == "zero_selector" else chain["two_law"]["nonzero_carrier_law"]
                        mw = law["MW_over_E_star"] * e_star
                        mz = law["MZ_over_E_star"] * e_star
                        if z_id == "z_stage3":
                            mz = mz / math.sqrt(1.0 + DELTA_RHO_STAGE3)
                        mw_off = mw - MW_TARGET
                        mz_off = mz - MZ_TARGET
                        landing = (
                            abs(mw_off) <= MW_SIGMA and abs(mz_off) <= MZ_SIGMA
                        )
                        entry = {
                            "id": entry_id,
                            "baseline": is_baseline,
                            "status": "evaluated",
                            "alpha_u": chain["alpha_u"],
                            "v_over_E_star": chain["v_over_E_star"],
                            "MW_GeV": mw,
                            "MZ_GeV": mz,
                            "MW_offset_GeV": mw_off,
                            "MZ_offset_GeV": mz_off,
                            "MW_pull_sigma": mw_off / MW_SIGMA,
                            "MZ_pull_sigma": mz_off / MZ_SIGMA,
                            "landing": landing,
                        }
                        entries.append(entry)
                        if landing:
                            landings.append(entry_id)

    evaluated = [e for e in entries if e["status"] == "evaluated"]
    failed = [e for e in entries if e["status"] == "chain_failed"]
    wall = time.time() - started

    result = {
        "artifact": "oph_ew_repair_preregistered_results",
        "date": "2026-07-14",
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "ledger_row": "CL-5",
        "menu": str(args.menu.name),
        "menu_size_total": len(entries),
        "targets_GeV": {
            "MW": [MW_TARGET, MW_SIGMA],
            "MZ": [MZ_TARGET, MZ_SIGMA],
        },
        "display_adapter": {
            "epsilon_Cs": clock["candidate_values_from_sources"]["epsilon_Cs"],
            "E_star_display_GeV": e_star,
        },
        "baseline_consistency": baseline_report,
        "failure_diagnoses": failure_diagnoses,
        "entries": entries,
        "summary": {
            "entries_evaluated": len(evaluated),
            "entries_chain_failed": len(failed),
            "landings": landings,
            "landing_count": len(landings),
            "verdict": (
                "at least one preregistered revision lands; its menu selection "
                "is counted and independent confirmation is the next frozen test"
                if landings
                else "no menu entry lands; the CL-5 suspect list is exhausted "
                "at one loop and the resolution moves to structure this menu "
                "does not contain"
            ),
        },
        "wall_time_seconds": wall,
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=False) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "output": str(args.output),
                "entries": len(entries),
                "chain_failures": len(failed),
                "landing_count": len(landings),
                "landings": landings,
                "wall_time_seconds": round(wall, 2),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
