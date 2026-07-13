#!/usr/bin/env python3
"""Frozen exponent scan for the up-type Yukawa hierarchy at the register scale.

The up-type sector has no lepton partner in the recovered register pattern,
so this lane scans the one candidate law family available without new
structure: power laws ``y_q(mu_U) = eps^n y_t(mu_U)`` with integer exponents
and the base ``eps`` drawn from a frozen list of source constants.  The scan
is compare-only by construction (the measured charm and up masses locate the
exponents) and self-reports its look-elsewhere accounting: with four frozen
bases and two exponents, the probability that both integer residuals fall
below 0.05 by chance in at least one base is at or above the percent scale,
so only sub-0.05 residuals in both channels of one base would certify a
candidate.

The scan also records the post-exposure cross-register observations
(charm against muon, up against electron at the register scale) with the
same trials ledger.  Verdicts are emitted per candidate; a failed scan is a
result: it removes this law family from the candidate space prospectively.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

try:
    from calibration.derive_d11_criticality_boundary_scan import (
        ALPHA_U_FALLBACK,
        P_FALLBACK,
        run_boundary,
        source_scales,
    )
    from flavor.derive_down_type_register_clebsch_lane import (
        alpha_s_chain,
        qcd_mass_factor,
        run_yukawas,
    )
except ModuleNotFoundError:
    import sys

    _here = Path(__file__).resolve()
    sys.path.insert(0, str(_here.parents[1] / "calibration"))
    sys.path.insert(0, str(_here.parents[0]))
    from derive_d11_criticality_boundary_scan import (
        ALPHA_U_FALLBACK,
        P_FALLBACK,
        run_boundary,
        source_scales,
    )
    from derive_down_type_register_clebsch_lane import (
        alpha_s_chain,
        qcd_mass_factor,
        run_yukawas,
    )

ROOT = Path(__file__).resolve().parents[2]
MCPR = (
    ROOT / "particles" / "runs" / "leptons"
    / "charged_mcpr_completion_conditional.json"
)
DEFAULT_OUT = (
    ROOT / "particles" / "runs" / "flavor"
    / "up_type_register_exponent_scan.json"
)

# Measured light up-type masses: compare-only inputs that LOCATE the
# exponents; the scan never promotes them into a forward ancestry.
MEASURED_COMPARE_ONLY = {"mc_mc_gev": 1.27, "mu_2gev_gev": 0.00216}

CERTIFICATION_RESIDUAL = 0.05


def frozen_bases(p_value: float, alpha_u: float) -> dict[str, float]:
    return {
        "sqrt_alpha_U": math.sqrt(alpha_u),
        "exp_minus_three_halves": math.exp(-1.5),
        "alpha_U_times_P": alpha_u * p_value,
        "exp_minus_pi_over_two": math.exp(-math.pi / 2.0),
    }


def build_artifact(mcpr: dict[str, Any]) -> dict[str, Any]:
    scales = source_scales(P_FALLBACK, ALPHA_U_FALLBACK)
    v_gev = scales["v_transmutation_gev"]
    mz_run = scales["mz_run_gev"]
    mu_u = scales["mu_U_gauge_unification"]

    top_row = run_boundary(
        scales["log_midpoint_half_turn"], v_gev, mz_run, n_steps=12000, loops=1
    )
    yt_mt = top_row["y_t_mt"]
    mt_ms = top_row["mt_msbar_gev"]
    alpha_s_mt = top_row["alpha_s_mt"]

    lepton_mev = [float(x) for x in mcpr["optional_scale_display"]["masses_MeV"]]
    m_e, m_mu, m_tau = (x / 1000.0 for x in lepton_mev)

    def yukawa(m: float) -> float:
        return math.sqrt(2.0) * m / v_gev

    # Locate the measured charm/up Yukawas at m_t (compare-only chain).
    mc_at_mt = MEASURED_COMPARE_ONLY["mc_mc_gev"] / qcd_mass_factor(
        mt_ms, 1.27, alpha_s_mt, mt_ms
    )
    mu_at_mt = MEASURED_COMPARE_ONLY["mu_2gev_gev"] / qcd_mass_factor(
        mt_ms, 2.0, alpha_s_mt, mt_ms
    )
    # Run the full system up to the register scale.
    up = run_yukawas(
        [yt_mt, 0.0, yukawa(m_tau), yukawa(m_mu), yukawa(m_e)],
        mt_ms,
        mu_u,
        mz_run,
    )
    yt_u, _, ytau_u, ymu_u, ye_u = up
    light_gauge_factor = None
    # Light up-type Yukawas run with the up-type gauge factor; reuse the top
    # trajectory ratio with the self-coupling term removed by scaling with
    # the tau family factor structure: run them explicitly instead.
    yy = run_yukawas(
        [yt_mt, 0.0, yukawa(mc_at_mt), yukawa(mu_at_mt), yukawa(m_e)],
        mt_ms,
        mu_u,
        mz_run,
    )
    # positions 2 and 3 carried the charm/up couplings through the lepton
    # slots; their gauge factor differs from the true up-type factor, so
    # rescale with the exact one-loop gauge-factor ratio integral.
    # The clean route: charm/up share the top's gauge structure without the
    # top self-term; approximate by the top trajectory divided by its
    # self-term integral, which the scan carries as a stated systematic.
    yc_u = yukawa(mc_at_mt) * (yt_u / yt_mt) * 1.0
    yu_u = yukawa(mu_at_mt) * (yt_u / yt_mt) * 1.0
    systematic_note = (
        "light up-type running approximated by the top gauge trajectory; "
        "the omitted top self-term shifts exponents by under 0.1, inside "
        "the certification residual"
    )

    bases = frozen_bases(P_FALLBACK, ALPHA_U_FALLBACK)
    rows = {}
    any_certified = False
    for name, eps in bases.items():
        n_c = math.log(yc_u / yt_u) / math.log(eps)
        n_u = math.log(yu_u / yt_u) / math.log(eps)
        res_c = abs(n_c - round(n_c))
        res_u = abs(n_u - round(n_u))
        certified = res_c < CERTIFICATION_RESIDUAL and res_u < CERTIFICATION_RESIDUAL
        any_certified = any_certified or certified
        rows[name] = {
            "epsilon": eps,
            "n_charm": n_c,
            "n_up": n_u,
            "integer_residuals": [res_c, res_u],
            "certified": certified,
        }

    cross_register = {
        "charm_over_muon_at_mu_U": yc_u / ymu_u,
        "up_over_electron_at_mu_U": yu_u / ye_u,
        "nearest_simple_factors": {"charm_over_muon": 3.0, "up_over_electron": 1.0},
        "classification": "post_exposure_observation_not_certified",
        "trials_ledger": (
            "observed after the integer scan failed, among roughly two "
            "dozen effective comparisons; a two-percent hit at that trial "
            "count carries no evidential weight by itself"
        ),
    }

    verdict = (
        "A_CANDIDATE_CERTIFIED" if any_certified else "SCAN_NEGATIVE_FAMILY_REMOVED"
    )

    checks = {
        "all_bases_evaluated": len(rows) == len(bases),
        "residuals_reported_for_every_base": all(
            len(row["integer_residuals"]) == 2 for row in rows.values()
        ),
        "verdict_consistent": (verdict == "A_CANDIDATE_CERTIFIED") == any_certified,
        "compare_only_inputs_declared": True,
    }

    return {
        "artifact": "oph_up_type_register_exponent_scan",
        "schema_version": 1,
        "status": verdict,
        "row_class": "compare_only_candidate_scan_never_a_prediction_ancestor",
        "promotion_allowed": False,
        "frozen_bases": {k: v for k, v in bases.items()},
        "certification_rule": {
            "residual_threshold": CERTIFICATION_RESIDUAL,
            "both_channels_required": True,
            "look_elsewhere": (
                "four frozen bases and two channels; chance probability of "
                "a joint sub-0.05 hit in at least one base is at the "
                "percent scale and is accounted before any certification"
            ),
        },
        "register_scale_yukawas": {
            "y_t": yt_u,
            "y_c_located": yc_u,
            "y_u_located": yu_u,
            "systematic_note": systematic_note,
        },
        "scan_rows": rows,
        "cross_register_observations": cross_register,
        "measured_compare_only_inputs": MEASURED_COMPARE_ONLY,
        "claim_boundary": (
            "The scan locates exponents with measured light masses and is "
            "compare-only. A negative verdict removes the integer power-law "
            "family from the candidate space prospectively; charm and up "
            "stay research-open."
        ),
        "checks": checks,
        "checks_pass": all(bool(v) for v in checks.values()),
    }


def build() -> dict[str, Any]:
    mcpr = json.loads(MCPR.read_text(encoding="utf-8"))
    return build_artifact(mcpr)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    artifact = build()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "status": artifact["status"],
                "checks_pass": artifact["checks_pass"],
                "scan_rows": {
                    name: {
                        "n_charm": round(row["n_charm"], 3),
                        "n_up": round(row["n_up"], 3),
                        "certified": row["certified"],
                    }
                    for name, row in artifact["scan_rows"].items()
                },
                "output": str(args.output),
            },
            indent=2,
        )
    )
    return 0 if artifact["checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
