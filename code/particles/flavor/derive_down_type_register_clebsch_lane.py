#!/usr/bin/env python3
"""Down-type quark masses from the charged leptons via register Clebsch factors.

The recovered quotient carries a unification register at ``mu_U`` where the
gauge couplings meet.  The historical candidate pairs the down-type quarks
with the charged leptons multiplet-wise and assigns the Clebsch factors
``(1, 1/3, 3)`` to ``(b/tau, s/mu, d/e)``:

    y_b(mu_U) = y_tau(mu_U),
    y_s(mu_U) = y_mu(mu_U)/3,
    y_d(mu_U) = 3 y_e(mu_U).

The candidate register origin of the factors is declared: ``3 = N_c`` and
``1/3`` is the invariant color measure of the transitive C3 action, the two
weights the corpus emits everywhere.  The selection artifact
(``clebsch_register_pairing_selection.json``) closes the pairing part
conditional on the #314 exterior matter contract and the unordered weight
set conditional on the frozen constraints F1/F2; the generation-order
assignment of the weights, ``GENERATION_REGISTER_ORDER``, remains open.  The
exhaustive six-assignment audit below is retrospective and target-informed;
it cannot discharge that premise.

The lane runs a five-Yukawa one-loop approximation (top, bottom, and three
charged leptons; strange and down share the computed bottom transport factor),
with top from the declared criticality boundary, applies the Clebsch boundary at
``mu_U``, runs down, and converts with flavor-banded QCD mass running.
The light ratio is protected by the common transport used here and is the
sharpest diagnostic of the route.  All six assignments are rejected by the
retrospective conservative FLAG gate.  ``sqrt(m_d/m_s)`` is retained only as
a Gatto-Sartori-Tonin texture diagnostic: the diagonal register ansatz itself
has a common eigenbasis and therefore gives the identity CKM matrix.  Measured
values appear only in the compare-only block.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any

try:
    from calibration.derive_d11_criticality_boundary_scan import (
        ALPHA_U_FALLBACK,
        P_FALLBACK,
        gauge_couplings,
        run_boundary,
        source_scales,
    )
except ModuleNotFoundError:
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "calibration"))
    from derive_d11_criticality_boundary_scan import (
        ALPHA_U_FALLBACK,
        P_FALLBACK,
        gauge_couplings,
        run_boundary,
        source_scales,
    )

ROOT = Path(__file__).resolve().parents[2]
MCPR = (
    ROOT / "particles" / "runs" / "leptons"
    / "charged_mcpr_completion_conditional.json"
)
DEFAULT_OUT = (
    ROOT / "particles" / "runs" / "flavor"
    / "down_type_register_clebsch_lane.json"
)
SELECTION = (
    ROOT / "particles" / "runs" / "flavor"
    / "clebsch_register_pairing_selection.json"
)
FLAG_FIXTURE = (
    ROOT / "particles" / "data"
    / "flag_2024_light_quark_ratio_fixture.json"
)
CALIBRATION_PRODUCER = (
    ROOT / "particles" / "calibration"
    / "derive_d11_criticality_boundary_scan.py"
)

CLEBSCH = {"b_over_tau": 1.0, "s_over_mu": 1.0 / 3.0, "d_over_e": 3.0}
ASSIGNMENT_SLOTS = ("b_over_tau", "s_over_mu", "d_over_e")
FACTOR_SET = (Fraction(1), Fraction(1, 3), Fraction(3))

# Compare-only references (never prediction ancestors).
COMPARE_ONLY = {
    "mb_mb_gev": 4.18,
    "ms_2gev_gev": 0.0935,
    "md_2gev_gev": 0.0047,
    "cabibbo": 0.2250,
}


def _yukawa_derivatives(y: list[float], mu: float, mz_run: float) -> list[float]:
    """One-loop SM Yukawa system [yt, yb, ytau, ymu, ye]."""

    g_y, g2, g3 = gauge_couplings(mu, mz_run)
    g1sq = (5.0 / 3.0) * g_y * g_y
    gauge_up = (17.0 / 20.0) * g1sq + 2.25 * g2 * g2 + 8.0 * g3 * g3
    gauge_down = 0.25 * g1sq + 2.25 * g2 * g2 + 8.0 * g3 * g3
    gauge_lepton = 2.25 * g1sq + 2.25 * g2 * g2
    yt, yb, ytau, ymu, ye = y
    trace = 3.0 * yt * yt + 3.0 * yb * yb + ytau * ytau + ymu * ymu + ye * ye
    kappa = 1.0 / (16.0 * math.pi**2)
    return [
        yt * kappa * (1.5 * yt * yt - 1.5 * yb * yb + trace - gauge_up),
        yb * kappa * (1.5 * yb * yb - 1.5 * yt * yt + trace - gauge_down),
        ytau * kappa * (1.5 * ytau * ytau + trace - gauge_lepton),
        ymu * kappa * (1.5 * ymu * ymu + trace - gauge_lepton),
        ye * kappa * (1.5 * ye * ye + trace - gauge_lepton),
    ]


def run_yukawas(
    y0: list[float], mu0: float, mu1: float, mz_run: float, n_steps: int = 12000
) -> list[float]:
    y = list(y0)
    t0, t1 = math.log(mu0), math.log(mu1)
    dt = (t1 - t0) / n_steps
    for i in range(n_steps):
        t = t0 + i * dt
        k1 = _yukawa_derivatives(y, math.exp(t), mz_run)
        y2 = [a + 0.5 * dt * b for a, b in zip(y, k1)]
        k2 = _yukawa_derivatives(y2, math.exp(t + 0.5 * dt), mz_run)
        y3 = [a + 0.5 * dt * b for a, b in zip(y, k2)]
        k3 = _yukawa_derivatives(y3, math.exp(t + 0.5 * dt), mz_run)
        y4 = [a + dt * b for a, b in zip(y, k3)]
        k4 = _yukawa_derivatives(y4, math.exp(t + dt), mz_run)
        y = [
            a + dt / 6.0 * (b1 + 2 * b2 + 2 * b3 + b4)
            for a, b1, b2, b3, b4 in zip(y, k1, k2, k3, k4)
        ]
    return y


def alpha_s_chain(mu: float, alpha_s_mt: float, mt: float) -> float:
    """One-loop flavor-banded strong coupling below m_t."""

    def run(a0: float, mu0: float, mu1: float, nf: int) -> float:
        b0 = 11.0 - 2.0 * nf / 3.0
        return 1.0 / (1.0 / a0 + b0 / (2.0 * math.pi) * math.log(mu1 / mu0))

    if mu >= 4.18:
        return run(alpha_s_mt, mt, mu, 5)
    a_mb = run(alpha_s_mt, mt, 4.18, 5)
    if mu >= 1.27:
        return run(a_mb, 4.18, mu, 4)
    a_mc = run(a_mb, 4.18, 1.27, 4)
    return run(a_mc, 1.27, mu, 3)


def qcd_mass_factor(mu0: float, mu1: float, alpha_s_mt: float, mt: float) -> float:
    """One-loop mass running factor across flavor bands from mu0 to mu1."""

    def gamma_exponent(nf: int) -> float:
        return 4.0 / (11.0 - 2.0 * nf / 3.0)

    factor = 1.0
    grid = [(4.18, 5), (1.27, 4)]
    current = mu0
    a_current = alpha_s_chain(current, alpha_s_mt, mt)
    for threshold, nf in grid:
        if mu1 < threshold < current:
            a_th = alpha_s_chain(threshold, alpha_s_mt, mt)
            factor *= (a_th / a_current) ** gamma_exponent(nf)
            current, a_current = threshold, a_th
    nf_final = 5 if mu1 >= 4.18 else (4 if mu1 >= 1.27 else 3)
    a_end = alpha_s_chain(mu1, alpha_s_mt, mt)
    factor *= (a_end / a_current) ** gamma_exponent(nf_final)
    return factor


def _load_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"expected JSON object: {path}")
    return payload


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _payload_sha256(payload: dict[str, Any]) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _factor_label(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else str(value)


def _declared_fraction(value: float | int | Fraction) -> Fraction:
    return Fraction(value).limit_denominator(3)


def _validate_selection(selection: dict[str, Any]) -> None:
    scan = selection.get("weight_set_scan", {})
    pairing = selection.get("pairing_theorem", {})
    if (
        selection.get("artifact")
        != "oph_clebsch_register_pairing_selection"
        or selection.get("promotion_allowed") is not False
        or scan.get("surviving_weight_set") != ["1/3", "1", "3"]
        or scan.get("unique_unordered_survivor") is not True
        or scan.get("order_assignment", {}).get("status") != "open"
        or pairing.get("classification")
        != "conditional_channel_compatibility_only"
        or pairing.get("register_relation_existence_proved") is not False
        or pairing.get("independent_yukawa_coefficients_equated") is not False
    ):
        raise ValueError(
            "Clebsch selection input does not preserve conditional channel "
            "compatibility, the unordered F1/F2 lemma, and open order"
        )


def _validate_flag_fixture(fixture: dict[str, Any]) -> None:
    boundary = fixture.get("claim_boundary", {})
    averages = fixture.get("averages", [])
    if (
        fixture.get("artifact")
        != "oph_flag_2024_light_quark_ratio_fixture"
        or fixture.get("schema")
        != "oph.flag_2024_light_quark_ratio_fixture.v1"
        or boundary.get("comparison_only") is not True
        or boundary.get("oph_fit_or_selection_input") is not False
        or boundary.get("oph_theory_uncertainty_supplied") is not False
        or boundary.get("significance_gate_preregistered") is not False
        or len(averages) != 2
    ):
        raise ValueError("FLAG fixture crossed its compare-only claim boundary")
    for row in averages:
        a = float(row["ms_over_mud"]["value"])
        sigma_a = float(row["ms_over_mud"]["standard_uncertainty"])
        u = float(row["mu_over_md"]["value"])
        sigma_u = float(row["mu_over_md"]["standard_uncertainty"])
        central = a * (1.0 + u) / 2.0
        independent = math.hypot(
            (1.0 + u) / 2.0 * sigma_a,
            a / 2.0 * sigma_u,
        )
        rho_plus_one = (
            (1.0 + u) / 2.0 * sigma_a + a / 2.0 * sigma_u
        )
        derived = row.get("derived_ms_over_md", {})
        if not (
            math.isclose(
                float(derived.get("value", "nan")),
                central,
                rel_tol=0.0,
                abs_tol=1e-13,
            )
            and math.isclose(
                float(
                    derived.get(
                        "independent_standard_uncertainty", "nan"
                    )
                ),
                independent,
                rel_tol=0.0,
                abs_tol=1e-13,
            )
            and math.isclose(
                float(
                    derived.get(
                        "rho_plus_one_standard_uncertainty", "nan"
                    )
                ),
                rho_plus_one,
                rel_tol=0.0,
                abs_tol=1e-13,
            )
        ):
            raise ValueError(
                f"FLAG derived ms/md fields drifted for Nf={row.get('nf')}"
            )


def _flag_compare(
    ms_over_md: float,
    fixture: dict[str, Any],
) -> list[dict[str, Any]]:
    rows = []
    threshold = float(
        fixture["derived_quantity"]["uncertainty_policy"][
            "conservative_rejection_threshold_sigma"
        ]
    )
    for reference in fixture["averages"]:
        derived = reference["derived_ms_over_md"]
        central = float(derived["value"])
        sigma_independent = float(
            derived["independent_standard_uncertainty"]
        )
        sigma_rho_plus_one = float(
            derived["rho_plus_one_standard_uncertainty"]
        )
        gap = abs(ms_over_md - central)
        conservative_gap_sigma = gap / sigma_rho_plus_one
        rows.append(
            {
                "nf": reference["nf"],
                "reference_ms_over_md": central,
                "independent_standard_uncertainty": sigma_independent,
                "rho_plus_one_standard_uncertainty": sigma_rho_plus_one,
                "absolute_relative_error": abs(ms_over_md / central - 1.0),
                "independent_gap_sigma": gap / sigma_independent,
                "rho_plus_one_gap_sigma": conservative_gap_sigma,
                "conservative_rejection_threshold_sigma": threshold,
                "conservative_rejection_triggered": (
                    conservative_gap_sigma >= threshold
                ),
            }
        )
    return rows


def build_artifact(
    mcpr: dict[str, Any],
    selection: dict[str, Any] | None = None,
    flag_fixture: dict[str, Any] | None = None,
    *,
    input_hashes: dict[str, str] | None = None,
    clebsch: dict[str, float] | None = None,
) -> dict[str, Any]:
    """Build the lane with injectable dependency artifacts for adversarial tests."""

    selection_injected = selection is not None
    flag_injected = flag_fixture is not None
    if selection is None:
        selection = _load_json(SELECTION)
    if flag_fixture is None:
        flag_fixture = _load_json(FLAG_FIXTURE)
    _validate_selection(selection)
    _validate_flag_fixture(flag_fixture)

    if input_hashes is None:
        input_hashes = {
            "mcpr": _payload_sha256(mcpr),
            "clebsch_selection": (
                _payload_sha256(selection)
                if selection_injected
                else _sha256(SELECTION)
            ),
            "flag_2024_fixture": (
                _payload_sha256(flag_fixture)
                if flag_injected
                else _sha256(FLAG_FIXTURE)
            ),
            "declared_boundary_calibration_producer": _sha256(
                CALIBRATION_PRODUCER
            ),
        }

    assumed_factors = dict(CLEBSCH if clebsch is None else clebsch)
    if set(assumed_factors) != set(ASSIGNMENT_SLOTS):
        raise ValueError("Clebsch assignment slots are incomplete")
    if sorted(
        _declared_fraction(value) for value in assumed_factors.values()
    ) != sorted(
        FACTOR_SET
    ):
        raise ValueError("Clebsch assignment is not a permutation of {1/3,1,3}")

    if mcpr.get("charged_reference_data_consumed") is not False:
        raise ValueError("MCPR input no longer preserves its runtime target boundary")

    shape_squares = [
        float(value)
        for value in mcpr["regular_C3_shape"]["root_squares_sorted"]
    ]
    if len(shape_squares) != 3 or min(shape_squares) <= 0:
        raise ValueError("MCPR dimensionless shape is not a positive triple")
    shape_e, shape_mu, _shape_tau = shape_squares

    scales = source_scales(P_FALLBACK, ALPHA_U_FALLBACK)
    v_gev = scales["v_transmutation_gev"]
    mz_run = scales["mz_run_gev"]
    mu_u = scales["mu_U_gauge_unification"]
    top_row = run_boundary(
        scales["log_midpoint_half_turn"],
        v_gev,
        mz_run,
        n_steps=12000,
        loops=1,
    )
    yt_mt = top_row["y_t_mt"]
    mt_ms = top_row["mt_msbar_gev"]
    alpha_s_mt = top_row["alpha_s_mt"]

    display_masses = mcpr["optional_scale_display"]["masses_MeV"]
    m_e, m_mu, m_tau = [float(value) / 1000.0 for value in display_masses]

    def yukawa(mass: float) -> float:
        return math.sqrt(2.0) * mass / v_gev

    up = run_yukawas(
        [yt_mt, 0.0, yukawa(m_tau), yukawa(m_mu), yukawa(m_e)],
        mt_ms,
        mu_u,
        mz_run,
    )
    yt_u, _, ytau_u, ymu_u, ye_u = up

    def evaluate_assignment(factors: dict[str, float]) -> dict[str, Any]:
        yb_u = factors["b_over_tau"] * ytau_u
        ys_u = factors["s_over_mu"] * ymu_u
        yd_u = factors["d_over_e"] * ye_u
        down = run_yukawas(
            [yt_u, yb_u, ytau_u, ymu_u, ye_u],
            mu_u,
            mt_ms,
            mz_run,
        )
        yb_mt = down[1]
        common_light_down_factor = yb_mt / yb_u
        ys_mt = ys_u * common_light_down_factor
        yd_mt = yd_u * common_light_down_factor
        mb_at_mt = yb_mt * v_gev / math.sqrt(2.0)
        ms_at_mt = ys_mt * v_gev / math.sqrt(2.0)
        md_at_mt = yd_mt * v_gev / math.sqrt(2.0)

        mb_mb = mb_at_mt
        for _ in range(25):
            mb_mb = mb_at_mt * qcd_mass_factor(
                mt_ms, mb_mb, alpha_s_mt, mt_ms
            )
        qcd_2gev = qcd_mass_factor(mt_ms, 2.0, alpha_s_mt, mt_ms)
        ms_2 = ms_at_mt * qcd_2gev
        md_2 = md_at_mt * qcd_2gev
        ratio = ms_2 / md_2
        cabibbo_gst = math.sqrt(md_2 / ms_2)

        exact_factor_ratio = (
            _declared_fraction(factors["s_over_mu"])
            / _declared_fraction(factors["d_over_e"])
        )
        ratio_only_shape = (
            float(exact_factor_ratio) * shape_mu / shape_e
        )
        register_ratio = (
            float(exact_factor_ratio) * ymu_u / ye_u
        )
        flag_rows = _flag_compare(ratio, flag_fixture)
        absolute_relatives = {
            "mb_relative": mb_mb / COMPARE_ONLY["mb_mb_gev"] - 1.0,
            "ms_relative": ms_2 / COMPARE_ONLY["ms_2gev_gev"] - 1.0,
            "md_relative": md_2 / COMPARE_ONLY["md_2gev_gev"] - 1.0,
        }
        predictions = {
            "mb_mb_gev": mb_mb,
            "ms_2gev_gev": ms_2,
            "md_2gev_gev": md_2,
            "ms_over_md": ratio,
            "cabibbo_gst_sqrt_md_over_ms": cabibbo_gst,
        }
        return {
            "assignment": dict(factors),
            "assignment_labels": {
                name: _factor_label(_declared_fraction(value))
                for name, value in factors.items()
            },
            "ratio_only_shape_core": {
                "ms_over_md_at_declared_clebsch_boundary_without_running":
                    ratio_only_shape,
                "source": "MCPR /regular_C3_shape/root_squares_sorted",
                "independent_of_optional_scale_display": True,
            },
            "boundary_values_at_mu_U": {
                "y_tau": ytau_u,
                "y_mu": ymu_u,
                "y_e": ye_u,
                "y_b": yb_u,
                "y_s": ys_u,
                "y_d": yd_u,
                "y_t": yt_u,
            },
            "common_transport": {
                "approximation": (
                    "five-Yukawa one-loop system [yt,yb,ytau,ymu,ye]; "
                    "ys and yd share the computed bottom transport factor"
                ),
                "common_light_down_factor": common_light_down_factor,
                "same_scale_qcd_factor": qcd_2gev,
                "generation_dependent_light_threshold_present": False,
                "factor_ratio_s_over_d": _factor_label(exact_factor_ratio),
                "register_scale_ratio_identity_rhs": register_ratio,
                "transported_ms_over_md": ratio,
                "relative_identity_gap": abs(ratio / register_ratio - 1.0),
            },
            "predictions": predictions,
            "absolute_output_diagnostics": {
                **absolute_relatives,
                "worst_absolute_mass_relative_error": max(
                    abs(value) for value in absolute_relatives.values()
                ),
                "depends_on_optional_scale_display": True,
                "depends_on_declared_boundary_calibration": True,
                "status": "calibration-dependent diagnostic, not ratio-only core",
            },
            "flag_2024_compare_only": flag_rows,
            "retrospective_metric": {
                "worst_absolute_mass_relative_error": max(
                    abs(value) for value in absolute_relatives.values()
                ),
                "light_ratio_absolute_relative_error_vs_flag_nf_2+1+1":
                    next(
                        row["absolute_relative_error"]
                        for row in flag_rows
                        if row["nf"] == "2+1+1"
                    ),
            },
            "arithmetic_checks": {
                "masses_positive": min(mb_mb, ms_2, md_2) > 0,
                "outputs_finite": all(
                    math.isfinite(value) for value in predictions.values()
                ),
                "common_transport_ratio_identity": (
                    abs(ratio / register_ratio - 1.0) < 1e-13
                ),
                "flag_comparisons_finite": all(
                    math.isfinite(row["rho_plus_one_gap_sigma"])
                    for row in flag_rows
                ),
            },
        }

    permutation_rows = []
    for permutation in itertools.permutations(FACTOR_SET):
        factors = {
            name: float(value)
            for name, value in zip(
                ASSIGNMENT_SLOTS, permutation, strict=True
            )
        }
        row = evaluate_assignment(factors)
        row["is_current_assumed_order"] = all(
            _declared_fraction(row["assignment"][name])
            == _declared_fraction(assumed_factors[name])
            for name in ASSIGNMENT_SLOTS
        )
        permutation_rows.append(row)

    if len(permutation_rows) != 6:
        raise ValueError("Clebsch permutation scan is incomplete")
    score = lambda row: (
        row["retrospective_metric"][
            "worst_absolute_mass_relative_error"
        ],
        row["retrospective_metric"][
            "light_ratio_absolute_relative_error_vs_flag_nf_2+1+1"
        ],
    )
    minimum_score = min(score(row) for row in permutation_rows)
    winners = [row for row in permutation_rows if score(row) == minimum_score]
    if len(winners) != 1:
        raise ValueError("retrospective discrepancy metric is not unique")
    for row in permutation_rows:
        row["retrospective_unique_least_discrepant"] = row is winners[0]
        row["conservative_flag_rejected_for_all_nf_rows"] = all(
            flag_row["conservative_rejection_triggered"]
            for flag_row in row["flag_2024_compare_only"]
        )

    assumed_row = next(
        row for row in permutation_rows if row["is_current_assumed_order"]
    )
    all_miss = all(
        row["conservative_flag_rejected_for_all_nf_rows"]
        for row in permutation_rows
    )
    arithmetic_checks = {
        "selection_artifact_validated": True,
        "flag_fixture_validated": True,
        "six_permutations_executed": len(permutation_rows) == 6,
        "all_permutation_arithmetic_checks_pass": all(
            all(row["arithmetic_checks"].values())
            for row in permutation_rows
        ),
        "current_assumed_order_present_once": sum(
            row["is_current_assumed_order"] for row in permutation_rows
        )
        == 1,
    }
    physical_premises = {
        "register_relation_existence_discharged": False,
        "independent_yukawa_coefficient_identification_discharged": False,
        "generation_register_order_selected": False,
        "third_generation_register_factor_selected": False,
        "mcpr_physical_family_attachment_discharged": False,
        "full_rg_threshold_matching_packet_supplied": False,
    }
    flag_nf211 = next(
        row
        for row in assumed_row["flag_2024_compare_only"]
        if row["nf"] == "2+1+1"
    )
    compare = {
        **assumed_row["absolute_output_diagnostics"],
        "ms_over_md_relative": (
            assumed_row["predictions"]["ms_over_md"]
            / flag_nf211["reference_ms_over_md"]
            - 1.0
        ),
        "cabibbo_relative": (
            assumed_row["predictions"]["cabibbo_gst_sqrt_md_over_ms"]
            / COMPARE_ONLY["cabibbo"]
            - 1.0
        ),
        "references": {
            **COMPARE_ONLY,
            "flag_ms_over_md_nf_2+1+1":
                flag_nf211["reference_ms_over_md"],
        },
        "flag_2024": assumed_row["flag_2024_compare_only"],
        "role": (
            "retrospective comparison only, outside prediction ancestry; "
            "no OPH theory uncertainty supplied"
        ),
    }
    lepton_register_ratio = ymu_u / ye_u
    permitted_factor_ratios = sorted(
        {
            numerator / denominator
            for numerator in FACTOR_SET
            for denominator in FACTOR_SET
            if numerator != denominator
        }
    )
    required_factor_ratio_flag = (
        Fraction(str(flag_nf211["reference_ms_over_md"]))
        / Fraction(str(lepton_register_ratio))
    )
    nearest_factor_ratio = min(
        permitted_factor_ratios,
        key=lambda candidate: abs(
            math.log(float(candidate / required_factor_ratio_flag))
        ),
    )
    gst_reference_ratio = 1.0 / COMPARE_ONLY["cabibbo"] ** 2
    quantization_audit = {
        "status": "retrospective_compare_only",
        "lepton_register_ratio": lepton_register_ratio,
        "permitted_distinct_assignment_factor_ratios": [
            _factor_label(value) for value in permitted_factor_ratios
        ],
        "excluded_values_from_prior_informal_note": ["1", "27", "1/27"],
        "required_factor_ratio_for_flag_nf_2+1+1": float(
            required_factor_ratio_flag
        ),
        "nearest_permitted_factor_ratio": _factor_label(
            nearest_factor_ratio
        ),
        "predicted_over_flag_nf_2+1+1": (
            assumed_row["predictions"]["ms_over_md"]
            / flag_nf211["reference_ms_over_md"]
        ),
        "gst_compare_only": {
            "cabibbo_reference": COMPARE_ONLY["cabibbo"],
            "implied_ms_over_md": gst_reference_ratio,
            "predicted_over_implied": (
                assumed_row["predictions"]["ms_over_md"]
                / gst_reference_ratio
            ),
            "status": (
                "sqrt(md/ms) is a texture diagnostic restatement, not a "
                "derived CKM angle"
            ),
        },
        "diagonal_ansatz_ckm_consequence": {
            "simultaneously_diagonal_mass_matrices": True,
            "ckm_matrix": "identity",
            "nonzero_cabibbo_angle_derived": False,
        },
    }
    arithmetic_checks["distinct_factor_ratio_set_is_exact"] = (
        permitted_factor_ratios
        == [Fraction(1, 9), Fraction(1, 3), Fraction(3), Fraction(9)]
    )

    return {
        "artifact": "oph_down_type_register_clebsch_lane",
        "schema_version": 2,
        "status": "CONDITIONAL_DECLARED_ROUTE_RETROSPECTIVELY_REJECTED",
        "row_class": "conditional_on_mcpr_leptons_and_open_register_relation",
        "promotion_allowed": False,
        "dependency_audit": {
            "direct_input_paths": [
                "particles/runs/leptons/charged_mcpr_completion_conditional.json",
                "particles/runs/flavor/clebsch_register_pairing_selection.json",
                "particles/data/flag_2024_light_quark_ratio_fixture.json",
                "particles/calibration/derive_d11_criticality_boundary_scan.py",
            ],
            "input_sha256": input_hashes,
            "flag_is_compare_only_not_prediction_ancestor": True,
        },
        "clebsch_boundary": {
            "factors": assumed_factors,
            "register_scale": "mu_U",
            "register_scale_status": "declared boundary scale",
            "current_order_status": (
                "assumed open-premise assignment; not source-selected"
            ),
            "conditional_origin": (
                "The matter receipt supplies channel compatibility only. "
                "F1/F2 give a conditional unordered-set lemma; existence of "
                "the cross-sector relation, coefficient identification, and "
                "generation order remain open."
            ),
            "selection_artifact":
                "runs/flavor/clebsch_register_pairing_selection.json",
        },
        "inputs": {
            "lepton_shape_source":
                "runs/leptons/charged_mcpr_completion_conditional.json",
            "top_trajectory": (
                "declared criticality boundary; calibration-dependent"
            ),
            "transport_approximation": (
                "five-Yukawa one-loop system plus common light-down "
                "transport and one-loop flavor-banded QCD mass running"
            ),
        },
        "ratio_only_core": {
            "independent_of_optional_scale_display": True,
            "source_json_pointer":
                "/regular_C3_shape/root_squares_sorted",
            "current_assumed_order_value": assumed_row[
                "ratio_only_shape_core"
            ][
                "ms_over_md_at_declared_clebsch_boundary_without_running"
            ],
        },
        "absolute_output_provenance": {
            "depends_on_mcpr_optional_scale_display": True,
            "depends_on_declared_boundary_calibration": True,
            "calibration_producer_sha256": input_hashes[
                "declared_boundary_calibration_producer"
            ],
            "status": (
                "absolute masses are calibration-dependent diagnostics, "
                "not source-only predictions"
            ),
        },
        "boundary_values_at_mu_U": assumed_row[
            "boundary_values_at_mu_U"
        ],
        "predictions": assumed_row["predictions"],
        "common_transport": assumed_row["common_transport"],
        "compare_only": compare,
        "quantization_audit": quantization_audit,
        "permutation_scan": {
            "factor_set": ["1/3", "1", "3"],
            "assignment_slots": list(ASSIGNMENT_SLOTS),
            "permutation_count": 6,
            "exhaustive": True,
            "rows": permutation_rows,
            "retrospective_metric": {
                "order": [
                    "minimize worst absolute-mass relative error",
                    "then minimize light-ratio relative error vs FLAG Nf=2+1+1",
                ],
                "target_informed": True,
                "preregistered": False,
                "physical_order_selected": False,
                "current_assumed_order_uniquely_least_discrepant":
                    winners[0]["is_current_assumed_order"],
            },
            "all_permutations_rejected_by_conservative_flag_gate":
                all_miss,
        },
        "retrospective_flag_rejection": {
            "prediction_preexisted_audit": True,
            "significance_gate_preregistered": False,
            "comparison_is_retrospective": True,
            "input_covariance_available": False,
            "oph_theory_uncertainty_supplied": False,
            "rho_plus_one_conservative_gate_used": True,
            "current_assumed_order_rejected_for_all_flag_nf_rows":
                assumed_row["conservative_flag_rejected_for_all_nf_rows"],
            "all_six_permutations_rejected": all_miss,
            "interpretation": (
                "The conservative FLAG comparison rejects this current "
                "declared-model route. It is not a prospective selection test "
                "and does not exclude a modified generation-dependent "
                "threshold lane."
            ),
        },
        "normalization_tension": {
            "statement": (
                "The current assumed order carries large retrospective "
                "absolute-mass residuals. They are reported, not accepted by "
                "a proximity criterion."
            ),
            "open_objects": [
                "REGISTER_RELATION_EXISTENCE",
                "INDEPENDENT_YUKAWA_COEFFICIENT_IDENTIFICATION",
                "GENERATION_REGISTER_ORDER",
                "THIRD_GENERATION_REGISTER_FACTOR",
                "FROZEN_RG_THRESHOLD_MATCHING_PACKET",
            ],
        },
        "physical_premises": physical_premises,
        "physical_premises_discharged": all(physical_premises.values()),
        "claim_boundary": (
            "This is a conditional, declared-model diagnostic. Arithmetic "
            "reproduction does not discharge any physical premise. The "
            "retrospective FLAG gate rejects all six declared assignments; "
            "the uniquely least-discrepant current order is target-informed "
            "and selects no physical generation order."
        ),
        "arithmetic_checks": arithmetic_checks,
        "arithmetic_checks_pass": all(arithmetic_checks.values()),
        "checks": arithmetic_checks,
        "checks_pass": all(arithmetic_checks.values()),
    }


def build() -> dict[str, Any]:
    mcpr = _load_json(MCPR)
    selection = _load_json(SELECTION)
    flag_fixture = _load_json(FLAG_FIXTURE)
    return build_artifact(
        mcpr,
        selection,
        flag_fixture,
        input_hashes={
            "mcpr": _sha256(MCPR),
            "clebsch_selection": _sha256(SELECTION),
            "flag_2024_fixture": _sha256(FLAG_FIXTURE),
            "declared_boundary_calibration_producer": _sha256(
                CALIBRATION_PRODUCER
            ),
        },
    )


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
                "physical_premises_discharged": artifact[
                    "physical_premises_discharged"
                ],
                "predictions": artifact["predictions"],
                "compare_only_relatives": {
                    k: round(v, 4)
                    for k, v in artifact["compare_only"].items()
                    if isinstance(v, float)
                },
                "output": str(args.output),
            },
            indent=2,
        )
    )
    return 0 if artifact["checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
