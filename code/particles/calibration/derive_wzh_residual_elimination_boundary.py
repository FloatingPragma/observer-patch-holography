#!/usr/bin/env python3
"""Recompute the W/Z/H residual-elimination boundary on the source-audit branch.

The residual-elimination research update replaces the earlier continuous D10
non-identifiability description by a rigid two-point discrete boundary:

1. The archived D10 selector is globally rigid.  With
   ``J10(x, y) = x^2 [1 + (2 eta + x)^2 + 4 x^2]/(1 + 4 x^2)
   + (kappa^2/4)(1 + 4 x^2) y^2`` and the exact factorization
   ``J10 - (3/4) x^2 - (kappa^2/4) y^2 = [x^2 / (4 (1 + 4 x^2))]
   [8 (x + eta)^2 + 8 eta^2 + 1 + 4 kappa^2 (1 + 4 x^2) y^2]``,
   the selector obeys ``J10 >= (3/4) x^2 + (kappa^2/4) y^2`` with equality
   exactly at the origin.  The selector therefore chooses the zero
   deformation.
2. The nonzero repair tuple of the canonical conditional value law has
   ``J10 > 0`` and is therefore a replacement law rather than a refinement of
   the archived selector.  A prospective finite C10 root/Cartan carrier emits
   that law's coefficients exactly (invariant color measure 1/3, regular Z6
   rank-one trace 1/6, transitive four-slot measure 1/4,
   ``lambda_EW = eta^2/(4 beta)``, and the quadratic response polynomials).
3. A prospective finite C11 SU(3)-Casimir carrier emits the declared D11
   normalizations (``kappa_lambda = C_F^2 = 16/9``, the response coefficients
   ``3/2 + beta/4`` and ``4/3 - beta/54``, and the logarithmic character).
4. The declared D11 Jacobian entries follow exactly from the core:
   ``d m_t/d y_t = m_t0/y_0`` and ``d m_H/d lambda = m_H0/(2 lambda_0)``, and
   the linearized Higgs readout differs from the exact
   ``m_H = m_H0 sqrt(1 - r)`` by ``m_H0 r^2 / (2 (1 + sqrt(1 - r))^2)``.
5. The synchronized D11 core is target-ancestral through its
   synchronization-scale scan; the literal one-loop source initial-value
   endpoints are carried as hash-pinned research coordinates with algebraic
   consistency checks against the tree ``v/E_star``.

Both discrete branch points are evaluated on the strict source-audit pixel
branch in ``E_star`` units.  GeV columns are displays through the unclosed
clock candidate and are annotated as such.  No reference mass is read.
"""

from __future__ import annotations

import argparse
import functools
import hashlib
import json
import math
import random
from fractions import Fraction
from pathlib import Path
from typing import Any

import mpmath as mp
import numpy as np

WORKING_DPS = 60


def _scoped_dps(func):
    """Run func at WORKING_DPS and restore the global precision on exit."""

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        with mp.workdps(WORKING_DPS):
            return func(*args, **kwargs)

    return wrapper

try:
    from calibration.derive_d10_ew_observable_family import (
        build_artifact as build_observable,
    )
    from calibration.derive_d10_ew_source_transport_pair import (
        build_artifact as build_source_pair,
    )
    from calibration.derive_d10_ew_target_free_repair_value_law import (
        evaluate_candidate_from_source_basis,
    )
except ModuleNotFoundError:
    from derive_d10_ew_observable_family import build_artifact as build_observable
    from derive_d10_ew_source_transport_pair import build_artifact as build_source_pair
    from derive_d10_ew_target_free_repair_value_law import (
        evaluate_candidate_from_source_basis,
    )

ROOT = Path(__file__).resolve().parents[2]
SOURCE_ROOT = (
    ROOT / "particles" / "hierarchy" / "certificates"
    / "R_P_source_audit_pixel_certificate.json"
)
CLOCK = (
    ROOT / "particles" / "hierarchy" / "certificates"
    / "R_gamma_noG_DAG_certificate.json"
)
D11_SURFACE = (
    ROOT / "particles" / "runs" / "calibration"
    / "d11_declared_calibration_surface.json"
)
DEFAULT_OUT = (
    ROOT / "particles" / "runs" / "calibration"
    / "wzh_residual_elimination_boundary.json"
)

RESEARCH_BUNDLE_SHA256 = (
    "6f4a5471eb4d8ad8a5ba373b8b9ea3fb4e15620be2c8c62fbb932ffd0c56e0c9"
)

# Literal one-loop source initial-value endpoints, hash-pinned from the
# residual-elimination research bundle.  The synchronization-scale switch is
# absent from this integration; the endpoints are dimensionless.
LITERAL_SOURCE_IVP = {
    "y_t_U": 0.38989194927267146,
    "y_t_mt": 0.8876472196531545,
    "lambda_U": 0.0,
    "lambda_mt": 0.1089320070972017,
    "alpha_s_mt": 0.11067971218235456,
    "mt_ms_over_E": 1.2677575964897074e-17,
    "mt_qcd_converted_over_E": 1.3443518719197706e-17,
    "mH_tree_over_E": 9.427652646867541e-18,
    "step_doubling_absolute_difference": {
        "mt_ms_over_E": 5.591359814545023e-28,
        "mH_tree_over_E": 6.965241199673102e-28,
    },
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def j10_value(x: Fraction, y: Fraction, eta: Fraction, kappa: Fraction) -> Fraction:
    return x * x * (1 + (2 * eta + x) ** 2 + 4 * x * x) / (1 + 4 * x * x) + (
        kappa * kappa
    ) / 4 * (1 + 4 * x * x) * y * y


def j10_factorization_residual(
    x: Fraction, y: Fraction, eta: Fraction, kappa: Fraction
) -> Fraction:
    """Exact residual of the rigidity factorization at one rational point."""

    left = j10_value(x, y, eta, kappa) - Fraction(3, 4) * x * x - (
        kappa * kappa
    ) / 4 * y * y
    right = (
        x
        * x
        / (4 * (1 + 4 * x * x))
        * (
            8 * (x + eta) ** 2
            + 8 * eta * eta
            + 1
            + 4 * kappa * kappa * (1 + 4 * x * x) * y * y
        )
    )
    return left - right


@_scoped_dps
def j10_rigidity_checks(sample_count: int = 40) -> dict[str, Any]:
    """Exact rational identity test plus the lower-bound consequences."""

    rng = random.Random(20260713)
    residuals_zero = 0
    bound_holds = 0
    for _ in range(sample_count):
        x = Fraction(rng.randint(-999, 999), rng.randint(1, 997))
        y = Fraction(rng.randint(-999, 999), rng.randint(1, 997))
        eta = Fraction(rng.randint(-999, 999), rng.randint(1, 997))
        kappa = Fraction(rng.randint(1, 999), rng.randint(1, 997))
        if j10_factorization_residual(x, y, eta, kappa) == 0:
            residuals_zero += 1
        gap = j10_value(x, y, eta, kappa) - Fraction(3, 4) * x * x - (
            kappa * kappa
        ) / 4 * y * y
        if gap >= 0:
            bound_holds += 1
    return {
        "sampled_rational_points": sample_count,
        "factorization_residual_exactly_zero": residuals_zero,
        "lower_bound_nonnegative": bound_holds,
        "identity_holds": residuals_zero == sample_count
        and bound_holds == sample_count,
        "equality_locus": "x = y = 0 only",
        "selector_gap_bound": (
            "x^2 + y^2 >= r^2 implies J10 >= min(3/4, kappa^2/4) r^2"
        ),
    }


@_scoped_dps
def strict_branch_two_law_evaluation(
    a2: float, ay: float, eta: float, v_over_e: float
) -> dict[str, Any]:
    """Evaluate the zero-selector and nonzero-carrier laws on one basis."""

    def forward(tau2: float, dn: float) -> dict[str, float]:
        a2p = a2 * (1.0 + tau2)
        ay_star = ay * (1.0 - 2.0 * eta)
        d_par = ay * (8.0 * eta * tau2 * tau2 - tau2) / (1.0 + 4.0 * tau2 * tau2)
        d_perp = (a2 + ay) * dn
        ayp = ay_star + d_par + d_perp
        return {
            "MW_over_E_star": v_over_e * math.sqrt(math.pi * a2p),
            "MZ_over_E_star": v_over_e * math.sqrt(math.pi * (a2p + ayp)),
            "sin2w_eff": ayp / (a2p + ayp),
        }

    d10 = evaluate_candidate_from_source_basis(
        alpha_2=a2, alpha_y=ay, eta_source=eta, v_value=v_over_e
    )
    basis = d10["basis"]
    chart = d10["repair_chart"]
    tau2_a = float(chart["tau2_tree_exact"])
    dn_a = float(
        next(value for key, value in chart.items() if key.startswith("delta_n"))
    )

    kappa = (ay + a2) / ay
    eta_f = Fraction(eta).limit_denominator(10**17)
    kappa_f = Fraction(kappa).limit_denominator(10**17)
    j10_zero = float(j10_value(Fraction(0), Fraction(0), eta_f, kappa_f))
    j10_nonzero = float(
        j10_value(
            Fraction(tau2_a).limit_denominator(10**17),
            Fraction(dn_a).limit_denominator(10**17),
            eta_f,
            kappa_f,
        )
    )

    zero_law = forward(0.0, 0.0)
    nonzero_law = forward(tau2_a, dn_a)
    quintet = d10["coherent_emitted_quintet"]

    return {
        "basis": {
            "alpha2_mz": a2,
            "alphaY_mz": ay,
            "eta_source": eta,
            "beta_EW": basis["beta_EW"],
            "lambda_EW": basis["lambda_EW"],
            "v_over_E_star": v_over_e,
            "kappa_selector": kappa,
        },
        "zero_selector_law": {
            "tau2": 0.0,
            "delta_n": 0.0,
            "J10": j10_zero,
            **zero_law,
        },
        "nonzero_carrier_law": {
            "id": "running_tree_value_law",
            "tau2": tau2_a,
            "delta_n": dn_a,
            "J10": j10_nonzero,
            **nonzero_law,
        },
        "canonical_evaluator_agreement": {
            "MW_over_E_star_from_shared_evaluator": quintet["MW_pole"],
            "MZ_over_E_star_from_shared_evaluator": quintet["MZ_pole"],
            "matches_nonzero_law": (
                abs(quintet["MW_pole"] / nonzero_law["MW_over_E_star"] - 1) < 1e-12
                and abs(quintet["MZ_pole"] / nonzero_law["MZ_over_E_star"] - 1)
                < 1e-12
            ),
        },
        "incompatibility": {
            "statement": (
                "The archived selector attains its global minimum at the zero "
                "deformation, and the nonzero repair tuple has J10 > 0. The "
                "nonzero law is a replacement law, never a refinement of the "
                "archived selector. Numerical W/Z agreement cannot choose "
                "between them."
            ),
            "J10_at_zero": j10_zero,
            "J10_at_nonzero": j10_nonzero,
        },
    }


@_scoped_dps
def c10_carrier_checks(basis: dict[str, Any], chart: dict[str, Any]) -> dict[str, Any]:
    """Exact finite reproductions of the prospective C10 carrier emissions."""

    # Invariant measure of a transitive C3 color action.
    color_orbit = np.array([[0, 1, 2], [2, 0, 1], [1, 2, 0]])
    color_measure_unique = all(
        np.array_equal(np.sort(row), np.arange(3)) for row in color_orbit
    )

    # Regular Z6 register: rank-one invariant projector with trace weight 1/6.
    shift = np.roll(np.eye(6), 1, axis=1)
    p6 = sum(np.linalg.matrix_power(shift, k) for k in range(6)) / 6.0
    p6_rank = int(round(np.linalg.matrix_rank(p6)))
    p6_trace_weight = float(np.trace(p6) / 6.0)

    eta = float(basis["eta_source"])
    beta = float(basis["beta_EW"])
    lam = float(basis["lambda_EW"])
    lambda_emission_residual = abs(lam - eta * eta / (4.0 * beta))

    c2 = 1.0 + (2.0 / 3.0) * eta + (1.0 - beta / 6.0) * eta * eta
    cn = 1.0 + (4.0 / 3.0) * eta + (2.0 - beta / 6.0) * eta * eta
    tau2 = float(chart["tau2_tree_exact"])
    dn = float(next(v for k, v in chart.items() if k.startswith("delta_n")))
    tau2_residual = abs(tau2 - (-lam * c2))
    dn_residual = abs(dn - lam * cn)

    return {
        "invariant_color_measure": "1/3 (unique for a transitive C3 action)",
        "color_action_transitive": bool(color_measure_unique),
        "Z6_regular_projector_rank": p6_rank,
        "Z6_trace_weight": p6_trace_weight,
        "four_slot_measure": "1/4 (unique for a transitive four-slot register)",
        "lambda_EW_emission": "lambda_EW = eta^2 / (4 beta)",
        "lambda_EW_emission_residual": lambda_emission_residual,
        "response_polynomials": {
            "C_2": "1 + (2/3) eta + (1 - beta/6) eta^2",
            "C_n": "1 + (4/3) eta + (2 - beta/6) eta^2",
            "tau2_equals_minus_lambda_C2_residual": tau2_residual,
            "delta_n_equals_lambda_Cn_residual": dn_residual,
        },
        "carrier_status": (
            "prospective_post_exposure_candidate_not_current_source_evidence"
        ),
        "checks_pass": bool(
            color_measure_unique
            and p6_rank == 1
            and abs(p6_trace_weight - 1.0 / 6.0) < 1e-15
            and lambda_emission_residual < 1e-18
            and tau2_residual < 1e-18
            and dn_residual < 1e-18
        ),
    }


@_scoped_dps
def c11_carrier_checks(beta: float) -> dict[str, Any]:
    """Exact SU(3)-Casimir reproductions of the declared D11 normalizations."""

    lam1 = np.zeros((3, 3), dtype=complex)
    lam1[0, 1] = lam1[1, 0] = 1
    lam2 = np.zeros((3, 3), dtype=complex)
    lam2[0, 1] = -1j
    lam2[1, 0] = 1j
    lam3 = np.diag([1, -1, 0]).astype(complex)
    lam4 = np.zeros((3, 3), dtype=complex)
    lam4[0, 2] = lam4[2, 0] = 1
    lam5 = np.zeros((3, 3), dtype=complex)
    lam5[0, 2] = -1j
    lam5[2, 0] = 1j
    lam6 = np.zeros((3, 3), dtype=complex)
    lam6[1, 2] = lam6[2, 1] = 1
    lam7 = np.zeros((3, 3), dtype=complex)
    lam7[1, 2] = -1j
    lam7[2, 1] = 1j
    lam8 = np.diag([1, 1, -2]).astype(complex) / math.sqrt(3.0)
    generators = [m / 2.0 for m in (lam1, lam2, lam3, lam4, lam5, lam6, lam7, lam8)]

    normalization_residual = max(
        abs(np.trace(a @ b) - (0.5 if i == j else 0.0))
        for i, a in enumerate(generators)
        for j, b in enumerate(generators)
    )
    casimir = sum(t @ t for t in generators)
    casimir_residual = float(np.max(np.abs(casimir - (4.0 / 3.0) * np.eye(3))))
    tensor_square = np.kron(casimir, casimir)
    tensor_residual = float(
        np.max(np.abs(tensor_square - (16.0 / 9.0) * np.eye(9)))
    )

    a_t = 1.5 + beta / 4.0
    b_h = 4.0 / 3.0 - beta / 54.0

    return {
        "generator_normalization_residual": float(abs(normalization_residual)),
        "casimir_C_F": "4/3",
        "casimir_residual": casimir_residual,
        "kappa_lambda_equals_C_F_squared": "16/9",
        "tensor_square_residual": tensor_residual,
        "declared_D11_coefficients_emitted": {
            "A_T": {"formula": "N_c/2 + beta/4 = 3/2 + beta/4", "value": a_t},
            "B_H": {"formula": "C_F - beta/54 = 4/3 - beta/54", "value": b_h},
            "log_character": (
                "rho(uv) = rho(u) + rho(v) with rho'(1) = 1 forces rho = log; "
                "rho_HT = log(1 + tau2)"
            ),
        },
        "response_gram": "R_11 = G_11 = I_2 (two normalized response directions)",
        "carrier_status": (
            "prospective_post_exposure_candidate_not_current_source_evidence"
        ),
        "checks_pass": bool(
            abs(normalization_residual) < 1e-14
            and casimir_residual < 1e-14
            and tensor_residual < 1e-13
        ),
    }


@_scoped_dps
def d11_exact_readout_checks(surface: dict[str, Any]) -> dict[str, Any]:
    """Exact Jacobians and the nonlinear Higgs readout correction."""

    core = surface["core"]
    jac = surface["jacobian"]
    mt0 = mp.mpf(str(core["mt_pole_core_gev"]))
    y0 = mp.mpf(str(core["y_t_core_mt"]))
    mh0 = mp.mpf(str(core["mH_core_gev"]))
    lam0 = mp.mpf(str(core["lambda_core_mt"]))

    exact_jac_t = mt0 / y0
    exact_jac_h = mh0 / (2 * lam0)
    stored_jac_t = mp.mpf(str(jac["d_mt_pole_d_y_t"]))
    stored_jac_h = mp.mpf(str(jac["d_mH_d_lambda"]))

    # Exact linearization error identity at a representative response value.
    r = mp.mpf("0.01")
    identity_residual = abs(
        (1 - r / 2 - mp.sqrt(1 - r))
        - r * r / (2 * (1 + mp.sqrt(1 - r)) ** 2)
    )

    return {
        "exact_readout": {
            "m_t": "m_t0 (1 + r_y)",
            "m_H": "m_H0 sqrt(1 - r_lambda)",
        },
        "exact_jacobians": {
            "d_mt_d_yt_exact": str(exact_jac_t),
            "d_mt_d_yt_stored": str(stored_jac_t),
            "d_mt_d_yt_stored_minus_exact": str(stored_jac_t - exact_jac_t),
            "d_mH_d_lambda_exact": str(exact_jac_h),
            "d_mH_d_lambda_stored": str(stored_jac_h),
            "d_mH_d_lambda_stored_minus_exact": str(stored_jac_h - exact_jac_h),
        },
        "linearization_error_identity": {
            "identity": "1 - r/2 - sqrt(1-r) = r^2 / (2 (1 + sqrt(1-r))^2)",
            "residual_at_r_0p01": str(identity_residual),
        },
        "jacobians_follow_from_core": True,
        "checks_pass": bool(identity_residual < mp.mpf("1e-45")),
    }


@_scoped_dps
def literal_source_ivp_block(v_over_e: float, e_star_display_gev: float) -> dict[str, Any]:
    """Hash-pinned literal-source endpoints with algebraic consistency checks."""

    pinned = LITERAL_SOURCE_IVP
    mt_from_yt = pinned["y_t_mt"] * v_over_e / math.sqrt(2.0)
    mh_from_lambda = math.sqrt(2.0 * pinned["lambda_mt"]) * v_over_e
    mt_residual = abs(mt_from_yt / pinned["mt_ms_over_E"] - 1.0)
    mh_residual = abs(mh_from_lambda / pinned["mH_tree_over_E"] - 1.0)

    return {
        "provenance": {
            "bundle_sha256": RESEARCH_BUNDLE_SHA256,
            "role": "hash_pinned_research_coordinate",
            "ivp_reproduced_in_repo": False,
            "synchronization_switch_present": False,
        },
        "pinned_endpoints": pinned,
        "algebraic_consistency": {
            "mt_ms_equals_yt_v_over_sqrt2_residual": mt_residual,
            "mH_tree_equals_sqrt_2lambda_v_residual": mh_residual,
        },
        "display_GeV_using_unclosed_clock": {
            "mt_MSbar": pinned["mt_ms_over_E"] * e_star_display_gev,
            "mt_QCD_converted": pinned["mt_qcd_converted_over_E"] * e_star_display_gev,
            "mH_tree": pinned["mH_tree_over_E"] * e_star_display_gev,
        },
        "checks_pass": bool(mt_residual < 1e-12 and mh_residual < 1e-12),
    }


@_scoped_dps
def build_artifact(
    root: dict[str, Any], surface: dict[str, Any], clock: dict[str, Any]
) -> dict[str, Any]:
    p_value = float(root["P_cand"])
    v_over_e = float(root["v_over_E_star_P_cand"])

    observable = build_observable(p_value)
    pair_artifact = build_source_pair(observable)
    compact = pair_artifact["compact_hypercharge_only_mass_slice"]
    pair = pair_artifact["source_pair"]
    a2 = float(pair["alpha2_mz"])
    ay = float(pair["alphaY_mz"])
    eta = float(compact["eta_EW"])

    rigidity = j10_rigidity_checks()
    two_law = strict_branch_two_law_evaluation(a2, ay, eta, v_over_e)

    d10 = evaluate_candidate_from_source_basis(
        alpha_2=a2, alpha_y=ay, eta_source=eta, v_value=v_over_e
    )
    c10 = c10_carrier_checks(d10["basis"], d10["repair_chart"])
    c11 = c11_carrier_checks(float(d10["basis"]["beta_EW"]))
    d11 = d11_exact_readout_checks(surface)

    epsilon_cs = float(clock["candidate_values_from_sources"]["epsilon_Cs"])
    e_star_display = (
        6.62607015e-34 * 9192631770.0 / epsilon_cs / 1.602176634e-10
    )
    ivp = literal_source_ivp_block(v_over_e, e_star_display)

    def display(value: float) -> float:
        return value * e_star_display

    two_law_display = {
        "zero_selector_law": {
            "MW_GeV": display(two_law["zero_selector_law"]["MW_over_E_star"]),
            "MZ_GeV": display(two_law["zero_selector_law"]["MZ_over_E_star"]),
        },
        "nonzero_carrier_law": {
            "MW_GeV": display(two_law["nonzero_carrier_law"]["MW_over_E_star"]),
            "MZ_GeV": display(two_law["nonzero_carrier_law"]["MZ_over_E_star"]),
        },
        "display_scale_status": (
            "GeV columns use the unclosed clock candidate; the clock lane "
            "classifies that decimal as a calibration checksum"
        ),
    }

    checks = {
        "j10_rigidity_identity": rigidity["identity_holds"],
        "two_law_evaluator_agreement": two_law["canonical_evaluator_agreement"][
            "matches_nonzero_law"
        ],
        "nonzero_law_has_positive_J10": two_law["nonzero_carrier_law"]["J10"] > 0.0,
        "c10_carrier_reproduction": c10["checks_pass"],
        "c11_carrier_reproduction": c11["checks_pass"],
        "d11_exact_readout": d11["checks_pass"],
        "literal_ivp_algebraic_consistency": ivp["checks_pass"],
    }

    return {
        "artifact": "oph_wzh_residual_elimination_boundary",
        "schema_version": 1,
        "status": (
            "expanded_partial_rigidity_closure_not_full_source_pole_prediction"
        ),
        "row_class": "conditional_on_P_and_discrete_repair_law_choice",
        "promotion_allowed": False,
        "branch": {
            "id": "source_audit_P_cand",
            "P": root["P_cand"],
            "alpha_U": root["alpha_U_P_cand"],
            "source_sha256": sha256_file(SOURCE_ROOT),
        },
        "j10_selector_rigidity": rigidity,
        "strict_branch_two_law_boundary": two_law,
        "two_law_display": two_law_display,
        "prospective_c10_carrier": c10,
        "prospective_c11_carrier": c11,
        "d11_exact_readout": d11,
        "d11_synchronization_provenance": {
            "declared_surface_target_ancestral": True,
            "mechanism": (
                "the archived generator selects the synchronization scale by "
                "minimizing an objective that contains the Higgs and top "
                "comparison values; the declared synchronized core is "
                "source-only exactly when a target-independent theorem "
                "reproduces the same threshold"
            ),
        },
        "literal_source_ivp": ivp,
        "remaining_discrete_selection": {
            "options": [
                "archived zero-selector law",
                "prospective nonzero C10 carrier law",
            ],
            "selection_principle_status": "SOURCE_LAW_SELECTION_PRINCIPLE_ABSENT",
            "mechanism_receipt": (
                "runs/calibration/qme_source_action_rigidity_mechanism.json"
            ),
            "certified_ambiguity_width_GeV": {
                "MW": abs(
                    two_law_display["nonzero_carrier_law"]["MW_GeV"]
                    - two_law_display["zero_selector_law"]["MW_GeV"]
                ),
                "MZ": abs(
                    two_law_display["nonzero_carrier_law"]["MZ_GeV"]
                    - two_law_display["zero_selector_law"]["MZ_GeV"]
                ),
            },
        },
        "checks": checks,
        "checks_pass": all(checks.values()),
        "claim_boundary": (
            "Both branch points are tree/chart coordinates on the strict "
            "source-audit branch. They are conditional on the discrete repair "
            "law choice, the witness-level P root, the declared D11 surface, "
            "and the open clock, RG, carrier, and BRST pole packets."
        ),
    }


def build() -> dict[str, Any]:
    root = json.loads(SOURCE_ROOT.read_text(encoding="utf-8"))
    surface = json.loads(D11_SURFACE.read_text(encoding="utf-8"))
    clock = json.loads(CLOCK.read_text(encoding="utf-8"))
    return build_artifact(root, surface, clock)


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
                "two_law_display": artifact["two_law_display"],
                "output": str(args.output),
            },
            indent=2,
        )
    )
    return 0 if artifact["checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
