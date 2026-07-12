#!/usr/bin/env python3
"""Evaluate the quarantined RSCC quark completion candidate.

The formulas are reconstructed from canonical repository artifacts and exact
structural integers.  They contain no quark comparison table and no five
legacy flavor-template decimals.  Formula and module-ledger discovery was
nevertheless informed by the known spectrum.  The pixel witness also has an
internal Stage-5 quark ancestor, while the optional D10 GeV display comes from
a different candidate/calibration branch.  This is therefore a dimensionless
completion diagnostic, not an OPH flavor theorem or physical postdiction.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import mpmath as mp


mp.mp.dps = 90

CODE_ROOT = Path(__file__).resolve().parents[2]
RUNS = CODE_ROOT / "particles" / "runs"
DEFAULT_OUTPUT = RUNS / "flavor" / "quark_rscc_completion_candidate.json"

PIXEL_CERTIFICATE = (
    CODE_ROOT
    / "particles"
    / "hierarchy"
    / "certificates"
    / "R_P_source_audit_pixel_certificate.json"
)
REPAIR_CERTIFICATE = (
    CODE_ROOT
    / "particles"
    / "hierarchy"
    / "certificates"
    / "R_m_rep_24_certificate.json"
)
D10_TRANSPORT = RUNS / "calibration" / "d10_ew_transport_kernel.json"
SINGLET_NO_GO = RUNS / "leptons" / "charged_12_24_singlet_no_go.json"
P_DERIVATION = CODE_ROOT / "P_derivation" / "paper_math.py"

CLAIM_CLASS = (
    "post_hoc_target_informed_rscc_module_incidence_and_cumulant_ansatz_"
    "not_physical_quark_postdiction"
)


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _m(value: object) -> mp.mpf:
    return mp.mpf(str(value))


def _fmt(value: mp.mpf, digits: int = 50) -> str:
    return mp.nstr(value, n=digits, strip_zeros=False)


def _center(values: list[mp.mpf]) -> list[mp.mpf]:
    mean = mp.fsum(values) / len(values)
    return [value - mean for value in values]


def load_repository_inputs() -> tuple[dict[str, object], dict[str, Any]]:
    pixel = _read_json(PIXEL_CERTIFICATE)
    repair = _read_json(REPAIR_CERTIFICATE)
    d10 = _read_json(D10_TRANSPORT)
    singlet_no_go = _read_json(SINGLET_NO_GO)
    p_source = P_DERIVATION.read_text(encoding="utf-8")

    oriented_slots = repair["representation_sector"]["oriented_support_dimension"]
    d10_core = d10["core_source_ref"]["generator_values_unrounded"]
    inputs: dict[str, object] = {
        "P": pixel["P_cand"],
        "alpha_U": pixel["alpha_U_P_cand"],
        "N_g": 3,
        "N_c": 3,
        "d_std_S3": 2,
        "order_S3": 6,
        "beta_EW": 4,
        "oriented_repair_slots": oriented_slots,
        "orientation_dimension": 2,
        "conditional_D10_v_gev": d10["coherent_output_quintet"]["v_report"],
    }

    dependency_paths = [
        PIXEL_CERTIFICATE,
        REPAIR_CERTIFICATE,
        D10_TRANSPORT,
        SINGLET_NO_GO,
        P_DERIVATION,
    ]
    provenance = {
        "dependency_files": {
            str(path.relative_to(CODE_ROOT)): {"sha256": _sha256(path)}
            for path in dependency_paths
        },
        "pixel_branch_status": pixel.get("status"),
        "pixel_branch_eligible_for_live_particle_prediction": False,
        "pixel_source_uses_internal_stage5_quark_model": (
            "quarks = self.diagonal_quark_masses(d10.v)" in p_source
            and 'mass_source="internal_stage5_continuation"' in p_source
        ),
        "repair_register_status": repair.get("status"),
        "repair_register_is_count_only_for_family_shape": True,
        "existing_singlet_no_go_artifact": singlet_no_go.get("artifact"),
        "d10_active_readout_selector_status": d10.get(
            "active_readout_selector_status"
        ),
        "d10_current_emission_status": d10.get("current_emission_status"),
        "d10_mixed_sources_detected": d10["coherence_witness"][
            "mixed_sources_detected"
        ],
        "d10_source_P": str(d10_core["p"]),
        "d10_source_alpha_U": str(d10_core["alpha_u"]),
        "rscc_P_matches_d10_P": str(pixel["P_cand"]) == str(d10_core["p"]),
        "rscc_alpha_U_matches_d10_alpha_U": str(pixel["alpha_U_P_cand"])
        == str(d10_core["alpha_u"]),
        "old_candidate_flavor_decimals_consumed": False,
        "quark_reference_values_consumed_at_runtime": False,
        "formula_and_module_ledger_discovery_target_informed": True,
        "public_prediction_allowed": False,
    }
    return inputs, provenance


def evaluate(
    *,
    include_d10_gev_display: bool = False,
    correction_mode: str = "full_rscc",
) -> dict[str, Any]:
    if correction_mode not in {"full_rscc", "lower_order_ablation"}:
        raise ValueError(
            "correction_mode must be 'full_rscc' or 'lower_order_ablation'"
        )
    inputs, provenance = load_repository_inputs()
    P = _m(inputs["P"])
    alpha = _m(inputs["alpha_U"])
    n_g = _m(inputs["N_g"])
    n_c = _m(inputs["N_c"])
    d_std = _m(inputs["d_std_S3"])
    order_s3 = _m(inputs["order_S3"])
    beta = _m(inputs["beta_EW"])
    repair = _m(inputs["oriented_repair_slots"])
    orientation = _m(inputs["orientation_dimension"])

    n_f = n_g + d_std
    f0 = n_f - 1
    if f0 != beta:
        raise ValueError("declared RSCC cross-lock dim(F0)=beta_EW failed")

    dimensions = {
        "F": int(n_f),
        "F0": int(f0),
        "C": int(n_c),
        "R_regular_S3": int(order_s3),
        "J_orientation": int(orientation),
        "O_oriented_repair": int(repair),
        "C_tensor_F": int(n_c * n_f),
        "M_mu_up=End(F)+F0": int(n_f * n_f + f0),
        "M_mu_down=O*R*C": int(repair * order_s3 * n_c),
        "M_a_up=J*(F+R)": int(orientation * (n_f + order_s3)),
        "M_a_down_linear=O+J*F0": int(repair + orientation * f0),
        "M_a_down_quadratic=O*(C+F0)*F": int(
            repair * (n_c + f0) * n_f
        ),
        "M_g_common=O*R*(C+F0)": int(repair * order_s3 * (n_c + f0)),
        "M_g_color=O*R*C": int(repair * order_s3 * n_c),
        "M_g_relabel=O*R*(F+R)": int(
            repair * order_s3 * (n_f + order_s3)
        ),
    }

    width = mp.pi * alpha
    scalar_fraction = 1 / n_f
    centered_fraction = f0 / n_f
    common_exposure = f0 / (n_c * n_f)

    tau_f = P / beta - width * scalar_fraction
    r = mp.exp(-3 * tau_f)
    rho = 3 / (2 + r)
    x = (r - 1) / (r + 1)

    mean_u = 3 * P + (n_f + common_exposure) * width
    mean_d = 2 * P + common_exposure * width
    a_u = 3 * P + (n_f + centered_fraction) * width
    a_d = (
        2 * P
        + (1 + 1 / dimensions["M_a_down_linear=O+J*F0"]) * width
    )
    delta_g = mp.mpf(0)
    if correction_mode == "full_rscc":
        mean_u += width**2 / dimensions["M_mu_up=End(F)+F0"]
        mean_d -= width**2 / dimensions["M_mu_down=O*R*C"]
        a_u += width**2 / dimensions["M_a_up=J*(F+R)"]
        a_d += (
            orientation
            * width**2
            / dimensions["M_a_down_quadratic=O*(C+F0)*F"]
        )
        delta_g = (
            P / dimensions["M_g_common=O*R*(C+F0)"]
            + width**2 / dimensions["M_g_color=O*R*C"]
            - width**2 / dimensions["M_g_relabel=O*R*(F+R)"]
        )
    g_ch_over_v = orientation * mp.exp(-2 * mp.pi + delta_g)
    v_over_e_star = P ** (-mp.mpf("0.5")) * mp.exp(
        -2 * mp.pi / (beta * alpha)
    )
    g_ch_over_e_star = g_ch_over_v * v_over_e_star

    linear = _center([-mp.mpf(1), x, mp.mpf(1)])
    quadratic = _center([mp.mpf(1), x**2, mp.mpf(1)])
    b_u_ray = a_u * (-rho * x + rho - x - 1) / (
        (1 + rho) * (x**2 - 1)
    )
    b_d_ray = a_d * (-rho * x - rho - x + 1) / (
        (1 + rho) * (x**2 - 1)
    )
    q_u = -width * rho / (orientation * n_f)
    q_d = -width / beta
    b_u = b_u_ray + q_u
    b_d = b_d_ray + q_d

    coefficient_a = 1 / (2 * (1 + rho - x**2))
    coefficient_b = 1 / (2 * (1 - x**2 - x**2 / (1 + rho)))
    sigma_seed = (mean_u + mean_d) / 2
    eta = (mean_u - mean_d) / 2
    g_u_over_e = g_ch_over_e_star * mp.exp(
        -(coefficient_a * sigma_seed - coefficient_b * eta)
    )
    g_d_over_e = g_ch_over_e_star * mp.exp(
        -(coefficient_a * sigma_seed + coefficient_b * eta)
    )

    centered_u = [
        a_u * linear[index] + b_u * quadratic[index] for index in range(3)
    ]
    centered_d = [
        a_d * linear[index] + b_d * quadratic[index] for index in range(3)
    ]
    up_over_e = [g_u_over_e * mp.exp(value) for value in centered_u]
    down_over_e = [g_d_over_e * mp.exp(value) for value in centered_d]

    tolerance = mp.mpf("1e-70")
    assert abs(mp.fsum(linear)) < tolerance
    assert abs(mp.fsum(quadratic)) < tolerance
    assert abs(mp.fsum(centered_u)) < tolerance
    assert abs(mp.fsum(centered_d)) < tolerance

    d10_display: dict[str, Any] | None = None
    if include_d10_gev_display:
        v_gev = _m(inputs["conditional_D10_v_gev"])
        e_star_gev = v_gev / v_over_e_star
        up_gev = [value * e_star_gev for value in up_over_e]
        down_gev = [value * e_star_gev for value in down_over_e]
        d10_display = {
            "status": (
                "compare_only_mixed_branch_D10_display; P/alpha and D10 scale "
                "do not share one certified source branch"
            ),
            "v_gev": _fmt(v_gev),
            "E_star_gev_inferred": _fmt(e_star_gev),
            "g_ch_gev": _fmt(g_ch_over_e_star * e_star_gev),
            "coordinates_gev": {
                "u": _fmt(up_gev[0]),
                "c": _fmt(up_gev[1]),
                "t": _fmt(up_gev[2]),
                "d": _fmt(down_gev[0]),
                "s": _fmt(down_gev[1]),
                "b": _fmt(down_gev[2]),
            },
        }

    return {
        "artifact": "oph_quark_rscc_completion_candidate_v1",
        "claim_class": CLAIM_CLASS,
        "promotion_allowed": False,
        "correction_mode": correction_mode,
        "repository_inputs": inputs,
        "provenance": provenance,
        "declared_module_dimensions": dimensions,
        "exact_arithmetic_status": (
            "exact_given_postulated_module_incidence_effect_ranks_signs_and_cumulant_law"
        ),
        "physical_hypotheses_not_derived": [
            "F=G+V_std is the unique minimal physical family-response carrier",
            "the eight composite module incidences are source-selected",
            "the declared projectors, ranks, orientations, and signs are physical",
            "full unitary isotropy and Gaussian two-cumulant truncation follow from MAR",
            "the count-only 24-slot register acquires a family non-singlet attachment",
            "the heat time, even responses, and affine mean law are source laws",
            "the residual functional is physically minimized",
        ],
        "derived_coordinates": {
            "w": _fmt(width),
            "tau_f": _fmt(tau_f),
            "r": _fmt(r),
            "rho": _fmt(rho),
            "x": _fmt(x),
            "mean_u": _fmt(mean_u),
            "mean_d": _fmt(mean_d),
            "a_u": _fmt(a_u),
            "a_d": _fmt(a_d),
            "b_u_ray": _fmt(b_u_ray),
            "b_d_ray": _fmt(b_d_ray),
            "q_u": _fmt(q_u),
            "q_d": _fmt(q_d),
            "b_u": _fmt(b_u),
            "b_d": _fmt(b_d),
            "A_ud": _fmt(coefficient_a),
            "B_ud": _fmt(coefficient_b),
            "g_ch_over_v": _fmt(g_ch_over_v),
            "v_over_E_star": _fmt(v_over_e_star),
            "g_ch_over_E_star": _fmt(g_ch_over_e_star),
        },
        "basis": {
            "L=ctr(-1,x,1)": [_fmt(value) for value in linear],
            "Q=ctr(1,x^2,1)": [_fmt(value) for value in quadratic],
        },
        "centered_log_coordinates": {
            "up": [_fmt(value) for value in centered_u],
            "down": [_fmt(value) for value in centered_d],
        },
        "dimensionless_output_coordinates_over_E_star": {
            "u": _fmt(up_over_e[0]),
            "c": _fmt(up_over_e[1]),
            "t": _fmt(up_over_e[2]),
            "d": _fmt(down_over_e[0]),
            "s": _fmt(down_over_e[1]),
            "b": _fmt(down_over_e[2]),
        },
        "conditional_D10_gev_display": d10_display,
        "blocking_conditions": [
            "derive a source-clean P/alpha root with no Stage-5 quark ancestry",
            "derive and uniquely select the RSCC carrier/module/effect ledger",
            "derive isotropic two-cumulant truncation and all signs under refinement",
            "supply a family non-singlet attachment beyond the invariant 24-slot count",
            "close a coherent absolute scale on the same source branch",
            "supply one common-scale Yukawa and RG/threshold/scheme packet",
            "freeze prospectively before an independent comparison",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--allow-retrospective-rscc", action="store_true")
    parser.add_argument("--include-d10-gev-display", action="store_true")
    parser.add_argument("--print-json", action="store_true")
    args = parser.parse_args()

    if not args.allow_retrospective_rscc:
        parser.error(
            "refusing to evaluate a target-informed RSCC ledger without "
            "--allow-retrospective-rscc"
        )
    artifact = evaluate(include_d10_gev_display=args.include_d10_gev_display)
    text = json.dumps(artifact, indent=2, sort_keys=True) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text, encoding="utf-8")
    if args.print_json:
        print(text, end="")
    else:
        print(f"saved: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
