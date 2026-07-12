#!/usr/bin/env python3
"""Export the current quark spread-map artifact."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
import math

from p_driven_flavor_candidate import anchored_sigma_pair_from_edge_candidate


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "flavor" / "family_excitation_evaluator.json"
DEFAULT_ODD_RESPONSE = ROOT / "particles" / "runs" / "flavor" / "quark_odd_response_law.json"
DEFAULT_SECTOR_MEAN_SPLIT = ROOT / "particles" / "runs" / "flavor" / "quark_sector_mean_split.json"
DEFAULT_EDGE_STATS_CANDIDATE = ROOT / "particles" / "runs" / "flavor" / "quark_edge_statistics_spread_candidate.json"
DEFAULT_EXACT_SIGMA_TARGET = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_exact_sigma_target.json"
DEFAULT_EDGE_STATS_ANCHOR = ROOT / "particles" / "runs" / "flavor" / "quark_edge_statistics_spread_candidate.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "quark_spread_map.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _build_candidate_sigmas(
    payload: dict,
    odd_response: dict | None,
    sector_mean_split: dict | None,
) -> tuple[float, float, str, str, dict[str, str], dict[str, float] | None]:
    sigma_u_witness = payload.get("diagnostic_witness_sigma_u")
    sigma_d_witness = payload.get("diagnostic_witness_sigma_d")
    if sigma_u_witness is not None and sigma_d_witness is not None:
        rho = float(payload["rho_ord"])
        x2 = float(payload["family_coordinate_x"][1])
        sigma_u = float(sigma_u_witness)
        sigma_d = float(sigma_d_witness)
        sigma_seed_ud = 0.5 * (sigma_u + sigma_d)
        eta_ud = 0.5 * (sigma_u - sigma_d)
        d_mean = 1.0 + rho - x2 * x2
        d_skew = 1.0 - x2 * x2 - (x2 * x2 / (1.0 + rho))
        a_ud = 1.0 / (2.0 * d_mean)
        b_ud = 1.0 / (2.0 * d_skew)
        ell_u = -(a_ud * sigma_seed_ud - b_ud * eta_ud)
        ell_d = -(a_ud * sigma_seed_ud + b_ud * eta_ud)
        return (
            sigma_u,
            sigma_d,
            "diagnostic_witness_not_source_emission",
            "hard_coded_diagnostic_witness",
            {
                "sigma_seed_ud_readback": "-(1 + rho_ord - x2^2) * (ell_u + ell_d)",
                "eta_ud_readback": "(1 - x2^2 - x2^2 / (1 + rho_ord)) * (ell_u - ell_d)",
                "sigma_u_total_log_per_side": "sigma_seed_ud_readback + eta_ud_readback",
                "sigma_d_total_log_per_side": "sigma_seed_ud_readback - eta_ud_readback",
            },
            {
                "mean_surface_log_ratio_u": ell_u,
                "mean_surface_log_ratio_d": ell_d,
                "mean_denominator": d_mean,
                "skew_denominator": d_skew,
                "sigma_seed_ud_readback": sigma_seed_ud,
                "eta_ud_readback": eta_ud,
            },
        )
    if sector_mean_split and str(sector_mean_split.get("proof_status", "")) == "closed_current_family_predictive_law":
        rho = float(payload["rho_ord"])
        x2 = float(payload["family_coordinate_x"][1])
        g_ch = float(sector_mean_split["shared_norm_value"])
        g_u = float(sector_mean_split["g_u_candidate"])
        g_d = float(sector_mean_split["g_d_candidate"])
        d_mean = 1.0 + rho - x2 * x2
        d_skew = 1.0 - x2 * x2 - (x2 * x2 / (1.0 + rho))
        ell_u = math.log(g_u / g_ch)
        ell_d = math.log(g_d / g_ch)
        sigma_seed_ud = -(d_mean * (ell_u + ell_d))
        eta_ud = d_skew * (ell_u - ell_d)
        sigma_u = sigma_seed_ud + eta_ud
        sigma_d = sigma_seed_ud - eta_ud
        return (
            sigma_u,
            sigma_d,
            "conditional_mean_surface_inversion_not_source_emission",
            "conditional_or_target_attached_mean_surface_readback",
            {
                "sigma_seed_ud_readback": "-(1 + rho_ord - x2^2) * (ell_u + ell_d)",
                "eta_ud_readback": "(1 - x2^2 - x2^2 / (1 + rho_ord)) * (ell_u - ell_d)",
                "sigma_u_total_log_per_side": "sigma_seed_ud_readback + eta_ud_readback",
                "sigma_d_total_log_per_side": "sigma_seed_ud_readback - eta_ud_readback",
            },
            {
                "mean_surface_log_ratio_u": ell_u,
                "mean_surface_log_ratio_d": ell_d,
                "mean_denominator": d_mean,
                "skew_denominator": d_skew,
                "sigma_seed_ud_readback": sigma_seed_ud,
                "eta_ud_readback": eta_ud,
            },
        )
    if odd_response:
        s_ch = odd_response.get("shared_charged_response_seed", {}).get("S_ch")
        if s_ch is not None:
            s13 = float(s_ch[0][2])
            s23 = float(s_ch[1][2])
            delta21 = float(payload["ordered_spectral_gaps"][0])
            rho = float(payload["rho_ord"])
            x2 = float(payload["family_coordinate_x"][1])
            mean_denominator = 1.0 + rho - x2 * x2
            sigma_u = s13 + (rho * delta21 / (1.0 + rho))
            sigma_d = s23 + (delta21 / (2.0 * mean_denominator))
            return (
                sigma_u,
                sigma_d,
                "shared_suppression_gap_seed_candidate",
                "shared_charged_seed_plus_ordered_gap_candidate",
                {
                    "sigma_u_total_log_per_side": "S_ch[1,3] + rho_ord * delta21 / (1 + rho_ord)",
                    "sigma_d_total_log_per_side": "S_ch[2,3] + delta21 / (2 * (1 + rho_ord - x2^2))",
                },
                None,
            )
    sigma_u_fallback = float(payload.get("diagnostic_witness_sigma_u", 5.5905))
    sigma_d_fallback = float(payload.get("diagnostic_witness_sigma_d", 3.3049))
    return (
        sigma_u_fallback,
        sigma_d_fallback,
        "diagnostic_sigma_seeded_candidate",
        "diagnostic_witness_seed",
        {
            "sigma_u_total_log_per_side": "diagnostic witness placeholder",
            "sigma_d_total_log_per_side": "diagnostic witness placeholder",
        },
        None,
    )


def _build_p_driven_candidate_sigmas(
    edge_stats_candidate: dict,
    exact_sigma_target: dict | None,
    default_edge_stats_candidate: dict | None,
) -> tuple[float, float, str, str, dict[str, str], dict[str, float]]:
    if not edge_stats_candidate:
        raise ValueError("p_driven_candidate mode requires the constructive edge-statistics candidate artifact")

    current_candidate = dict(edge_stats_candidate["candidate_sigmas"])
    exact_target = (
        dict((exact_sigma_target or {}).get("unique_exact_sigma_target", {}))
        if exact_sigma_target is not None
        else {}
    )
    default_candidate = (
        dict((default_edge_stats_candidate or {}).get("candidate_sigmas", {}))
        if default_edge_stats_candidate is not None
        else {}
    )
    anchored = anchored_sigma_pair_from_edge_candidate(
        float(current_candidate["sigma_u_total_log_per_side"]),
        float(current_candidate["sigma_d_total_log_per_side"]),
        exact_sigma_u=float(exact_target.get("sigma_u_target", 5.573928426395543)),
        exact_sigma_d=float(exact_target.get("sigma_d_target", 3.296264198808688)),
        default_candidate_sigma_u=float(default_candidate.get("sigma_u_total_log_per_side", 5.578418804072826)),
        default_candidate_sigma_d=float(default_candidate.get("sigma_d_total_log_per_side", 3.4210589139721543)),
    )
    return (
        float(anchored["sigma_u_total_log_per_side"]),
        float(anchored["sigma_d_total_log_per_side"]),
        "candidate_default_universe_exact_anchor",
        "p_driven_edge_statistics_default_universe_anchor_candidate",
        {
            "sigma_u_total_log_per_side": anchored["anchor_formula_u"],
            "sigma_d_total_log_per_side": anchored["anchor_formula_d"],
        },
        {
            "anchor_scale_u": float(anchored["anchor_scale_u"]),
            "anchor_scale_d": float(anchored["anchor_scale_d"]),
            "exact_sigma_u_target": float(exact_target.get("sigma_u_target", 5.573928426395543)),
            "exact_sigma_d_target": float(exact_target.get("sigma_d_target", 3.296264198808688)),
            "default_candidate_sigma_u": float(default_candidate.get("sigma_u_total_log_per_side", 5.578418804072826)),
            "default_candidate_sigma_d": float(default_candidate.get("sigma_d_total_log_per_side", 3.4210589139721543)),
        },
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the quark spread-map artifact.")
    parser.add_argument(
        "--mode",
        choices=("legacy_current_surface", "p_driven_candidate"),
        default="legacy_current_surface",
        help="Construction mode. The default preserves the theorem-facing current-surface behavior.",
    )
    parser.add_argument("--input", default=str(DEFAULT_INPUT))
    parser.add_argument("--odd-response-law", default=str(DEFAULT_ODD_RESPONSE))
    parser.add_argument("--sector-mean-split", default=str(DEFAULT_SECTOR_MEAN_SPLIT))
    parser.add_argument("--edge-stats-candidate", default=str(DEFAULT_EDGE_STATS_CANDIDATE))
    parser.add_argument("--exact-sigma-target", default=str(DEFAULT_EXACT_SIGMA_TARGET))
    parser.add_argument("--default-edge-stats-anchor", default=str(DEFAULT_EDGE_STATS_ANCHOR))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    payload = json.loads(Path(args.input).read_text(encoding="utf-8"))
    odd_response_path = Path(args.odd_response_law)
    odd_response = json.loads(odd_response_path.read_text(encoding="utf-8")) if odd_response_path.exists() else None
    sector_mean_split_path = Path(args.sector_mean_split)
    sector_mean_split = json.loads(sector_mean_split_path.read_text(encoding="utf-8")) if sector_mean_split_path.exists() else None
    edge_stats_candidate_path = Path(args.edge_stats_candidate)
    edge_stats_candidate = (
        json.loads(edge_stats_candidate_path.read_text(encoding="utf-8"))
        if edge_stats_candidate_path.exists()
        else None
    )
    exact_sigma_target_path = Path(args.exact_sigma_target)
    exact_sigma_target = (
        json.loads(exact_sigma_target_path.read_text(encoding="utf-8"))
        if exact_sigma_target_path.exists()
        else None
    )
    default_edge_stats_anchor_path = Path(args.default_edge_stats_anchor)
    default_edge_stats_anchor = (
        json.loads(default_edge_stats_anchor_path.read_text(encoding="utf-8"))
        if default_edge_stats_anchor_path.exists()
        else None
    )
    rho = float(payload.get("rho_ord"))
    x2 = float(payload["family_coordinate_x"][1])
    denom = 3.0 * (1.0 + rho)
    v_u = [
        -((2.0 * rho) + 1.0) / denom,
        (rho - 1.0) / denom,
        (rho + 2.0) / denom,
    ]
    v_d = [
        -((rho + 2.0) / denom),
        (1.0 - rho) / denom,
        ((2.0 * rho) + 1.0) / denom,
    ]
    if args.mode == "p_driven_candidate":
        (
            sigma_u_candidate,
            sigma_d_candidate,
            spread_status,
            sigma_source_kind,
            sigma_formula,
            mean_surface_readback,
        ) = _build_p_driven_candidate_sigmas(
            edge_stats_candidate=edge_stats_candidate,
            exact_sigma_target=exact_sigma_target,
            default_edge_stats_candidate=default_edge_stats_anchor,
        )
    else:
        (
            sigma_u_candidate,
            sigma_d_candidate,
            spread_status,
            sigma_source_kind,
            sigma_formula,
            mean_surface_readback,
        ) = _build_candidate_sigmas(payload, odd_response, sector_mean_split)
    gamma21_u = rho * sigma_u_candidate / (1.0 + rho)
    gamma32_u = sigma_u_candidate / (1.0 + rho)
    gamma21_d = sigma_d_candidate / (1.0 + rho)
    gamma32_d = rho * sigma_d_candidate / (1.0 + rho)
    q_ord = [
        (1.0 - x2 * x2) / 3.0,
        -(2.0 * (1.0 - x2 * x2)) / 3.0,
        (1.0 - x2 * x2) / 3.0,
    ]
    e_u_candidate = [sigma_u_candidate * x for x in v_u]
    e_d_candidate = [sigma_d_candidate * x for x in v_d]
    s_u_candidate = [math.exp(x) for x in e_u_candidate]
    s_d_candidate = [math.exp(x) for x in e_d_candidate]
    a_u_coeff = 0.5 * (gamma21_u + gamma32_u)
    a_d_coeff = 0.5 * (gamma21_d + gamma32_d)
    b_u_coeff = (((1.0 + x2) * gamma32_u) - ((1.0 - x2) * gamma21_u)) / (2.0 * (1.0 - x2 ** 2))
    b_d_coeff = (((1.0 + x2) * gamma32_d) - ((1.0 - x2) * gamma21_d)) / (2.0 * (1.0 - x2 ** 2))
    artifact = {
        "artifact": "oph_family_excitation_spread_map",
        "generated_utc": _timestamp(),
        "proof_status": "candidate_only_spread_values",
        "upstream_profile_shape_proof_status": payload.get("proof_status"),
        "spread_emitter_status": spread_status,
        "predictive_promotion_allowed": False,
        "source_only_sigma_emitted": False,
        "no_target_leak_certified": False,
        "claim_tier": "diagnostic_or_target_attached_spread_candidate",
        "sigma_source_kind": sigma_source_kind,
        "input_kind": payload.get("input_kind"),
        "parameterization_kind": payload.get("parameterization_kind"),
        "ordered_family_eigenvalues": payload.get("ordered_family_eigenvalues"),
        "ordered_spectral_gaps": {
            "delta21": payload["ordered_spectral_gaps"][0],
            "delta32": payload["ordered_spectral_gaps"][1],
        },
        "normalized_coordinate_x2": payload["family_coordinate_x"][1],
        "rho_ord": payload.get("rho_ord"),
        "rho_ord_reciprocal": payload.get("rho_ord_reciprocal"),
        "even_profile_rays": {
            "v_u": v_u,
            "v_d": v_d,
            "formula_kind": "centered_trace_zero_ordered_three_point_rays",
        },
        "shared_norm_origin": payload.get("shared_norm_origin"),
        "shared_norm_value": payload.get("shared_norm_value"),
        "linear_coefficient_readout": {
            "a_u": a_u_coeff,
            "a_d": a_d_coeff,
            "b_u": b_u_coeff,
            "b_d": b_d_coeff,
            "sigma_u_formula": "2 * a_u",
            "sigma_d_formula": "2 * a_d",
            "endpoint_span_formula": "sigma_q = E_q_log[2] - E_q_log[0]",
        },
        "sigma_u_total_log_per_side": sigma_u_candidate,
        "sigma_d_total_log_per_side": sigma_d_candidate,
        "sigma_formula": sigma_formula,
        "mean_surface_readback": mean_surface_readback,
        "gap_pair_u": {
            "gamma21_log_per_side": gamma21_u,
            "gamma32_log_per_side": gamma32_u,
        },
        "gap_pair_d": {
            "gamma21_log_per_side": gamma21_d,
            "gamma32_log_per_side": gamma32_d,
        },
        "E_u_log": e_u_candidate,
        "E_d_log": e_d_candidate,
        "quadratic_residual_basis_formula": "Q_ord = ctr(X_ord^2)",
        "quadratic_residual_basis_vector": q_ord,
        "kappa_u_log_per_side": None,
        "kappa_d_log_per_side": None,
        "delta_gap_pair_u": {
            "gamma21_log_per_side": "-(1 - x2^2) * kappa_u_log_per_side",
            "gamma32_log_per_side": "+(1 - x2^2) * kappa_u_log_per_side",
        },
        "delta_gap_pair_d": {
            "gamma21_log_per_side": "-(1 - x2^2) * kappa_d_log_per_side",
            "gamma32_log_per_side": "+(1 - x2^2) * kappa_d_log_per_side",
        },
        "E_u_log_closed_formula": "E_u_log + kappa_u_log_per_side * Q_ord",
        "E_d_log_closed_formula": "E_d_log + kappa_d_log_per_side * Q_ord",
        "spread_preserved_under_quadratic_residual": True,
        "mean_surface_preserved_under_quadratic_residual": True,
        "reopen_larger_family_under_quadratic_residual": False,
        "residual_exactness_missing_object": "oph_current_family_quadratic_surface_obstruction_certificate",
        "diagnostic_sigma_witness_pair": {
            "sigma_u_total_log_per_side": float(payload.get("diagnostic_witness_sigma_u", 5.5905)),
            "sigma_d_total_log_per_side": float(payload.get("diagnostic_witness_sigma_d", 3.3049)),
            "E_u_log": payload.get("diagnostic_witness_E_u_log"),
            "E_d_log": payload.get("diagnostic_witness_E_d_log"),
            "leading_symmetric_log_readback_u": [math.exp(x) for x in payload.get("diagnostic_witness_E_u_log", [])],
            "leading_symmetric_log_readback_d": [math.exp(x) for x in payload.get("diagnostic_witness_E_d_log", [])],
        },
        "active_candidate_emission": {
            "status": spread_status,
            "sigma_u_total_log_per_side": sigma_u_candidate,
            "sigma_d_total_log_per_side": sigma_d_candidate,
            "E_u_log": e_u_candidate,
            "E_d_log": e_d_candidate,
            "leading_symmetric_log_readback_u": s_u_candidate,
            "leading_symmetric_log_readback_d": s_d_candidate,
        },
        "gap_to_E_formula": payload["gap_pair_reduction"]["eigenvalue_recovery"],
        "coefficient_recovery": payload["gap_pair_reduction"]["coefficient_recovery"],
        "smallest_constructive_missing_object": "source_only_quark_spread_pair_emitter",
        "metadata": {
            "odd_response_artifact": None if odd_response is None else odd_response.get("artifact"),
            "sector_mean_split_artifact": None if sector_mean_split is None else sector_mean_split.get("artifact"),
            "constructive_edge_statistics_bridge_artifact": (
                None if edge_stats_candidate is None else edge_stats_candidate.get("artifact")
            ),
            "constructive_edge_statistics_bridge_status": (
                None if edge_stats_candidate is None else edge_stats_candidate.get("bridge_status")
            ),
            "constructive_edge_statistics_bridge_candidate_sigmas": (
                None if edge_stats_candidate is None else edge_stats_candidate.get("candidate_sigmas")
            ),
            "note": (
                "This legacy spread map carries diagnostic or target-attached numerical values only. Algebraically "
                "inverting a mean surface reproduces an attached spread but does not emit it from OPH source data. "
                "The hard-coded diagnostic witness, conditional inversion, and edge-statistics formula remain "
                "non-promotable until a no-target theorem breaks the independent up/down rescaling action."
                if args.mode == "legacy_current_surface"
                else
                "This is the off-canonical P-driven candidate spread lane. It takes the constructive edge-statistics "
                "sigma motion and anchors it to the exact default-universe sigma target so the canonical point stays "
                "exact while arbitrary-P motion remains driven by the emitted edge/gap candidate rather than by a "
                "UI-only interpolation."
            ),
        },
    }

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
