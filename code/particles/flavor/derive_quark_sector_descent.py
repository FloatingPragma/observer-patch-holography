#!/usr/bin/env python3
"""Assemble the active quark descent from closed mean data plus open odd data.

Chain role: combine the closed common-sector pieces with the remaining odd
source lane to build explicit up/down Yukawa tensors.

Mathematics: factorized diagonal/excitation decomposition, left/right diagonal
log shifts, and complex matrix reconstruction from suppression and phase data.

OPH-derived inputs: the flavor response law, shared quark normalization,
sector-mean split, spread map, and diagonal gap/B-odd evaluator artifacts.

Output: the source-side quark descent artifact together with the exact residual
objects that still block certification.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import pathlib
import sys
from datetime import datetime, timezone
from typing import Any

import numpy as np


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_OBSERVABLE = ROOT / "particles" / "runs" / "flavor" / "flavor_observable_artifact.json"
DEFAULT_CHARGED_BUDGET = ROOT / "particles" / "runs" / "flavor" / "charged_budget_transport.json"
DEFAULT_TENSORS = ROOT / "particles" / "runs" / "flavor" / "suppression_phase_tensors.json"
DEFAULT_RESPONSE_LAW = ROOT / "particles" / "runs" / "flavor" / "quark_odd_response_law.json"
DEFAULT_SHARED_NORM_BINDING = ROOT / "particles" / "runs" / "flavor" / "quark_shared_absolute_norm_binding.json"
DEFAULT_SECTOR_MEAN_SPLIT = ROOT / "particles" / "runs" / "flavor" / "quark_sector_mean_split.json"
DEFAULT_SPREAD_MAP = ROOT / "particles" / "runs" / "flavor" / "quark_spread_map.json"
DEFAULT_DIAGONAL_GAP_SHIFT = ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_gap_shift_map.json"
DEFAULT_DIAGONAL_GAP_SHIFT_SCALAR_EVALUATOR = ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_gap_shift_scalar_evaluator.json"
DEFAULT_DIAGONAL_GAP_SHIFT_TAU_MAP = ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_gap_shift_tau_map.json"
DEFAULT_DIAGONAL_GAP_SHIFT_EMITTER = ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_gap_shift_emitter.json"
DEFAULT_B_ODD_SOURCE_SCALAR_EVALUATOR = ROOT / "particles" / "runs" / "flavor" / "quark_diagonal_B_odd_source_scalar_evaluator.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "quark_sector_descent.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _as_diag(values: np.ndarray) -> list[float]:
    return [float(x) for x in values.tolist()]


def _yukawa_from_factorized(
    g_scalar: float,
    diag_left: np.ndarray,
    suppression: np.ndarray,
    phase: np.ndarray,
    diag_right: np.ndarray,
) -> np.ndarray:
    return (
        float(g_scalar)
        * np.diag(np.exp(diag_left))
        @ (np.exp(-suppression) * np.exp(1j * phase))
        @ np.diag(np.exp(diag_right))
    )


def build_artifact(
    response_law: dict[str, Any],
    shared_norm_binding: dict[str, Any] | None = None,
    sector_mean_split: dict[str, Any] | None = None,
    spread_map: dict[str, Any] | None = None,
    diagonal_gap_shift_map: dict[str, Any] | None = None,
    diagonal_gap_shift_scalar_evaluator: dict[str, Any] | None = None,
    diagonal_gap_shift_tau_map: dict[str, Any] | None = None,
    diagonal_gap_shift_emitter: dict[str, Any] | None = None,
    b_odd_source_scalar_evaluator: dict[str, Any] | None = None,
) -> dict[str, Any]:
    odd_splitter = np.asarray(response_law["Xi_q_real"], dtype=float) + 1j * np.asarray(response_law["Xi_q_imag"], dtype=float)
    seed = dict(response_law.get("shared_charged_response_seed", {}))
    s_ch = np.asarray(seed["S_ch"], dtype=float)
    phi_ch = np.asarray(seed["Phi_ch"], dtype=float)
    delta_s = np.asarray(response_law["Delta_S_q"], dtype=float)
    delta_phi = np.asarray(response_law["Delta_Phi_q"], dtype=float)
    delta_logd_left = np.asarray(response_law["Delta_logD_left_q"], dtype=float)
    delta_logd_right = np.asarray(response_law["Delta_logD_right_q"], dtype=float)
    delta_logg_q = float(response_law["delta_logg_q"])

    s_u = s_ch + delta_s
    s_d = s_ch - delta_s
    phi_u = phi_ch + delta_phi
    phi_d = phi_ch - delta_phi
    logd_ch_left = np.asarray(seed.get("logD_ch_left", [0.0, 0.0, 0.0]), dtype=float)
    logd_ch_right = np.asarray(seed.get("logD_ch_right", [0.0, 0.0, 0.0]), dtype=float)
    logd_u_left = logd_ch_left + delta_logd_left
    logd_d_left = logd_ch_left - delta_logd_left
    logd_u_right = logd_ch_right + delta_logd_right
    logd_d_right = logd_ch_right - delta_logd_right

    shared_norm_binding = shared_norm_binding or {}
    if sector_mean_split and sector_mean_split.get("g_u_candidate") is not None and sector_mean_split.get("g_d_candidate") is not None:
        g_ch = float(sector_mean_split.get("shared_norm_value", (shared_norm_binding or {}).get("g_ch", seed.get("g_ch_candidate", 1.0))))
        g_u = float(sector_mean_split["g_u_candidate"])
        g_d = float(sector_mean_split["g_d_candidate"])
        shared_norm_origin = sector_mean_split.get("artifact", "oph_quark_sector_mean_split")
        shared_norm_status = sector_mean_split.get("proof_status", "current_family_candidate_only")
    elif shared_norm_binding.get("g_ch") is not None:
        g_ch = float(shared_norm_binding["g_ch"])
        g_u = float(shared_norm_binding.get("g_u", g_ch))
        g_d = float(shared_norm_binding.get("g_d", g_ch))
        shared_norm_origin = shared_norm_binding.get("artifact", "quark_shared_absolute_norm_binding")
        shared_norm_status = shared_norm_binding.get("proof_status", "current_family_writeback_complete")
    else:
        g_ch = float(seed.get("g_ch_candidate", 1.0))
        g_u = g_ch * float(np.exp(delta_logg_q))
        g_d = g_ch * float(np.exp(-delta_logg_q))
        shared_norm_origin = "charged_dirac_scalarization_restriction"
        shared_norm_status = "numerically_pinned_candidate"

    y_u = _yukawa_from_factorized(g_u, logd_u_left, s_u, phi_u, logd_u_right)
    y_d = _yukawa_from_factorized(g_d, logd_d_left, s_d, phi_d, logd_d_right)
    h_u = y_u @ y_u.conj().T
    h_d = y_d @ y_d.conj().T
    commutator = h_u @ h_d - h_d @ h_u
    noncentrality = float(np.linalg.norm(commutator, ord="fro"))

    sector_distinctness = float(
        max(
            np.max(np.abs(s_u - s_d)),
            np.max(np.abs(phi_u - phi_d)),
            np.max(np.abs(logd_u_left - logd_d_left)),
        )
    )
    local_shape_closed = (
        bool(response_law.get("coefficient_free", False))
        and response_law.get("lift_constants") in (None, {})
        and str(response_law.get("delta_logg_q_status", "")) == "closed_zero_corollary"
        and str(response_law.get("lift_parameterization_kind", "")) == "single_kappa_signed_projector_ray"
    )
    spread_map = spread_map or {}
    diagonal_gap_shift_map = diagonal_gap_shift_map or {}
    diagonal_gap_shift_scalar_evaluator = diagonal_gap_shift_scalar_evaluator or {}
    diagonal_gap_shift_tau_map = diagonal_gap_shift_tau_map or {}
    diagonal_gap_shift_emitter = diagonal_gap_shift_emitter or {}
    b_odd_source_scalar_evaluator = b_odd_source_scalar_evaluator or {}
    sigma_u = spread_map.get("sigma_u_total_log_per_side")
    sigma_d = spread_map.get("sigma_d_total_log_per_side")
    even_excitation_proof_status = "missing"
    e_u_log = None
    e_d_log = None
    tau_u = diagonal_gap_shift_tau_map.get("tau_u_log_per_side")
    tau_d = diagonal_gap_shift_tau_map.get("tau_d_log_per_side")
    tau_source = "tau_map" if tau_u is not None and tau_d is not None else "none"
    if tau_source == "none":
        tau_u = diagonal_gap_shift_emitter.get("tau_u_log_per_side")
        tau_d = diagonal_gap_shift_emitter.get("tau_d_log_per_side")
        tau_source = "emitter" if tau_u is not None and tau_d is not None else "none"
    if tau_source == "none":
        tau_u = diagonal_gap_shift_map.get("tau_u_log_per_side")
        tau_d = diagonal_gap_shift_map.get("tau_d_log_per_side")
        tau_source = "map" if tau_u is not None and tau_d is not None else "none"
    if spread_map.get("E_u_log") is not None and spread_map.get("E_d_log") is not None:
        e_u_log = [float(value) for value in spread_map["E_u_log"]]
        e_d_log = [float(value) for value in spread_map["E_d_log"]]
        even_excitation_proof_status = str(spread_map.get("spread_emitter_status", "candidate_spread_values_populated"))
    elif sigma_u is not None and sigma_d is not None:
        v_u = np.asarray(spread_map["even_profile_rays"]["v_u"], dtype=float)
        v_d = np.asarray(spread_map["even_profile_rays"]["v_d"], dtype=float)
        e_u_log = (float(sigma_u) * v_u).tolist()
        e_d_log = (float(sigma_d) * v_d).tolist()
        even_excitation_proof_status = "candidate_spread_values_populated"
    if e_u_log is not None and e_d_log is not None and tau_u is not None and tau_d is not None:
        b_ord = np.asarray(
            diagonal_gap_shift_emitter.get("B_ord", diagonal_gap_shift_map.get("B_ord", [-1.0, 0.0, 1.0])),
            dtype=float,
        )
        e_u_log = (np.asarray(e_u_log, dtype=float) + float(tau_u) * b_ord).tolist()
        e_d_log = (np.asarray(e_d_log, dtype=float) + float(tau_d) * b_ord).tolist()
        diagonal_status = str(
            diagonal_gap_shift_tau_map.get("proof_status", "candidate_only")
            if tau_source == "tau_map"
            else diagonal_gap_shift_emitter.get("proof_status", "candidate_only")
            if tau_source == "emitter"
            else diagonal_gap_shift_map.get("proof_status", "candidate_only")
        )
        even_excitation_proof_status = "closed" if diagonal_status in {"closed", "exact"} else "diagonal_gap_shift_candidate_populated"
    b_odd_source_values_open = (
        b_odd_source_scalar_evaluator.get("artifact") == "oph_quark_diagonal_B_odd_source_scalar_evaluator"
        and (
            b_odd_source_scalar_evaluator.get("J_B_source_u") is None
            or b_odd_source_scalar_evaluator.get("J_B_source_d") is None
        )
    )
    tau_pair_open = tau_u is None or tau_d is None
    if tau_pair_open and even_excitation_proof_status == "closed":
        even_excitation_proof_status = "pure_B_source_readback_open"
    elif b_odd_source_values_open and even_excitation_proof_status == "closed":
        even_excitation_proof_status = "pure_B_odd_source_values_open"
    exact_missing_object = None
    if tau_pair_open and (
        diagonal_gap_shift_tau_map.get("artifact") == "oph_family_excitation_diagonal_gap_shift_tau_map"
        or diagonal_gap_shift_scalar_evaluator.get("artifact") == "oph_family_excitation_diagonal_gap_shift_scalar_evaluator"
    ):
        exact_missing_object = diagonal_gap_shift_tau_map.get(
            "smallest_constructive_missing_object",
            diagonal_gap_shift_scalar_evaluator.get(
                "smallest_constructive_missing_object",
                "tau_u_log_per_side_and_tau_d_log_per_side",
            ),
        )
    elif b_odd_source_values_open:
        exact_missing_object = b_odd_source_scalar_evaluator.get(
            "smallest_constructive_missing_object",
            "J_B_source_u_and_J_B_source_d",
        )
    elif not str(even_excitation_proof_status).startswith("current_family_exact") and even_excitation_proof_status != "closed":
        exact_missing_object = (
            diagonal_gap_shift_scalar_evaluator.get(
                "smallest_constructive_missing_object",
                "oph_family_excitation_diagonal_common_gap_shift_source_emission",
            )
            if tau_u is None
            and tau_d is None
            and (
                diagonal_gap_shift_scalar_evaluator.get("artifact") == "oph_family_excitation_diagonal_gap_shift_scalar_evaluator"
                or diagonal_gap_shift_tau_map.get("artifact") == "oph_family_excitation_diagonal_gap_shift_tau_map"
                or diagonal_gap_shift_emitter.get("artifact") == "oph_family_excitation_diagonal_gap_shift_emitter"
                or diagonal_gap_shift_map.get("artifact") == "oph_family_excitation_diagonal_gap_shift_map"
            )
            else "oph_family_excitation_spread_map"
        )

    return {
        "artifact": "oph_quark_sector_descent",
        "generated_utc": _timestamp(),
        "channel_splitter_kind": "odd_projector_resolved_splitter",
        "proof_status": "odd_splitter_candidate",
        "uses_full_projector_algebra": True,
        "quark_descent_proof_status": (
            "closed"
            if local_shape_closed and not bool(response_law.get("degenerate_placeholder_fallback_used", False))
            else "open"
        ),
        "response_map_id": response_law.get("response_map_id"),
        "coefficient_free": bool(response_law.get("coefficient_free", False)),
        "orientation_fixed": bool(response_law.get("orientation_fixed", False)),
        "delta_logg_q_status": response_law.get("delta_logg_q_status"),
        "shared_scalarization_law_status": response_law.get("shared_scalarization_law_status"),
        "degenerate_placeholder_fallback_used": bool(response_law.get("degenerate_placeholder_fallback_used", False)),
        "projector_split_coefficients": list(response_law.get("projector_split_coefficients", [])),
        "Xi_q_real": np.real(odd_splitter).tolist(),
        "Xi_q_imag": np.imag(odd_splitter).tolist(),
        "Delta_S_q": delta_s.tolist(),
        "Delta_Phi_q": delta_phi.tolist(),
        "Delta_logD_left_q": _as_diag(delta_logd_left),
        "Delta_logD_right_q": _as_diag(delta_logd_right),
        "delta_logg_q": float(delta_logg_q),
        "S_u": s_u.tolist(),
        "S_d": s_d.tolist(),
        "Phi_u": phi_u.tolist(),
        "Phi_d": phi_d.tolist(),
        "D_u_left_log": _as_diag(logd_u_left),
        "D_u_right_log": _as_diag(logd_u_right),
        "D_d_left_log": _as_diag(logd_d_left),
        "D_d_right_log": _as_diag(logd_d_right),
        "g_ch": g_ch,
        "g_u": g_u,
        "g_d": g_d,
        "u_raw_channel_norm_candidate": g_u,
        "d_raw_channel_norm_candidate": g_d,
        "shared_norm_origin": shared_norm_origin,
        "shared_norm_status": shared_norm_status,
        "sector_mean_split_artifact": None if sector_mean_split is None else sector_mean_split.get("artifact"),
        "sector_mean_split_status": None if sector_mean_split is None else sector_mean_split.get("proof_status"),
        "even_excitation_theorem_candidate": "oph_family_excitation_evaluator",
        "even_excitation_gap_map_candidate": "oph_family_excitation_gap_map",
        "even_excitation_spread_map_candidate": "oph_family_excitation_spread_map",
        "even_excitation_diagonal_gap_shift_candidate": "oph_family_excitation_diagonal_gap_shift_map",
        "even_excitation_input_kind": "ordered_branch_generator_spectral_package",
        "even_evaluator_family": "trace_zero_quadratic_on_ordered_branch_generator",
        "ordered_branch_coordinate_formula": "x_a = 2*(lambda_a-lambda_1)/(lambda_3-lambda_1) - 1",
        "even_evaluator_centering_mode": "trace_zero",
        "exact_missing_object": exact_missing_object,
        "hierarchy_driver_status": (
            "candidate_spread_values_populated"
            if e_u_log is not None and e_d_log is not None
            else "ratio_closed_spread_only_missing"
        ),
        "even_excitation_source_artifact": spread_map.get("artifact"),
        "even_excitation_source_status": spread_map.get("spread_emitter_status"),
        "diagonal_gap_shift_artifact": diagonal_gap_shift_map.get("artifact"),
        "diagonal_gap_shift_status": diagonal_gap_shift_map.get("proof_status"),
        "diagonal_gap_shift_scalar_evaluator_artifact": diagonal_gap_shift_scalar_evaluator.get("artifact"),
        "diagonal_gap_shift_scalar_evaluator_status": diagonal_gap_shift_scalar_evaluator.get("proof_status"),
        "diagonal_gap_shift_tau_map_artifact": diagonal_gap_shift_tau_map.get("artifact"),
        "diagonal_gap_shift_tau_map_status": diagonal_gap_shift_tau_map.get("proof_status"),
        "diagonal_gap_shift_emitter_artifact": diagonal_gap_shift_emitter.get("artifact"),
        "diagonal_gap_shift_emitter_status": diagonal_gap_shift_emitter.get("proof_status"),
        "b_odd_source_scalar_evaluator_artifact": b_odd_source_scalar_evaluator.get("artifact"),
        "b_odd_source_scalar_evaluator_status": b_odd_source_scalar_evaluator.get("proof_status"),
        "J_B_source_u": b_odd_source_scalar_evaluator.get("J_B_source_u"),
        "J_B_source_d": b_odd_source_scalar_evaluator.get("J_B_source_d"),
        "tau_u_log_per_side": tau_u,
        "tau_d_log_per_side": tau_d,
        "even_excitation_proof_status": even_excitation_proof_status,
        "E_u_log": e_u_log,
        "E_d_log": e_d_log,
        "budget_neutrality_certificate": {
            "status": "closed",
            "delta_logg_q": float(delta_logg_q),
            "u_plus_d_odd_part_sum": 0.0,
        },
        "sector_distinctness_witness": {
            "status": "closed" if sector_distinctness > 1.0e-12 else "failed",
            "value": sector_distinctness,
        },
        "noncentrality_witness": {
            "status": "closed" if noncentrality > 1.0e-18 else "failed",
            "fro_norm": noncentrality,
        },
        "dense_entrywise_amplitude_used": False,
        "metadata": {
            "shared_charged_seed_kind": seed.get("seed_kind"),
            "note": "Candidate odd quark-sector descent built from an explicit shared charged response seed plus the separate odd-response-law artifact. The local shape map is explicit and coefficient-free on the current single-kappa ray. The reduced even family is fixed as a trace-zero quadratic on the ordered branch coordinate; when a current-family spread pair is available, the centered even-log profiles are read back directly from the exported spread map. If the next diagonal gap-shift family is later populated, the builder can lift those logs by B_ord without reopening the old larger branch. A compact current-family sector-mean split may refine g_u and g_d without changing those rays. No reference-fit exact-readout artifact is consumed here.",
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Derive a candidate odd quark-sector splitter.")
    parser.add_argument("--observable", default=str(DEFAULT_OBSERVABLE), help="Input flavor observable JSON path.")
    parser.add_argument("--charged-budget", default=str(DEFAULT_CHARGED_BUDGET), help="Input charged-budget JSON path.")
    parser.add_argument("--tensors", default=str(DEFAULT_TENSORS), help="Input suppression/phase tensors JSON path.")
    parser.add_argument("--response-law", default=str(DEFAULT_RESPONSE_LAW), help="Input quark odd-response-law JSON path.")
    parser.add_argument("--shared-norm-binding", default=str(DEFAULT_SHARED_NORM_BINDING), help="Optional shared quark norm binding JSON path.")
    parser.add_argument("--sector-mean-split", default=str(DEFAULT_SECTOR_MEAN_SPLIT), help="Optional quark sector-mean split JSON path.")
    parser.add_argument("--spread-map", default=str(DEFAULT_SPREAD_MAP), help="Optional quark spread-map JSON path.")
    parser.add_argument("--diagonal-gap-shift", default=str(DEFAULT_DIAGONAL_GAP_SHIFT), help="Optional quark diagonal gap-shift JSON path.")
    parser.add_argument("--diagonal-gap-shift-scalar-evaluator", default=str(DEFAULT_DIAGONAL_GAP_SHIFT_SCALAR_EVALUATOR), help="Optional quark diagonal gap-shift scalar evaluator JSON path.")
    parser.add_argument("--diagonal-gap-shift-tau-map", default=str(DEFAULT_DIAGONAL_GAP_SHIFT_TAU_MAP), help="Optional quark diagonal gap-shift tau-map JSON path.")
    parser.add_argument("--diagonal-gap-shift-emitter", default=str(DEFAULT_DIAGONAL_GAP_SHIFT_EMITTER), help="Optional quark diagonal gap-shift emitter JSON path.")
    parser.add_argument("--b-odd-source-scalar-evaluator", default=str(DEFAULT_B_ODD_SOURCE_SCALAR_EVALUATOR), help="Optional quark pure-B odd source scalar evaluator JSON path.")
    parser.add_argument("--output", default=str(DEFAULT_OUT), help="Output quark-sector-descent JSON path.")
    args = parser.parse_args()

    response_path = pathlib.Path(args.response_law)
    if response_path.exists():
        response_law = json.loads(response_path.read_text(encoding="utf-8"))
    else:
        module_path = ROOT / "particles" / "flavor" / "derive_quark_odd_response_law.py"
        spec = importlib.util.spec_from_file_location("derive_quark_odd_response_law", module_path)
        if spec is None or spec.loader is None:
            raise SystemExit(f"unable to load {module_path}")
        module = importlib.util.module_from_spec(spec)
        sys.modules.setdefault("derive_quark_odd_response_law", module)
        spec.loader.exec_module(module)
        build_response_law = module.build_artifact

        observable = json.loads(pathlib.Path(args.observable).read_text(encoding="utf-8"))
        charged_budget = json.loads(pathlib.Path(args.charged_budget).read_text(encoding="utf-8"))
        tensors = json.loads(pathlib.Path(args.tensors).read_text(encoding="utf-8"))
        response_law = build_response_law(observable, charged_budget, tensors)
    binding_path = pathlib.Path(args.shared_norm_binding)
    shared_norm_binding = json.loads(binding_path.read_text(encoding="utf-8")) if binding_path.exists() else None
    sector_mean_split_path = pathlib.Path(args.sector_mean_split)
    sector_mean_split = json.loads(sector_mean_split_path.read_text(encoding="utf-8")) if sector_mean_split_path.exists() else None
    spread_map_path = pathlib.Path(args.spread_map)
    spread_map = json.loads(spread_map_path.read_text(encoding="utf-8")) if spread_map_path.exists() else None
    diagonal_gap_shift_path = pathlib.Path(args.diagonal_gap_shift)
    diagonal_gap_shift_map = json.loads(diagonal_gap_shift_path.read_text(encoding="utf-8")) if diagonal_gap_shift_path.exists() else None
    diagonal_gap_shift_scalar_evaluator_path = pathlib.Path(args.diagonal_gap_shift_scalar_evaluator)
    diagonal_gap_shift_scalar_evaluator = (
        json.loads(diagonal_gap_shift_scalar_evaluator_path.read_text(encoding="utf-8"))
        if diagonal_gap_shift_scalar_evaluator_path.exists()
        else None
    )
    diagonal_gap_shift_tau_map_path = pathlib.Path(args.diagonal_gap_shift_tau_map)
    diagonal_gap_shift_tau_map = (
        json.loads(diagonal_gap_shift_tau_map_path.read_text(encoding="utf-8"))
        if diagonal_gap_shift_tau_map_path.exists()
        else None
    )
    diagonal_gap_shift_emitter_path = pathlib.Path(args.diagonal_gap_shift_emitter)
    diagonal_gap_shift_emitter = (
        json.loads(diagonal_gap_shift_emitter_path.read_text(encoding="utf-8"))
        if diagonal_gap_shift_emitter_path.exists()
        else None
    )
    b_odd_source_scalar_evaluator_path = pathlib.Path(args.b_odd_source_scalar_evaluator)
    b_odd_source_scalar_evaluator = (
        json.loads(b_odd_source_scalar_evaluator_path.read_text(encoding="utf-8"))
        if b_odd_source_scalar_evaluator_path.exists()
        else None
    )
    artifact = build_artifact(
        response_law,
        shared_norm_binding,
        sector_mean_split,
        spread_map,
        diagonal_gap_shift_map,
        diagonal_gap_shift_scalar_evaluator,
        diagonal_gap_shift_tau_map,
        diagonal_gap_shift_emitter,
        b_odd_source_scalar_evaluator,
    )

    out_path = pathlib.Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
