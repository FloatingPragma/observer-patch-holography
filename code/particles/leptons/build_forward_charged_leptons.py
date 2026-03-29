#!/usr/bin/env python3
"""Attach the charged absolute scale and expose the forward lepton candidate.

Chain role: combine the charged shape artifact with the local scale/norm lane
and compute the singular-value masses and bases.

Mathematics: matrix assembly, SVD, and explicit bookkeeping for whether the
missing hierarchy/support scalar has actually been emitted.

OPH-derived inputs: the charged shape artifact, channel norm, and shared/local
scale binding from the active `/particles` lepton lane.

Output: the forward charged artifact together with blocker metadata that says
whether the candidate is promotable or still structurally incomplete.
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib
from datetime import datetime, timezone
from typing import Any

import numpy as np


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "leptons" / "lepton_log_spectrum_readout.json"
DEFAULT_CHANNEL_NORM = ROOT / "particles" / "runs" / "leptons" / "lepton_channel_norm.json"
DEFAULT_SHARED_SCALE_BINDING = ROOT / "particles" / "runs" / "leptons" / "lepton_shared_absolute_scale_binding.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "leptons" / "forward_charged_leptons.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _encode_complex_matrix(matrix: np.ndarray) -> dict[str, Any]:
    return {
        "real": np.real(matrix).tolist(),
        "imag": np.imag(matrix).tolist(),
    }


def _svd_data(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    u_left, singular_values, v_right_h = np.linalg.svd(matrix)
    order = np.argsort(singular_values)
    singular_values = singular_values[order]
    u_left = u_left[:, order]
    v_right = v_right_h.conj().T[:, order]
    return singular_values, u_left, v_right


def build_artifact(
    payload: dict[str, Any],
    channel_norm_payload: dict[str, Any] | None,
    shared_scale_binding_payload: dict[str, Any] | None,
) -> dict[str, Any]:
    shape_shift_missing = False
    hierarchy_split_missing = False
    if payload.get("artifact") != "oph_lepton_log_spectrum_readout":
        raise ValueError("Input must be the active lepton_log_spectrum_readout artifact")
    if "gap_pair_e" not in payload or "E_e_log_centered" not in payload:
        raise ValueError("Input must provide the active ordered-gap charged payload")

    hierarchy_mode = payload.get("hierarchy_mode")
    eta_e_split = payload.get("eta_e_split_log_per_side")
    if hierarchy_mode == "two_scalar_ordered_gap_pair_candidate":
        sigma_e = float(payload["sigma_e_total_log_per_side"])
        x2 = float(payload["ordered_family_coordinate"][1])
        if eta_e_split is None:
            hierarchy_split_missing = True
            eta_e_value = float(payload["eta_e_rigid_fallback"])
        else:
            eta_e_value = float(eta_e_split)
        gamma21 = (((1.0 + x2) * sigma_e) - eta_e_value) / 2.0
        gamma32 = (((1.0 - x2) * sigma_e) + eta_e_value) / 2.0
    else:
        gap_pair = dict(payload["gap_pair_e"])
        gamma21 = float(gap_pair["gamma21_log_per_side"])
        gamma32 = float(gap_pair["gamma32_log_per_side"])

    e_e_log_centered = np.asarray(
        [
            -((2.0 * gamma21) + gamma32) / 3.0,
            (gamma21 - gamma32) / 3.0,
            (gamma21 + (2.0 * gamma32)) / 3.0,
        ],
        dtype=float,
    )
    if e_e_log_centered.shape != (3,):
        raise ValueError("E_e_log_centered must be length-3")
    shape_log_shift = payload.get("shape_log_shift_e")
    if shape_log_shift is None and not hierarchy_split_missing:
        shape_shift_missing = True
        e_e_shape = e_e_log_centered
    else:
        e_e_shape = e_e_log_centered if shape_log_shift is None else e_e_log_centered + float(shape_log_shift)
    y_e_shape = np.diag(np.exp(e_e_shape))
    singular_values_shape, u_e_left, u_e_right = _svd_data(y_e_shape)
    shape_closed = bool(payload.get("shape_closed", False))
    channel_norm_payload = channel_norm_payload or {}
    shared_scale_binding_payload = shared_scale_binding_payload or {}
    channel_norm_closed = bool(channel_norm_payload.get("channel_norm_closed", False))
    proof_status = str(channel_norm_payload.get("proof_status", "open"))
    closure_route = channel_norm_payload.get("closure_route")
    smallest_missing_object = payload.get("smallest_constructive_missing_object")
    g_e = channel_norm_payload.get("g_e")
    binding_g_e = shared_scale_binding_payload.get("g_e")
    if binding_g_e is None and shared_scale_binding_payload.get("mu_e_absolute_log_candidate") is not None:
        binding_g_e = math.exp(float(shared_scale_binding_payload["mu_e_absolute_log_candidate"]))
    binding_g_e_effective_candidate = shared_scale_binding_payload.get("g_e_effective_candidate")
    effective_g_e = g_e if (channel_norm_closed and g_e is not None) else binding_g_e
    y_e = None
    singular_values_abs = None
    if effective_g_e is not None:
        g_e_float = float(effective_g_e)
        y_e = g_e_float * y_e_shape
        singular_values_abs = [g_e_float * float(item) for item in singular_values_shape.tolist()]

    if hierarchy_split_missing:
        closure_state = "hierarchy_split_missing"
    elif shape_shift_missing:
        closure_state = "shape_shift_missing"
    elif smallest_missing_object == "eta_source_support_extension_log_per_side":
        closure_state = "support_extension_missing"
    elif shape_closed and channel_norm_closed and y_e is not None:
        closure_state = "absolute_scale_closed"
    elif shape_closed and binding_g_e is not None:
        closure_state = "shared_writeback_current_family"
    elif shape_closed and proof_status == "shared_budget_only":
        closure_state = "shared_budget_only"
    elif shape_closed:
        closure_state = "ratio_closed"
    else:
        closure_state = "open"

    return {
        "artifact": "oph_forward_charged_leptons",
        "generated_utc": _timestamp(),
        "source_mode": "ordered_ratio_readout",
        "labels": payload.get("labels"),
        "Y_e_shape": _encode_complex_matrix(y_e_shape),
        "singular_values_shape": [float(x) for x in singular_values_shape.tolist()],
        "Y_e": None if y_e is None else _encode_complex_matrix(y_e),
        "singular_values_abs": singular_values_abs,
        "E_e_log_centered": payload.get("E_e_log_centered"),
        "hierarchy_mode": payload.get("hierarchy_mode"),
        "eta_e_split_log_per_side": payload.get("eta_e_split_log_per_side"),
        "eta_e_rigid_fallback": payload.get("eta_e_rigid_fallback"),
        "shape_log_shift_e": payload.get("shape_log_shift_e"),
        "shape_shift_missing": shape_shift_missing,
        "hierarchy_split_missing": hierarchy_split_missing,
        "U_e_left": _encode_complex_matrix(u_e_left),
        "U_e_right": _encode_complex_matrix(u_e_right),
        "shape_closed": shape_closed,
        "channel_norm": {
            "channel_norm_closed": channel_norm_closed,
            "g_e_candidate": channel_norm_payload.get("g_e_candidate"),
            "g_e": g_e,
            "g_e_by_refinement": channel_norm_payload.get("g_e_by_refinement"),
            "proof_status": proof_status,
            "closure_route": closure_route,
            "shared_scale_binding_artifact": shared_scale_binding_payload.get("artifact"),
            "shared_scale_binding_g_e": binding_g_e,
            "shared_scale_binding_g_e_effective_candidate": binding_g_e_effective_candidate,
            "shared_scale_binding_candidate_status": shared_scale_binding_payload.get("candidate_status"),
            "shared_scale_binding_mu_common_log": shared_scale_binding_payload.get("mu_common_log"),
            "shared_scale_binding_gamma_min_common": shared_scale_binding_payload.get("gamma_min_common"),
            "shared_scale_binding_status": shared_scale_binding_payload.get("proof_status"),
            "scale_scope": channel_norm_payload.get("scale_scope"),
            "shared_budget_key": channel_norm_payload.get("shared_budget_key"),
            "shared_budget_share_e": channel_norm_payload.get("shared_budget_share_e"),
            "beta_e_by_refinement": channel_norm_payload.get("beta_e_by_refinement"),
            "charged_budget_total_by_refinement": channel_norm_payload.get("charged_budget_total_by_refinement"),
        },
        "basis_contract": {
            "labels": payload.get("labels"),
            "orientation_preserved": bool(payload.get("projector_basis_provenance", {}).get("orientation_preserved", False)),
        },
        "closure_state": closure_state,
        "theorem_status": closure_state,
        "absolute_scale_closed": closure_state == "absolute_scale_closed",
        "metadata": {
            **dict(payload.get("metadata", {})),
            "channel_norm_artifact": channel_norm_payload.get("artifact"),
            "shared_scale_binding_artifact": shared_scale_binding_payload.get("artifact"),
            "closure_route": closure_route,
            "input_artifact": payload.get("artifact"),
            "shape_contract": (
                "two_scalar_ordered_gap_pair_plus_common_shift"
                if payload.get("hierarchy_mode") == "two_scalar_ordered_gap_pair_candidate"
                else "ordered_gap_centered_plus_common_shift"
            ),
            "active_hierarchy_blocker": (
                "eta_e_split_log_per_side"
                if hierarchy_split_missing
                else smallest_missing_object
            ),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build blind charged-lepton shape and absolute artifacts.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT), help="Input JSON path.")
    parser.add_argument(
        "--channel-norm",
        default=str(DEFAULT_CHANNEL_NORM),
        help="Optional channel-norm JSON path.",
    )
    parser.add_argument(
        "--shared-scale-binding",
        default=str(DEFAULT_SHARED_SCALE_BINDING),
        help="Optional shared absolute-scale binding JSON path.",
    )
    parser.add_argument("--output", default=str(DEFAULT_OUT), help="Output JSON path.")
    args = parser.parse_args()

    input_path = pathlib.Path(args.input)
    payload = json.loads(input_path.read_text(encoding="utf-8"))
    channel_norm_path = pathlib.Path(args.channel_norm)
    channel_norm_payload = json.loads(channel_norm_path.read_text(encoding="utf-8")) if channel_norm_path.exists() else None
    binding_path = pathlib.Path(args.shared_scale_binding)
    shared_scale_binding_payload = json.loads(binding_path.read_text(encoding="utf-8")) if binding_path.exists() else None
    artifact = build_artifact(payload, channel_norm_payload, shared_scale_binding_payload)

    out_path = pathlib.Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
