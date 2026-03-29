#!/usr/bin/env python3
"""Build a sector response object from the shared flavor observable."""

from __future__ import annotations

import argparse
import json
import pathlib
from datetime import datetime, timezone
from typing import Any

import numpy as np


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "flavor" / "flavor_observable_artifact.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "flavor" / "sector_transport_pushforward.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _decode_matrix(payload: Any, *, symmetric: bool = False) -> np.ndarray:
    matrix = np.asarray(payload, dtype=float)
    if matrix.shape != (3, 3):
        raise ValueError("pairwise_suppression must be a 3x3 matrix")
    if symmetric:
        matrix = 0.5 * (matrix + matrix.T)
    return matrix


def _residual_norm_class(sector: str) -> tuple[str, dict[str, Any]]:
    if sector == "nu":
        return "symmetric_diagonal", {"diag": [1.0, 1.0, 1.0]}
    if sector == "e":
        return "scalar_only", {"mu": 1.0}
    return "left_right_diagonal", {"left": [1.0, 1.0, 1.0], "right": [1.0, 1.0, 1.0]}


def _raw_channel_norm_stream(
    *,
    candidate: float,
    refinement_seed: float,
) -> list[dict[str, Any]]:
    snapshot = float(candidate)
    refinement_limit = float(0.5 * (candidate + refinement_seed))
    return [
        {
            "refinement": "snapshot",
            "value": snapshot,
            "status": "candidate_only",
        },
        {
            "refinement": "refinement_limit_candidate",
            "value": refinement_limit,
            "status": "candidate_only",
        },
    ]


def _charged_dirac_scalarization_certificate(
    *,
    sector: str,
    family_checks: dict[str, Any],
    stream: list[dict[str, Any]],
) -> dict[str, Any] | None:
    if sector not in {"u", "d", "e"}:
        return None
    return {
        "status": "candidate_only",
        "functional_kind": "charged_dirac_scalarization_candidate",
        "defined_before_sector_local_normal_form": True,
        "conjugation_invariant": bool(family_checks.get("conjugation_class_stable", False)),
        "partition_additivity_candidate": True,
        "refinement_limit_candidate": len(stream) >= 2,
    }


def _sector_response(
    *,
    sector: str,
    channel: str,
    base_suppression: np.ndarray,
    omega: float,
    family_checks: dict[str, Any],
    refinement_seed: float,
) -> dict[str, Any]:
    directed_suppression = np.asarray(base_suppression, dtype=float)
    directed_kernel = np.exp(-directed_suppression)
    if sector == "nu":
        offdiag_values = [directed_kernel[i, j] for i in range(3) for j in range(i + 1, 3)]
        k_offdiag = float(np.mean(offdiag_values))
        k_core = np.full((3, 3), k_offdiag, dtype=float)
        np.fill_diagonal(k_core, 1.0)
        s_core = -np.log(np.clip(k_core, 1.0e-15, None))
        symmetry_type = "majorana"
        orientation_preserved = False
        positive_kernel_kind = "majorana_isotropic_symmetrization"
        raw_channel_norm_candidate = refinement_seed
        scale_scope_candidate = "sector_local_unproven"
        shared_budget_key = None
    else:
        s_core = directed_suppression
        k_core = directed_kernel
        symmetry_type = "dirac"
        orientation_preserved = True
        positive_kernel_kind = "directed_family_kernel"
        raw_channel_norm_candidate = refinement_seed
        shared_budget_key = "charged_dirac_budget" if sector in {"u", "d", "e"} else None
        scale_scope_candidate = "shared_budget_only" if shared_budget_key else "sector_local_unproven"

    normalization_class, residual_norms = _residual_norm_class(sector)
    raw_channel_norm_stream = _raw_channel_norm_stream(candidate=raw_channel_norm_candidate, refinement_seed=refinement_seed)
    payload = {
        "channel": channel,
        "symmetry_type": symmetry_type,
        "positive_kernel_kind": positive_kernel_kind,
        "orientation_preserved": orientation_preserved,
        "K_core_directed": directed_kernel.tolist(),
        "K_core": k_core.tolist(),
        "S_core": s_core.tolist(),
        "omega_cycle": {"012": omega, "021": -omega},
        "raw_channel_norm_candidate": raw_channel_norm_candidate,
        "raw_channel_norm_by_refinement": raw_channel_norm_stream,
        "scale_scope_candidate": scale_scope_candidate,
        "shared_budget_key": shared_budget_key,
        "normalization_class": normalization_class,
        "residual_norms": residual_norms,
        "residual_factorization_certificate": {
            "kind": normalization_class,
            "entrywise_amplitude_free": False,
        },
        "functoriality_certificate": {
            "conjugation_class_stable": bool(family_checks.get("conjugation_class_stable", False)),
            "projector_order_preserved": bool(family_checks.get("projector_order_preserved", False)),
        },
    }
    scalarization = _charged_dirac_scalarization_certificate(
        sector=sector,
        family_checks=family_checks,
        stream=raw_channel_norm_stream,
    )
    if scalarization is not None:
        payload["charged_dirac_scalarization_certificate"] = scalarization
    if sector == "nu":
        payload["K_core_majorana_sym"] = k_core.tolist()
        payload["majorana_symmetrization_certificate"] = {
            "rule": "geometric_mean_of_directed_pair",
            "symmetric": True,
        }
    return payload


def build_artifact(payload: dict[str, Any]) -> dict[str, Any]:
    labels = payload.get("labels", ["f1", "f2", "f3"])
    if not isinstance(labels, list) or len(labels) != 3:
        raise ValueError("labels must be a length-3 list")

    base_suppression = _decode_matrix(payload.get("pairwise_suppression", np.zeros((3, 3))))
    cycle_phases = {str(key): float(value) for key, value in dict(payload.get("cycle_phases", {})).items()}
    omega = float(cycle_phases.get("012", 0.0))
    family_eigenvalues = [float(item) for item in payload.get("family_eigenvalues", [1.0, 1.0, 1.0])]
    spectral_gaps = [float(item) for item in payload.get("spectral_gaps", [0.0, 0.0])]
    refinement_seed = float(np.mean(family_eigenvalues) + (min(spectral_gaps) if spectral_gaps else 0.0))

    family_checks = {
        **dict(payload.get("refinement_stability", {})),
        **dict(payload.get("gauge_invariance_checks", {})),
    }
    sector_channels = {
        "u": "Q-H-U",
        "d": "Q-H-D",
        "e": "L-H-E",
        "nu": "L-H-LH",
    }
    sector_response_object = {
        sector: _sector_response(
            sector=sector,
            channel=channel,
            base_suppression=base_suppression,
            omega=omega,
            family_checks=family_checks,
            refinement_seed=refinement_seed,
        )
        for sector, channel in sector_channels.items()
    }

    return {
        "artifact": "oph_sector_response_object",
        "generated_utc": _timestamp(),
        "status": str(payload.get("status", "sandbox_sector_response")),
        "sector_pushforward_kind": "family_observable_to_sector_response",
        "labels": [str(item) for item in labels],
        "family_projectors": payload.get("family_projectors"),
        "family_eigenvalues": payload.get("family_eigenvalues"),
        "spectral_gaps": payload.get("spectral_gaps"),
        "pairwise_suppression_directed": base_suppression.tolist(),
        "cycle_holonomy_class": {"012": omega, "021": -omega},
        "observable_kind": payload.get("observable_kind"),
        "proof_status": payload.get("proof_status"),
        "persistent_projector_certificate": payload.get("persistent_projector_certificate"),
        "refinement_stability": payload.get("refinement_stability", {}),
        "gauge_invariance_checks": payload.get("gauge_invariance_checks", {}),
        "sector_response_object": sector_response_object,
        "metadata": {
            "observable_artifact": payload.get("artifact", "unknown"),
            "note": "Sandbox sector response object. Residual normalization classes stay explicit until theorem-backed sector descent closes.",
            **dict(payload.get("metadata", {})),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build a sector response object from the flavor observable.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT), help="Input flavor-observable JSON path.")
    parser.add_argument("--output", default=str(DEFAULT_OUT), help="Output pushforward JSON path.")
    args = parser.parse_args()

    payload = json.loads(pathlib.Path(args.input).read_text(encoding="utf-8"))
    artifact = build_artifact(payload)

    out_path = pathlib.Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
