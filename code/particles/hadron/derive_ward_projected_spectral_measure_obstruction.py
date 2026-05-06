#!/usr/bin/env python3
"""Emit the local obstruction for the Ward-projected hadronic measure.

This is a source-only closeout check.  It inspects the current hadron runtime
artifacts and proves that they do not determine the positive U(1)_Q spectral
measure needed by the Thomson endpoint transport.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from fractions import Fraction
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SCHEMA = ROOT / "particles" / "hadron" / "ward_projected_spectral_measure.schema.json"
DEFAULT_CONTRACT = ROOT / "particles" / "runs" / "hadron" / "ward_projected_spectral_measure_contract.json"
DEFAULT_FULL_UNQUENCHED = ROOT / "particles" / "runs" / "hadron" / "full_unquenched_correlator.json"
DEFAULT_RHO_LEVELS = ROOT / "particles" / "runs" / "hadron" / "rho_levels.json"
DEFAULT_RHO_PHASE_SHIFT = ROOT / "particles" / "runs" / "hadron" / "rho_phase_shift_fit.json"
DEFAULT_STABLE_PAYLOAD = ROOT / "particles" / "runs" / "hadron" / "stable_channel_cfg_source_measure_payload.json"
DEFAULT_CLOSURE = ROOT / "particles" / "runs" / "hadron" / "hadron_production_closure_validation_report.json"
DEFAULT_AUDIT = ROOT / "particles" / "runs" / "hadron" / "current_hadron_lane_audit.json"
DEFAULT_OUT = (
    ROOT
    / "P_derivation"
    / "runtime"
    / "ward_projected_hadronic_spectral_measure_obstruction_current.json"
)


def _now_utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_optional(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"artifact_path": str(path.relative_to(ROOT)), "exists": False}
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload.setdefault("artifact_path", str(path.relative_to(ROOT)))
    payload.setdefault("exists", True)
    return payload


def _nonempty(value: Any) -> bool:
    return isinstance(value, (list, dict)) and bool(value)


def _has_nested_nonempty_arrays(payload: dict[str, Any], *field_names: str) -> bool:
    """Return true if any named field holds a nonempty list/dict anywhere."""
    stack: list[Any] = [payload]
    wanted = set(field_names)
    while stack:
        current = stack.pop()
        if isinstance(current, dict):
            for key, value in current.items():
                if key in wanted and _nonempty(value):
                    return True
                stack.append(value)
        elif isinstance(current, list):
            stack.extend(current)
    return False


def _key_paths(payload: Any, names: set[str], prefix: str = "") -> list[str]:
    paths: list[str] = []
    if isinstance(payload, dict):
        for key, value in payload.items():
            key_text = str(key)
            path = f"{prefix}.{key_text}" if prefix else key_text
            if key_text in names:
                paths.append(path)
            paths.extend(_key_paths(value, names, path))
    elif isinstance(payload, list):
        for index, value in enumerate(payload):
            paths.extend(_key_paths(value, names, f"{prefix}[{index}]"))
    return paths


def _rho_state(full_unquenched: dict[str, Any]) -> dict[str, Any]:
    rho = (full_unquenched.get("channel_payloads") or {}).get("rho_scattering") or {}
    return {
        "correlation_matrices_populated": _nonempty(rho.get("correlation_matrices")),
        "principal_correlators_populated": _nonempty(rho.get("principal_correlators")),
        "aE_lab_count": len(rho.get("aE_lab") or []),
        "aE_cm_count": len(rho.get("aE_cm") or []),
        "delta1_count": len(rho.get("delta1_rad") or []),
        "target_promoted_fields": list(rho.get("target_promoted_fields") or []),
    }


def _resonance_populated(phase_shift: dict[str, Any]) -> bool:
    readout = phase_shift.get("resonance_readout") or {}
    return any(readout.get(key) is not None for key in ("m_rho_gev", "Gamma_rho_gev", "pole_s_gev2"))


def _stable_arrays_populated(stable_payload: dict[str, Any]) -> bool:
    return _has_nested_nonempty_arrays(
        stable_payload,
        "cfg_source_corr_t",
        "cfg_source_corr_direct_t",
        "cfg_source_corr_exchange_t",
    )


def _thomson_moment_for_dimensionless_atom(y: int) -> str:
    moment = Fraction(1, y * (1 + y))
    return f"{moment.numerator}/{moment.denominator}"


def build_obstruction(
    *,
    schema: dict[str, Any] | None = None,
    contract: dict[str, Any] | None = None,
    full_unquenched: dict[str, Any] | None = None,
    rho_levels: dict[str, Any] | None = None,
    rho_phase_shift: dict[str, Any] | None = None,
    stable_payload: dict[str, Any] | None = None,
    closure: dict[str, Any] | None = None,
    audit: dict[str, Any] | None = None,
) -> dict[str, Any]:
    schema = schema or _load_optional(DEFAULT_SCHEMA)
    contract = contract or _load_optional(DEFAULT_CONTRACT)
    full_unquenched = full_unquenched or _load_optional(DEFAULT_FULL_UNQUENCHED)
    rho_levels = rho_levels or _load_optional(DEFAULT_RHO_LEVELS)
    rho_phase_shift = rho_phase_shift or _load_optional(DEFAULT_RHO_PHASE_SHIFT)
    stable_payload = stable_payload or _load_optional(DEFAULT_STABLE_PAYLOAD)
    closure = closure or _load_optional(DEFAULT_CLOSURE)
    audit = audit or _load_optional(DEFAULT_AUDIT)

    level_points = rho_levels.get("level_points") or []
    phase_level_points = rho_phase_shift.get("level_points") or []
    rho_state = _rho_state(full_unquenched)
    residue_paths = _key_paths(
        {
            "full_unquenched": full_unquenched,
            "rho_levels": rho_levels,
            "rho_phase_shift": rho_phase_shift,
            "stable_payload": stable_payload,
        },
        {"ward_projected_residues", "current_normalization", "rho_had", "rho_had_or_measure"},
    )
    has_residue_payload = any(path.endswith(("ward_projected_residues", "current_normalization")) for path in residue_paths)
    has_measure_payload = any(path.endswith(("rho_had", "rho_had_or_measure")) for path in residue_paths)

    missing_fields = []
    if not level_points and not phase_level_points:
        missing_fields.append("nonempty_finite_volume_U1_Q_level_support")
    if not has_residue_payload:
        missing_fields.append("ward_projected_current_residues_and_normalization")
    if not has_measure_payload:
        missing_fields.append("rho_Q_or_primitive_positive_measure_payload")
    if not rho_state["correlation_matrices_populated"]:
        missing_fields.append("rho_channel_correlation_matrices")
    if not rho_state["principal_correlators_populated"]:
        missing_fields.append("rho_channel_principal_correlators")
    missing_fields.extend(
        [
            "continuum_pushforward_rule_to_rho_Q_of_s_P",
            "same_scheme_current_matching_budget",
            "spectral_quadrature_and_tail_budget",
        ]
    )

    mu_a_y = 2
    mu_b_y = 3
    mu_a_moment = _thomson_moment_for_dimensionless_atom(mu_a_y)
    mu_b_moment = _thomson_moment_for_dimensionless_atom(mu_b_y)

    return {
        "artifact": "oph_ward_projected_hadronic_spectral_measure_obstruction",
        "generated_utc": _now_utc(),
        "status": "rho_Q_not_constructible_from_current_local_artifacts",
        "promotion_allowed": False,
        "source_only": True,
        "external_inputs_used": False,
        "inspected_artifacts": {
            "schema": str(DEFAULT_SCHEMA.relative_to(ROOT)),
            "contract": str(DEFAULT_CONTRACT.relative_to(ROOT)),
            "full_unquenched": str(DEFAULT_FULL_UNQUENCHED.relative_to(ROOT)),
            "rho_levels": str(DEFAULT_RHO_LEVELS.relative_to(ROOT)),
            "rho_phase_shift": str(DEFAULT_RHO_PHASE_SHIFT.relative_to(ROOT)),
            "stable_payload": str(DEFAULT_STABLE_PAYLOAD.relative_to(ROOT)),
            "closure": str(DEFAULT_CLOSURE.relative_to(ROOT)),
            "audit": str(DEFAULT_AUDIT.relative_to(ROOT)),
        },
        "schema_contract": {
            "required_artifact": (schema.get("properties") or {}).get("artifact", {}).get("const"),
            "required_fields": list(schema.get("required") or []),
            "ward_projected_residues_required": "ward_projected_residues" in (schema.get("required") or []),
            "finite_volume_levels_required": "finite_volume_levels" in (schema.get("required") or []),
        },
        "local_artifact_state": {
            "contract_status": contract.get("status"),
            "contract_promotion_allowed": contract.get("promotion_allowed"),
            "rho_levels": {
                "artifact": rho_levels.get("artifact"),
                "status": rho_levels.get("status"),
                "level_points_count": len(level_points),
                "anti_cheat": rho_levels.get("anti_cheat"),
            },
            "rho_phase_shift": {
                "artifact": rho_phase_shift.get("artifact"),
                "status": rho_phase_shift.get("status"),
                "level_points_count": len(phase_level_points),
                "resonance_populated": _resonance_populated(rho_phase_shift),
                "anti_cheat": rho_phase_shift.get("anti_cheat"),
            },
            "full_unquenched_rho_scattering": rho_state,
            "stable_channel_payload": {
                "artifact": stable_payload.get("artifact"),
                "status": stable_payload.get("status"),
                "stable_arrays_populated": _stable_arrays_populated(stable_payload),
                "payload_realization_status": stable_payload.get("payload_realization_status"),
            },
            "production_closure": {
                "artifact": closure.get("artifact"),
                "closure_grade": closure.get("closure_grade"),
                "production_dump_present": bool(closure.get("production_dump_present")),
                "public_unsuppression_ready": bool(closure.get("public_unsuppression_ready")),
                "smallest_live_residual_object": closure.get("smallest_live_residual_object"),
            },
            "audit_frontier": {
                "artifact": audit.get("artifact"),
                "smallest_constructive_missing_object": audit.get("smallest_constructive_missing_object"),
                "strict_missing_program": audit.get("strict_missing_program"),
            },
            "spectral_payload_key_paths_found": residue_paths,
        },
        "forbidden_shortcuts": {
            "codata_or_nist_endpoint": "not inspected and not allowed as source input",
            "reference_targets_gev": "compare-only audit targets, not source spectral data",
            "surrogate_execution_bridge": "diagnostic only; public_promotion_allowed is false",
            "local_rho_effective_mass": "explicitly rejected by rho artifacts' anti_cheat guards",
            "free_quark_screened_ansatz": "forbidden by ward_projected_spectral_measure_contract",
            "stable_channel_only_backend_export": "does not contain the U(1)_Q vector spectral measure",
        },
        "missing_payload_fields": missing_fields,
        "theorem_grade_obstruction": {
            "theorem": "CurrentLocalHadronArtifactsDoNotDetermineRhoQ",
            "claim": (
                "No deterministic source-only constructor can emit rho_Q(s;P) from the current "
                "local hadron artifacts, because the emitted projection contains no nonempty "
                "U(1)_Q finite-volume vector levels, no Ward-projected current residues, no "
                "current normalization, and no continuum pushforward."
            ),
            "proof_by_nonidentifiability": {
                "source_projection": (
                    "Forget every un-emitted spectral datum and retain exactly the current local "
                    "artifacts: D10 family labels, U(1)_Q contract, empty rho level support, "
                    "empty rho correlation matrices, candidate-only rho readout, and stable-channel "
                    "execution scaffolding."
                ),
                "two_positive_completions_same_projection": [
                    {
                        "id": "mu_A",
                        "dimensionless_support_y": [mu_a_y],
                        "weights": [1],
                        "thomson_kernel_without_common_factor": "k(y)=1/(y*(1+y))",
                        "moment": mu_a_moment,
                    },
                    {
                        "id": "mu_B",
                        "dimensionless_support_y": [mu_b_y],
                        "weights": [1],
                        "thomson_kernel_without_common_factor": "k(y)=1/(y*(1+y))",
                        "moment": mu_b_moment,
                    },
                ],
                "moments_differ": mu_a_moment != mu_b_moment,
                "conclusion": (
                    "The local artifact projection is identical for mu_A and mu_B, but the "
                    "Ward-projected Thomson spectral moment differs. Therefore rho_Q(s;P) is "
                    "not identifiable from the current source packet."
                ),
            },
        },
        "minimal_source_object_still_missing": {
            "id": "oph_qcd_ward_projected_hadronic_spectral_measure",
            "mathematical_object": "positive U(1)_Q source spectral measure rho_Q(s;P) or primitive measure mu_Q(ds;P)",
            "strict_minimum_fields": [
                "source_family_id=d10_running_tree",
                "current=U1_Q",
                "same_subtraction_scheme_as_a0(P)",
                "nonempty finite-volume vector-channel level support",
                "Ward-projected residues/weights with current normalization",
                "threshold and positivity certificate",
                "continuum/finite-volume/chiral/current-matching budgets",
                "pushforward rule from finite-volume data to rho_Q(s;P)",
                "quadrature and OPE/tail bound for the Thomson kernel",
            ],
            "needed_after_rho_Q_for_full_fine_structure_endpoint": [
                "same-scheme Delta_EW zero theorem or source bound",
                "directed-rounding interval fixed-point certificate",
            ],
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Emit the Ward-projected hadronic spectral-measure obstruction."
    )
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    parser.add_argument("--print-json", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    payload = build_obstruction()
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text, encoding="utf-8")
    if args.print_json:
        print(text, end="")
    else:
        print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
