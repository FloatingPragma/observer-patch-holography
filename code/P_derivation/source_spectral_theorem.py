#!/usr/bin/env python3
"""Emit the Ward-projected source-spectral theorem artifact.

The theorem closes the formal reduction from a source-emitted electromagnetic
spectral measure to the Thomson endpoint. It does not manufacture the
spectral measure from comparison data or from the residual scalar.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from fractions import Fraction
import json
from pathlib import Path
from typing import Any

from thomson_spectral_transport import (
    THEOREM_GRADE_REQUIRED_FIELDS,
    blocked_missing_source_transport,
    validate_source_transport_payload,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "P_derivation" / "runtime" / "source_spectral_theorem_current.json"
DEFAULT_HADRON_CONTRACT = (
    ROOT / "particles" / "runs" / "hadron" / "ward_projected_spectral_measure_contract.json"
)
DEFAULT_HADRON_CLOSURE = (
    ROOT / "particles" / "runs" / "hadron" / "hadron_production_closure_validation_report.json"
)
DEFAULT_RHO_LEVELS = ROOT / "particles" / "runs" / "hadron" / "rho_levels.json"
DEFAULT_SCHEMA = ROOT / "particles" / "hadron" / "ward_projected_spectral_measure.schema.json"


def _now_utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_optional(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"artifact_path": path.relative_to(ROOT).as_posix(), "exists": False}
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload.setdefault("artifact_path", path.relative_to(ROOT).as_posix())
    payload.setdefault("exists", True)
    return payload


def _bool_from_path(payload: dict[str, Any], *keys: str) -> bool:
    cursor: Any = payload
    for key in keys:
        if not isinstance(cursor, dict):
            return False
        cursor = cursor.get(key)
    return bool(cursor)


def _levels_populated(rho_levels: dict[str, Any]) -> bool:
    level_points = rho_levels.get("level_points")
    return isinstance(level_points, list) and len(level_points) > 0


def _fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def _dimensionless_thomson_kernel(y: Fraction) -> Fraction:
    return Fraction(1, 1) / (y * (y + 1))


def _constructive_nonidentifiability_witness() -> dict[str, Any]:
    """Return two current-corpus-compatible measures with unequal moments."""
    weight = Fraction(1, 1)
    support_a = Fraction(2, 1)
    support_b = Fraction(3, 1)
    moment_a = weight * _dimensionless_thomson_kernel(support_a)
    moment_b = weight * _dimensionless_thomson_kernel(support_b)
    moment_difference = moment_a - moment_b

    common_projection = {
        "source_family_id": "d10_running_tree",
        "current": "U1_Q",
        "same_anchor_a0": True,
        "same_mZ": True,
        "same_lepton_kernel": True,
        "same_currently_emitted_hadron_fields": [],
        "same_total_positive_weight": _fraction_text(weight),
        "source_measure_payload_populated": False,
    }

    return {
        "status": "constructive_counterexample_emitted",
        "logical_form": (
            "There exist two positive Ward-projected spectral-measure completions whose projection "
            "to every invariant emitted by the current corpus is identical, while their Thomson "
            "dispersion moments differ."
        ),
        "normalization": "Use y=s/mZ(P)^2 and omit the common positive factor 1/(3*pi*mZ(P)^2).",
        "kernel": "k(y)=1/(y*(1+y))",
        "common_current_corpus_projection": common_projection,
        "measure_a": {
            "description": "one positive atom at y=2",
            "atoms": [{"support_y": _fraction_text(support_a), "weight": _fraction_text(weight)}],
            "positive": True,
        },
        "measure_b": {
            "description": "one positive atom at y=3",
            "atoms": [{"support_y": _fraction_text(support_b), "weight": _fraction_text(weight)}],
            "positive": True,
        },
        "dimensionless_thomson_moment_a": _fraction_text(moment_a),
        "dimensionless_thomson_moment_b": _fraction_text(moment_b),
        "dimensionless_thomson_moment_difference": _fraction_text(moment_difference),
        "moments_equal": moment_a == moment_b,
        "external_inputs_used": False,
        "why_this_blocks_exact_alpha": (
            "The current corpus emits no finite-volume vector levels, Ward-projected residues, "
            "continuum pushforward, or same-scheme source remainder. The projection that the corpus "
            "can see is therefore identical for these measures, but the endpoint functional is not."
        ),
    }


def _current_corpus_obstruction(
    *,
    hadron_contract: dict[str, Any],
    hadron_closure: dict[str, Any],
    rho_levels: dict[str, Any],
) -> dict[str, Any]:
    return {
        "hadron_contract_status": hadron_contract.get("status", "missing"),
        "production_dump_present": _bool_from_path(hadron_closure, "production_dump_present"),
        "public_unsuppression_ready": _bool_from_path(hadron_closure, "public_unsuppression_ready"),
        "finite_volume_levels_populated": _levels_populated(rho_levels),
        "ward_projected_residues_populated": False,
        "source_measure_payload_populated": False,
        "smallest_live_residual_object": hadron_closure.get(
            "smallest_live_residual_object",
            "oph_qcd_ward_projected_hadronic_spectral_measure",
        ),
    }


def build_source_spectral_theorem(
    *,
    source_transport_payload: dict[str, Any] | None = None,
    hadron_contract: dict[str, Any] | None = None,
    hadron_closure: dict[str, Any] | None = None,
    rho_levels: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Build the theorem artifact and evaluate any supplied source payload."""
    hadron_contract = hadron_contract or _load_optional(DEFAULT_HADRON_CONTRACT)
    hadron_closure = hadron_closure or _load_optional(DEFAULT_HADRON_CLOSURE)
    rho_levels = rho_levels or _load_optional(DEFAULT_RHO_LEVELS)

    validation = (
        validate_source_transport_payload(source_transport_payload).to_json()
        if source_transport_payload is not None
        else blocked_missing_source_transport()
    )
    source_payload_supplied = source_transport_payload is not None
    source_payload_accepted = bool(validation["promotion_allowed"])
    obstruction = _current_corpus_obstruction(
        hadron_contract=hadron_contract,
        hadron_closure=hadron_closure,
        rho_levels=rho_levels,
    )

    if source_payload_accepted:
        status = "source_spectral_payload_contract_satisfied"
    else:
        status = "source_spectral_reduction_theorem_emitted_measure_payload_absent"

    return {
        "artifact": "oph_ward_projected_source_spectral_theorem",
        "generated_utc": _now_utc(),
        "github_issue": 235,
        "theorem_id": "WardProjectedHadronicSpectralEmission_Q",
        "status": status,
        "promotion_allowed": source_payload_accepted,
        "external_inputs_used": False,
        "source_payload_supplied": source_payload_supplied,
        "source_payload_validation": validation,
        "source_only_guard": {
            "codata_allowed": False,
            "measured_alpha_allowed": False,
            "residual_fit_allowed": False,
            "free_quark_screen_promotable": False,
            "accepted_artifact": "oph_source_ward_projected_thomson_transport",
        },
        "closed_reduction_theorem": {
            "statement": (
                "For a fixed D10 source family and the Ward-projected charge current Q=T3+Y, "
                "a positive source-emitted U(1)_Q spectral measure plus a same-scheme electroweak "
                "remainder determines the Thomson endpoint by the subtracted dispersion functional."
            ),
            "ward_projected_correlator": (
                "Pi_Q^{mu nu}(q;P)=(q^mu q^nu-q^2 eta^{mu nu}) Pi_Q(q^2;P)"
            ),
            "hadronic_transport": (
                "Delta_had(P)=mZ(P)^2/(3*pi) * integral rho_Q(s;P)/(s*(s+mZ(P)^2)) ds"
            ),
            "endpoint_map": (
                "A_T(P)=a0(P)+Delta_lep_src(P)+Delta_had_src(P)+Delta_EW_src(P)"
            ),
            "fixed_point_map": "P=phi+sqrt(pi)/A_T(P)",
        },
        "accepted_source_payload_contract": {
            "schema": DEFAULT_SCHEMA.relative_to(ROOT).as_posix(),
            "required_artifact": "oph_qcd_ward_projected_hadronic_spectral_measure",
            "required_transport_wrapper": "oph_source_ward_projected_thomson_transport",
            "required_fields": list(THEOREM_GRADE_REQUIRED_FIELDS),
        },
        "current_corpus_obstruction": obstruction,
        "nonidentifiability_corollary": {
            "status": "constructive_no_external_input_no_go_closed_for_current_source_packet",
            "statement": (
                "The current D10 invariant packet fixes the mZ anchor, lepton kernel, naive quark "
                "kernel, and first-order screen. It contains no finite-volume vector levels, "
                "current residues, or continuum pushforward for rho_Q(s;P). Distinct positive "
                "spectral measures can share that packet and give different Thomson moments."
            ),
            "constructive_witness": _constructive_nonidentifiability_witness(),
            "fitted_scalars_rejected": [
                "missing_source_transport_delta_alpha_inv",
                "required_screening_factor",
                "residual_second_order_coefficient",
                "c_Q",
                "S_required",
            ],
        },
        "conclusion": {
            "source_spectral_reduction_closed": True,
            "source_transport_payload_accepted": source_payload_accepted,
            "exact_alpha_emitted_from_current_corpus": False,
            "exact_alpha_promotion_allowed_for_supplied_payload": source_payload_accepted,
            "reason": (
                "The formal Ward-projected spectral reduction is emitted. The numerical fine-structure "
                "endpoint requires a populated source spectral measure and same-scheme remainder."
            ),
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Emit the Ward-projected source-spectral theorem.")
    parser.add_argument("--source-transport-payload", help="Optional source transport JSON payload.")
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    parser.add_argument("--print-json", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    source_payload = None
    if args.source_transport_payload:
        source_payload = json.loads(Path(args.source_transport_payload).read_text(encoding="utf-8"))
    payload = build_source_spectral_theorem(source_transport_payload=source_payload)
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
