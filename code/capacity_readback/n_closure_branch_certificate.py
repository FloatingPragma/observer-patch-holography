#!/usr/bin/env python3
"""Fail-closed arithmetic packet for the two retrospective N-closure branches.

The packet evaluates two different reserve readings on the same conditional
screen/electroweak baseline:

* finite presence: ``N = N0 * (1 - P/24)``;
* Poisson/projective limit: ``N = N0 * exp(-P/24)``.

The finite edge-center receipt types ``1 - P/24`` as the finite one-step
presence value.  It types ``exp(-P/24)`` as a refinement limit, not as the
same finite record.  The second branch therefore carries a separate,
undischarged mean-count or projective-limit carrier premise.

Both rows are retrospective and unselected.  This packet supplies arithmetic,
source pins, and type guards only.  It does not identify either row with the
global capacity, a horizon record count, or a prediction.  In particular,
self-reference forces an equality only after both sides have been proved to
read the same typed quantity.  That common-load and physical-attachment work
remains open here.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from mpmath import mp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
OUTPUT_PATH = HERE / "manifests" / "n_closure_branch_certificate.json"

P_SOURCE_PATH = (
    ROOT
    / "code"
    / "particles"
    / "hierarchy"
    / "certificates"
    / "R_P_source_audit_pixel_certificate.json"
)
FINITE_SEMANTICS_PATH = (
    ROOT / "code" / "cosmology" / "manifests" / "edge_center_clock_certificate.json"
)
MENU_PATH = HERE / "manifests" / "capacity_semantics_menu_reference.json"

SCHEMA = "oph.n_closure_branch_certificate.v1"
ARTIFACT = "oph_n_closure_branch_certificate"
STATUS = "RETROSPECTIVE_UNSELECTED_CONDITIONAL_BRANCH_MENU"
ISSUE = 648
PRECISION = 100

SOURCE_STATUS = "source_audit_branch_witness_not_full_endpoint_proof"
FINITE_SCHEMA = "oph.edge_center_clock_certificate.v3"
FINITE_STATUS = "conditional_edge_center_arithmetic_with_open_source_and_clock_gates"
MENU_SCHEMA = "oph.capacity_semantics_menu_certificate.v1"

FINITE_BRANCH = "finite_presence"
POISSON_BRANCH = "poisson_projective_limit"
BRANCH_IDS = (FINITE_BRANCH, POISSON_BRANCH)


class CertificateError(ValueError):
    """Fail-closed validation error with a stable code."""

    def __init__(self, code: str, message: str) -> None:
        super().__init__(f"{code}: {message}")
        self.code = code
        self.message = message


def require(condition: bool, code: str, message: str) -> None:
    if not condition:
        raise CertificateError(code, message)


def canonical_json_bytes(value: Any) -> bytes:
    return (
        json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + "\n"
    ).encode("utf-8")


def _load_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise CertificateError("SOURCE_READ", f"cannot read {path}: {exc}") from exc
    require(isinstance(value, dict), "SOURCE_TYPE", f"{path} is not a JSON object")
    return value


def _source_pin(path: Path) -> dict[str, Any]:
    try:
        raw = path.read_bytes()
    except OSError as exc:
        raise CertificateError("SOURCE_READ", f"cannot read {path}: {exc}") from exc
    return {
        "path": path.relative_to(ROOT).as_posix(),
        "sha256": "sha256:" + hashlib.sha256(raw).hexdigest(),
        "bytes": len(raw),
    }


def _number(value: Any, label: str) -> Any:
    require(isinstance(value, str), "SOURCE_TYPE", f"{label} must be a decimal string")
    try:
        result = mp.mpf(value)
    except (TypeError, ValueError) as exc:
        raise CertificateError("SOURCE_NUMBER", f"{label} is not numeric") from exc
    require(mp.isfinite(result), "SOURCE_NUMBER", f"{label} is not finite")
    return result


def _decimal(value: Any, digits: int = 80) -> str:
    return mp.nstr(value, digits, strip_zeros=False)


def _load_sources() -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    p_source = _load_json(P_SOURCE_PATH)
    finite = _load_json(FINITE_SEMANTICS_PATH)
    menu = _load_json(MENU_PATH)

    require(
        p_source.get("artifact") == "R_P_source_audit_pixel_certificate",
        "SOURCE_IDENTITY",
        "P/alpha_U source artifact drift",
    )
    require(
        p_source.get("status") == SOURCE_STATUS,
        "SOURCE_STATUS",
        "P/alpha_U source is not the declared unpromoted branch witness",
    )
    require(
        finite.get("schema") == FINITE_SCHEMA and finite.get("status") == FINITE_STATUS,
        "SOURCE_STATUS",
        "finite survival receipt status drift",
    )
    require(
        finite.get("generator", {}).get("operational_clock_bound") is False,
        "TYPE_CONFUSION",
        "finite survival receipt unexpectedly promotes an operational clock",
    )
    presence_statement = (
        finite.get("antecedents", {}).get("presence_reading", {}).get("statement")
    )
    require(
        isinstance(presence_statement, str)
        and "finite one-step survival is the presence value 1 - P/24" in presence_statement
        and "e^(-P/24) is the depth limit" in presence_statement,
        "TYPE_CONFUSION",
        "finite and projective-limit semantics are not separated by the source receipt",
    )
    require(
        menu.get("schema") == MENU_SCHEMA,
        "SOURCE_IDENTITY",
        "capacity semantics menu schema drift",
    )
    require(
        menu.get("campaign_verdict", {}).get("source_only_fixed_point_selector")
        == "conditional_open_interface",
        "SOURCE_STATUS",
        "capacity semantics menu unexpectedly selects a source-only branch",
    )
    reserve_axes = [
        row
        for row in menu.get("semantic_axes", [])
        if isinstance(row, dict) and row.get("axis") == "reserve_semantics"
    ]
    require(len(reserve_axes) == 1, "SOURCE_TYPE", "reserve-semantics axis is missing")
    alternatives = reserve_axes[0].get("alternatives")
    require(
        isinstance(alternatives, list)
        and any("exp(-P/24)" in str(row) for row in alternatives)
        and any("1 - P/24" in str(row) for row in alternatives),
        "TYPE_CONFUSION",
        "source menu does not retain both reserve readings",
    )
    return p_source, finite, menu


def build() -> dict[str, Any]:
    """Build the deterministic retrospective branch packet."""
    mp.dps = PRECISION
    p_source, finite, _menu = _load_sources()

    P = _number(p_source.get("P_cand"), "P_cand")
    alpha_u = _number(p_source.get("alpha_U_P_cand"), "alpha_U_P_cand")
    require(P > 0, "DOMAIN", "P must be positive")
    require(alpha_u > 0, "DOMAIN", "alpha_U must be positive")

    reserve = P / 24
    finite_factor = 1 - reserve
    poisson_factor = mp.exp(-reserve)
    require(finite_factor > 0, "DOMAIN", "finite presence factor must be positive")
    require(
        finite_factor < poisson_factor < 1,
        "BRANCH_SEPARATION",
        "finite-presence and Poisson factors are not strictly separated",
    )

    log_n0_over_pi = 6 * mp.pi / (P * alpha_u)
    n0 = mp.pi * mp.exp(log_n0_over_pi)
    finite_n = n0 * finite_factor
    poisson_n = n0 * poisson_factor

    sources = [
        _source_pin(P_SOURCE_PATH),
        _source_pin(FINITE_SEMANTICS_PATH),
        _source_pin(MENU_PATH),
    ]

    shared_open_gates = {
        "inherited_screen_electroweak_bridge_premises": "open",
        "common_load_identification": "open",
        "global_capacity_attachment": "open",
        "horizon_record_identity": "open",
        "physical_repair_law_attachment": "open",
        "one_class_reserve_attachment": "open",
        "scalar_weighted_reserve_receipt": "open",
    }
    branches = [
        {
            "branch_id": FINITE_BRANCH,
            "status": "conditional_unselected_retrospective",
            "selected": False,
            "source_type": "finite_one_step_presence",
            "formula": "N = N0 * (1 - P/24)",
            "factor": _decimal(finite_factor),
            "N": _decimal(finite_n),
            "log_N_over_pi": _decimal(mp.log(finite_n / mp.pi)),
            "mean_count_or_projective_limit_carrier_required": False,
            "source_semantics": (
                "the finite edge-center receipt types 1 - P/24 as the "
                "finite one-step presence value"
            ),
            "open_gates": shared_open_gates,
        },
        {
            "branch_id": POISSON_BRANCH,
            "status": "conditional_unselected_retrospective",
            "selected": False,
            "source_type": "poisson_or_projective_limit",
            "formula": "N = N0 * exp(-P/24)",
            "factor": _decimal(poisson_factor),
            "N": _decimal(poisson_n),
            "log_N_over_pi": _decimal(mp.log(poisson_n / mp.pi)),
            "mean_count_or_projective_limit_carrier_required": True,
            "source_semantics": (
                "the finite edge-center receipt supplies exp(-P/24) only as "
                "a refinement limit; consuming it as a capacity factor "
                "requires a separate mean-count or projective-limit carrier"
            ),
            "additional_open_gate": {
                "mean_count_or_projective_limit_carrier": "open"
            },
            "open_gates": shared_open_gates,
        },
    ]

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "artifact": ARTIFACT,
        "issue": ISSUE,
        "status": STATUS,
        "scope": {
            "retrospective": True,
            "target_blind_forecast": False,
            "branch_selected": False,
            "comparison_data_consumed": False,
            "global_capacity_derived": False,
            "horizon_attachment_derived": False,
            "prediction_promoted": False,
        },
        "self_reference_boundary": {
            "principle": (
                "self-reference forces equality after the simulator-side and "
                "simulated-side readings are proved to be the same typed quantity"
            ),
            "same_typed_quantity_identified": False,
            "source_return_map_physically_attached": False,
            "uniqueness_from_self_reference_alone": False,
        },
        "source_pins": sources,
        "source_inputs": {
            "P": _decimal(P),
            "alpha_U": _decimal(alpha_u),
            "source_artifact_status": p_source["status"],
            "source_artifact_scope": (
                "source-audit branch witness, not a full endpoint proof"
            ),
        },
        "conditional_baseline": {
            "formula": "N0 = pi * exp(6*pi/(P*alpha_U))",
            "log_N0_over_pi": _decimal(log_n0_over_pi),
            "N0": _decimal(n0),
            "physical_status": (
                "conditional screen/electroweak bridge; inherited premises "
                "and common-load identification remain open"
            ),
        },
        "reserve_coordinate": {
            "P_over_24": _decimal(reserve),
            "finite_factor_strictly_below_poisson_factor": True,
        },
        "branches": branches,
        "promotion_controls": {
            "global_capacity": False,
            "horizon_record_count": False,
            "physical_cosmological_constant": False,
            "prediction": False,
            "prospective_freeze": False,
        },
        "verdict": (
            "ARITHMETIC_REPRODUCED_TWO_RETROSPECTIVE_BRANCHES_UNSELECTED_"
            "PHYSICAL_ATTACHMENTS_OPEN"
        ),
    }
    validate(payload, verify_sources=True)
    return payload


def _branch_map(payload: Mapping[str, Any]) -> dict[str, Mapping[str, Any]]:
    rows = payload.get("branches")
    require(isinstance(rows, list), "PAYLOAD_TYPE", "branches must be a list")
    require(
        all(isinstance(row, Mapping) for row in rows),
        "PAYLOAD_TYPE",
        "every branch must be an object",
    )
    result = {str(row.get("branch_id")): row for row in rows}
    require(
        tuple(row.get("branch_id") for row in rows) == BRANCH_IDS
        and set(result) == set(BRANCH_IDS),
        "BRANCH_SET",
        "branch order or identity drift",
    )
    return result


def validate(payload: Mapping[str, Any], *, verify_sources: bool = True) -> None:
    """Validate status, source pins, arithmetic, and semantic separation."""
    mp.dps = PRECISION
    require(isinstance(payload, Mapping), "PAYLOAD_TYPE", "payload is not an object")
    require(payload.get("schema") == SCHEMA, "SCHEMA", "schema drift")
    require(payload.get("artifact") == ARTIFACT, "SCHEMA", "artifact drift")
    require(payload.get("issue") == ISSUE, "SCHEMA", "issue drift")
    require(payload.get("status") == STATUS, "STATUS", "status drift")

    scope = payload.get("scope")
    require(isinstance(scope, Mapping), "PAYLOAD_TYPE", "scope is missing")
    expected_scope = {
        "retrospective": True,
        "target_blind_forecast": False,
        "branch_selected": False,
        "comparison_data_consumed": False,
        "global_capacity_derived": False,
        "horizon_attachment_derived": False,
        "prediction_promoted": False,
    }
    require(dict(scope) == expected_scope, "PROMOTION", "scope or promotion drift")

    boundary = payload.get("self_reference_boundary")
    require(isinstance(boundary, Mapping), "PAYLOAD_TYPE", "self-reference boundary missing")
    for key in (
        "same_typed_quantity_identified",
        "source_return_map_physically_attached",
        "uniqueness_from_self_reference_alone",
    ):
        require(boundary.get(key) is False, "PROMOTION", f"{key} was promoted")

    pins = payload.get("source_pins")
    require(isinstance(pins, list) and len(pins) == 3, "SOURCE_PINS", "source pins incomplete")
    paths = [P_SOURCE_PATH, FINITE_SEMANTICS_PATH, MENU_PATH]
    expected_pins = [_source_pin(path) for path in paths]
    require(pins == expected_pins, "SOURCE_PINS", "source pins do not match consumed bytes")
    if verify_sources:
        _load_sources()

    source_inputs = payload.get("source_inputs")
    require(isinstance(source_inputs, Mapping), "PAYLOAD_TYPE", "source inputs missing")
    require(
        source_inputs.get("source_artifact_status") == SOURCE_STATUS,
        "SOURCE_STATUS",
        "source input promotion drift",
    )
    P = _number(source_inputs.get("P"), "source_inputs.P")
    alpha_u = _number(source_inputs.get("alpha_U"), "source_inputs.alpha_U")
    require(P > 0 and alpha_u > 0, "DOMAIN", "P and alpha_U must be positive")

    baseline = payload.get("conditional_baseline")
    require(isinstance(baseline, Mapping), "PAYLOAD_TYPE", "baseline missing")
    reserve_record = payload.get("reserve_coordinate")
    require(isinstance(reserve_record, Mapping), "PAYLOAD_TYPE", "reserve coordinate missing")
    branches = _branch_map(payload)

    reserve = P / 24
    finite_factor = 1 - reserve
    poisson_factor = mp.exp(-reserve)
    log_n0 = 6 * mp.pi / (P * alpha_u)
    n0 = mp.pi * mp.exp(log_n0)
    expected_numbers = {
        ("reserve_coordinate", "P_over_24"): _decimal(reserve),
        ("conditional_baseline", "log_N0_over_pi"): _decimal(log_n0),
        ("conditional_baseline", "N0"): _decimal(n0),
        (FINITE_BRANCH, "factor"): _decimal(finite_factor),
        (FINITE_BRANCH, "N"): _decimal(n0 * finite_factor),
        (FINITE_BRANCH, "log_N_over_pi"): _decimal(mp.log(n0 * finite_factor / mp.pi)),
        (POISSON_BRANCH, "factor"): _decimal(poisson_factor),
        (POISSON_BRANCH, "N"): _decimal(n0 * poisson_factor),
        (POISSON_BRANCH, "log_N_over_pi"): _decimal(mp.log(n0 * poisson_factor / mp.pi)),
    }
    for (owner, key), expected in expected_numbers.items():
        if owner == "reserve_coordinate":
            actual = reserve_record.get(key)
        elif owner == "conditional_baseline":
            actual = baseline.get(key)
        else:
            actual = branches[owner].get(key)
        require(actual == expected, "ARITHMETIC", f"{owner}.{key} drift")

    require(
        reserve_record.get("finite_factor_strictly_below_poisson_factor") is True
        and finite_factor < poisson_factor < 1,
        "BRANCH_SEPARATION",
        "branch-separation control failed",
    )
    finite = branches[FINITE_BRANCH]
    poisson = branches[POISSON_BRANCH]
    for branch_id, row in branches.items():
        require(row.get("selected") is False, "SELECTION", f"{branch_id} was selected")
        require(
            row.get("status") == "conditional_unselected_retrospective",
            "SELECTION",
            f"{branch_id} status drift",
        )
        gates = row.get("open_gates")
        require(
            isinstance(gates, Mapping)
            and gates
            and all(value == "open" for value in gates.values()),
            "PROMOTION",
            f"{branch_id} physical gate was discharged",
        )
    require(
        finite.get("source_type") == "finite_one_step_presence"
        and finite.get("mean_count_or_projective_limit_carrier_required") is False,
        "TYPE_CONFUSION",
        "finite-presence row changed type",
    )
    require(
        poisson.get("source_type") == "poisson_or_projective_limit"
        and poisson.get("mean_count_or_projective_limit_carrier_required") is True
        and poisson.get("additional_open_gate")
        == {"mean_count_or_projective_limit_carrier": "open"},
        "TYPE_CONFUSION",
        "Poisson row lost its separate carrier premise",
    )

    controls = payload.get("promotion_controls")
    require(
        isinstance(controls, Mapping)
        and controls
        and all(value is False for value in controls.values()),
        "PROMOTION",
        "a physical or prediction promotion was enabled",
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="validate the on-disk output and source pins without rewriting it",
    )
    args = parser.parse_args()

    if args.validate_only:
        payload = _load_json(args.output)
        validate(payload, verify_sources=True)
        print("N_CLOSURE_BRANCH_CERTIFICATE_VALID")
        return 0

    payload = build()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(canonical_json_bytes(payload))
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
