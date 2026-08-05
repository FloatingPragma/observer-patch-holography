#!/usr/bin/env python3
"""Build the evidence-derived H0 source-only closure preflight.

The inventory is reconstructed from independent registered sources.  In
particular, the capacity inventory comes from the seventeen-row semantics
menu rather than from an order copied into this module's contract.  Historical
and target-exposed hierarchy packets are included with explicit dispositions.

This certificate reports progress over registered packets.  It does not
exhaust future source laws, open a comparison, attach a laboratory quantity,
or close issue 708.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from collections import Counter
from decimal import Decimal
from pathlib import Path
from typing import Any, Mapping

import rc_load_source_only


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent

CONTRACT = HERE / "data/source_only_closure_contract_v3.json"
RC_PROVENANCE = HERE / "data/rc_load_workspace_provenance_v1.json"
RC_SOURCE = HERE / "rc_load_source_only.py"
RC_OUTPUT = HERE / "outputs/rc_load_source_only_certificate.json"
OUTPUT = HERE / "outputs/source_only_closure_preflight.json"
PRODUCER_SOURCE = HERE / "source_only_closure_preflight.py"
INDEPENDENT_VERIFIER = HERE / "verify_source_only_closure_preflight_independent.py"

P_INTERVAL = ROOT / "code/P_derivation/runtime/p_interval_contraction_certificate_2026-07-14.json"
P_TRUNK = ROOT / "code/P_derivation/runtime/p_closure_trunk_current.json"
CAP_K = ROOT / "code/capacity_readback/runtime/F_candidate_capK_certificates.json"
CAP_P = ROOT / "code/capacity_readback/runtime/F_candidate_capP_certificates.json"
CAP_L = ROOT / "code/capacity_readback/runtime/F_candidate_capL_certificates.json"
COUPLED = ROOT / "code/capacity_readback/runtime/F_candidate_coupled_certificates.json"
DIRECT_N = ROOT / "code/capacity_readback/runtime/direct_n_closure_verdict.json"
N_BRANCH = ROOT / "code/capacity_readback/manifests/n_closure_branch_certificate.json"
NAMED_N = ROOT / "code/capacity_readback/runtime/named_law_n_closure_verdict.json"
SEMANTICS_MENU = ROOT / "code/capacity_readback/manifests/capacity_semantics_menu_reference.json"
ARCHIVE_190 = ROOT / "code/capacity_readback/runtime/F_construction_comparison_2026-07-14.json"
EXPOSURE = ROOT / "code/particles/forecast_contract/data/candidate_inventory_v2.json"
FORECAST_STATE = ROOT / "code/particles/forecast_contract/outputs/forecast_contract_state.json"

HIERARCHY_DIR = ROOT / "code/particles/hierarchy/certificates"
HIER_P_PUBLIC = HIERARCHY_DIR / "R_P_public_pixel_certificate.json"
HIER_P_SOURCE = HIERARCHY_DIR / "R_P_source_audit_pixel_certificate.json"
HIER_PN_JOINT = HIERARCHY_DIR / "R_PN_joint_fixed_point_certificate_report.json"
HIER_N_TICK = HIERARCHY_DIR / "R_N_global_repair_tick_certificate.json"
HIER_EW_CAPACITY = HIERARCHY_DIR / "R_EW_global_capacity_certificate.json"

SCHEMA = "oph.source_only_closure_preflight.v3"
CONTRACT_SCHEMA = "oph.source_only_closure_preflight_contract.v3"
ISSUE = 708
P_MODES = (
    "thomson_structured_running",
    "thomson_structured_running_plus_gauge_width",
)
CAPACITY_FAMILY_COUNTS = {"CAP-K": 4, "CAP-P": 6, "CAP-L": 5, "CAP-B": 1, "coupled": 1}
QUALIFYING_ROW_KINDS = {"candidate_map", "candidate_family"}
QUALIFYING_EXPOSURES = {"source_only_precomparison", "sealed_prospective"}
P_SELECTION_ALLOWLIST = {
    "ELIGIBLE_SOURCE_VISIBLE",
}
RC_EQUATION = "ln(N/pi) * (6*pi/(P*alpha_U) - ln(N/pi)) = 6*pi"
COMMON_LOAD_FORMULA = "N0 = pi * exp(6*pi/(P*alpha_U))"
RESERVE_FORMULAS = {
    "finite_presence": "N = N0 * (1 - P/24)",
    "poisson_projective_limit": "N = N0 * exp(-P/24)",
}


def canonical_json_bytes(value: Any) -> bytes:
    return (
        json.dumps(value, ensure_ascii=True, separators=(",", ":"), sort_keys=True)
        + "\n"
    ).encode("ascii")


def tagged(raw: bytes) -> str:
    return "sha256:" + hashlib.sha256(raw).hexdigest()


def load(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    require(isinstance(value, dict), f"JSON root is not an object: {path}")
    return value


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def pointer(artifact: Path, json_pointer: str, observed: Any) -> dict[str, Any]:
    return {
        "artifact": artifact.relative_to(ROOT).as_posix(),
        "json_pointer": json_pointer,
        "observed": observed,
    }


def pin(path: Path, raw: bytes | None = None) -> dict[str, Any]:
    payload = path.read_bytes() if raw is None else raw
    return {
        "path": path.relative_to(ROOT).as_posix(),
        "bytes": len(payload),
        "sha256": tagged(payload),
    }


def gate(attained: bool, evidence: list[dict[str, Any]]) -> dict[str, Any]:
    require(bool(evidence), "gate evidence cannot be empty")
    return {
        "attained": bool(attained),
        "state": "ATTAINED" if attained else "OPEN",
        "evidence": evidence,
    }


def _run_validator(label: str, command: list[str], *, cwd: Path = ROOT) -> dict[str, str]:
    result = subprocess.run(
        command,
        cwd=cwd,
        text=True,
        capture_output=True,
        check=False,
    )
    require(
        result.returncode == 0,
        f"parent validator failed ({label}): {result.stdout}{result.stderr}",
    )
    return {"validator": label, "status": "PASS"}


def validate_parent_artifacts() -> list[dict[str, str]]:
    """Replay cheap parent validators before issuing the H0 projection."""

    py = sys.executable
    return [
        _run_validator(
            "capacity_semantics_menu",
            [py, str(ROOT / "code/capacity_readback/capacity_semantics_menu_certificate.py"), "--verify"],
        ),
        _run_validator(
            "direct_n_closure",
            [py, str(ROOT / "code/capacity_readback/direct_n_closure_verdict.py"), "--verify"],
        ),
        _run_validator(
            "bounded_packet_lift_independent",
            [py, str(ROOT / "code/capacity_readback/verify_complete_packet_lift_independent.py")],
        ),
        _run_validator(
            "n_closure_branch",
            [py, str(ROOT / "code/capacity_readback/n_closure_branch_certificate.py"), "--validate-only"],
        ),
        _run_validator(
            "named_law_n_closure",
            [py, str(ROOT / "code/capacity_readback/named_law_n_closure_verdict.py"), "--verify"],
        ),
        _run_validator(
            "forecast_contract_state",
            [
                py,
                str(ROOT / "code/particles/forecast_contract/forecast_contract.py"),
                "--verify",
            ],
        ),
        _run_validator(
            "hierarchy_joint_fixed_point",
            [
                py,
                str(ROOT / "code/particles/hierarchy/validators/validate_joint_fixed_point_certificate.py"),
                str(HIER_PN_JOINT),
            ],
        ),
        _run_validator(
            "hierarchy_global_repair_tick",
            [
                py,
                str(ROOT / "code/particles/hierarchy/validators/validate_global_repair_tick_certificate.py"),
                str(HIER_N_TICK),
            ],
        ),
        _run_validator(
            "hierarchy_ew_capacity",
            [
                py,
                str(ROOT / "code/particles/hierarchy/validators/validate_issue_344_exact_capacity.py"),
                str(HIER_EW_CAPACITY),
            ],
        ),
    ]


def _load_contract() -> dict[str, Any]:
    value = load(CONTRACT)
    require(value.get("schema") == CONTRACT_SCHEMA, "contract schema drift")
    require(value.get("issue") == ISSUE, "contract issue drift")
    boundary = value.get("scope_boundary", {})
    require(boundary.get("future_source_laws_exhausted") is False, "future-law exhaustion enabled")
    require(boundary.get("issue_closure_authorized") is False, "issue closure enabled")
    require(
        boundary.get("physical_or_laboratory_attachment_required") is False,
        "downstream attachment entered H0",
    )
    required = value.get("source_only_required_gates", [])
    completion = value.get("preflight_completion_gates", [])
    completion_state = value.get("preflight_completion_state", {})
    require(len(required) == len(set(required)) == 7, "source-only gate inventory drift")
    require(len(completion) == len(set(completion)) == 5, "completion gate inventory drift")
    require(list(completion_state) == completion, "completion gate state order drift")
    require(
        all(
            isinstance(completion_state[name], dict)
            and isinstance(completion_state[name].get("attained"), bool)
            and isinstance(completion_state[name].get("reason"), str)
            and bool(completion_state[name]["reason"])
            for name in completion
        ),
        "completion gate state malformed",
    )
    require(
        [completion_state[name]["attained"] for name in completion]
        == [False, False, False, False, True],
        "completion gate evidence changed without an executable gate upgrade",
    )
    require(
        not set(required).intersection(value.get("explicitly_downstream_gates", [])),
        "downstream gate entered source-only qualification",
    )
    decision = value.get("decision_rule", {})
    require(decision.get("frozen") is True, "decision rule is not frozen")
    require(
        decision.get("archive_enumeration_is_not_a_source_return_map") is True,
        "archive/return-map boundary missing",
    )
    return value


def _row(
    *,
    candidate_id: str,
    registry_source: str,
    row_kind: str,
    quantity: str,
    family: str,
    mathematical_status: str,
    disposition: str,
    exposure_class: str,
    solver_target_payload_read: bool,
    exact_fixed_point_progress: bool,
    gates: dict[str, dict[str, Any]],
    details: dict[str, Any],
    required: list[str],
) -> dict[str, Any]:
    require(set(gates) == set(required), f"row gate drift: {candidate_id}")
    attained = sum(gates[name]["attained"] is True for name in required)
    return {
        "candidate_id": candidate_id,
        "registry_source": registry_source,
        "row_kind": row_kind,
        "quantity": quantity,
        "family": family,
        "mathematical_status": mathematical_status,
        "disposition": disposition,
        "exposure_class": exposure_class,
        "solver_target_payload_read": bool(solver_target_payload_read),
        "exact_fixed_point_progress": bool(exact_fixed_point_progress),
        "source_only_gates": gates,
        "attained_gate_count": attained,
        "row_gate_complete": attained == len(required),
        "qualifies_source_only": False,
        "open_source_only_gates": [
            name for name in required if gates[name]["attained"] is not True
        ],
        "details": details,
    }


def _finalize_qualification(
    rows: list[dict[str, Any]], required: list[str], completion_complete: bool
) -> list[str]:
    qualifying: list[str] = []
    for row in rows:
        typed_exposure = row["exposure_class"] in QUALIFYING_EXPOSURES
        eligible = (
            row["row_kind"] in QUALIFYING_ROW_KINDS
            and row["solver_target_payload_read"] is False
            and typed_exposure
            and row["row_gate_complete"] is True
            and completion_complete
        )
        row["qualifies_source_only"] = eligible
        if eligible:
            qualifying.append(row["candidate_id"])
    return qualifying


def derive_progress_status(
    rows: list[dict[str, Any]],
    required_gates: list[str],
    completion_gates: Mapping[str, Mapping[str, Any]],
) -> str:
    completion_complete = all(
        completion_gates[name]["attained"] is True for name in completion_gates
    )
    eligible_gate_complete = any(
        row.get("row_kind") in QUALIFYING_ROW_KINDS
        and row.get("solver_target_payload_read") is False
        and row.get("exposure_class") in QUALIFYING_EXPOSURES
        and all(
            row.get("source_only_gates", {}).get(name, {}).get("attained") is True
            for name in required_gates
        )
        for row in rows
    )
    if eligible_gate_complete and completion_complete:
        return "SOURCE_ONLY_CANDIDATE_PRESENT"
    if eligible_gate_complete:
        return "OPEN_GATE_COMPLETE_CANDIDATE_BLOCKED_BY_PREFLIGHT_COMPLETENESS"
    if any(row.get("exact_fixed_point_progress") is True for row in rows):
        return "OPEN_REGISTERED_CANDIDATES_FAIL_SOURCE_ONLY_GATES"
    return "OPEN_NO_FIXED_POINT_EVIDENCE"


def _p_rows(
    p_interval: dict[str, Any],
    p_trunk: dict[str, Any],
    exposure: dict[str, Any],
    required: list[str],
) -> list[dict[str, Any]]:
    require(
        p_interval.get("claim_status")
        == "interval_contraction_certificate_for_declared_closure_map",
        "P interval status drift",
    )
    require(tuple(p_interval.get("modes", {})) == P_MODES, "P mode registry drift")
    consumer = p_trunk.get("consumer_policy", {})
    inventory = exposure.get("candidates", {}).get("x_p_certificate_alpha", {})
    eligibility = inventory.get("eligibility")
    freeze = inventory.get("freeze_evidence")
    freeze_exposure = freeze.get("exposure_class") if isinstance(freeze, dict) else None
    selected = bool(
        eligibility in P_SELECTION_ALLOWLIST
        and isinstance(freeze, dict)
        and freeze.get("frozen_before_comparison") is True
        and freeze_exposure in QUALIFYING_EXPOSURES
        and isinstance(freeze.get("content_sha256"), str)
        and str(freeze["content_sha256"]).startswith("sha256:")
        and len(str(freeze["content_sha256"])) == 71
    )
    same_quantity = inventory.get("same_quantity_constructed") is True
    source_complete = bool(
        p_interval.get("promotion_allowed") is True
        and consumer.get("may_feed_live_particle_predictions") is True
        and consumer.get("default_thomson_endpoint_allowed") is True
    )
    rows: list[dict[str, Any]] = []
    for mode_name in P_MODES:
        mode = p_interval["modes"][mode_name]
        banach = mode["banach"]
        enclosure = mode["certified_enclosure"]
        declared_coordinate = str(mode.get("space", "")).startswith(
            "alpha (fine-structure coupling)"
        )
        existence = banach.get("existence") is True
        uniqueness = banach.get("uniqueness_in_interval") is True
        stability = bool(
            banach.get("contraction") is True
            and Decimal(str(banach.get("lipschitz_bound"))) < 1
        )
        gates = {
            "declared_dimensionless_coordinate": gate(
                declared_coordinate,
                [pointer(P_INTERVAL, f"/modes/{mode_name}/space", mode.get("space"))],
            ),
            "target_independent_candidate_selection": gate(
                selected,
                [
                    pointer(EXPOSURE, "/candidates/x_p_certificate_alpha/eligibility", eligibility),
                    pointer(EXPOSURE, "/candidates/x_p_certificate_alpha/freeze_evidence", freeze),
                ],
            ),
            "same_quantity_constructed": gate(
                same_quantity,
                [pointer(EXPOSURE, "/candidates/x_p_certificate_alpha/same_quantity_constructed", inventory.get("same_quantity_constructed"))],
            ),
            "source_return_map_complete": gate(
                source_complete,
                [
                    pointer(P_INTERVAL, "/promotion_allowed", p_interval.get("promotion_allowed")),
                    pointer(P_TRUNK, "/consumer_policy/may_feed_live_particle_predictions", consumer.get("may_feed_live_particle_predictions")),
                    pointer(P_TRUNK, "/consumer_policy/default_thomson_endpoint_allowed", consumer.get("default_thomson_endpoint_allowed")),
                ],
            ),
            "existence_certified": gate(
                existence,
                [pointer(P_INTERVAL, f"/modes/{mode_name}/banach/existence", banach.get("existence"))],
            ),
            "uniqueness_certified": gate(
                uniqueness,
                [pointer(P_INTERVAL, f"/modes/{mode_name}/banach/uniqueness_in_interval", banach.get("uniqueness_in_interval"))],
            ),
            "stability_certified": gate(
                stability,
                [pointer(P_INTERVAL, f"/modes/{mode_name}/banach/lipschitz_bound", banach.get("lipschitz_bound"))],
            ),
        }
        rows.append(
            _row(
                candidate_id=f"p.{mode_name}",
                registry_source=P_INTERVAL.relative_to(ROOT).as_posix(),
                row_kind="candidate_map",
                quantity="P",
                family="P_INTERVAL_MAP",
                mathematical_status="INTERVAL_FIXED_POINT_CERTIFIED",
                disposition="CONDITIONAL_ARITHMETIC_FIXED_POINT__SELECTION_IDENTITY_AND_RETURN_OPEN",
                exposure_class=(
                    str(freeze_exposure)
                    if selected
                    else "registered_source_audit__not_frozen_precomparison"
                ),
                solver_target_payload_read=bool(consumer.get("hidden_external_alpha_allowed")),
                exact_fixed_point_progress=existence and uniqueness and stability,
                gates=gates,
                details={
                    "mode": mode_name,
                    "certified_enclosure": {
                        "P": enclosure.get("P"),
                        "alpha_inverse": enclosure.get("alpha_inv"),
                    },
                },
                required=required,
            )
        )
    require(
        Decimal(rows[1]["details"]["certified_enclosure"]["P"]["hi"])
        < Decimal(rows[0]["details"]["certified_enclosure"]["P"]["lo"]),
        "P fixed-point intervals are not disjoint",
    )
    return rows


def _menu_crosswalk(
    menu: dict[str, Any],
    cap_k: dict[str, Any],
    cap_p: dict[str, Any],
    cap_l: dict[str, Any],
    coupled: dict[str, Any],
    required: list[str],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    menu_rows = menu.get("menu_rows", [])
    require(menu.get("menu_row_count") == len(menu_rows) == 17, "capacity menu row count drift")
    family_counts = Counter(row.get("family") for row in menu_rows)
    require(dict(family_counts) == CAPACITY_FAMILY_COUNTS, "capacity menu family count drift")
    require(
        menu.get("campaign_verdict", {}).get("semantics_enumeration")
        == "declared_menu_complete_for_executed_families",
        "capacity menu completeness status drift",
    )
    menu_ids = [str(row.get("row_id")) for row in menu_rows]
    require(len(menu_ids) == len(set(menu_ids)), "duplicate capacity menu row")

    cap_k_by_id = {row["branch"]: row for row in cap_k.get("rows", [])}
    cap_p_by_id = {row["branch"]: row for row in cap_p.get("rows", [])}
    cap_l_groups: dict[str, list[dict[str, Any]]] = {}
    for row in cap_l.get("rows", []):
        prefix = ".".join(str(row["branch"]).split(".")[:2])
        cap_l_groups.setdefault(prefix, []).append(row)
    require(len(cap_k_by_id) == 4, "CAP-K runtime row count drift")
    require(len(cap_p_by_id) == 6, "CAP-P runtime row count drift")
    require(set(cap_l_groups) == {f"capL.R{i}" for i in range(1, 6)}, "CAP-L group drift")
    require(all(len(rows) == 36 for rows in cap_l_groups.values()), "CAP-L sublattice count drift")

    selector = menu["campaign_verdict"]["source_only_fixed_point_selector"]
    no_target_import = menu["campaign_verdict"].get("no_target_import")
    selector_attained = bool(
        selector in {"source_only_precomparison", "sealed_prospective"}
        and no_target_import == "verified"
    )
    rows: list[dict[str, Any]] = []
    for index, menu_row in enumerate(menu_rows):
        row_id = str(menu_row["row_id"])
        family = str(menu_row["family"])
        channel = menu_row["semantics"]["channel"]
        declared_coordinate = bool(channel)
        same_quantity = menu["campaign_verdict"].get("same_quantity_constructed") == "verified"
        source_return = False
        solver_target = False
        exact_progress = False
        existence = False
        uniqueness = False
        stability = False
        details: dict[str, Any]
        exposure_class = "declared_before_comparison__retrospective_comparison_recorded"

        if family == "CAP-K":
            source = cap_k_by_id[row_id]
            stability = Decimal(source["s"]["hi"]) < 1
            status = source["status"]
            disposition = "EXACT_ROW_NO_POSITIVE_FIXED_POINT__BOUNDED_FAMILY_ONLY"
            details = {"menu_row_id": row_id, "map": source["map"], "s": source["s"]}
            runtime_path = CAP_K
            runtime_gate_pointer = f"/rows/{list(cap_k_by_id).index(row_id)}"
        elif family == "CAP-P":
            source = cap_p_by_id[row_id]
            status = source["status"]
            fixed = status == "fixed_point_certified"
            existence = fixed
            uniqueness = fixed and source.get("fixed_point_exact") is not None
            if fixed:
                stability = source.get("contraction_certificate", {}).get("banach_pass") is True
            else:
                stability = Decimal(source["s"]["hi"]) < 1
            exact_progress = bool(existence and uniqueness and stability)
            disposition = (
                "EXACT_NORMALIZATION_FIXED_POINT_N_EQUALS_PI__SOURCE_SELECTION_IDENTITY_AND_RETURN_OPEN"
                if fixed
                else "EXACT_ROW_NO_POSITIVE_FIXED_POINT__BOUNDED_FAMILY_ONLY"
            )
            details = {
                "menu_row_id": row_id,
                "map": source["map"],
                "s": source["s"],
                "fixed_point_exact": source.get("fixed_point_exact"),
                "fixed_point_enclosure": source.get("fixed_point", {}).get("enclosure"),
            }
            runtime_path = CAP_P
            runtime_gate_pointer = f"/rows/{list(cap_p_by_id).index(row_id)}"
        elif family == "CAP-L":
            group = cap_l_groups[row_id]
            counts = Counter(row["status"] for row in group)
            recorded_counts = menu_row.get("recorded_status_counts")
            require(dict(counts) == recorded_counts, f"CAP-L status count drift: {row_id}")
            existence = counts.get("fixed_point_certified", 0) > 0
            uniqueness = False
            stability = not any(
                counts.get(name, 0) > 0
                for name in ("fixed_point_unstable_rejected", "rejected_no_contraction")
            )
            exact_progress = existence
            status = "MIXED_36_ROW_SUBLATTICE"
            disposition = "ARCHIVED_MIXED_SUBLATTICE__NO_SOURCE_SELECTOR"
            details = {
                "menu_row_id": row_id,
                "sublattice_row_count": len(group),
                "status_counts": dict(sorted(counts.items())),
            }
            runtime_path = SEMANTICS_MENU
            runtime_gate_pointer = f"/menu_rows/{index}/recorded_status_counts"
        elif family == "CAP-B":
            require(row_id == "capB.bridge_constant", "CAP-B identity drift")
            status = menu_row["executed_verdict"]
            disposition = "BARRED_TARGET_BRIDGE_PRE_EVALUATION__NEVER_EXECUTED"
            exposure_class = "barred_target_bridge_pre_evaluation"
            details = {"menu_row_id": row_id, "executed_verdict": status}
            runtime_path = SEMANTICS_MENU
            runtime_gate_pointer = f"/menu_rows/{index}/executed_verdict"
        elif family == "coupled":
            require(row_id == "coupled.cp1_cp2_cp3", "coupled identity drift")
            load_certificate = coupled["certificate"]["load_coordinate"]
            contraction = load_certificate["contraction_certificate"]
            existence = load_certificate["fixed_point"]["located"] is True
            uniqueness = load_certificate["exact_solution_check"]["x_ew_inside_certified_box"] is True
            stability = contraction["banach_pass"] is True
            exact_progress = existence and uniqueness and stability
            status = coupled["status"]
            disposition = "CONDITIONAL_FIXED_POINT__CP1_CP2_CP3_OPEN__TARGET_EXPOSED"
            exposure_class = "target_exposed_theorem_coupled"
            solver_target = True
            details = {
                "menu_row_id": row_id,
                "lambda": load_certificate["lambda"],
                "fixed_point_enclosure": load_certificate["fixed_point"]["enclosure"],
                "cl7_status": coupled["cl7_status"],
                "evaluation_cone_contains_bridge_expression": coupled["blindness"]["cone_contains_cl3_bridge_expression"],
            }
            runtime_path = COUPLED
            runtime_gate_pointer = "/certificate/load_coordinate"
        else:  # pragma: no cover - guarded by family count equality
            raise ValueError(f"unknown capacity menu family: {family}")

        gates = {
            "declared_dimensionless_coordinate": gate(
                declared_coordinate,
                [pointer(SEMANTICS_MENU, f"/menu_rows/{index}/semantics/channel", channel)],
            ),
            "target_independent_candidate_selection": gate(
                selector_attained,
                [
                    pointer(SEMANTICS_MENU, "/campaign_verdict/source_only_fixed_point_selector", selector),
                    pointer(SEMANTICS_MENU, "/campaign_verdict/no_target_import", no_target_import),
                    pointer(SEMANTICS_MENU, f"/menu_rows/{index}/executed_verdict", menu_row["executed_verdict"]),
                ],
            ),
            "same_quantity_constructed": gate(
                same_quantity,
                [pointer(SEMANTICS_MENU, "/campaign_verdict/same_quantity_constructed", menu["campaign_verdict"].get("same_quantity_constructed"))],
            ),
            "source_return_map_complete": gate(
                source_return,
                [
                    pointer(
                        SEMANTICS_MENU,
                        "/campaign_verdict/horizon_area_and_ew_bridge",
                        menu["campaign_verdict"].get("horizon_area_and_ew_bridge"),
                    )
                ],
            ),
            "existence_certified": gate(
                existence,
                [pointer(runtime_path, runtime_gate_pointer, _resolve_pointer(load(runtime_path), runtime_gate_pointer))],
            ),
            "uniqueness_certified": gate(
                uniqueness,
                [pointer(runtime_path, runtime_gate_pointer, _resolve_pointer(load(runtime_path), runtime_gate_pointer))],
            ),
            "stability_certified": gate(
                stability,
                [pointer(runtime_path, runtime_gate_pointer, _resolve_pointer(load(runtime_path), runtime_gate_pointer))],
            ),
        }
        row_kind = "candidate_family" if family == "CAP-L" else "candidate_map"
        if family == "CAP-B":
            row_kind = "registered_control"
        rows.append(
            _row(
                candidate_id=f"n.menu.{row_id}",
                registry_source=SEMANTICS_MENU.relative_to(ROOT).as_posix(),
                row_kind=row_kind,
                quantity="N_CAPACITY_MENU",
                family=family,
                mathematical_status=status,
                disposition=disposition,
                exposure_class=exposure_class,
                solver_target_payload_read=solver_target,
                exact_fixed_point_progress=exact_progress,
                gates=gates,
                details=details,
                required=required,
            )
        )

    crosswalk = {
        "authoritative_menu_path": SEMANTICS_MENU.relative_to(ROOT).as_posix(),
        "authoritative_menu_row_count": len(menu_rows),
        "authoritative_menu_row_ids": menu_ids,
        "family_counts": dict(family_counts),
        "runtime_counts": {"CAP-K": 4, "CAP-P": 6, "CAP-L": 180, "CAP-B": 0, "coupled": 1},
        "inventory_order_derived_from_authoritative_sources": True,
        "contract_candidate_order_used": False,
    }
    return rows, crosswalk


def _resolve_pointer(document: Any, json_pointer: str) -> Any:
    current = document
    for raw in json_pointer.lstrip("/").split("/"):
        token = raw.replace("~1", "/").replace("~0", "~")
        if isinstance(current, list):
            current = current[int(token)]
        elif isinstance(current, dict):
            if token not in current:
                return None
            current = current[token]
        else:
            return None
    return current


def _direct_n_row(direct: dict[str, Any], required: list[str]) -> dict[str, Any]:
    fixed = direct["fixed_cutoff_result"]
    generated = direct["bounded_generation_register_result"]
    controls = direct["source_controls"]
    declared = (
        direct["typed_coordinates"]["input_coordinate"] == "N_in=log(D)"
        and direct["typed_coordinates"]["output_coordinate"] == "N_out=log(M0)"
    )
    selected = bool(
        controls.get("target_clean") is True
        and generated.get("unique_source_zero_entailed") is True
        and direct.get("freeze_evidence", {}).get("frozen_before_comparison") is True
    )
    source_complete = bool(
        generated.get("universal_all_rung_membership_proved") is True
        and generated.get("executable_lean_membership_bridge_proved") is True
    )
    gates = {
        "declared_dimensionless_coordinate": gate(declared, [pointer(DIRECT_N, "/typed_coordinates", direct["typed_coordinates"])]),
        "target_independent_candidate_selection": gate(
            selected,
            [
                pointer(DIRECT_N, "/source_controls", controls),
                pointer(
                    DIRECT_N,
                    "/source_controls/constructor_reads_desired_capacity",
                    controls.get("constructor_reads_desired_capacity"),
                ),
                pointer(
                    DIRECT_N,
                    "/bounded_generation_register_result/unique_source_zero_entailed",
                    generated.get("unique_source_zero_entailed"),
                ),
                pointer(DIRECT_N, "/freeze_evidence", direct.get("freeze_evidence")),
            ],
        ),
        "same_quantity_constructed": gate(direct.get("same_quantity_constructed") is True, [pointer(DIRECT_N, "/same_quantity_constructed", direct.get("same_quantity_constructed"))]),
        "source_return_map_complete": gate(source_complete, [pointer(DIRECT_N, "/bounded_generation_register_result/universal_all_rung_membership_proved", generated.get("universal_all_rung_membership_proved")), pointer(DIRECT_N, "/bounded_generation_register_result/executable_lean_membership_bridge_proved", generated.get("executable_lean_membership_bridge_proved"))]),
        "existence_certified": gate(False, [pointer(DIRECT_N, "/bounded_generation_register_result/unique_source_zero_entailed", generated.get("unique_source_zero_entailed"))]),
        "uniqueness_certified": gate(False, [pointer(DIRECT_N, "/bounded_generation_register_result/unique_source_zero_entailed", generated.get("unique_source_zero_entailed"))]),
        "stability_certified": gate(False, [pointer(DIRECT_N, "/stability_certificate", direct.get("stability_certificate"))]),
    }
    return _row(
        candidate_id="n.direct_correctable_record",
        registry_source=DIRECT_N.relative_to(ROOT).as_posix(),
        row_kind="fixed_cutoff_control",
        quantity="N_DIRECT",
        family="DIRECT_CORRECTABLE_RECORD",
        mathematical_status=direct["status"],
        disposition="FIXED_CUTOFF_CONTROL__COSMIC_SOURCE_ANTECEDENT_INCOMPLETE",
        exposure_class="target_clean_bounded_control__not_a_candidate_map",
        solver_target_payload_read=not bool(controls.get("target_clean")),
        exact_fixed_point_progress=False,
        gates=gates,
        details={
            "fixed_cutoff_control": True,
            "D": fixed.get("D"),
            "M0": fixed.get("M0"),
            "cosmic_value_selected": fixed.get("cosmic_value_selected"),
        },
        required=required,
    )


def _common_load_rows(branch: dict[str, Any], named: dict[str, Any], required: list[str]) -> list[dict[str, Any]]:
    scope = branch["scope"]
    identity = branch["self_reference_boundary"]
    baseline = branch["conditional_baseline"]
    source_inputs = branch["source_inputs"]
    common_source = bool(
        scope.get("global_capacity_derived") is True
        and source_inputs.get("source_artifact_status") == "source_complete_endpoint_proof"
    )
    common_same = identity.get("same_typed_quantity_identified") is True
    common_selected = bool(
        scope.get("target_blind_forecast") is True
        and scope.get("retrospective") is False
        and scope.get("branch_selected") is True
    )
    common_unique = identity.get("unique_source_fixed_point_proved") is True
    common_stable = named.get("stability_certified") is True
    require(baseline.get("formula") == COMMON_LOAD_FORMULA, "common-load formula drift")

    def gates_for(formula_pointer: str, formula: Any, selected: bool) -> dict[str, dict[str, Any]]:
        return {
            "declared_dimensionless_coordinate": gate(bool(formula), [pointer(N_BRANCH, formula_pointer, formula)]),
            "target_independent_candidate_selection": gate(bool(common_selected and selected), [pointer(N_BRANCH, "/scope", scope)]),
            "same_quantity_constructed": gate(common_same, [pointer(N_BRANCH, "/self_reference_boundary/same_typed_quantity_identified", identity.get("same_typed_quantity_identified"))]),
            "source_return_map_complete": gate(common_source, [pointer(N_BRANCH, "/scope/global_capacity_derived", scope.get("global_capacity_derived")), pointer(N_BRANCH, "/source_inputs/source_artifact_status", source_inputs.get("source_artifact_status"))]),
            "existence_certified": gate(False, [pointer(N_BRANCH, "/status", branch.get("status"))]),
            "uniqueness_certified": gate(common_unique, [pointer(N_BRANCH, "/self_reference_boundary/unique_source_fixed_point_proved", identity.get("unique_source_fixed_point_proved"))]),
            "stability_certified": gate(common_stable, [pointer(NAMED_N, "/stability_certified", named.get("stability_certified"))]),
        }

    rows = [
        _row(
            candidate_id="n.common_load_baseline",
            registry_source=N_BRANCH.relative_to(ROOT).as_posix(),
            row_kind="candidate_map",
            quantity="N_COMMON_LOAD",
            family="COMMON_LOAD",
            mathematical_status="EXACT_CONDITIONAL_BASELINE_ARITHMETIC",
            disposition="RETROSPECTIVE_CONDITIONAL_ARITHMETIC__NO_SOURCE_SELECTION_OR_RETURN",
            exposure_class="exposed_retrospective",
            solver_target_payload_read=scope.get("comparison_data_consumed") is True,
            exact_fixed_point_progress=False,
            gates=gates_for("/conditional_baseline/formula", baseline.get("formula"), scope.get("branch_selected") is True),
            details={"formula": baseline.get("formula"), "evaluated_N": baseline.get("N0")},
            required=required,
        )
    ]
    branch_rows = branch.get("branches", [])
    require([row.get("branch_id") for row in branch_rows] == ["finite_presence", "poisson_projective_limit"], "reserve branch grammar drift")
    for index, source in enumerate(branch_rows):
        branch_id = source["branch_id"]
        require(
            source.get("formula") == RESERVE_FORMULAS[branch_id],
            f"reserve formula drift: {branch_id}",
        )
        rows.append(
            _row(
                candidate_id=f"n.reserve_{branch_id}",
                registry_source=N_BRANCH.relative_to(ROOT).as_posix(),
                row_kind="candidate_map",
                quantity="N_COMMON_LOAD",
                family="RESERVE_BRANCH",
                mathematical_status=source["status"],
                disposition="RETROSPECTIVE_UNSELECTED_CONDITIONAL_ARITHMETIC",
                exposure_class="exposed_retrospective",
                solver_target_payload_read=scope.get("comparison_data_consumed") is True,
                exact_fixed_point_progress=False,
                gates=gates_for(f"/branches/{index}/formula", source.get("formula"), source.get("selected") is True),
                details={
                    "branch_id": branch_id,
                    "formula": source.get("formula"),
                    "evaluated_N": source.get("N"),
                    "selected": source.get("selected"),
                    "mean_count_or_projective_limit_carrier_required": source.get("mean_count_or_projective_limit_carrier_required"),
                },
                required=required,
            )
        )
    return rows


def _rc_load_row(rc: dict[str, Any], required: list[str]) -> dict[str, Any]:
    validation = rc["validation"]
    scope = rc["scope"]
    exposure = rc["exposure"]
    existence = validation["residual_encloses_zero"] is True
    uniqueness = validation["unique_fixed_point_in_domain"] is True
    stability = validation["contraction"] is True
    require(rc.get("equation") == RC_EQUATION, "RC-LOAD equation drift")
    selected = bool(
        exposure.get("class") in QUALIFYING_EXPOSURES
        and exposure.get("promotion_from_landing") is False
        and rc.get("freeze_evidence", {}).get("frozen_before_comparison") is True
    )
    gates = {
        "declared_dimensionless_coordinate": gate(rc["equation"].startswith("ln(N/pi)"), [pointer(RC_OUTPUT, "/equation", rc["equation"])]),
        "target_independent_candidate_selection": gate(selected, [pointer(RC_OUTPUT, "/exposure", exposure), pointer(RC_OUTPUT, "/freeze_evidence", rc.get("freeze_evidence"))]),
        "same_quantity_constructed": gate(scope["same_quantity_identity_proved_here"] is True, [pointer(RC_OUTPUT, "/scope/same_quantity_identity_proved_here", scope["same_quantity_identity_proved_here"])]),
        "source_return_map_complete": gate(scope["pricing_law_proved_here"] is True, [pointer(RC_OUTPUT, "/scope/pricing_law_proved_here", scope["pricing_law_proved_here"])]),
        "existence_certified": gate(existence, [pointer(RC_OUTPUT, "/validation/residual_encloses_zero", validation["residual_encloses_zero"])]),
        "uniqueness_certified": gate(uniqueness, [pointer(RC_OUTPUT, "/validation/unique_fixed_point_in_domain", validation["unique_fixed_point_in_domain"])]),
        "stability_certified": gate(stability, [pointer(RC_OUTPUT, "/validation/contraction", validation["contraction"])]),
    }
    return _row(
        candidate_id="n.rc_load",
        registry_source=RC_OUTPUT.relative_to(ROOT).as_posix(),
        row_kind="candidate_map",
        quantity="N_RC_LOAD",
        family="RC_LOAD",
        mathematical_status=rc["status"],
        disposition="CONDITIONAL_FIXED_POINT__EXPOSED_CAMPAIGN__PRICING_AND_IDENTITY_OPEN",
        exposure_class="exposed_retrospective_upstream_campaign__target_free_replay",
        solver_target_payload_read=scope["observational_comparison_imported"] is True,
        exact_fixed_point_progress=existence and uniqueness and stability,
        gates=gates,
        details={
            "equation": rc["equation"],
            "budget_x": rc["budget_x"],
            "fixed_point_X": rc["fixed_point"]["X"],
        },
        required=required,
    )


def _archive_row(
    archive: dict[str, Any], menu_crosswalk: dict[str, Any], required: list[str]
) -> dict[str, Any]:
    archive_rows = archive.get("rows", [])
    branch_ids = {str(row["branch"]) for row in archive_rows}
    expected_archive_ids = set(load(CAP_K).get("rows", [])[i]["branch"] for i in range(4))
    expected_archive_ids |= {row["branch"] for row in load(CAP_P).get("rows", [])}
    expected_archive_ids |= {row["branch"] for row in load(CAP_L).get("rows", [])}
    require(branch_ids == expected_archive_ids, "190-row archive/runtime crosswalk drift")
    require(archive.get("total_rows") == len(archive_rows) == 190, "archive row count drift")
    counts = Counter(row["status"] for row in archive_rows)
    certified = archive.get("certified_rows")
    gates = {
        "declared_dimensionless_coordinate": gate(bool(archive.get("certified_fixed_point_range_nats")), [pointer(ARCHIVE_190, "/certified_fixed_point_range_nats", archive.get("certified_fixed_point_range_nats"))]),
        "target_independent_candidate_selection": gate(False, [pointer(ARCHIVE_190, "/landing_criterion", archive.get("landing_criterion"))]),
        "same_quantity_constructed": gate(False, [pointer(ARCHIVE_190, "/verdict", archive.get("verdict"))]),
        "source_return_map_complete": gate(
            False,
            [
                pointer(
                    SEMANTICS_MENU,
                    "/campaign_verdict/semantics_enumeration",
                    "declared_menu_complete_for_executed_families",
                ),
                pointer(
                    CONTRACT,
                    "/decision_rule/archive_enumeration_is_not_a_source_return_map",
                    True,
                ),
            ],
        ),
        "existence_certified": gate(isinstance(certified, int) and certified > 0, [pointer(ARCHIVE_190, "/certified_rows", certified)]),
        "uniqueness_certified": gate(False, [pointer(ARCHIVE_190, "/certified_rows", certified)]),
        "stability_certified": gate(
            False,
            [
                pointer(ARCHIVE_190, "/total_rows", archive.get("total_rows")),
                pointer(ARCHIVE_190, "/certified_rows", certified),
            ],
        ),
    }
    return _row(
        candidate_id="n.capacity_map_archive_190",
        registry_source=ARCHIVE_190.relative_to(ROOT).as_posix(),
        row_kind="archive_aggregate",
        quantity="N_ARCHIVED_MAP_LATTICE",
        family="CAP_L_P_K_ARCHIVE",
        mathematical_status="DECLARED_EXECUTED_FAMILY_ARCHIVE",
        disposition="AGGREGATE_METADATA_ONLY__ENUMERATION_IS_NOT_RETURN_MAP",
        exposure_class="declared_before_comparison__retrospective_comparison_archive",
        solver_target_payload_read=True,
        exact_fixed_point_progress=isinstance(certified, int) and certified > 0,
        gates=gates,
        details={
            "total_rows": len(archive_rows),
            "certified_rows": certified,
            "landed_rows": archive.get("landed_rows"),
            "status_counts": dict(sorted(counts.items())),
            "covered_runtime_families": ["CAP-K", "CAP-P", "CAP-L"],
            "omitted_registered_menu_rows": ["capB.bridge_constant", "coupled.cp1_cp2_cp3"],
            "source_return_map_constructed": False,
            "authoritative_menu_row_count": menu_crosswalk["authoritative_menu_row_count"],
        },
        required=required,
    )


def _hierarchy_rows(required: list[str]) -> list[dict[str, Any]]:
    p_public = load(HIER_P_PUBLIC)
    p_source = load(HIER_P_SOURCE)
    joint = load(HIER_PN_JOINT)
    tick = load(HIER_N_TICK)
    ew = load(HIER_EW_CAPACITY)

    def empty_gates(path: Path, status_pointer: str, status: Any, *, declared: bool = True) -> dict[str, dict[str, Any]]:
        return {
            "declared_dimensionless_coordinate": gate(declared, [pointer(path, status_pointer, status)]),
            "target_independent_candidate_selection": gate(False, [pointer(path, status_pointer, status)]),
            "same_quantity_constructed": gate(False, [pointer(path, status_pointer, status)]),
            "source_return_map_complete": gate(False, [pointer(path, status_pointer, status)]),
            "existence_certified": gate(False, [pointer(path, status_pointer, status)]),
            "uniqueness_certified": gate(False, [pointer(path, status_pointer, status)]),
            "stability_certified": gate(False, [pointer(path, status_pointer, status)]),
        }

    rows = [
        _row(
            candidate_id="hierarchy.p_public_endpoint",
            registry_source=HIER_P_PUBLIC.relative_to(ROOT).as_posix(),
            row_kind="packet_disposition",
            quantity="P",
            family="HIERARCHY_PACKET",
            mathematical_status=p_public["status"],
            disposition="TARGET_EXPOSED_PUBLIC_ENDPOINT_CONTROL__INELIGIBLE_SOURCE_ONLY",
            exposure_class="target_exposed_public_endpoint_locator",
            solver_target_payload_read=True,
            exact_fixed_point_progress=False,
            gates=empty_gates(HIER_P_PUBLIC, "/status", p_public["status"]),
            details={"artifact": p_public["artifact"], "public_endpoint_convention": p_public["public_endpoint_convention"], "not_supplied_here": p_public["proof_rule"]["not_supplied_here"]},
            required=required,
        ),
        _row(
            candidate_id="hierarchy.p_source_audit_witness",
            registry_source=HIER_P_SOURCE.relative_to(ROOT).as_posix(),
            row_kind="packet_disposition",
            quantity="P",
            family="HIERARCHY_PACKET",
            mathematical_status=p_source["status"],
            disposition="SOURCE_AUDIT_WITNESS__NOT_DISTINCT_CANDIDATE__ENDPOINT_PROOF_OPEN",
            exposure_class="source_audit_witness__not_frozen_candidate",
            solver_target_payload_read=False,
            exact_fixed_point_progress=False,
            gates=empty_gates(HIER_P_SOURCE, "/status", p_source["status"]),
            details={"artifact": p_source["artifact"], "P_cand": p_source["P_cand"], "needed_for_public_endpoint_completion": p_source["needed_for_public_endpoint_completion"]},
            required=required,
        ),
        _row(
            candidate_id="hierarchy.pn_product_contraction_adapter",
            registry_source=HIER_PN_JOINT.relative_to(ROOT).as_posix(),
            row_kind="packet_disposition",
            quantity="P_N_PRODUCT",
            family="HIERARCHY_PACKET",
            mathematical_status=joint["status"],
            disposition="CONDITIONAL_PRODUCT_ADAPTER__SOURCE_N_COMPONENT_ABSENT__TARGET_EXPOSED_DIAGNOSTIC",
            exposure_class="target_exposed_numeric_demo_and_circular_backsolve",
            solver_target_payload_read=True,
            exact_fixed_point_progress=False,
            gates=empty_gates(HIER_PN_JOINT, "/status", joint["status"]),
            details={"issue": joint["issue"], "product_condition": joint["product_contraction_certificate"], "coupled_boundary": joint["coupled_contraction_certificate"], "N_display_warning": joint["N_display_warning"], "N_backsolved_warning": joint["N_backsolved_warning"]},
            required=required,
        ),
        _row(
            candidate_id="hierarchy.n_global_repair_tick_adapter",
            registry_source=HIER_N_TICK.relative_to(ROOT).as_posix(),
            row_kind="packet_disposition",
            quantity="N_DEPENDENT_TICK",
            family="HIERARCHY_PACKET",
            mathematical_status=tick["status"],
            disposition="CONDITIONAL_TICK_IDENTITY__NOT_AN_N_PRODUCER__TARGET_EXPOSED_DISPLAY",
            exposure_class="downstream_conditional__contains_exposed_rounded_display",
            solver_target_payload_read=True,
            exact_fixed_point_progress=False,
            gates=empty_gates(HIER_N_TICK, "/status", tick["status"]),
            details={"artifact": tick["artifact"], "normalization": tick["normalization"], "declared_not_derived": tick["claim_boundary"]["declared_not_derived"], "numeric_display": tick["numeric_display_for_rounded_capacity"]},
            required=required,
        ),
    ]

    ew_contraction = ew["contraction_certificate"]
    ew_gates = {
        "declared_dimensionless_coordinate": gate(True, [pointer(HIER_EW_CAPACITY, "/definitions/exact_log_capacity", ew["definitions"]["exact_log_capacity"])]),
        "target_independent_candidate_selection": gate(False, [pointer(HIER_EW_CAPACITY, "/source_values/P_star_branch_locator", ew["source_values"]["P_star_branch_locator"])]),
        "same_quantity_constructed": gate(False, [pointer(HIER_EW_CAPACITY, "/branch_selection", ew["branch_selection"])]),
        "source_return_map_complete": gate(False, [pointer(HIER_EW_CAPACITY, "/claim_boundary/scope", ew["claim_boundary"]["scope"])]),
        "existence_certified": gate(ew["accepted"] is True, [pointer(HIER_EW_CAPACITY, "/accepted", ew["accepted"])]),
        "uniqueness_certified": gate(ew_contraction["banach_unique_fixed_point"] is True, [pointer(HIER_EW_CAPACITY, "/contraction_certificate/banach_unique_fixed_point", ew_contraction["banach_unique_fixed_point"])]),
        "stability_certified": gate(Decimal(ew_contraction["lipschitz_constant"]) < 1, [pointer(HIER_EW_CAPACITY, "/contraction_certificate/lipschitz_constant", ew_contraction["lipschitz_constant"])]),
    }
    rows.append(
        _row(
            candidate_id="hierarchy.n_ew_bridge_defined_capacity",
            registry_source=HIER_EW_CAPACITY.relative_to(ROOT).as_posix(),
            row_kind="candidate_map",
            quantity="N_EW_BRIDGE_DEFINED",
            family="HIERARCHY_PACKET",
            mathematical_status=ew["status"],
            disposition="TARGET_EXPOSED_BRIDGE_DEFINED_FIXED_POINT__NOT_INDEPENDENT_CAPACITY",
            exposure_class="target_exposed_public_endpoint_branch_locator",
            solver_target_payload_read=True,
            exact_fixed_point_progress=True,
            gates=ew_gates,
            details={"certificate_id": ew["certificate_id"], "target_relation": ew["target_relation"], "lambda": ew["source_values"]["lambda"], "fixed_point": ew["exact_capacity_fixed_point"]["N_CRC_EW"], "rounded_capacity_diagnostic": ew["rounded_capacity_diagnostic"]},
            required=required,
        )
    )
    return rows


def _completion_gates(
    contract: dict[str, Any], forecast_state: dict[str, Any], menu: dict[str, Any]
) -> dict[str, dict[str, Any]]:
    declared = contract["preflight_completion_state"]

    def declared_gate(name: str, extra: list[dict[str, Any]] | None = None) -> dict[str, Any]:
        state = declared[name]
        evidence = [
            pointer(CONTRACT, f"/preflight_completion_state/{name}/attained", state["attained"]),
            pointer(CONTRACT, f"/preflight_completion_state/{name}/reason", state["reason"]),
        ]
        if extra:
            evidence.extend(extra)
        return gate(state["attained"], evidence)

    return {
        "input_catalogue_frozen": declared_gate(
            "input_catalogue_frozen",
            [
                pointer(
                    FORECAST_STATE,
                    "/inventory_completeness/status",
                    forecast_state["inventory_completeness"]["status"],
                ),
                pointer(
                    FORECAST_STATE,
                    "/controls/candidate_inventory_exhaustiveness",
                    forecast_state["controls"]["candidate_inventory_exhaustiveness"],
                ),
            ],
        ),
        "nuisance_directions_frozen": declared_gate("nuisance_directions_frozen"),
        "units_and_schemes_frozen": declared_gate("units_and_schemes_frozen"),
        "candidate_grammar_frozen": declared_gate(
            "candidate_grammar_frozen",
            [
                pointer(
                    FORECAST_STATE,
                    "/selection_rule/status",
                    forecast_state["selection_rule"]["status"],
                ),
                pointer(
                    SEMANTICS_MENU,
                    "/campaign_verdict/semantics_enumeration",
                    menu["campaign_verdict"]["semantics_enumeration"],
                ),
                pointer(
                    CONTRACT,
                    "/scope_boundary/future_source_laws_exhausted",
                    contract["scope_boundary"]["future_source_laws_exhausted"],
                ),
            ],
        ),
        "decision_rule_frozen": declared_gate(
            "decision_rule_frozen",
            [pointer(CONTRACT, "/decision_rule/frozen", contract["decision_rule"]["frozen"])],
        ),
    }


def build(*, validate_parents: bool = True) -> dict[str, Any]:
    contract = _load_contract()
    required = contract["source_only_required_gates"]
    completion_names = contract["preflight_completion_gates"]
    parent_validation = validate_parent_artifacts() if validate_parents else []

    p_interval = load(P_INTERVAL)
    p_trunk = load(P_TRUNK)
    cap_k = load(CAP_K)
    cap_p = load(CAP_P)
    cap_l = load(CAP_L)
    coupled = load(COUPLED)
    direct = load(DIRECT_N)
    branch = load(N_BRANCH)
    named = load(NAMED_N)
    menu = load(SEMANTICS_MENU)
    archive = load(ARCHIVE_190)
    exposure = load(EXPOSURE)
    forecast_state = load(FORECAST_STATE)
    rc = rc_load_source_only.build()
    rc_bytes = canonical_json_bytes(rc)

    rows = _p_rows(p_interval, p_trunk, exposure, required)
    menu_rows, menu_crosswalk = _menu_crosswalk(menu, cap_k, cap_p, cap_l, coupled, required)
    rows.extend(menu_rows)
    rows.append(_direct_n_row(direct, required))
    rows.extend(_common_load_rows(branch, named, required))
    rows.append(_rc_load_row(rc, required))
    rows.append(_archive_row(archive, menu_crosswalk, required))
    rows.extend(_hierarchy_rows(required))
    ids = [row["candidate_id"] for row in rows]
    require(len(ids) == len(set(ids)) == 30, "derived H0 inventory count or identity drift")

    completion = _completion_gates(contract, forecast_state, menu)
    require(list(completion) == completion_names, "completion gate order drift")
    completion_complete = all(row["attained"] is True for row in completion.values())
    qualifying = _finalize_qualification(rows, required, completion_complete)
    status = derive_progress_status(rows, required, completion)
    strongest_count = max(row["attained_gate_count"] for row in rows)

    hierarchy_ids = [row["candidate_id"] for row in rows if row["family"] == "HIERARCHY_PACKET"]
    value: dict[str, Any] = {
        "schema": SCHEMA,
        "artifact": "oph_source_only_dimensionless_closure_preflight",
        "issue": ISSUE,
        "status": status,
        "scope": {
            "registered_packets_only": True,
            "future_source_laws_exhausted": False,
            "issue_closure_claimed": False,
            "physical_or_laboratory_attachment_used_as_gate": False,
            "new_observational_comparison_opened": False,
            "prediction_or_postdiction_promoted": False,
        },
        "inventory_crosswalk": {
            **menu_crosswalk,
            "archive": {
                "executed_row_count": archive["total_rows"],
                "covered_families": ["CAP-K", "CAP-P", "CAP-L"],
                "does_not_cover_menu_rows": ["capB.bridge_constant", "coupled.cp1_cp2_cp3"],
                "is_not_an_independent_candidate": True,
                "is_not_a_source_return_map": True,
            },
            "hierarchy_packet_count": len(hierarchy_ids),
            "hierarchy_packet_ids": hierarchy_ids,
            "derived_inventory_row_count": len(rows),
        },
        "preflight_completion": {
            "required_gates": completion_names,
            "gates": completion,
            "all_complete": completion_complete,
            "candidate_present_forbidden_while_incomplete": True,
        },
        "qualification": {
            "required_source_only_gates": required,
            "qualifying_exposure_classes": sorted(QUALIFYING_EXPOSURES),
            "explicitly_downstream_gates": contract["explicitly_downstream_gates"],
            "qualifying_candidate_ids": qualifying,
            "qualifying_candidate_count": len(qualifying),
        },
        "progress": {
            "registered_row_count": len(rows),
            "authoritative_capacity_menu_row_count": len(menu_rows),
            "rows_with_exact_fixed_point_progress": sum(row["exact_fixed_point_progress"] is True for row in rows),
            "strongest_attained_gate_count": strongest_count,
            "strongest_progress_candidate_ids": [row["candidate_id"] for row in rows if row["attained_gate_count"] == strongest_count],
            "status_is_evidence_derived": True,
            "direct_D_equals_M0_control_counted_as_fixed_point_progress": False,
        },
        "candidate_rows": rows,
        "parent_validation": parent_validation,
        "parent_pins": [
            pin(CONTRACT),
            pin(RC_PROVENANCE),
            pin(PRODUCER_SOURCE),
            pin(INDEPENDENT_VERIFIER),
            pin(RC_SOURCE),
            pin(RC_OUTPUT, rc_bytes),
            pin(P_INTERVAL),
            pin(P_TRUNK),
            pin(CAP_K),
            pin(CAP_P),
            pin(CAP_L),
            pin(COUPLED),
            pin(DIRECT_N),
            pin(N_BRANCH),
            pin(NAMED_N),
            pin(SEMANTICS_MENU),
            pin(ARCHIVE_190),
            pin(EXPOSURE),
            pin(FORECAST_STATE),
            pin(HIER_P_PUBLIC),
            pin(HIER_P_SOURCE),
            pin(HIER_PN_JOINT),
            pin(HIER_N_TICK),
            pin(HIER_EW_CAPACITY),
        ],
    }
    value["certificate_sha256"] = tagged(canonical_json_bytes(value))
    return value


def write_outputs() -> None:
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    rc_value = rc_load_source_only.build()
    RC_OUTPUT.write_bytes(canonical_json_bytes(rc_value))
    OUTPUT.write_bytes(canonical_json_bytes(build(validate_parents=True)))


def verify(path: Path = OUTPUT) -> None:
    expected_rc = rc_load_source_only.build()
    require(load(RC_OUTPUT) == expected_rc, "RC-LOAD source-only replay is stale")
    require(RC_OUTPUT.read_bytes() == canonical_json_bytes(expected_rc), "RC-LOAD replay is not canonical")
    expected = build(validate_parents=True)
    actual = load(path)
    require(actual == expected, "source-only closure preflight is stale")
    require(path.read_bytes() == canonical_json_bytes(expected), "source-only closure preflight is not canonical")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    try:
        if args.write:
            write_outputs()
        if args.verify:
            verify()
            print("SOURCE_ONLY_CLOSURE_PREFLIGHT_VALID")
        if not args.write and not args.verify:
            print(canonical_json_bytes(build(validate_parents=True)).decode("ascii"), end="")
    except (KeyError, TypeError, ValueError, OSError, json.JSONDecodeError) as error:
        print(f"SOURCE_ONLY_CLOSURE_PREFLIGHT_INVALID: {error}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
