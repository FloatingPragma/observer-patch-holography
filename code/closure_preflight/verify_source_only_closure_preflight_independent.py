#!/usr/bin/env python3
"""Independent verifier for the H0 source-only closure preflight.

The verifier imports neither H0 producer.  It reconstructs the thirty-row
inventory from the registered P packets, the authoritative seventeen-row
capacity menu, the direct/common/RC controls, and the five hierarchy packets.
It also recomputes the RC interval theorem and the common-load arithmetic.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from decimal import Decimal
from pathlib import Path, PurePosixPath
from typing import Any, Mapping

from mpmath import iv, mp, mpf


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent

CERTIFICATE = HERE / "outputs/source_only_closure_preflight.json"
CONTRACT = HERE / "data/source_only_closure_contract_v3.json"
RC_PROVENANCE = HERE / "data/rc_load_workspace_provenance_v1.json"
RC_SOURCE = HERE / "rc_load_source_only.py"
RC_OUTPUT = HERE / "outputs/rc_load_source_only_certificate.json"
PRODUCER = HERE / "source_only_closure_preflight.py"
VERIFIER = HERE / "verify_source_only_closure_preflight_independent.py"

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

HIERARCHY = ROOT / "code/particles/hierarchy/certificates"
HIER_P_PUBLIC = HIERARCHY / "R_P_public_pixel_certificate.json"
HIER_P_SOURCE = HIERARCHY / "R_P_source_audit_pixel_certificate.json"
HIER_PN_JOINT = HIERARCHY / "R_PN_joint_fixed_point_certificate_report.json"
HIER_N_TICK = HIERARCHY / "R_N_global_repair_tick_certificate.json"
HIER_EW_CAPACITY = HIERARCHY / "R_EW_global_capacity_certificate.json"

SCHEMA = "oph.source_only_closure_preflight.v3"
CONTRACT_SCHEMA = "oph.source_only_closure_preflight_contract.v3"
ISSUE = 708
P_MODES = (
    "thomson_structured_running",
    "thomson_structured_running_plus_gauge_width",
)
REQUIRED_GATES = (
    "declared_dimensionless_coordinate",
    "target_independent_candidate_selection",
    "same_quantity_constructed",
    "source_return_map_complete",
    "existence_certified",
    "uniqueness_certified",
    "stability_certified",
)
COMPLETION_GATES = (
    "input_catalogue_frozen",
    "nuisance_directions_frozen",
    "units_and_schemes_frozen",
    "candidate_grammar_frozen",
    "decision_rule_frozen",
)
DOWNSTREAM_GATES = (
    "laboratory_attachment",
    "cosmological_attachment",
    "absolute_scale_binding",
)
QUALIFYING_ROW_KINDS = {"candidate_map", "candidate_family"}
QUALIFYING_EXPOSURES = {"source_only_precomparison", "sealed_prospective"}
P_ELIGIBILITY_ALLOWLIST = {"ELIGIBLE_SOURCE_VISIBLE"}

CAPACITY_IDS = (
    "capK.s_poisson",
    "capK.s_presence",
    "capK.s_nat_share",
    "capK.s_edge_share",
    "capP.s_poisson_port",
    "capP.s_presence_port",
    "capP.s_poisson_pair",
    "capP.s_presence_pair",
    "capP.add_slot",
    "capP.add_port",
    "capL.R1",
    "capL.R2",
    "capL.R3",
    "capL.R4",
    "capL.R5",
    "capB.bridge_constant",
    "coupled.cp1_cp2_cp3",
)
CAPACITY_FAMILIES = {
    "CAP-K": 4,
    "CAP-P": 6,
    "CAP-L": 5,
    "CAP-B": 1,
    "coupled": 1,
}
ROW_IDS = (
    *(f"p.{mode}" for mode in P_MODES),
    *(f"n.menu.{row_id}" for row_id in CAPACITY_IDS),
    "n.direct_correctable_record",
    "n.common_load_baseline",
    "n.reserve_finite_presence",
    "n.reserve_poisson_projective_limit",
    "n.rc_load",
    "n.capacity_map_archive_190",
    "hierarchy.p_public_endpoint",
    "hierarchy.p_source_audit_witness",
    "hierarchy.pn_product_contraction_adapter",
    "hierarchy.n_global_repair_tick_adapter",
    "hierarchy.n_ew_bridge_defined_capacity",
)

EXPECTED_PARENT_PATHS = (
    "code/closure_preflight/data/source_only_closure_contract_v3.json",
    "code/closure_preflight/data/rc_load_workspace_provenance_v1.json",
    "code/closure_preflight/source_only_closure_preflight.py",
    "code/closure_preflight/verify_source_only_closure_preflight_independent.py",
    "code/closure_preflight/rc_load_source_only.py",
    "code/closure_preflight/outputs/rc_load_source_only_certificate.json",
    "code/P_derivation/runtime/p_interval_contraction_certificate_2026-07-14.json",
    "code/P_derivation/runtime/p_closure_trunk_current.json",
    "code/capacity_readback/runtime/F_candidate_capK_certificates.json",
    "code/capacity_readback/runtime/F_candidate_capP_certificates.json",
    "code/capacity_readback/runtime/F_candidate_capL_certificates.json",
    "code/capacity_readback/runtime/F_candidate_coupled_certificates.json",
    "code/capacity_readback/runtime/direct_n_closure_verdict.json",
    "code/capacity_readback/manifests/n_closure_branch_certificate.json",
    "code/capacity_readback/runtime/named_law_n_closure_verdict.json",
    "code/capacity_readback/manifests/capacity_semantics_menu_reference.json",
    "code/capacity_readback/runtime/F_construction_comparison_2026-07-14.json",
    "code/particles/forecast_contract/data/candidate_inventory_v2.json",
    "code/particles/forecast_contract/outputs/forecast_contract_state.json",
    "code/particles/hierarchy/certificates/R_P_public_pixel_certificate.json",
    "code/particles/hierarchy/certificates/R_P_source_audit_pixel_certificate.json",
    "code/particles/hierarchy/certificates/R_PN_joint_fixed_point_certificate_report.json",
    "code/particles/hierarchy/certificates/R_N_global_repair_tick_certificate.json",
    "code/particles/hierarchy/certificates/R_EW_global_capacity_certificate.json",
)
EXPECTED_PARENT_VALIDATORS = (
    "capacity_semantics_menu",
    "direct_n_closure",
    "bounded_packet_lift_independent",
    "n_closure_branch",
    "named_law_n_closure",
    "forecast_contract_state",
    "hierarchy_joint_fixed_point",
    "hierarchy_global_repair_tick",
    "hierarchy_ew_capacity",
)

RC_EQUATION = "ln(N/pi) * (6*pi/(P*alpha_U) - ln(N/pi)) = 6*pi"
COMMON_LOAD_FORMULA = "N0 = pi * exp(6*pi/(P*alpha_U))"
RESERVE_FORMULAS = {
    "finite_presence": "N = N0 * (1 - P/24)",
    "poisson_projective_limit": "N = N0 * exp(-P/24)",
}


class VerificationError(ValueError):
    """Fail-closed verification error."""


def check(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def load(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    check(isinstance(value, dict), f"JSON root is not an object: {path}")
    return value


def canonical_json_bytes(value: Any) -> bytes:
    return (
        json.dumps(value, ensure_ascii=True, separators=(",", ":"), sort_keys=True)
        + "\n"
    ).encode("ascii")


def tagged(raw: bytes) -> str:
    return "sha256:" + hashlib.sha256(raw).hexdigest()


def tagged_digest(value: Any) -> bool:
    if not isinstance(value, str) or not value.startswith("sha256:"):
        return False
    digest = value[7:]
    return len(digest) == 64 and all(character in "0123456789abcdef" for character in digest)


def relative(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def pointer(path: Path, json_pointer: str, observed: Any) -> dict[str, Any]:
    return {"artifact": relative(path), "json_pointer": json_pointer, "observed": observed}


def gate(attained: bool, evidence: list[dict[str, Any]]) -> dict[str, Any]:
    check(bool(evidence), "expected gate has no evidence")
    return {
        "attained": bool(attained),
        "state": "ATTAINED" if attained else "OPEN",
        "evidence": evidence,
    }


def row(
    *,
    candidate_id: str,
    registry_source: Path,
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
) -> dict[str, Any]:
    check(tuple(gates) == REQUIRED_GATES, f"expected gate order drift: {candidate_id}")
    attained = sum(gates[name]["attained"] is True for name in REQUIRED_GATES)
    return {
        "candidate_id": candidate_id,
        "registry_source": relative(registry_source),
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
        "row_gate_complete": attained == len(REQUIRED_GATES),
        "qualifies_source_only": False,
        "open_source_only_gates": [
            name for name in REQUIRED_GATES if gates[name]["attained"] is not True
        ],
        "details": details,
    }


def verify_self_hash(value: Mapping[str, Any], field: str) -> None:
    work = dict(value)
    recorded = work.pop(field, None)
    check(tagged_digest(recorded), f"{field} is not a tagged digest")
    check(recorded == tagged(canonical_json_bytes(work)), f"{field} mismatch")


def verify_contract(contract: dict[str, Any]) -> None:
    check(contract.get("schema") == CONTRACT_SCHEMA, "contract schema drift")
    check(contract.get("issue") == ISSUE, "contract issue drift")
    check(tuple(contract.get("source_only_required_gates", [])) == REQUIRED_GATES, "required gate drift")
    check(tuple(contract.get("preflight_completion_gates", [])) == COMPLETION_GATES, "completion gate drift")
    check(tuple(contract.get("explicitly_downstream_gates", [])) == DOWNSTREAM_GATES, "downstream gate drift")
    check(
        contract.get("scope_boundary")
        == {
            "future_source_laws_exhausted": False,
            "issue_closure_authorized": False,
            "physical_or_laboratory_attachment_required": False,
        },
        "contract scope boundary drift",
    )
    declared = contract.get("preflight_completion_state")
    check(isinstance(declared, dict) and tuple(declared) == COMPLETION_GATES, "completion declaration drift")
    for name in COMPLETION_GATES:
        state = declared[name]
        check(isinstance(state, dict), f"completion declaration malformed: {name}")
        check(isinstance(state.get("attained"), bool), f"completion attained is not boolean: {name}")
        check(isinstance(state.get("reason"), str) and state["reason"].strip(), f"completion reason missing: {name}")
    check([declared[name]["attained"] for name in COMPLETION_GATES] == [False, False, False, False, True], "completion declaration promoted")
    check(
        contract.get("decision_rule")
        == {
            "frozen": True,
            "candidate_present_requires_every_row_gate": True,
            "candidate_present_requires_every_preflight_completion_gate": True,
            "archive_enumeration_is_not_a_source_return_map": True,
            "fixed_cutoff_controls_are_not_cosmic_fixed_points": True,
        },
        "decision rule drift",
    )
    check(not set(REQUIRED_GATES).intersection(DOWNSTREAM_GATES), "downstream gate entered H0")


def verify_forecast_state(state: dict[str, Any]) -> None:
    check(state.get("schema") == "oph.forecast_contract_state.v2", "forecast-state schema drift")
    check(state.get("issue") == 639, "forecast-state issue drift")
    verify_self_hash(state, "state_sha256")
    check(state.get("eligible_candidates") == [], "forecast state unexpectedly has an eligible candidate")
    completeness = state.get("inventory_completeness", {})
    check(completeness.get("status") == "PROVISIONAL_NOT_EXHAUSTIVE", "forecast inventory was treated as exhaustive")
    check(completeness.get("closure_criterion_frozen") is False, "forecast inventory criterion unexpectedly frozen")
    check(state.get("controls", {}).get("candidate_inventory_exhaustiveness") == "NOT_ESTABLISHED", "forecast inventory exhaustiveness drift")
    check(state.get("selection_rule", {}).get("status") == "DRAFT_UNFROZEN", "forecast selection rule unexpectedly frozen")


def verify_parent_pins(pins: Any) -> None:
    check(isinstance(pins, list), "parent pins are not a list")
    paths = [entry.get("path") for entry in pins if isinstance(entry, dict)]
    check(tuple(paths) == EXPECTED_PARENT_PATHS, "parent pin order/inventory mismatch")
    check(len(paths) == len(set(paths)), "duplicate parent pin")
    for entry in pins:
        path_text = entry["path"]
        pure = PurePosixPath(path_text)
        check(path_text == pure.as_posix() and not pure.is_absolute() and ".." not in pure.parts, f"unsafe parent path: {path_text}")
        raw = (ROOT / path_text).read_bytes()
        check(entry.get("bytes") == len(raw), f"parent byte drift: {path_text}")
        check(entry.get("sha256") == tagged(raw), f"parent hash drift: {path_text}")


def resolve_pointer(document: Any, json_pointer: str) -> Any:
    check(isinstance(json_pointer, str) and json_pointer.startswith("/"), "malformed JSON pointer")
    check("*" not in json_pointer, f"wildcard JSON pointer rejected: {json_pointer}")
    current = document
    for raw_token in json_pointer[1:].split("/"):
        token = raw_token.replace("~1", "/").replace("~0", "~")
        if isinstance(current, list):
            check(token.isdigit(), f"nonnumeric list pointer token: {json_pointer}")
            index = int(token)
            check(0 <= index < len(current), f"list pointer out of range: {json_pointer}")
            current = current[index]
        elif isinstance(current, dict):
            if token not in current:
                return None
            current = current[token]
        else:
            return None
    return current


def verify_all_evidence(rows: list[dict[str, Any]], completion: dict[str, Any]) -> None:
    cache: dict[str, Any] = {}
    evidence_items: list[dict[str, Any]] = []
    for candidate in rows:
        for name in REQUIRED_GATES:
            evidence_items.extend(candidate["source_only_gates"][name]["evidence"])
    for name in COMPLETION_GATES:
        evidence_items.extend(completion["gates"][name]["evidence"])
    for item in evidence_items:
        check(isinstance(item, dict) and set(item) == {"artifact", "json_pointer", "observed"}, "evidence shape drift")
        artifact = item["artifact"]
        pointer_text = item["json_pointer"]
        check(isinstance(artifact, str), "evidence artifact is not a string")
        pure = PurePosixPath(artifact)
        check(artifact == pure.as_posix() and not pure.is_absolute() and ".." not in pure.parts, f"unsafe evidence path: {artifact}")
        check("*" not in pointer_text, f"wildcard JSON pointer rejected: {pointer_text}")
        path = ROOT / artifact
        check(path.is_file(), f"evidence artifact missing: {artifact}")
        document = cache.setdefault(artifact, load(path))
        check(resolve_pointer(document, pointer_text) == item["observed"], f"evidence mismatch: {artifact}{pointer_text}")


def scalar_endpoint(value: Any) -> str:
    text = str(value)
    if text.startswith("[") and text.endswith("]"):
        left, right = text[1:-1].split(",", 1)
        check(left.strip() == right.strip(), "interval endpoint is not scalar")
        return left.strip()
    return text


def verify_rc(rc: dict[str, Any], p_interval: dict[str, Any], provenance: dict[str, Any]) -> None:
    check(rc.get("schema") == "oph.rc_load_source_only_replay.v1", "RC schema drift")
    verify_self_hash(rc, "certificate_sha256")
    check(rc.get("equation") == RC_EQUATION, "RC equation drift")
    check(provenance.get("certificate", {}).get("equation") == RC_EQUATION, "RC provenance equation drift")
    check(provenance.get("certificate", {}).get("verdict") == "PASS", "RC provenance did not pass")
    check(provenance.get("certificate", {}).get("pricing_law_proved_here") is False, "RC pricing law was silently promoted")
    check(provenance.get("certificate", {}).get("typed_identity_proved_here") is False, "RC identity was silently promoted")
    check(provenance.get("external_files_required_for_replay") is False, "RC replay needs an external file")
    expected_scope = {
        "pricing_law_proved_here": False,
        "same_quantity_identity_proved_here": False,
        "observational_comparison_imported": False,
        "physical_or_laboratory_attachment_evaluated": False,
        "workspace_root_checkout_required": False,
    }
    check(rc.get("scope") == expected_scope, "RC no-target scope drift")
    check(
        rc.get("exposure")
        == {
            "class": provenance["campaign"]["exposure_class"],
            "promotion_from_landing": provenance["campaign"]["promotion_from_landing"],
        },
        "RC exposure drift",
    )
    check(rc["exposure"]["class"] == "exposed_retrospective", "RC exposure was retyped")
    check(rc["exposure"]["promotion_from_landing"] is False, "RC landing was promoted")

    mode = p_interval["modes"]["thomson_structured_running"]
    p_box = mode["certified_enclosure"]["P"]
    au_box = mode["interval_diagnostics"]["alpha_u"]
    check(
        rc.get("inputs")
        == {
            "P": {"lo": p_box["lo"], "hi": p_box["hi"]},
            "alpha_U": {"lo": au_box["lo"], "hi": au_box["hi"]},
            "P_mode": "thomson_structured_running",
        },
        "RC source inputs drift",
    )
    mp.dps = 80
    iv.dps = 80
    p_value = iv.mpf([p_box["lo"], p_box["hi"]])
    alpha_u = iv.mpf([au_box["lo"], au_box["hi"]])
    budget = 6 * iv.pi / (p_value * alpha_u)
    discriminant = budget**2 - 24 * iv.pi
    root = (budget + iv.sqrt(discriminant)) / 2
    micro = (budget - iv.sqrt(discriminant)) / 2
    domain = iv.mpf([mpf(10), mpf(10)])
    image_at_lo = budget - 6 * iv.pi / domain
    lipschitz = (6 * iv.pi / domain**2).b
    local_lipschitz = (6 * iv.pi / root**2).b
    residual = root * (budget - root) - 6 * iv.pi
    enclosure = lambda value: {"lo": scalar_endpoint(value.a), "hi": scalar_endpoint(value.b)}
    checks = {
        "positive_discriminant": bool(discriminant.a > 0),
        "maps_into_domain": bool(image_at_lo.a > 10 and budget.b > image_at_lo.b),
        "contraction": bool(lipschitz < 1),
        "unique_fixed_point_in_domain": bool(image_at_lo.a > 10 and budget.b > image_at_lo.b and lipschitz < 1),
        "residual_encloses_zero": bool(residual.a <= 0 <= residual.b),
        "micro_root_outside_domain": bool(micro.b < 10),
    }
    check(rc.get("budget_x") == enclosure(budget), "RC budget mismatch")
    check(rc.get("fixed_point", {}).get("X") == enclosure(root), "RC fixed point mismatch")
    check(rc.get("fixed_point", {}).get("residual") == enclosure(residual), "RC residual mismatch")
    check(rc.get("excluded_micro_root", {}).get("X") == enclosure(micro), "RC micro-root mismatch")
    check(rc.get("validation") == checks and all(checks.values()), "RC validation mismatch")
    check(rc.get("banach") == {
        "domain_lo": "10",
        "domain_hi": scalar_endpoint(budget.b),
        "maps_into_domain": checks["maps_into_domain"],
        "sup_lipschitz": scalar_endpoint(lipschitz),
        "local_lipschitz_at_fixed_point": scalar_endpoint(local_lipschitz),
        "unique_fixed_point_in_domain": checks["unique_fixed_point_in_domain"],
    }, "RC Banach payload mismatch")
    check(rc.get("status") == "CONDITIONAL_FIXED_POINT_CERTIFIED" and rc.get("verdict") == "PASS", "RC did not pass")
    expected_rc_parents = []
    for source in (P_INTERVAL, RC_PROVENANCE):
        raw = source.read_bytes()
        expected_rc_parents.append({"path": relative(source), "bytes": len(raw), "sha256": tagged(raw)})
    check(rc.get("parent_pins") == expected_rc_parents, "RC parent pins drift")
    check(rc.get("upstream_workspace_provenance") == {
        "source_sha256": provenance["source"]["sha256"],
        "certificate_sha256": provenance["certificate"]["sha256"],
        "campaign_sha256": provenance["campaign"]["sha256"],
    }, "RC upstream provenance drift")


def expected_p_rows(p_interval: dict[str, Any], p_trunk: dict[str, Any], exposure: dict[str, Any]) -> list[dict[str, Any]]:
    check(p_interval.get("claim_status") == "interval_contraction_certificate_for_declared_closure_map", "P certificate status drift")
    check(tuple(p_interval.get("modes", {})) == P_MODES, "P mode inventory drift")
    consumer = p_trunk.get("consumer_policy", {})
    inventory = exposure.get("candidates", {}).get("x_p_certificate_alpha", {})
    eligibility = inventory.get("eligibility")
    freeze = inventory.get("freeze_evidence")
    freeze_exposure = freeze.get("exposure_class") if isinstance(freeze, dict) else None
    selected = bool(
        eligibility in P_ELIGIBILITY_ALLOWLIST
        and isinstance(freeze, dict)
        and freeze.get("frozen_before_comparison") is True
        and freeze_exposure in QUALIFYING_EXPOSURES
        and tagged_digest(freeze.get("content_sha256"))
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
        existence = banach.get("existence") is True
        uniqueness = banach.get("uniqueness_in_interval") is True
        stability = banach.get("contraction") is True and Decimal(str(banach.get("lipschitz_bound"))) < 1
        gates = {
            "declared_dimensionless_coordinate": gate(str(mode.get("space", "")).startswith("alpha (fine-structure coupling)"), [pointer(P_INTERVAL, f"/modes/{mode_name}/space", mode.get("space"))]),
            "target_independent_candidate_selection": gate(selected, [pointer(EXPOSURE, "/candidates/x_p_certificate_alpha/eligibility", eligibility), pointer(EXPOSURE, "/candidates/x_p_certificate_alpha/freeze_evidence", freeze)]),
            "same_quantity_constructed": gate(same_quantity, [pointer(EXPOSURE, "/candidates/x_p_certificate_alpha/same_quantity_constructed", inventory.get("same_quantity_constructed"))]),
            "source_return_map_complete": gate(source_complete, [pointer(P_INTERVAL, "/promotion_allowed", p_interval.get("promotion_allowed")), pointer(P_TRUNK, "/consumer_policy/may_feed_live_particle_predictions", consumer.get("may_feed_live_particle_predictions")), pointer(P_TRUNK, "/consumer_policy/default_thomson_endpoint_allowed", consumer.get("default_thomson_endpoint_allowed"))]),
            "existence_certified": gate(existence, [pointer(P_INTERVAL, f"/modes/{mode_name}/banach/existence", banach.get("existence"))]),
            "uniqueness_certified": gate(uniqueness, [pointer(P_INTERVAL, f"/modes/{mode_name}/banach/uniqueness_in_interval", banach.get("uniqueness_in_interval"))]),
            "stability_certified": gate(stability, [pointer(P_INTERVAL, f"/modes/{mode_name}/banach/lipschitz_bound", banach.get("lipschitz_bound"))]),
        }
        enclosure = mode["certified_enclosure"]
        rows.append(row(
            candidate_id=f"p.{mode_name}", registry_source=P_INTERVAL, row_kind="candidate_map", quantity="P", family="P_INTERVAL_MAP",
            mathematical_status="INTERVAL_FIXED_POINT_CERTIFIED",
            disposition="CONDITIONAL_ARITHMETIC_FIXED_POINT__SELECTION_IDENTITY_AND_RETURN_OPEN",
            exposure_class=str(freeze_exposure) if selected else "registered_source_audit__not_frozen_precomparison",
            solver_target_payload_read=bool(consumer.get("hidden_external_alpha_allowed")),
            exact_fixed_point_progress=existence and uniqueness and stability,
            gates=gates,
            details={"mode": mode_name, "certified_enclosure": {"P": enclosure.get("P"), "alpha_inverse": enclosure.get("alpha_inv")}},
        ))
    check(Decimal(rows[1]["details"]["certified_enclosure"]["P"]["hi"]) < Decimal(rows[0]["details"]["certified_enclosure"]["P"]["lo"]), "P intervals overlap")
    return rows


def expected_menu_rows(menu: dict[str, Any], cap_k: dict[str, Any], cap_p: dict[str, Any], cap_l: dict[str, Any], coupled: dict[str, Any]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    menu_rows = menu.get("menu_rows", [])
    check(menu.get("schema") == "oph.capacity_semantics_menu_certificate.v1", "capacity menu schema drift")
    check(menu.get("menu_row_count") == len(menu_rows) == len(CAPACITY_IDS), "capacity menu count drift")
    check(tuple(str(item.get("row_id")) for item in menu_rows) == CAPACITY_IDS, "capacity menu identity/order drift")
    family_counts = Counter(item.get("family") for item in menu_rows)
    check(dict(family_counts) == CAPACITY_FAMILIES, "capacity menu family counts drift")
    check(menu.get("campaign_verdict", {}).get("semantics_enumeration") == "declared_menu_complete_for_executed_families", "capacity menu status drift")

    cap_k_rows = cap_k.get("rows", [])
    cap_p_rows = cap_p.get("rows", [])
    cap_l_rows = cap_l.get("rows", [])
    check([item.get("branch") for item in cap_k_rows] == list(CAPACITY_IDS[:4]), "CAP-K runtime identity drift")
    check([item.get("branch") for item in cap_p_rows] == list(CAPACITY_IDS[4:10]), "CAP-P runtime identity/count drift")
    check(len(cap_p_rows) == 6, "CAP-P must have exactly six runtime rows")
    cap_l_groups: dict[str, list[dict[str, Any]]] = {f"capL.R{i}": [] for i in range(1, 6)}
    for source in cap_l_rows:
        prefix = ".".join(str(source.get("branch")).split(".")[:2])
        check(prefix in cap_l_groups, f"unexpected CAP-L prefix: {prefix}")
        cap_l_groups[prefix].append(source)
    check(len(cap_l_rows) == 180 and all(len(group) == 36 for group in cap_l_groups.values()), "CAP-L 5x36 lattice drift")

    selector = menu["campaign_verdict"]["source_only_fixed_point_selector"]
    no_target_import = menu["campaign_verdict"].get("no_target_import")
    selector_attained = selector in QUALIFYING_EXPOSURES and no_target_import == "verified"
    same_quantity = menu["campaign_verdict"].get("same_quantity_constructed") == "verified"
    result: list[dict[str, Any]] = []
    for index, menu_row in enumerate(menu_rows):
        row_id = menu_row["row_id"]
        family = menu_row["family"]
        channel = menu_row["semantics"]["channel"]
        source_return = False
        solver_target = False
        exact_progress = False
        existence = False
        uniqueness = False
        stability = False
        exposure_class = "declared_before_comparison__retrospective_comparison_recorded"
        if family == "CAP-K":
            source = cap_k_rows[index]
            check(source["branch"] == row_id and source["status"] == menu_row["executed_verdict"], f"CAP-K crosswalk drift: {row_id}")
            stability = Decimal(source["s"]["hi"]) < 1
            status = source["status"]
            disposition = "EXACT_ROW_NO_POSITIVE_FIXED_POINT__BOUNDED_FAMILY_ONLY"
            details = {"menu_row_id": row_id, "map": source["map"], "s": source["s"]}
            runtime_path, runtime_pointer = CAP_K, f"/rows/{index}"
        elif family == "CAP-P":
            source_index = index - 4
            source = cap_p_rows[source_index]
            check(source["branch"] == row_id, f"CAP-P crosswalk drift: {row_id}")
            status = source["status"]
            fixed = status == "fixed_point_certified"
            existence = fixed
            uniqueness = fixed and source.get("fixed_point_exact") == "N = pi (the unique solution of (N/pi)^s = N/pi with s != 1)"
            stability = source.get("contraction_certificate", {}).get("banach_pass") is True if fixed else Decimal(source["s"]["hi"]) < 1
            exact_progress = bool(existence and uniqueness and stability)
            check(menu_row["executed_verdict"] == ("excluded" if fixed else "no_positive_fixed_point"), f"CAP-P menu verdict drift: {row_id}")
            disposition = "EXACT_NORMALIZATION_FIXED_POINT_N_EQUALS_PI__SOURCE_SELECTION_IDENTITY_AND_RETURN_OPEN" if fixed else "EXACT_ROW_NO_POSITIVE_FIXED_POINT__BOUNDED_FAMILY_ONLY"
            details = {"menu_row_id": row_id, "map": source["map"], "s": source["s"], "fixed_point_exact": source.get("fixed_point_exact"), "fixed_point_enclosure": source.get("fixed_point", {}).get("enclosure")}
            runtime_path, runtime_pointer = CAP_P, f"/rows/{source_index}"
        elif family == "CAP-L":
            group = cap_l_groups[row_id]
            counts = Counter(item["status"] for item in group)
            check(dict(counts) == menu_row.get("recorded_status_counts"), f"CAP-L count drift: {row_id}")
            existence = counts.get("fixed_point_certified", 0) > 0
            uniqueness = False
            stability = not any(counts.get(name, 0) > 0 for name in ("fixed_point_unstable_rejected", "rejected_no_contraction"))
            exact_progress = existence
            status = "MIXED_36_ROW_SUBLATTICE"
            disposition = "ARCHIVED_MIXED_SUBLATTICE__NO_SOURCE_SELECTOR"
            details = {"menu_row_id": row_id, "sublattice_row_count": len(group), "status_counts": dict(sorted(counts.items()))}
            runtime_path, runtime_pointer = SEMANTICS_MENU, f"/menu_rows/{index}/recorded_status_counts"
        elif family == "CAP-B":
            check(row_id == "capB.bridge_constant" and menu_row["executed_verdict"] == "excluded_pre_evaluation", "CAP-B crosswalk drift")
            status = menu_row["executed_verdict"]
            disposition = "BARRED_TARGET_BRIDGE_PRE_EVALUATION__NEVER_EXECUTED"
            exposure_class = "barred_target_bridge_pre_evaluation"
            details = {"menu_row_id": row_id, "executed_verdict": status}
            runtime_path, runtime_pointer = SEMANTICS_MENU, f"/menu_rows/{index}/executed_verdict"
        elif family == "coupled":
            check(row_id == "coupled.cp1_cp2_cp3" and menu_row["executed_verdict"] == "conditional_open", "coupled crosswalk drift")
            source = coupled["certificate"]["load_coordinate"]
            contraction = source["contraction_certificate"]
            existence = source["fixed_point"]["located"] is True
            uniqueness = source["exact_solution_check"]["x_ew_inside_certified_box"] is True
            stability = contraction["banach_pass"] is True
            exact_progress = existence and uniqueness and stability
            status = coupled["status"]
            disposition = "CONDITIONAL_FIXED_POINT__CP1_CP2_CP3_OPEN__TARGET_EXPOSED"
            exposure_class = "target_exposed_theorem_coupled"
            solver_target = True
            details = {"menu_row_id": row_id, "lambda": source["lambda"], "fixed_point_enclosure": source["fixed_point"]["enclosure"], "cl7_status": coupled["cl7_status"], "evaluation_cone_contains_bridge_expression": coupled["blindness"]["cone_contains_cl3_bridge_expression"]}
            runtime_path, runtime_pointer = COUPLED, "/certificate/load_coordinate"
        else:
            raise VerificationError(f"unexpected menu family: {family}")
        runtime_value = resolve_pointer(load(runtime_path), runtime_pointer)
        gates = {
            "declared_dimensionless_coordinate": gate(bool(channel), [pointer(SEMANTICS_MENU, f"/menu_rows/{index}/semantics/channel", channel)]),
            "target_independent_candidate_selection": gate(selector_attained, [pointer(SEMANTICS_MENU, "/campaign_verdict/source_only_fixed_point_selector", selector), pointer(SEMANTICS_MENU, "/campaign_verdict/no_target_import", no_target_import), pointer(SEMANTICS_MENU, f"/menu_rows/{index}/executed_verdict", menu_row["executed_verdict"])]),
            "same_quantity_constructed": gate(same_quantity, [pointer(SEMANTICS_MENU, "/campaign_verdict/same_quantity_constructed", menu["campaign_verdict"].get("same_quantity_constructed"))]),
            "source_return_map_complete": gate(source_return, [pointer(SEMANTICS_MENU, "/campaign_verdict/horizon_area_and_ew_bridge", menu["campaign_verdict"].get("horizon_area_and_ew_bridge"))]),
            "existence_certified": gate(existence, [pointer(runtime_path, runtime_pointer, runtime_value)]),
            "uniqueness_certified": gate(uniqueness, [pointer(runtime_path, runtime_pointer, runtime_value)]),
            "stability_certified": gate(stability, [pointer(runtime_path, runtime_pointer, runtime_value)]),
        }
        result.append(row(
            candidate_id=f"n.menu.{row_id}", registry_source=SEMANTICS_MENU,
            row_kind="registered_control" if family == "CAP-B" else ("candidate_family" if family == "CAP-L" else "candidate_map"),
            quantity="N_CAPACITY_MENU", family=family, mathematical_status=status,
            disposition=disposition, exposure_class=exposure_class,
            solver_target_payload_read=solver_target, exact_fixed_point_progress=exact_progress,
            gates=gates, details=details,
        ))
    crosswalk = {
        "authoritative_menu_path": relative(SEMANTICS_MENU),
        "authoritative_menu_row_count": 17,
        "authoritative_menu_row_ids": list(CAPACITY_IDS),
        "family_counts": dict(family_counts),
        "runtime_counts": {"CAP-K": 4, "CAP-P": 6, "CAP-L": 180, "CAP-B": 0, "coupled": 1},
        "inventory_order_derived_from_authoritative_sources": True,
        "contract_candidate_order_used": False,
    }
    return result, crosswalk


def expected_direct_row(direct: dict[str, Any]) -> dict[str, Any]:
    fixed = direct["fixed_cutoff_result"]
    generated = direct["bounded_generation_register_result"]
    controls = direct["source_controls"]
    check(fixed == {"D": 24, "M0": 24, "cosmic_value_selected": False, "status": "SOURCE_DERIVED_FIXED_CUTOFF_PACKET"}, "direct D=M0 control drift")
    declared = direct["typed_coordinates"]["input_coordinate"] == "N_in=log(D)" and direct["typed_coordinates"]["output_coordinate"] == "N_out=log(M0)"
    selected = bool(controls.get("target_clean") is True and generated.get("unique_source_zero_entailed") is True and direct.get("freeze_evidence", {}).get("frozen_before_comparison") is True)
    source_complete = bool(generated.get("universal_all_rung_membership_proved") is True and generated.get("executable_lean_membership_bridge_proved") is True)
    gates = {
        "declared_dimensionless_coordinate": gate(declared, [pointer(DIRECT_N, "/typed_coordinates", direct["typed_coordinates"])]),
        "target_independent_candidate_selection": gate(selected, [pointer(DIRECT_N, "/source_controls", controls), pointer(DIRECT_N, "/source_controls/constructor_reads_desired_capacity", controls.get("constructor_reads_desired_capacity")), pointer(DIRECT_N, "/bounded_generation_register_result/unique_source_zero_entailed", generated.get("unique_source_zero_entailed")), pointer(DIRECT_N, "/freeze_evidence", direct.get("freeze_evidence"))]),
        "same_quantity_constructed": gate(direct.get("same_quantity_constructed") is True, [pointer(DIRECT_N, "/same_quantity_constructed", direct.get("same_quantity_constructed"))]),
        "source_return_map_complete": gate(source_complete, [pointer(DIRECT_N, "/bounded_generation_register_result/universal_all_rung_membership_proved", generated.get("universal_all_rung_membership_proved")), pointer(DIRECT_N, "/bounded_generation_register_result/executable_lean_membership_bridge_proved", generated.get("executable_lean_membership_bridge_proved"))]),
        "existence_certified": gate(False, [pointer(DIRECT_N, "/bounded_generation_register_result/unique_source_zero_entailed", generated.get("unique_source_zero_entailed"))]),
        "uniqueness_certified": gate(False, [pointer(DIRECT_N, "/bounded_generation_register_result/unique_source_zero_entailed", generated.get("unique_source_zero_entailed"))]),
        "stability_certified": gate(False, [pointer(DIRECT_N, "/stability_certificate", direct.get("stability_certificate"))]),
    }
    return row(
        candidate_id="n.direct_correctable_record", registry_source=DIRECT_N, row_kind="fixed_cutoff_control", quantity="N_DIRECT", family="DIRECT_CORRECTABLE_RECORD",
        mathematical_status=direct["status"], disposition="FIXED_CUTOFF_CONTROL__COSMIC_SOURCE_ANTECEDENT_INCOMPLETE",
        exposure_class="target_clean_bounded_control__not_a_candidate_map", solver_target_payload_read=not bool(controls.get("target_clean")),
        exact_fixed_point_progress=False, gates=gates,
        details={"fixed_cutoff_control": True, "D": fixed["D"], "M0": fixed["M0"], "cosmic_value_selected": fixed["cosmic_value_selected"]},
    )


def decimal80(value: Any) -> str:
    return mp.nstr(value, 80, strip_zeros=False)


def verify_common_arithmetic(branch: dict[str, Any]) -> None:
    check(branch.get("conditional_baseline", {}).get("formula") == COMMON_LOAD_FORMULA, "common-load formula drift")
    branch_rows = branch.get("branches", [])
    check([item.get("branch_id") for item in branch_rows] == list(RESERVE_FORMULAS), "reserve branch identity/order drift")
    for item in branch_rows:
        check(item.get("formula") == RESERVE_FORMULAS[item["branch_id"]], f"reserve formula drift: {item['branch_id']}")
    mp.dps = 100
    p_value = mp.mpf(branch["source_inputs"]["P"])
    alpha_u = mp.mpf(branch["source_inputs"]["alpha_U"])
    reserve = p_value / 24
    log_n0 = 6 * mp.pi / (p_value * alpha_u)
    n0 = mp.pi * mp.exp(log_n0)
    expected = {
        "baseline_log": decimal80(log_n0),
        "baseline_n": decimal80(n0),
        "finite_factor": decimal80(1 - reserve),
        "finite_n": decimal80(n0 * (1 - reserve)),
        "finite_log": decimal80(mp.log(n0 * (1 - reserve) / mp.pi)),
        "poisson_factor": decimal80(mp.exp(-reserve)),
        "poisson_n": decimal80(n0 * mp.exp(-reserve)),
        "poisson_log": decimal80(mp.log(n0 * mp.exp(-reserve) / mp.pi)),
    }
    baseline = branch["conditional_baseline"]
    finite, poisson = branch_rows
    check(baseline["log_N0_over_pi"] == expected["baseline_log"] and baseline["N0"] == expected["baseline_n"], "common-load arithmetic mismatch")
    check((finite["factor"], finite["N"], finite["log_N_over_pi"]) == (expected["finite_factor"], expected["finite_n"], expected["finite_log"]), "finite-reserve arithmetic mismatch")
    check((poisson["factor"], poisson["N"], poisson["log_N_over_pi"]) == (expected["poisson_factor"], expected["poisson_n"], expected["poisson_log"]), "Poisson-reserve arithmetic mismatch")
    check(finite["mean_count_or_projective_limit_carrier_required"] is False and poisson["mean_count_or_projective_limit_carrier_required"] is True, "reserve carrier typing drift")


def expected_common_rows(branch: dict[str, Any], named: dict[str, Any]) -> list[dict[str, Any]]:
    verify_common_arithmetic(branch)
    scope = branch["scope"]
    identity = branch["self_reference_boundary"]
    source_inputs = branch["source_inputs"]
    common_source = scope.get("global_capacity_derived") is True and source_inputs.get("source_artifact_status") == "source_complete_endpoint_proof"
    common_same = identity.get("same_typed_quantity_identified") is True
    common_selected = scope.get("target_blind_forecast") is True and scope.get("retrospective") is False and scope.get("branch_selected") is True
    common_unique = identity.get("unique_source_fixed_point_proved") is True
    common_stable = named.get("stability_certified") is True

    def gates(formula_pointer: str, formula: Any, selected: bool) -> dict[str, dict[str, Any]]:
        return {
            "declared_dimensionless_coordinate": gate(bool(formula), [pointer(N_BRANCH, formula_pointer, formula)]),
            "target_independent_candidate_selection": gate(bool(common_selected and selected), [pointer(N_BRANCH, "/scope", scope)]),
            "same_quantity_constructed": gate(common_same, [pointer(N_BRANCH, "/self_reference_boundary/same_typed_quantity_identified", identity.get("same_typed_quantity_identified"))]),
            "source_return_map_complete": gate(common_source, [pointer(N_BRANCH, "/scope/global_capacity_derived", scope.get("global_capacity_derived")), pointer(N_BRANCH, "/source_inputs/source_artifact_status", source_inputs.get("source_artifact_status"))]),
            "existence_certified": gate(False, [pointer(N_BRANCH, "/status", branch.get("status"))]),
            "uniqueness_certified": gate(common_unique, [pointer(N_BRANCH, "/self_reference_boundary/unique_source_fixed_point_proved", identity.get("unique_source_fixed_point_proved"))]),
            "stability_certified": gate(common_stable, [pointer(NAMED_N, "/stability_certified", named.get("stability_certified"))]),
        }

    baseline = branch["conditional_baseline"]
    result = [row(
        candidate_id="n.common_load_baseline", registry_source=N_BRANCH, row_kind="candidate_map", quantity="N_COMMON_LOAD", family="COMMON_LOAD",
        mathematical_status="EXACT_CONDITIONAL_BASELINE_ARITHMETIC", disposition="RETROSPECTIVE_CONDITIONAL_ARITHMETIC__NO_SOURCE_SELECTION_OR_RETURN",
        exposure_class="exposed_retrospective", solver_target_payload_read=scope.get("comparison_data_consumed") is True,
        exact_fixed_point_progress=False, gates=gates("/conditional_baseline/formula", baseline["formula"], scope.get("branch_selected") is True),
        details={"formula": baseline["formula"], "evaluated_N": baseline["N0"]},
    )]
    for index, source in enumerate(branch["branches"]):
        branch_id = source["branch_id"]
        result.append(row(
            candidate_id=f"n.reserve_{branch_id}", registry_source=N_BRANCH, row_kind="candidate_map", quantity="N_COMMON_LOAD", family="RESERVE_BRANCH",
            mathematical_status=source["status"], disposition="RETROSPECTIVE_UNSELECTED_CONDITIONAL_ARITHMETIC",
            exposure_class="exposed_retrospective", solver_target_payload_read=scope.get("comparison_data_consumed") is True,
            exact_fixed_point_progress=False, gates=gates(f"/branches/{index}/formula", source["formula"], source["selected"] is True),
            details={"branch_id": branch_id, "formula": source["formula"], "evaluated_N": source["N"], "selected": source["selected"], "mean_count_or_projective_limit_carrier_required": source["mean_count_or_projective_limit_carrier_required"]},
        ))
    return result


def expected_rc_row(rc: dict[str, Any]) -> dict[str, Any]:
    validation = rc["validation"]
    scope = rc["scope"]
    exposure = rc["exposure"]
    existence = validation["residual_encloses_zero"] is True
    uniqueness = validation["unique_fixed_point_in_domain"] is True
    stability = validation["contraction"] is True
    selected = bool(exposure.get("class") in QUALIFYING_EXPOSURES and exposure.get("promotion_from_landing") is False and rc.get("freeze_evidence", {}).get("frozen_before_comparison") is True)
    gates = {
        "declared_dimensionless_coordinate": gate(rc["equation"] == RC_EQUATION, [pointer(RC_OUTPUT, "/equation", rc["equation"])]),
        "target_independent_candidate_selection": gate(selected, [pointer(RC_OUTPUT, "/exposure", exposure), pointer(RC_OUTPUT, "/freeze_evidence", rc.get("freeze_evidence"))]),
        "same_quantity_constructed": gate(scope["same_quantity_identity_proved_here"] is True, [pointer(RC_OUTPUT, "/scope/same_quantity_identity_proved_here", scope["same_quantity_identity_proved_here"])]),
        "source_return_map_complete": gate(scope["pricing_law_proved_here"] is True, [pointer(RC_OUTPUT, "/scope/pricing_law_proved_here", scope["pricing_law_proved_here"])]),
        "existence_certified": gate(existence, [pointer(RC_OUTPUT, "/validation/residual_encloses_zero", validation["residual_encloses_zero"])]),
        "uniqueness_certified": gate(uniqueness, [pointer(RC_OUTPUT, "/validation/unique_fixed_point_in_domain", validation["unique_fixed_point_in_domain"])]),
        "stability_certified": gate(stability, [pointer(RC_OUTPUT, "/validation/contraction", validation["contraction"])]),
    }
    return row(
        candidate_id="n.rc_load", registry_source=RC_OUTPUT, row_kind="candidate_map", quantity="N_RC_LOAD", family="RC_LOAD",
        mathematical_status=rc["status"], disposition="CONDITIONAL_FIXED_POINT__EXPOSED_CAMPAIGN__PRICING_AND_IDENTITY_OPEN",
        exposure_class="exposed_retrospective_upstream_campaign__target_free_replay",
        solver_target_payload_read=scope["observational_comparison_imported"] is True,
        exact_fixed_point_progress=existence and uniqueness and stability, gates=gates,
        details={"equation": rc["equation"], "budget_x": rc["budget_x"], "fixed_point_X": rc["fixed_point"]["X"]},
    )


def expected_archive_row(archive: dict[str, Any], cap_k: dict[str, Any], cap_p: dict[str, Any], cap_l: dict[str, Any]) -> dict[str, Any]:
    archive_rows = archive.get("rows", [])
    expected_ids = {item["branch"] for item in cap_k["rows"] + cap_p["rows"] + cap_l["rows"]}
    check(len(expected_ids) == 190, "runtime archive source identity count drift")
    check({str(item.get("branch")) for item in archive_rows} == expected_ids, "190-row archive membership drift")
    check(archive.get("total_rows") == len(archive_rows) == 190, "190-row archive size drift")
    counts = Counter(item["status"] for item in archive_rows)
    certified = archive.get("certified_rows")
    check(certified == counts.get("fixed_point_certified", 0), "archive certified count drift")
    gates = {
        "declared_dimensionless_coordinate": gate(bool(archive.get("certified_fixed_point_range_nats")), [pointer(ARCHIVE_190, "/certified_fixed_point_range_nats", archive.get("certified_fixed_point_range_nats"))]),
        "target_independent_candidate_selection": gate(False, [pointer(ARCHIVE_190, "/landing_criterion", archive.get("landing_criterion"))]),
        "same_quantity_constructed": gate(False, [pointer(ARCHIVE_190, "/verdict", archive.get("verdict"))]),
        "source_return_map_complete": gate(False, [pointer(SEMANTICS_MENU, "/campaign_verdict/semantics_enumeration", "declared_menu_complete_for_executed_families"), pointer(CONTRACT, "/decision_rule/archive_enumeration_is_not_a_source_return_map", True)]),
        "existence_certified": gate(isinstance(certified, int) and certified > 0, [pointer(ARCHIVE_190, "/certified_rows", certified)]),
        "uniqueness_certified": gate(False, [pointer(ARCHIVE_190, "/certified_rows", certified)]),
        "stability_certified": gate(False, [pointer(ARCHIVE_190, "/total_rows", archive.get("total_rows")), pointer(ARCHIVE_190, "/certified_rows", certified)]),
    }
    return row(
        candidate_id="n.capacity_map_archive_190", registry_source=ARCHIVE_190, row_kind="archive_aggregate", quantity="N_ARCHIVED_MAP_LATTICE", family="CAP_L_P_K_ARCHIVE",
        mathematical_status="DECLARED_EXECUTED_FAMILY_ARCHIVE", disposition="AGGREGATE_METADATA_ONLY__ENUMERATION_IS_NOT_RETURN_MAP",
        exposure_class="declared_before_comparison__retrospective_comparison_archive", solver_target_payload_read=True,
        exact_fixed_point_progress=isinstance(certified, int) and certified > 0, gates=gates,
        details={"total_rows": 190, "certified_rows": certified, "landed_rows": archive.get("landed_rows"), "status_counts": dict(sorted(counts.items())), "covered_runtime_families": ["CAP-K", "CAP-P", "CAP-L"], "omitted_registered_menu_rows": ["capB.bridge_constant", "coupled.cp1_cp2_cp3"], "source_return_map_constructed": False, "authoritative_menu_row_count": 17},
    )


def expected_hierarchy_rows() -> list[dict[str, Any]]:
    p_public = load(HIER_P_PUBLIC)
    p_source = load(HIER_P_SOURCE)
    joint = load(HIER_PN_JOINT)
    tick = load(HIER_N_TICK)
    ew = load(HIER_EW_CAPACITY)

    def empty(path: Path, status: Any, declared: bool = True) -> dict[str, dict[str, Any]]:
        evidence = [pointer(path, "/status", status)]
        return {
            "declared_dimensionless_coordinate": gate(declared, evidence),
            "target_independent_candidate_selection": gate(False, evidence),
            "same_quantity_constructed": gate(False, evidence),
            "source_return_map_complete": gate(False, evidence),
            "existence_certified": gate(False, evidence),
            "uniqueness_certified": gate(False, evidence),
            "stability_certified": gate(False, evidence),
        }

    check(p_public.get("status") == "conditional_public_endpoint_certificate", "hierarchy public-P status drift")
    check(p_source.get("status") == "source_audit_branch_witness_not_full_endpoint_proof", "hierarchy source-P status drift")
    check(joint.get("status") == "closed_product_branch_theorem_with_explicit_coupled_branch_boundary", "hierarchy product status drift")
    check("CIRCULAR_DIAGNOSTIC_ONLY" in joint.get("N_backsolved_warning", ""), "hierarchy product circularity warning missing")
    check(tick.get("status") == "closed_global_repair_tick_theorem_with_derived_round_count", "hierarchy tick status drift")
    check(ew.get("status") == "closed_bridge_refined_global_capacity_fixed_point_certificate", "hierarchy EW status drift")
    ew_scope = ew.get("claim_boundary", {}).get("scope", "")
    check(
        isinstance(ew_scope, str)
        and "restricted to the source-side closed-form fixed point" in ew_scope
        and "not an exact bridge certificate" in ew_scope,
        "hierarchy EW scope drift",
    )

    result = [
        row(candidate_id="hierarchy.p_public_endpoint", registry_source=HIER_P_PUBLIC, row_kind="packet_disposition", quantity="P", family="HIERARCHY_PACKET", mathematical_status=p_public["status"], disposition="TARGET_EXPOSED_PUBLIC_ENDPOINT_CONTROL__INELIGIBLE_SOURCE_ONLY", exposure_class="target_exposed_public_endpoint_locator", solver_target_payload_read=True, exact_fixed_point_progress=False, gates=empty(HIER_P_PUBLIC, p_public["status"]), details={"artifact": p_public["artifact"], "public_endpoint_convention": p_public["public_endpoint_convention"], "not_supplied_here": p_public["proof_rule"]["not_supplied_here"]}),
        row(candidate_id="hierarchy.p_source_audit_witness", registry_source=HIER_P_SOURCE, row_kind="packet_disposition", quantity="P", family="HIERARCHY_PACKET", mathematical_status=p_source["status"], disposition="SOURCE_AUDIT_WITNESS__NOT_DISTINCT_CANDIDATE__ENDPOINT_PROOF_OPEN", exposure_class="source_audit_witness__not_frozen_candidate", solver_target_payload_read=False, exact_fixed_point_progress=False, gates=empty(HIER_P_SOURCE, p_source["status"]), details={"artifact": p_source["artifact"], "P_cand": p_source["P_cand"], "needed_for_public_endpoint_completion": p_source["needed_for_public_endpoint_completion"]}),
        row(candidate_id="hierarchy.pn_product_contraction_adapter", registry_source=HIER_PN_JOINT, row_kind="packet_disposition", quantity="P_N_PRODUCT", family="HIERARCHY_PACKET", mathematical_status=joint["status"], disposition="CONDITIONAL_PRODUCT_ADAPTER__SOURCE_N_COMPONENT_ABSENT__TARGET_EXPOSED_DIAGNOSTIC", exposure_class="target_exposed_numeric_demo_and_circular_backsolve", solver_target_payload_read=True, exact_fixed_point_progress=False, gates=empty(HIER_PN_JOINT, joint["status"]), details={"issue": joint["issue"], "product_condition": joint["product_contraction_certificate"], "coupled_boundary": joint["coupled_contraction_certificate"], "N_display_warning": joint["N_display_warning"], "N_backsolved_warning": joint["N_backsolved_warning"]}),
        row(candidate_id="hierarchy.n_global_repair_tick_adapter", registry_source=HIER_N_TICK, row_kind="packet_disposition", quantity="N_DEPENDENT_TICK", family="HIERARCHY_PACKET", mathematical_status=tick["status"], disposition="CONDITIONAL_TICK_IDENTITY__NOT_AN_N_PRODUCER__TARGET_EXPOSED_DISPLAY", exposure_class="downstream_conditional__contains_exposed_rounded_display", solver_target_payload_read=True, exact_fixed_point_progress=False, gates=empty(HIER_N_TICK, tick["status"]), details={"artifact": tick["artifact"], "normalization": tick["normalization"], "declared_not_derived": tick["claim_boundary"]["declared_not_derived"], "numeric_display": tick["numeric_display_for_rounded_capacity"]}),
    ]
    contraction = ew["contraction_certificate"]
    ew_gates = {
        "declared_dimensionless_coordinate": gate(True, [pointer(HIER_EW_CAPACITY, "/definitions/exact_log_capacity", ew["definitions"]["exact_log_capacity"])]),
        "target_independent_candidate_selection": gate(False, [pointer(HIER_EW_CAPACITY, "/source_values/P_star_branch_locator", ew["source_values"]["P_star_branch_locator"])]),
        "same_quantity_constructed": gate(False, [pointer(HIER_EW_CAPACITY, "/branch_selection", ew["branch_selection"])]),
        "source_return_map_complete": gate(False, [pointer(HIER_EW_CAPACITY, "/claim_boundary/scope", ew["claim_boundary"]["scope"])]),
        "existence_certified": gate(ew["accepted"] is True, [pointer(HIER_EW_CAPACITY, "/accepted", ew["accepted"])]),
        "uniqueness_certified": gate(contraction["banach_unique_fixed_point"] is True, [pointer(HIER_EW_CAPACITY, "/contraction_certificate/banach_unique_fixed_point", contraction["banach_unique_fixed_point"])]),
        "stability_certified": gate(Decimal(contraction["lipschitz_constant"]) < 1, [pointer(HIER_EW_CAPACITY, "/contraction_certificate/lipschitz_constant", contraction["lipschitz_constant"])]),
    }
    result.append(row(candidate_id="hierarchy.n_ew_bridge_defined_capacity", registry_source=HIER_EW_CAPACITY, row_kind="candidate_map", quantity="N_EW_BRIDGE_DEFINED", family="HIERARCHY_PACKET", mathematical_status=ew["status"], disposition="TARGET_EXPOSED_BRIDGE_DEFINED_FIXED_POINT__NOT_INDEPENDENT_CAPACITY", exposure_class="target_exposed_public_endpoint_branch_locator", solver_target_payload_read=True, exact_fixed_point_progress=True, gates=ew_gates, details={"certificate_id": ew["certificate_id"], "target_relation": ew["target_relation"], "lambda": ew["source_values"]["lambda"], "fixed_point": ew["exact_capacity_fixed_point"]["N_CRC_EW"], "rounded_capacity_diagnostic": ew["rounded_capacity_diagnostic"]}))
    return result


def expected_completion(contract: dict[str, Any], forecast: dict[str, Any], menu: dict[str, Any]) -> dict[str, Any]:
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

    gates = {
        "input_catalogue_frozen": declared_gate("input_catalogue_frozen", [pointer(FORECAST_STATE, "/inventory_completeness/status", forecast["inventory_completeness"]["status"]), pointer(FORECAST_STATE, "/controls/candidate_inventory_exhaustiveness", forecast["controls"]["candidate_inventory_exhaustiveness"])]),
        "nuisance_directions_frozen": declared_gate("nuisance_directions_frozen"),
        "units_and_schemes_frozen": declared_gate("units_and_schemes_frozen"),
        "candidate_grammar_frozen": declared_gate("candidate_grammar_frozen", [pointer(FORECAST_STATE, "/selection_rule/status", forecast["selection_rule"]["status"]), pointer(SEMANTICS_MENU, "/campaign_verdict/semantics_enumeration", menu["campaign_verdict"]["semantics_enumeration"]), pointer(CONTRACT, "/scope_boundary/future_source_laws_exhausted", contract["scope_boundary"]["future_source_laws_exhausted"])]),
        "decision_rule_frozen": declared_gate("decision_rule_frozen", [pointer(CONTRACT, "/decision_rule/frozen", contract["decision_rule"]["frozen"])]),
    }
    complete = all(gates[name]["attained"] is True for name in COMPLETION_GATES)
    return {"required_gates": list(COMPLETION_GATES), "gates": gates, "all_complete": complete, "candidate_present_forbidden_while_incomplete": True}


def finalize(rows: list[dict[str, Any]], completion_complete: bool) -> list[str]:
    qualifying: list[str] = []
    for candidate in rows:
        expected_eligible = bool(
            candidate["row_kind"] in QUALIFYING_ROW_KINDS
            and candidate["solver_target_payload_read"] is False
            and candidate["exposure_class"] in QUALIFYING_EXPOSURES
            and candidate["row_gate_complete"] is True
            and completion_complete
        )
        candidate["qualifies_source_only"] = expected_eligible
        if expected_eligible:
            qualifying.append(candidate["candidate_id"])
    return qualifying


def derived_status(rows: list[dict[str, Any]], completion_complete: bool) -> str:
    gate_complete = any(
        candidate["row_kind"] in QUALIFYING_ROW_KINDS
        and candidate["solver_target_payload_read"] is False
        and candidate["exposure_class"] in QUALIFYING_EXPOSURES
        and candidate["row_gate_complete"] is True
        for candidate in rows
    )
    if gate_complete and completion_complete:
        return "SOURCE_ONLY_CANDIDATE_PRESENT"
    if gate_complete:
        return "OPEN_GATE_COMPLETE_CANDIDATE_BLOCKED_BY_PREFLIGHT_COMPLETENESS"
    if any(candidate["exact_fixed_point_progress"] is True for candidate in rows):
        return "OPEN_REGISTERED_CANDIDATES_FAIL_SOURCE_ONLY_GATES"
    return "OPEN_NO_FIXED_POINT_EVIDENCE"


def verify(path: Path) -> None:
    raw = path.read_bytes()
    certificate = json.loads(raw)
    check(isinstance(certificate, dict), "certificate root is not an object")
    check(raw == canonical_json_bytes(certificate), "certificate is not canonical")
    check(set(certificate) == {"schema", "artifact", "issue", "status", "scope", "inventory_crosswalk", "preflight_completion", "qualification", "progress", "candidate_rows", "parent_validation", "parent_pins", "certificate_sha256"}, "top-level field inventory drift")
    check(certificate.get("schema") == SCHEMA and certificate.get("artifact") == "oph_source_only_dimensionless_closure_preflight", "certificate identity drift")
    check(certificate.get("issue") == ISSUE, "certificate issue drift")
    verify_self_hash(certificate, "certificate_sha256")

    contract = load(CONTRACT)
    provenance = load(RC_PROVENANCE)
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
    forecast = load(FORECAST_STATE)
    rc = load(RC_OUTPUT)

    verify_contract(contract)
    verify_forecast_state(forecast)
    verify_rc(rc, p_interval, provenance)
    verify_parent_pins(certificate.get("parent_pins"))
    check(certificate.get("parent_validation") == [{"validator": name, "status": "PASS"} for name in EXPECTED_PARENT_VALIDATORS], "parent validation inventory/status drift")

    rows = expected_p_rows(p_interval, p_trunk, exposure)
    menu_rows, menu_crosswalk = expected_menu_rows(menu, cap_k, cap_p, cap_l, coupled)
    rows.extend(menu_rows)
    rows.append(expected_direct_row(direct))
    rows.extend(expected_common_rows(branch, named))
    rows.append(expected_rc_row(rc))
    rows.append(expected_archive_row(archive, cap_k, cap_p, cap_l))
    rows.extend(expected_hierarchy_rows())
    check(tuple(candidate["candidate_id"] for candidate in rows) == ROW_IDS and len(rows) == 30, "independent row inventory drift")

    completion = expected_completion(contract, forecast, menu)
    qualifying = finalize(rows, completion["all_complete"])
    check(certificate.get("candidate_rows") == rows, "candidate row semantics/gates/details drift")
    check(certificate.get("preflight_completion") == completion, "preflight completion drift")
    verify_all_evidence(rows, completion)

    expected_scope = {
        "registered_packets_only": True,
        "future_source_laws_exhausted": False,
        "issue_closure_claimed": False,
        "physical_or_laboratory_attachment_used_as_gate": False,
        "new_observational_comparison_opened": False,
        "prediction_or_postdiction_promoted": False,
    }
    check(certificate.get("scope") == expected_scope, "scope/promotion drift")
    check(certificate["scope"]["new_observational_comparison_opened"] is False, "H0 opened an observational comparison")

    hierarchy_ids = [candidate["candidate_id"] for candidate in rows if candidate["family"] == "HIERARCHY_PACKET"]
    expected_crosswalk = {
        **menu_crosswalk,
        "archive": {
            "executed_row_count": archive["total_rows"],
            "covered_families": ["CAP-K", "CAP-P", "CAP-L"],
            "does_not_cover_menu_rows": ["capB.bridge_constant", "coupled.cp1_cp2_cp3"],
            "is_not_an_independent_candidate": True,
            "is_not_a_source_return_map": True,
        },
        "hierarchy_packet_count": 5,
        "hierarchy_packet_ids": hierarchy_ids,
        "derived_inventory_row_count": 30,
    }
    check(certificate.get("inventory_crosswalk") == expected_crosswalk, "inventory crosswalk drift")
    check(certificate["inventory_crosswalk"]["archive"]["is_not_a_source_return_map"] is True, "archive promoted to a return map")
    check(certificate["candidate_rows"][24]["details"]["omitted_registered_menu_rows"] == ["capB.bridge_constant", "coupled.cp1_cp2_cp3"], "archive omission disclosure drift")

    expected_qualification = {
        "required_source_only_gates": list(REQUIRED_GATES),
        "qualifying_exposure_classes": sorted(QUALIFYING_EXPOSURES),
        "explicitly_downstream_gates": list(DOWNSTREAM_GATES),
        "qualifying_candidate_ids": qualifying,
        "qualifying_candidate_count": len(qualifying),
    }
    check(certificate.get("qualification") == expected_qualification, "qualification drift")
    check(not completion["all_complete"] and not qualifying, "candidate present while preflight is incomplete")

    strongest = max(candidate["attained_gate_count"] for candidate in rows)
    expected_progress = {
        "registered_row_count": 30,
        "authoritative_capacity_menu_row_count": 17,
        "rows_with_exact_fixed_point_progress": sum(candidate["exact_fixed_point_progress"] is True for candidate in rows),
        "strongest_attained_gate_count": strongest,
        "strongest_progress_candidate_ids": [candidate["candidate_id"] for candidate in rows if candidate["attained_gate_count"] == strongest],
        "status_is_evidence_derived": True,
        "direct_D_equals_M0_control_counted_as_fixed_point_progress": False,
    }
    check(certificate.get("progress") == expected_progress, "progress summary drift")
    status = derived_status(rows, completion["all_complete"])
    check(certificate.get("status") == status, "derived status drift")
    check(not (certificate.get("status") == "SOURCE_ONLY_CANDIDATE_PRESENT" and not completion["all_complete"]), "candidate-present status bypasses completion")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--certificate", type=Path, default=CERTIFICATE)
    args = parser.parse_args()
    try:
        verify(args.certificate.resolve())
    except (KeyError, IndexError, TypeError, ValueError, OSError, json.JSONDecodeError) as error:
        print(f"SOURCE_ONLY_CLOSURE_PREFLIGHT_INVALID: {error}")
        return 1
    print("SOURCE_ONLY_CLOSURE_PREFLIGHT_INDEPENDENT_VALID")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
