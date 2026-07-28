#!/usr/bin/env python3
"""Independent verifier for the bounded #630 Higgs/Yukawa frontier.

This checker imports no producer code.  It validates both schemas, resolves
and re-hashes every parent, repeats the target firewall, reclassifies the
conditional operator space, verifies all completion witnesses, and rejects
coefficient or physical-promotion mutations.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any, Iterable

from jsonschema import Draft202012Validator


PACKAGE = Path(__file__).resolve().parent
REPO_ROOT = PACKAGE.parents[3]
DEFAULT_INPUT = PACKAGE / "outputs" / "higgs_yukawa_source_frontier.json"
SCHEMA_PATH = PACKAGE / "schemas" / "higgs_yukawa_source_frontier_v1.schema.json"
POLICY_PATH = PACKAGE / "data" / "higgs_yukawa_source_policy_v1.json"
POLICY_SCHEMA_PATH = PACKAGE / "schemas" / "higgs_yukawa_source_policy_v1.schema.json"

EXPECTED_ROLES = {
    "typed_realized_event_domain",
    "finite_current_algebra",
    "finite_chiral_matter",
    "finite_global_form",
    "conditional_family_attachment",
    "scalar_multiplicity_boundary",
    "scalar_chain_boundary",
    "higgs_top_fiber_boundary",
}

ALWAYS_FORBIDDEN_OUTPUT_KEYS = {
    "coefficient_value",
    "coefficient_values",
    "selected_value",
    "selected_values",
    "matrix_entries",
    "physical_mass",
    "physical_masses",
    "prediction",
}


def fail(code: str, message: str) -> None:
    raise SystemExit(f"{code}: {message}")


def require(condition: bool, code: str, message: str) -> None:
    if not condition:
        fail(code, message)


def canonical_json_bytes(value: Any) -> bytes:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("utf-8")


def canonical_sha256(value: Any) -> str:
    return hashlib.sha256(canonical_json_bytes(value)).hexdigest()


def walk_keys(value: Any) -> Iterable[str]:
    if isinstance(value, dict):
        for key, child in value.items():
            yield str(key)
            yield from walk_keys(child)
    elif isinstance(value, list):
        for child in value:
            yield from walk_keys(child)


def validate_schema(payload: dict[str, Any], schema_path: Path, code: str) -> None:
    schema = json.loads(schema_path.read_text(encoding="utf-8"))
    Draft202012Validator.check_schema(schema)
    errors = sorted(
        Draft202012Validator(schema).iter_errors(payload),
        key=lambda error: [str(part) for part in error.path],
    )
    require(not errors, code, errors[0].message if errors else "schema validation failed")


def resolve_pin(pin: dict[str, Any]) -> Any:
    relative = pin["path"]
    path = (REPO_ROOT / relative).resolve()
    require(
        path.is_relative_to(REPO_ROOT.resolve()),
        "PATH_TRAVERSAL",
        f"path escapes repository: {relative}",
    )
    require(path.is_file(), "PIN_MISSING", f"pinned file missing: {relative}")
    raw = path.read_bytes()
    require(len(raw) == pin["bytes"], "PIN_SIZE", f"byte size changed: {relative}")
    require(
        hashlib.sha256(raw).hexdigest() == pin["byte_sha256"],
        "PIN_HASH",
        f"byte hash changed: {relative}",
    )
    if pin["kind"] == "text":
        return raw.decode("utf-8")
    parsed = json.loads(raw)
    require(
        canonical_sha256(parsed) == pin["canonical_json_sha256"],
        "PIN_CANONICAL_HASH",
        f"canonical JSON hash changed: {relative}",
    )
    return parsed


def check_policy_and_sources(
    payload: dict[str, Any],
) -> tuple[dict[str, Any], dict[str, Any]]:
    policy_pin = payload["source_inputs"]["policy"]
    require(
        policy_pin["path"] == POLICY_PATH.relative_to(REPO_ROOT).as_posix(),
        "POLICY_PATH",
        "frontier does not pin the canonical policy",
    )
    policy = resolve_pin(policy_pin)
    validate_schema(policy, POLICY_SCHEMA_PATH, "POLICY_SCHEMA")
    require(
        policy == json.loads(POLICY_PATH.read_text(encoding="utf-8")),
        "POLICY_DRIFT",
        "resolved policy differs from canonical policy",
    )
    policy_roles = [row["role"] for row in policy["source_inputs"]]
    require(len(policy_roles) == len(set(policy_roles)), "POLICY_ROLE_DUPLICATE", "policy roles duplicate")
    require(set(policy_roles) == EXPECTED_ROLES, "POLICY_ROLES", "policy role set changed")
    policy_by_role = {row["role"]: row for row in policy["source_inputs"]}
    topology = policy["issue_topology"]
    require(
        topology["semantic_dependencies"] == [314, 566, 567, 569, 634],
        "POLICY_DEPENDENCIES",
        "semantic dependency set changed",
    )
    require(
        topology["open_blocking_dependencies"] == [569, 634]
        and 503 not in topology["open_blocking_dependencies"],
        "POLICY_OPEN_DEPENDENCIES",
        "#503 was promoted from partial-receipt ancestry to a closure dependency",
    )
    require(
        topology["non_gating_partial_receipt_ancestry"]
        == [
            {
                "issue": 503,
                "path": "code/geometry/runs/realized_event_receipt_report.json",
                "closure_required": False,
                "scope": "completed finite E1/E2/E4 screen-sheet receipts only; E3 bulk depth remains open and is not used",
            }
        ],
        "POLICY_NON_GATING_ANCESTRY",
        "#503 partial-receipt ancestry changed",
    )
    require(
        payload["issue_topology"] == topology,
        "ISSUE_TOPOLOGY",
        "artifact issue topology differs from the source policy",
    )

    rows = payload["source_inputs"]["artifacts"]
    require({row["role"] for row in rows} == EXPECTED_ROLES, "SOURCE_ROLES", "source role set changed")
    forbidden_paths = set(policy["target_firewall"]["forbidden_source_paths"])
    allowed_roots = tuple(policy["target_firewall"]["allowed_source_roots"])
    resolved: dict[str, Any] = {}
    for row in rows:
        role = row["role"]
        spec = policy_by_role[role]
        require(row["path"] not in forbidden_paths, "TARGET_PATH", f"forbidden path entered: {row['path']}")
        require(
            any(row["path"] == root or row["path"].startswith(root + "/") for root in allowed_roots),
            "PATH_ALLOWLIST",
            f"path outside allowlist: {row['path']}",
        )
        require(row["path"] == spec["path"], "SOURCE_PATH", f"path changed for {role}")
        require(row["kind"] == spec["kind"], "SOURCE_KIND", f"kind changed for {role}")
        require(row["status"] == spec["status"], "SOURCE_STATUS", f"status changed for {role}")
        resolved[role] = resolve_pin(row)

    expected_bundle = canonical_sha256(
        {
            "policy": policy_pin,
            "artifacts": rows,
        }
    )
    require(
        payload["source_inputs"]["input_bundle_digest"] == expected_bundle,
        "INPUT_DIGEST",
        "input bundle digest mismatch",
    )

    forbidden_keys = {key.casefold() for key in policy["target_firewall"]["forbidden_structured_keys"]}
    hits: list[str] = []
    for role, source in resolved.items():
        if not isinstance(source, (dict, list)):
            continue
        for key in walk_keys(source):
            if key.casefold() in forbidden_keys:
                hits.append(f"{role}:{key}")
    require(not hits, "TARGET_KEY", f"forbidden target keys entered source ancestry: {hits}")
    firewall = payload["target_firewall"]
    require(firewall["target_paths_present"] is False, "TARGET_FLAG", "target path flag changed")
    require(
        firewall["forbidden_structured_keys_present"] is False,
        "TARGET_KEY_FLAG",
        "target key flag changed",
    )
    require(
        firewall["physical_reference_store_consumed"] is False,
        "REFERENCE_STORE",
        "physical reference store entered source ancestry",
    )
    require(
        firewall["external_593_packet_consumed"] is False,
        "EXTERNAL_593",
        "#593 calculation packet entered source ancestry",
    )
    return policy, resolved


def check_parent_semantics(resolved: dict[str, Any]) -> tuple[int, list[list[str]]]:
    event = resolved["typed_realized_event_domain"]
    require(event.get("issue") == 503, "EVENT_ISSUE", "wrong realized-event issue")
    receipts = event.get("receipts_witnessed", {})
    require(
        receipts.get("e1_screen_population") is True
        and receipts.get("e2_certified_separation") is True
        and receipts.get("e4_moebius_cocycle") is True,
        "EVENT_RECEIPTS",
        "required finite event receipts changed",
    )
    require(receipts.get("e3_rank_four_bulk_depth") is False, "EVENT_PROMOTION", "bulk-depth gate changed")

    current = resolved["finite_current_algebra"]
    require(
        current.get("issue") == 566
        and current.get("schema") == "oph.port_current_inner_receipt.v5",
        "CURRENT_PARENT",
        "finite current parent changed",
    )
    require(
        current.get("source_definedness", {}).get("algebraic_construction_verified") is True,
        "CURRENT_SOURCE",
        "finite current source gate changed",
    )

    matter = resolved["finite_chiral_matter"]
    require(
        matter.get("issue") == 314
        and matter.get("schema") == "oph.super_tannakian_matter_receipt.v5",
        "MATTER_PARENT",
        "finite matter parent changed",
    )
    scalar = matter["scalar_and_channel_selection"]
    require(scalar["admissible_scalar_charges"] == [3, -3], "SCALAR_CHARGE", "scalar charge pair changed")
    require("not derived" in scalar["scalar_content_status"], "SCALAR_PROMOTION", "scalar status changed")
    channels = sorted(tuple(row) for row in scalar["derived_channels_for_declared_representative"])
    expected = sorted(
        [
            ("Q", "S", "u_c"),
            ("Q", "Sbar", "d_c"),
            ("L", "Sbar", "e_c"),
        ]
    )
    require(channels == expected, "CHANNELS", "invariant channel set changed")
    require(matter["yukawa_sector"]["invariant_sector_dimension"] == 3, "CHANNEL_DIM", "channel dimension changed")

    global_form = resolved["finite_global_form"]
    require(
        global_form.get("schema") == "oph.axis_center_descent_manifest.v4"
        and global_form.get("hypercharge_convention") == "q = 6Y",
        "GLOBAL_PARENT",
        "global-form parent changed",
    )

    family = resolved["conditional_family_attachment"]
    require(family.get("issue") == 569, "FAMILY_ISSUE", "wrong family issue")
    require(
        family.get("selection", {}).get("minimizer") == "3"
        and family.get("selection", {}).get("strict") is True,
        "FAMILY_SELECTION",
        "rank-three selection changed",
    )
    require(
        family.get("named_interface", {}).get("class") == "conditional_open_interface",
        "FAMILY_PROMOTION",
        "family physical interface changed",
    )
    require(family.get("generation", {}).get("weyl_state_count") == 15, "FAMILY_RANK", "family rank changed")

    multiplicity = resolved["scalar_multiplicity_boundary"]
    require(multiplicity.get("issues") == [616, 617], "MULTIPLICITY_ISSUES", "boundary issue set changed")
    scalar_boundary = multiplicity["scalar_response_multiplicity"]
    require(
        scalar_boundary["verdict"]["scalar_existence"] == "not_source_determined",
        "SCALAR_BOUNDARY",
        "scalar boundary changed",
    )
    require(
        [row["name"] for row in scalar_boundary["countermodel_battery"]["configurations"]]
        == ["n0_no_scalar", "n2_duplicate_identical_charge", "n2_one_inert"],
        "SCALAR_COUNTERMODELS",
        "scalar countermodel set changed",
    )

    chain = resolved["scalar_chain_boundary"]
    require(
        chain.get("issue") == 623
        and chain.get("bounded_exit") == "limited_subchain_audit_physical_interface_open",
        "CHAIN_BOUNDARY",
        "limited scalar-chain boundary changed",
    )
    require(
        chain.get("discriminating_interface", {}).get("interface_class")
        == "conditional_open_interface",
        "CHAIN_PROMOTION",
        "physical scalar interface changed",
    )

    no_go = resolved["higgs_top_fiber_boundary"]
    for line in (
        "linear_a2 = u * u",
        "born_a2 = u",
        "expected_higgs = 2 * u * (1 - u)",
        "expected_top = u * (1 - u) / 2",
    ):
        require(line in no_go, "HIGGS_TOP_BOUNDARY", f"missing exact boundary line: {line}")
    return int(family["selection"]["minimizer"]), [list(row) for row in expected]


def check_coefficient_space(
    payload: dict[str, Any],
    family_dimension: int,
    expected_channels: list[list[str]],
) -> None:
    space = payload["conditional_coefficient_space"]
    require(
        space["status"] == "exact_conditional_operator_basis__all_physical_coefficients_open",
        "CLASSIFICATION_STATUS",
        "coefficient classification status changed",
    )
    require(family_dimension == 3, "CLASSIFICATION_FAMILY", "family dimension changed")
    require(space["scalar_kinetic"]["real_basis_dimension"] == 1, "KINETIC_DIM", "kinetic basis changed")
    require(space["scalar_kinetic"]["source_selection"] == "not_emitted", "KINETIC_SELECTION", "kinetic coefficient promoted")
    require(space["scalar_potential"]["real_basis_dimension"] == 2, "POTENTIAL_DIM", "potential basis changed")
    require(
        space["scalar_potential"]["operator_basis"] == ["S^dagger S", "(S^dagger S)^2"],
        "POTENTIAL_BASIS",
        "potential operator basis changed",
    )
    require(space["scalar_potential"]["source_selection"] == "not_emitted", "POTENTIAL_SELECTION", "potential promoted")
    require(
        space["scalar_potential"]["stable_branch_selection"] == "not_emitted",
        "VACUUM_SELECTION",
        "stable branch promoted",
    )
    yukawa = space["yukawa_sector"]
    require(yukawa["channel_types"] == expected_channels, "YUKAWA_CHANNELS", "channel ordering or content changed")
    per_channel = family_dimension * family_dimension
    total = len(expected_channels) * per_channel
    require(per_channel == 9 and total == 27, "YUKAWA_ARITHMETIC", "classification arithmetic changed")
    require(yukawa["complex_dimension_per_channel"] == per_channel, "YUKAWA_MATRIX_DIM", "matrix dimension mismatch")
    require(yukawa["total_complex_dimension"] == total, "YUKAWA_TOTAL_DIM", "total dimension mismatch")
    require(yukawa["source_selection"] == "not_emitted", "YUKAWA_SELECTION", "Yukawa matrices promoted")
    require(
        space["physical_coefficient_assignments_emitted"] is False,
        "COEFFICIENT_PROMOTION",
        "coefficient assignment flag changed",
    )


def check_witnesses(
    payload: dict[str, Any],
    policy: dict[str, Any],
) -> None:
    witnesses = payload["nonidentifiability_witnesses"]
    expected_parent_digest = canonical_sha256(
        {
            "source_pins": payload["source_inputs"]["artifacts"],
            "conditional_premises": policy["conditional_classification_premises"],
        }
    )
    require(
        {row["parent_projection_digest"] for row in witnesses.values()} == {expected_parent_digest},
        "PARENT_PROJECTION",
        "completion witnesses do not share the exact finite parent projection",
    )
    for name, row in witnesses.items():
        require(row["physical_promotion_allowed"] is False, "WITNESS_PROMOTION", f"{name} was promoted")
        require(
            row["coefficient_assignments_emitted"] is False,
            "WITNESS_COEFFICIENT",
            f"{name} emits coefficient assignments",
        )

    scalar = witnesses["scalar_content"]
    require(scalar["completion_A"] != scalar["completion_B"], "SCALAR_WITNESS", "scalar completions collapsed")
    require(scalar["common_finite_parent"] is True, "SCALAR_PARENT", "scalar parents differ")
    require(
        scalar["exact_relation"] == "rank_scalar(A) != rank_scalar(B)",
        "SCALAR_RELATION",
        "scalar separation relation changed",
    )

    yukawa = witnesses["yukawa_matrices"]
    require(yukawa["completion_A"] != yukawa["completion_B"], "YUKAWA_WITNESS", "Yukawa completions collapsed")
    require(yukawa["witness_domain"] == "Mat_3(C)^3", "YUKAWA_DOMAIN", "Yukawa domain changed")
    require(yukawa["exact_relation"] == "Y_A != Y_B", "YUKAWA_RELATION", "Yukawa separation changed")

    fj = witnesses["v_chart_to_v_F"]
    require(fj["completion_A"] != fj["completion_B"], "FJ_WITNESS", "coordinate completions collapsed")
    require(
        fj["witness_domain"] == "invertible_local_coordinate_one_jets",
        "FJ_DOMAIN",
        "coordinate witness domain changed",
    )
    require(fj["exact_relation"] == "jet(Phi_A) != jet(Phi_B)", "FJ_RELATION", "coordinate relation changed")

    issue_521 = witnesses["issue_521_exact_separation"]
    require(issue_521["source_coordinate_domain"] == "0 < u < 1", "SEPARATION_DOMAIN", "#521 domain changed")
    require(
        issue_521["squared_lift_relations"]
        == {
            "linear": "a_linear^2 = u^2",
            "Born": "a_Born^2 = u",
        },
        "SEPARATION_LIFTS",
        "#521 lift relations changed",
    )
    require(
        issue_521["exact_separation"]
        == {
            "Higgs": "r_H(linear)-r_H(Born)=2u(1-u)>0",
            "top": "r_t(Born)-r_t(linear)=u(1-u)/2>0",
        },
        "SEPARATION",
        "#521 exact separation changed",
    )
    # Independent exact polynomial reconstruction in ascending powers of u.
    linear_higgs = (Fraction(2), Fraction(0), Fraction(-2))
    born_higgs = (Fraction(2), Fraction(-2), Fraction(0))
    linear_top = (Fraction(0), Fraction(0), Fraction(1, 2))
    born_top = (Fraction(0), Fraction(1, 2), Fraction(0))
    require(
        tuple(a - b for a, b in zip(linear_higgs, born_higgs))
        == (Fraction(0), Fraction(2), Fraction(-2)),
        "SEPARATION_HIGGS_ALGEBRA",
        "#521 Higgs polynomial reconstruction failed",
    )
    require(
        tuple(a - b for a, b in zip(born_top, linear_top))
        == (Fraction(0), Fraction(1, 2), Fraction(-1, 2)),
        "SEPARATION_TOP_ALGEBRA",
        "#521 top polynomial reconstruction failed",
    )


def check(payload: dict[str, Any]) -> None:
    output_hits = sorted(
        {key for key in walk_keys(payload) if key.casefold() in ALWAYS_FORBIDDEN_OUTPUT_KEYS}
    )
    require(not output_hits, "COEFFICIENT_VALUE", f"forbidden coefficient/value keys found: {output_hits}")
    require(payload.get("promotion_allowed") is False, "PROMOTION", "frontier cannot promote")
    require(
        payload.get("physical_source_action_emitted") is False,
        "SOURCE_ACTION_PROMOTION",
        "positive source action is not emitted",
    )
    require(
        payload.get("coefficient_assignments_emitted") is False,
        "COEFFICIENT_PROMOTION",
        "coefficient assignments are not emitted",
    )
    witnesses = payload.get("nonidentifiability_witnesses", {})
    scalar = witnesses.get("scalar_content", {})
    require(
        scalar.get("completion_A") != scalar.get("completion_B"),
        "SCALAR_WITNESS",
        "scalar completions collapsed",
    )
    yukawa = witnesses.get("yukawa_matrices", {})
    require(
        yukawa.get("completion_A") != yukawa.get("completion_B"),
        "YUKAWA_WITNESS",
        "Yukawa completions collapsed",
    )
    fj = witnesses.get("v_chart_to_v_F", {})
    require(
        fj.get("completion_A") != fj.get("completion_B"),
        "FJ_WITNESS",
        "coordinate completions collapsed",
    )
    issue_521 = witnesses.get("issue_521_exact_separation", {})
    require(
        issue_521.get("exact_separation")
        == {
            "Higgs": "r_H(linear)-r_H(Born)=2u(1-u)>0",
            "top": "r_t(Born)-r_t(linear)=u(1-u)/2>0",
        },
        "SEPARATION",
        "#521 exact separation changed",
    )
    validate_schema(payload, SCHEMA_PATH, "SCHEMA")
    require(payload["issue"] == 630, "ISSUE", "wrong issue")
    require(
        payload["status"] == "BOUNDED_NONPROMOTING_FRONTIER__POSITIVE_SOURCE_ACTION_OPEN",
        "STATUS",
        "frontier status changed",
    )

    policy, resolved = check_policy_and_sources(payload)
    family_dimension, expected_channels = check_parent_semantics(resolved)
    check_coefficient_space(payload, family_dimension, expected_channels)
    check_witnesses(payload, policy)

    for name, row in payload["positive_source_objects"].items():
        require(row["status"] == "not_emitted", "POSITIVE_OBJECT", f"{name} was emitted")
        require(bool(row["reason"]), "POSITIVE_REASON", f"{name} lacks an explicit blocker")

    statuses = [row["status"] for row in payload["acceptance_map"]]
    require(statuses.count("open") == 6, "ACCEPTANCE_OPEN", "six positive acceptance criteria must remain open")
    require(statuses.count("frontier_only") == 1, "ACCEPTANCE_FRONTIER", "frontier replay row changed")
    require(
        {row["issue"] for row in payload["consumer_effect"]} == {32, 34, 594},
        "CONSUMERS",
        "consumer set changed",
    )
    require(
        "closure of issue 630" in payload["strictly_not_claimed"],
        "NONCLAIM",
        "issue closure nonclaim is absent",
    )

    subject = dict(payload)
    recorded_digest = subject.pop("subject_digest")
    require(
        recorded_digest == canonical_sha256(subject),
        "SUBJECT_DIGEST",
        "subject digest mismatch",
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    args = parser.parse_args()
    payload = json.loads(args.input.read_text(encoding="utf-8"))
    check(payload)
    print("higgs_yukawa_source_frontier: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
