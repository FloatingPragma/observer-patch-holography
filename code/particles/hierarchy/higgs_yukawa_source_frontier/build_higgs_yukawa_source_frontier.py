#!/usr/bin/env python3
"""Build the exact, non-promoting Higgs/Yukawa source frontier for #630.

The producer resolves a closed list of finite positive and boundary parents.
It classifies the coefficient space of a conditional one-doublet,
three-family action and records exact completion-fiber witnesses.  It emits no
scalar coefficient, Yukawa entry, vacuum coordinate, FJ map, mass, or physical
promotion.
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
POLICY_PATH = PACKAGE / "data" / "higgs_yukawa_source_policy_v1.json"
POLICY_SCHEMA_PATH = PACKAGE / "schemas" / "higgs_yukawa_source_policy_v1.schema.json"
DEFAULT_OUTPUT = PACKAGE / "outputs" / "higgs_yukawa_source_frontier.json"
SCHEMA = "oph.higgs_yukawa_source_frontier.v1"

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


class FrontierError(ValueError):
    """Fail-closed producer error with a stable machine-readable code."""

    def __init__(self, code: str, message: str) -> None:
        super().__init__(f"{code}: {message}")
        self.code = code


def require(condition: bool, code: str, message: str) -> None:
    if not condition:
        raise FrontierError(code, message)


def canonical_json_bytes(value: Any) -> bytes:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("utf-8")


def canonical_sha256(value: Any) -> str:
    return hashlib.sha256(canonical_json_bytes(value)).hexdigest()


def _keys(value: Any) -> Iterable[str]:
    if isinstance(value, dict):
        for key, child in value.items():
            yield str(key)
            yield from _keys(child)
    elif isinstance(value, list):
        for child in value:
            yield from _keys(child)


def validate_policy(policy: dict[str, Any]) -> None:
    schema = json.loads(POLICY_SCHEMA_PATH.read_text(encoding="utf-8"))
    Draft202012Validator.check_schema(schema)
    errors = sorted(
        Draft202012Validator(schema).iter_errors(policy),
        key=lambda error: [str(part) for part in error.path],
    )
    require(not errors, "POLICY_SCHEMA", errors[0].message if errors else "invalid policy")
    roles = [row["role"] for row in policy["source_inputs"]]
    require(len(roles) == len(set(roles)), "POLICY_ROLE_DUPLICATE", "source roles must be unique")
    require(set(roles) == EXPECTED_ROLES, "POLICY_ROLES", "source role set changed")
    require(
        {row.get("issue") for row in policy["consumers"]} == {32, 34, 594},
        "POLICY_CONSUMERS",
        "consumer issue set changed",
    )
    topology = policy["issue_topology"]
    require(
        topology["semantic_dependencies"] == [314, 566, 567, 569, 634]
        and topology["open_blocking_dependencies"] == [636, 637, 638],
        "POLICY_DEPENDENCIES",
        "semantic or open blocking dependency set changed",
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
        "#503 must remain completed partial-receipt ancestry rather than a closure dependency",
    )


def _repo_path(relative: str) -> Path:
    path = (REPO_ROOT / relative).resolve()
    require(
        path.is_relative_to(REPO_ROOT.resolve()),
        "PATH_TRAVERSAL",
        f"path escapes repository: {relative}",
    )
    return path


def json_pin(
    path: Path,
    *,
    role: str | None = None,
    status: str | None = None,
) -> dict[str, Any]:
    raw = path.read_bytes()
    parsed = json.loads(raw)
    pin: dict[str, Any] = {
        "path": path.relative_to(REPO_ROOT).as_posix(),
        "kind": "json",
        "bytes": len(raw),
        "byte_sha256": hashlib.sha256(raw).hexdigest(),
        "canonical_json_sha256": canonical_sha256(parsed),
    }
    if role is not None:
        pin["role"] = role
    if status is not None:
        pin["status"] = status
    return pin


def text_pin(path: Path, *, role: str, status: str) -> dict[str, Any]:
    raw = path.read_bytes()
    return {
        "path": path.relative_to(REPO_ROOT).as_posix(),
        "kind": "text",
        "bytes": len(raw),
        "byte_sha256": hashlib.sha256(raw).hexdigest(),
        "role": role,
        "status": status,
    }


def resolve_sources(
    policy: dict[str, Any],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    allowed_roots = tuple(policy["target_firewall"]["allowed_source_roots"])
    forbidden_paths = set(policy["target_firewall"]["forbidden_source_paths"])
    rows: list[dict[str, Any]] = []
    payloads: dict[str, Any] = {}
    for entry in policy["source_inputs"]:
        relative = entry["path"]
        require(relative not in forbidden_paths, "TARGET_PATH", f"forbidden input path: {relative}")
        require(
            any(relative == root or relative.startswith(root + "/") for root in allowed_roots),
            "PATH_ALLOWLIST",
            f"path lies outside the closed allowlist: {relative}",
        )
        path = _repo_path(relative)
        require(path.is_file(), "SOURCE_MISSING", f"source input is missing: {relative}")
        if entry["kind"] == "json":
            payload = json.loads(path.read_text(encoding="utf-8"))
            row = json_pin(path, role=entry["role"], status=entry["status"])
        else:
            payload = path.read_text(encoding="utf-8")
            row = text_pin(path, role=entry["role"], status=entry["status"])
        rows.append(row)
        payloads[entry["role"]] = payload
    require(len(rows) == 8, "SOURCE_COUNT", "exactly eight finite source and boundary parents are required")
    return rows, payloads


def target_firewall(payloads: dict[str, Any], policy: dict[str, Any]) -> None:
    forbidden = {key.casefold() for key in policy["target_firewall"]["forbidden_structured_keys"]}
    hits: list[str] = []
    for role, payload in payloads.items():
        if not isinstance(payload, (dict, list)):
            continue
        for key in _keys(payload):
            if key.casefold() in forbidden:
                hits.append(f"{role}:{key}")
    require(not hits, "TARGET_KEY", f"forbidden structured target keys entered source ancestry: {hits}")


def validate_parent_semantics(payloads: dict[str, Any]) -> dict[str, Any]:
    event = payloads["typed_realized_event_domain"]
    require(event.get("issue") == 503, "EVENT_ISSUE", "realized event parent must be issue #503")
    require(
        event.get("receipts_witnessed", {}).get("e1_screen_population") is True
        and event.get("receipts_witnessed", {}).get("e2_certified_separation") is True
        and event.get("receipts_witnessed", {}).get("e4_moebius_cocycle") is True,
        "EVENT_RECEIPTS",
        "required finite event receipts are absent",
    )
    require(
        event.get("receipts_witnessed", {}).get("e3_rank_four_bulk_depth") is False,
        "EVENT_OVERPROMOTION",
        "the selected parent must retain the open bulk-depth boundary",
    )

    current = payloads["finite_current_algebra"]
    require(current.get("issue") == 566, "CURRENT_ISSUE", "current parent must be issue #566")
    require(
        current.get("schema") == "oph.port_current_inner_receipt.v5",
        "CURRENT_SCHEMA",
        "unexpected current receipt schema",
    )
    require(
        current.get("source_definedness", {}).get("algebraic_construction_verified") is True,
        "CURRENT_SOURCE",
        "finite current construction is not verified",
    )
    require(
        current.get("conditional_algebraic_gate", {}).get("passed") is True
        and current.get("physical_source_gate", {}).get("passed") is False
        and current.get("physical_source_gate", {}).get("response_model_source_bound")
        is False,
        "CURRENT_SCOPE",
        "finite current parent must remain a declared conditional fixture",
    )

    matter = payloads["finite_chiral_matter"]
    require(matter.get("issue") == 314, "MATTER_ISSUE", "matter parent must be issue #314")
    require(
        matter.get("schema") == "oph.super_tannakian_matter_receipt.v5",
        "MATTER_SCHEMA",
        "unexpected matter receipt schema",
    )
    require(
        matter.get("conditional_algebraic_gate", {}).get("passed") is True
        and matter.get("physical_source_gate", {}).get("passed") is False
        and matter.get("physical_source_gate", {}).get("matter_lift_source_bound")
        is False,
        "MATTER_SCOPE",
        "finite matter parent must remain a declared conditional fixture",
    )
    scalar = matter["scalar_and_channel_selection"]
    require(
        scalar["admissible_scalar_charges"] == [3, -3],
        "SCALAR_CHARGE",
        "conditional scalar charge pair changed",
    )
    require(
        "not derived" in scalar["scalar_content_status"],
        "SCALAR_SELECTION",
        "scalar existence must remain underived",
    )
    channels = sorted(tuple(row) for row in scalar["derived_channels_for_declared_representative"])
    expected_channels = sorted(
        [
            ("Q", "S", "u_c"),
            ("Q", "Sbar", "d_c"),
            ("L", "Sbar", "e_c"),
        ]
    )
    require(channels == expected_channels, "YUKAWA_CHANNELS", "conditional invariant channel set changed")
    require(
        matter["yukawa_sector"]["invariant_sector_dimension"] == 3,
        "YUKAWA_DIMENSION",
        "invariant channel-space dimension changed",
    )

    global_form = payloads["finite_global_form"]
    require(
        global_form.get("schema") == "oph.axis_center_descent_manifest.v4",
        "GLOBAL_SCHEMA",
        "unexpected global-form schema",
    )
    require(
        global_form.get("hypercharge_convention") == "q = 6Y",
        "GLOBAL_CHARGE",
        "q=6Y convention changed",
    )
    require(
        "do not select a physical global form"
        in str(global_form.get("description", "")).lower(),
        "GLOBAL_SCOPE",
        "global-form parent must retain its physical nonselection boundary",
    )

    family = payloads["conditional_family_attachment"]
    require(family.get("issue") == 569, "FAMILY_ISSUE", "family parent must be issue #569")
    require(
        family.get("selection", {}).get("minimizer") == "3"
        and family.get("selection", {}).get("strict") is True,
        "FAMILY_SELECTION",
        "rank-three screen selection changed",
    )
    require(
        family.get("named_interface", {}).get("class") == "conditional_open_interface",
        "FAMILY_PROMOTION",
        "family parent must retain the physical attachment interface",
    )
    family_scope = family.get("conditional_structural_scope", {})
    require(
        family_scope.get("physical_current_source_bound") is False
        and family_scope.get("physical_matter_lift_source_bound") is False
        and family_scope.get("physical_global_form_selected") is False,
        "FAMILY_SCOPE",
        "family parent must retain its conditional structural scope",
    )
    require(
        family.get("generation", {}).get("weyl_state_count") == 15,
        "FAMILY_RANK",
        "per-family matter rank changed",
    )

    multiplicity = payloads["scalar_multiplicity_boundary"]
    require(
        multiplicity.get("issues") == [616, 617],
        "MULTIPLICITY_ISSUES",
        "multiplicity boundary issue set changed",
    )
    scalar_boundary = multiplicity["scalar_response_multiplicity"]
    require(
        scalar_boundary["verdict"]["scalar_existence"] == "not_source_determined",
        "SCALAR_BOUNDARY",
        "scalar existence boundary was removed",
    )
    names = [row["name"] for row in scalar_boundary["countermodel_battery"]["configurations"]]
    require(
        names == ["n0_no_scalar", "n2_duplicate_identical_charge", "n2_one_inert"],
        "SCALAR_COUNTERMODELS",
        "registered scalar countermodels changed",
    )

    chain_boundary = payloads["scalar_chain_boundary"]
    require(chain_boundary.get("issue") == 623, "CHAIN_ISSUE", "chain boundary must be issue #623")
    require(
        chain_boundary.get("bounded_exit") == "limited_subchain_audit_physical_interface_open",
        "CHAIN_SCOPE",
        "limited scalar-chain boundary changed",
    )
    require(
        chain_boundary.get("discriminating_interface", {}).get("interface_class")
        == "conditional_open_interface",
        "CHAIN_PROMOTION",
        "physical scalar discriminator interface must remain open",
    )

    no_go = payloads["higgs_top_fiber_boundary"]
    required_lines = (
        "linear_a2 = u * u",
        "born_a2 = u",
        "expected_higgs = 2 * u * (1 - u)",
        "expected_top = u * (1 - u) / 2",
    )
    require(
        all(line in no_go for line in required_lines),
        "HIGGS_TOP_BOUNDARY",
        "#521 exact separation implementation changed",
    )
    return {
        "family_dimension": int(family["selection"]["minimizer"]),
        "conditional_channels": [list(row) for row in expected_channels],
    }


def coefficient_space(context: dict[str, Any]) -> dict[str, Any]:
    family_dimension = context["family_dimension"]
    channel_count = len(context["conditional_channels"])
    per_channel_complex_dimension = family_dimension * family_dimension
    total_complex_dimension = channel_count * per_channel_complex_dimension
    require(family_dimension == 3, "CLASSIFICATION_FAMILY", "conditional family dimension must be three")
    require(channel_count == 3, "CLASSIFICATION_CHANNELS", "three invariant channel types are required")
    require(per_channel_complex_dimension == 9, "CLASSIFICATION_MATRIX", "family matrix dimension changed")
    require(total_complex_dimension == 27, "CLASSIFICATION_TOTAL", "total Yukawa dimension changed")
    return {
        "status": "exact_conditional_operator_basis__all_physical_coefficients_open",
        "classification_scope": "scalar and Yukawa operators for one declared complex weak doublet on the conditional rank-three family branch, modulo field-independent constants and total derivatives",
        "premises": [
            "one_complex_weak_doublet",
            "rank_three_family_space",
            "four_dimensional_power_counting_through_degree_four",
            "positive_local_scalar_kinetic_form",
        ],
        "scalar_kinetic": {
            "operator_basis": ["(D S)^dagger (D S)"],
            "real_basis_dimension": 1,
            "coefficient_domain": "positive real ray",
            "source_selection": "not_emitted",
        },
        "scalar_potential": {
            "operator_basis": ["S^dagger S", "(S^dagger S)^2"],
            "real_basis_dimension": 2,
            "coefficient_domain": "real plane with a stability subcone",
            "source_selection": "not_emitted",
            "stable_branch_selection": "not_emitted",
        },
        "yukawa_sector": {
            "channel_types": context["conditional_channels"],
            "invariant_channel_type_count": channel_count,
            "family_space": "C^3",
            "coefficient_space_per_channel": "Mat_3(C)",
            "complex_dimension_per_channel": per_channel_complex_dimension,
            "complete_coefficient_space": "Mat_3(C)^3",
            "total_complex_dimension": total_complex_dimension,
            "source_selection": "not_emitted",
        },
        "product_space": "R_positive x R^2 x Mat_3(C)^3",
        "normalization_quotient_status": "not_formed_without_a_source_kinetic_and_field_redefinition_rule",
        "flavor_basis_quotient_status": "not_formed_without_source-selected family bases and equivalence maps",
        "physical_coefficient_assignments_emitted": False,
    }


def nonidentifiability_witnesses(
    source_rows: list[dict[str, Any]],
    policy: dict[str, Any],
) -> dict[str, Any]:
    parent_projection = {
        "source_pins": source_rows,
        "conditional_premises": policy["conditional_classification_premises"],
    }
    digest = canonical_sha256(parent_projection)
    common = {
        "parent_projection_digest": digest,
        "physical_promotion_allowed": False,
        "coefficient_assignments_emitted": False,
    }
    return {
        "scalar_content": {
            **common,
            "completion_A": "registered_empty_scalar_extension",
            "completion_B": "declared_one_weak_doublet_extension",
            "common_finite_parent": True,
            "distinguishing_projection": "scalar_carrier_multiplicity",
            "exact_relation": "rank_scalar(A) != rank_scalar(B)",
            "status": "current_finite_parent_not_constant_on_scalar_content_fiber",
            "scope": "grammar and declared-completion level; no physical scalar is selected",
        },
        "yukawa_matrices": {
            **common,
            "completion_A": "formal_family_matrix_point_Y_A",
            "completion_B": "formal_family_matrix_point_Y_B",
            "common_finite_parent": True,
            "witness_domain": "Mat_3(C)^3",
            "exact_relation": "Y_A != Y_B",
            "existence_reason": "the classified complex coefficient space has positive dimension",
            "distinguishing_projection": "family_matrix_coordinate_projection",
            "status": "current_finite_parent_not_constant_on_yukawa_coefficient_fiber",
            "scope": "formal coefficient completion only; no matrix entries are emitted",
        },
        "v_chart_to_v_F": {
            **common,
            "completion_A": "formal_invertible_coordinate_jet_Phi_A",
            "completion_B": "formal_invertible_coordinate_jet_Phi_B",
            "common_finite_parent": True,
            "witness_domain": "invertible_local_coordinate_one_jets",
            "exact_relation": "jet(Phi_A) != jet(Phi_B)",
            "existence_reason": "the finite parents impose no scheme, scale, Jacobian, or finite-shift equation",
            "distinguishing_projection": "coordinate_one_jet",
            "status": "current_finite_parent_not_constant_on_FJ_coordinate_map_fiber",
            "scope": "formal coordinate completion only; no FJ map is emitted",
        },
        "issue_521_exact_separation": {
            **common,
            "source_coordinate_domain": "0 < u < 1",
            "completion_A": "linear_probability_lift",
            "completion_B": "Born_amplitude_lift",
            "squared_lift_relations": {
                "linear": "a_linear^2 = u^2",
                "Born": "a_Born^2 = u",
            },
            "formal_readouts": {
                "Higgs": "r_H = 2(1-a^2)",
                "top": "r_t = a^2/2",
            },
            "exact_separation": {
                "Higgs": "r_H(linear)-r_H(Born)=2u(1-u)>0",
                "top": "r_t(Born)-r_t(linear)=u(1-u)/2>0",
            },
            "status": "exact_current_reduct_nonidentifiability",
            "scope": "dimensionless formal readout coefficients only; no physical source action or pole value",
        },
    }


def not_emitted_objects() -> dict[str, Any]:
    reasons = {
        "physical_scalar_carrier_and_multiplicity": "finite parents retain scalar-content countermodels",
        "canonical_scalar_kinetic_normalization": "no source normalization functional is present",
        "scalar_potential_coefficients": "operator basis is classified but coefficients are not selected",
        "complete_Yu_Yd_Ye_matrices": "channel types do not select family-space matrix entries",
        "stable_vacuum_branch": "potential coefficients and stability selector are absent",
        "v_chart_to_v_F_map": "scheme, scale, Jacobian, and finite-shift data are absent",
        "source_uncertainty_packet": "no positive source action has been emitted",
    }
    return {
        name: {
            "status": "not_emitted",
            "reason": reason,
        }
        for name, reason in reasons.items()
    }


def validate_issue_521_exact_algebra() -> None:
    """Recheck the formal polynomial separation without choosing a value of u."""

    # Coefficients are stored in ascending powers of u.
    linear_higgs = (Fraction(2), Fraction(0), Fraction(-2))
    born_higgs = (Fraction(2), Fraction(-2), Fraction(0))
    linear_top = (Fraction(0), Fraction(0), Fraction(1, 2))
    born_top = (Fraction(0), Fraction(1, 2), Fraction(0))
    higgs_difference = tuple(a - b for a, b in zip(linear_higgs, born_higgs))
    top_difference = tuple(a - b for a, b in zip(born_top, linear_top))
    require(
        higgs_difference == (Fraction(0), Fraction(2), Fraction(-2)),
        "HIGGS_SEPARATION_ALGEBRA",
        "#521 Higgs polynomial identity failed",
    )
    require(
        top_difference == (Fraction(0), Fraction(1, 2), Fraction(-1, 2)),
        "TOP_SEPARATION_ALGEBRA",
        "#521 top polynomial identity failed",
    )


def build_artifact(
    policy: dict[str, Any],
    *,
    policy_path: Path = POLICY_PATH,
) -> dict[str, Any]:
    validate_policy(policy)
    require(policy_path.resolve() == POLICY_PATH.resolve(), "POLICY_PATH", "only the canonical policy may be pinned")
    source_rows, payloads = resolve_sources(policy)
    target_firewall(payloads, policy)
    context = validate_parent_semantics(payloads)
    validate_issue_521_exact_algebra()
    policy_pin = json_pin(policy_path)
    witnesses = nonidentifiability_witnesses(source_rows, policy)
    artifact: dict[str, Any] = {
        "schema": SCHEMA,
        "issue": 630,
        "status": "BOUNDED_NONPROMOTING_FRONTIER__POSITIVE_SOURCE_ACTION_OPEN",
        "promotion_allowed": False,
        "physical_source_action_emitted": False,
        "coefficient_assignments_emitted": False,
        "issue_topology": policy["issue_topology"],
        "source_inputs": {
            "policy": policy_pin,
            "artifacts": source_rows,
            "input_bundle_digest": canonical_sha256(
                {
                    "policy": policy_pin,
                    "artifacts": source_rows,
                }
            ),
        },
        "target_firewall": {
            "status": "closed_allowlist_and_structured_key_scan",
            "target_paths_present": False,
            "forbidden_structured_keys_present": False,
            "physical_reference_store_consumed": False,
            "external_593_packet_consumed": False,
        },
        "conditional_coefficient_space": coefficient_space(context),
        "nonidentifiability_witnesses": witnesses,
        "positive_source_objects": not_emitted_objects(),
        "acceptance_map": [
            {
                "criterion": "all physical coefficients emitted before target mounting",
                "status": "open",
                "evidence": "coefficient assignments are not emitted",
            },
            {
                "criterion": "physical scalar existence and multiplicity selected",
                "status": "open",
                "evidence": "scalar-content witness pair survives the finite parent projection",
            },
            {
                "criterion": "complete Yu Yd Ye matrices emitted",
                "status": "open",
                "evidence": "three channel types classify Mat_3(C)^3 but do not choose a point",
            },
            {
                "criterion": "kinetic normalization potential stability and vacuum selected on one action",
                "status": "open",
                "evidence": "the conditional operator basis has no source coefficient assignment",
            },
            {
                "criterion": "symbolic v_chart to v_F map with scheme scale Jacobian and uncertainty",
                "status": "open",
                "evidence": "two formal coordinate one-jets survive",
            },
            {
                "criterion": "existing scalar and Higgs/top countermodels distinguished by positive source data",
                "status": "open",
                "evidence": "countermodels are retained and explicitly separated",
            },
            {
                "criterion": "independent replay rejects target and promotion mutations",
                "status": "frontier_only",
                "evidence": "independent checker and mutation tests guard this bounded packet",
            },
        ],
        "consumer_effect": [
            {
                "issue": row["issue"],
                "status": row["status_after_this_frontier"],
            }
            for row in policy["consumers"]
        ],
        "strictly_not_claimed": [
            "physical scalar existence or multiplicity",
            "canonical scalar kinetic normalization",
            "scalar potential or stable vacuum",
            "any Yukawa matrix entry",
            "a v_chart to v_F coordinate map",
            "a Higgs top W or Z source value",
            "closure of issue 630",
        ],
    }
    forbidden_output_keys = {
        key.casefold() for key in policy["target_firewall"]["forbidden_output_keys"]
    }
    output_hits = sorted({key for key in _keys(artifact) if key.casefold() in forbidden_output_keys})
    require(not output_hits, "OUTPUT_VALUE_KEY", f"forbidden output keys present: {output_hits}")
    artifact["subject_digest"] = canonical_sha256(artifact)
    return artifact


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    policy = json.loads(POLICY_PATH.read_text(encoding="utf-8"))
    artifact = build_artifact(policy)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
