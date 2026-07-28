#!/usr/bin/env python3
"""Independently resolve the route-neutral issue-#633 source-clock frontier.

This checker imports no producer code. It verifies the closed schema, policy
and diagnostic byte pins, strict dependency semantics, exact rational interval
inversion, provenance DAG, target firewall, and the absence of every physical
output.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from fractions import Fraction
from pathlib import Path, PurePosixPath
from typing import Any, Mapping

from jsonschema import Draft202012Validator


HERE = Path(__file__).resolve().parent
DEFAULT_REPO_ROOT = HERE.parents[3]
DEFAULT_FRONTIER = HERE / "outputs" / "source_clock_frontier.json"
SCHEMA_PATH = HERE / "schemas" / "source_clock_frontier_v1.schema.json"
POLICY_REL = (
    "code/particles/hierarchy/source_clock_frontier/"
    "data/source_clock_policy_v1.json"
)

EXPECTED_SCHEMA = "oph.source_clock.frontier.v1"
EXPECTED_STATUS = "NONPROMOTING_ROUTE_NEUTRAL_CONTRACT__PHYSICAL_CLOCK_ATTACHMENT_OPEN"
EXPECTED_HARD_DEPENDENCIES = [634]
EXPECTED_OPTIONAL_ROUTE_OWNERS = [32, 34, 317, 318, 425, 522, 545, 546, 569, 633]
EXPECTED_DOWNSTREAM_ONLY_ISSUES = [334]
EXPECTED_CANDIDATE_ROUTE_ID = "cesium_133_hyperfine"
EXPECTED_COMPONENT_IDS = [
    "R_U",
    "R_alpha",
    "R_e_abs",
    "R_QCD_nuc_133Cs",
    "R_atom_133Cs",
]
EXPECTED_DIAGNOSTIC_ROLES = [
    "legacy_gravity_checksum_skeleton",
    "synthetic_feshbach_fixture",
]
EXPECTED_GATE_IDS = [
    "route_neutral_dimensionless_clock_observable",
    "admissible_physical_clock_route_and_process_receipt",
    "target_clean_clock_to_si_attachment",
]

H_SI = Fraction(662607015, 10**42)
NU_CS_HZ = Fraction(9192631770, 1)
GEV_J = Fraction(1602176634, 10**19)
C_SI = Fraction(299792458, 1)
SYNTHETIC_EPSILON_LO = Fraction(1, 3)
SYNTHETIC_EPSILON_HI = Fraction(1, 2)


class FrontierVerificationError(ValueError):
    """Raised when independent source-clock resolution fails closed."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise FrontierVerificationError(message)


def canonical_json_bytes(value: Any) -> bytes:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
    ).encode("utf-8")


def canonical_sha256(value: Any) -> str:
    return hashlib.sha256(canonical_json_bytes(value)).hexdigest()


def byte_sha256(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def load_json(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise FrontierVerificationError(f"cannot read JSON {path}: {exc}") from exc


def safe_repo_file(repo_root: Path, relative_path: str) -> Path:
    require(isinstance(relative_path, str) and bool(relative_path), "empty repository path")
    pure = PurePosixPath(relative_path)
    require(not pure.is_absolute(), f"absolute path is forbidden: {relative_path}")
    require(".." not in pure.parts, f"path traversal is forbidden: {relative_path}")
    require("\\" not in relative_path, f"non-POSIX path is forbidden: {relative_path}")
    require(":" not in pure.parts[0], f"drive-qualified path is forbidden: {relative_path}")
    root = repo_root.resolve()
    cursor = root
    for part in pure.parts:
        cursor = cursor / part
        require(not cursor.is_symlink(), f"symlinked input is forbidden: {relative_path}")
    resolved = cursor.resolve()
    try:
        resolved.relative_to(root)
    except ValueError as exc:
        raise FrontierVerificationError(
            f"path resolves outside repository: {relative_path}"
        ) from exc
    require(resolved.is_file(), f"referenced file is missing: {relative_path}")
    return resolved


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def parse_fraction(value: str, label: str) -> Fraction:
    require(
        isinstance(value, str) and re.fullmatch(r"-?[0-9]+/[0-9]+", value) is not None,
        f"{label} is not a canonical rational string",
    )
    try:
        parsed = Fraction(value)
    except (ValueError, ZeroDivisionError) as exc:
        raise FrontierVerificationError(f"invalid rational at {label}") from exc
    require(fraction_text(parsed) == value, f"{label} is not reduced canonically")
    return parsed


def validate_schema(frontier: Mapping[str, Any]) -> None:
    schema = load_json(SCHEMA_PATH)
    errors = sorted(
        Draft202012Validator(schema).iter_errors(frontier),
        key=lambda error: list(error.absolute_path),
    )
    if errors:
        details = "; ".join(
            f"{'/'.join(map(str, error.absolute_path)) or '$'}: {error.message}"
            for error in errors
        )
        raise FrontierVerificationError(f"frontier schema validation failed: {details}")


def verify_policy(frontier: Mapping[str, Any], repo_root: Path) -> dict[str, Any]:
    pin = frontier["policy_pin"]
    require(pin["path"] == POLICY_REL, "policy path drifted")
    path = safe_repo_file(repo_root, pin["path"])
    raw = path.read_bytes()
    policy = load_json(path)
    require(isinstance(policy, dict), "policy must be an object")
    require(pin["bytes"] == len(raw), "policy byte count mismatch")
    require(pin["byte_sha256"] == byte_sha256(raw), "policy byte digest mismatch")
    require(
        pin["canonical_json_sha256"] == canonical_sha256(policy),
        "policy canonical digest mismatch",
    )
    require(
        policy.get("schema") == "oph.source_clock.frontier_policy.v1",
        "policy schema mismatch",
    )
    require(policy.get("issue") == 633, "policy issue mismatch")
    require(policy.get("status") == EXPECTED_STATUS, "policy status mismatch")
    require(
        frontier["dependency_semantics"] == policy["dependency_semantics"],
        "dependency semantics do not match the pinned policy",
    )
    require(
        frontier["candidate_routes"] == policy["candidate_routes"],
        "candidate routes do not match the pinned policy",
    )
    require(frontier["open_gates"] == policy["open_gates"], "open gates do not match policy")
    require(
        frontier["downstream_consumers"] == policy["downstream_consumers"],
        "downstream consumer contract does not match policy",
    )
    require(
        frontier["target_firewall"]["forbidden_source_tokens"]
        == policy["forbidden_source_tokens"],
        "forbidden-token policy drifted",
    )
    return policy


def verify_dependency_semantics(frontier: Mapping[str, Any]) -> None:
    semantics = frontier["dependency_semantics"]
    require(
        semantics["hard_issue_dependencies"] == EXPECTED_HARD_DEPENDENCIES,
        "hard issue dependency set or order drifted",
    )
    require(semantics["alternative_routes_allowed"] is True, "alternative routes disabled")
    require(
        semantics["optional_route_owner_issues"] == EXPECTED_OPTIONAL_ROUTE_OWNERS,
        "optional route owner issue set or order drifted",
    )
    require(
        semantics["downstream_only_issues"] == EXPECTED_DOWNSTREAM_ONLY_ISSUES,
        "downstream-only issue set drifted",
    )
    require(
        not set(EXPECTED_OPTIONAL_ROUTE_OWNERS).intersection(EXPECTED_HARD_DEPENDENCIES),
        "optional route owner was promoted to a hard dependency",
    )
    require(
        not set(EXPECTED_DOWNSTREAM_ONLY_ISSUES).intersection(EXPECTED_HARD_DEPENDENCIES),
        "downstream issue was promoted to a hard dependency",
    )


def verify_candidate_routes(frontier: Mapping[str, Any]) -> None:
    routes = frontier["candidate_routes"]
    require(len(routes) == 1, "candidate route count drifted")
    route = routes[0]
    require(route["route_id"] == EXPECTED_CANDIDATE_ROUTE_ID, "candidate route id drifted")
    require(
        route["status"] == "optional_incomplete_candidate_not_hard_dependency",
        "candidate route status drifted",
    )
    require(route["required_for_issue_633"] is False, "cesium route became required")
    require(
        route["owner_issues"] == EXPECTED_OPTIONAL_ROUTE_OWNERS,
        "candidate route owner issue set or order drifted",
    )
    components = route["component_contracts"]
    require(
        [row["component_id"] for row in components] == EXPECTED_COMPONENT_IDS,
        "component set or order drifted",
    )
    expected_fields = {
        "source_map",
        "positive_interval",
        "scheme",
        "refinement_rule",
        "provenance_dag",
        "target_clean_receipt",
    }
    for row in components:
        require(set(row["required_fields"]) == expected_fields, "component fields incomplete")
        require(
            row["current_status"] != "physical_source_emitted",
            "a physical source component was forged",
        )
        require(
            set(row["owner_issues"]).issubset(set(EXPECTED_OPTIONAL_ROUTE_OWNERS)),
            "candidate component owner escapes the optional-route owner set",
        )
    require(frontier["source_payloads"] == [], "source payloads must remain absent")


def verify_diagnostics(
    frontier: Mapping[str, Any],
    policy: Mapping[str, Any],
    repo_root: Path,
) -> None:
    pins = frontier["diagnostic_reference_pins"]
    specs = policy["diagnostic_references"]
    require([pin["role"] for pin in pins] == EXPECTED_DIAGNOSTIC_ROLES, "diagnostic roles drifted")
    require(len(pins) == len(specs), "diagnostic reference count drifted")
    for pin, spec in zip(pins, specs, strict=True):
        require(pin["role"] == spec["role"], "diagnostic role mismatch")
        require(pin["path"] == spec["path"], "diagnostic path mismatch")
        require(
            pin["expected_artifact"] == spec["expected_artifact"],
            "diagnostic artifact contract mismatch",
        )
        require(
            pin["expected_status"] == spec["expected_status"],
            "diagnostic status contract mismatch",
        )
        require(pin["source_ancestry_allowed"] is False, "diagnostic source ancestry enabled")
        require(spec["source_ancestry_allowed"] is False, "policy enables diagnostic ancestry")
        path = safe_repo_file(repo_root, pin["path"])
        raw = path.read_bytes()
        payload = load_json(path)
        require(isinstance(payload, dict), "diagnostic payload must be an object")
        require(pin["bytes"] == len(raw), f"diagnostic byte count mismatch: {pin['path']}")
        require(
            pin["byte_sha256"] == byte_sha256(raw),
            f"diagnostic byte digest mismatch: {pin['path']}",
        )
        require(
            pin["canonical_json_sha256"] == canonical_sha256(payload),
            f"diagnostic canonical digest mismatch: {pin['path']}",
        )
        require(payload.get("artifact") == pin["expected_artifact"], "diagnostic artifact drift")
        require(payload.get("status") == pin["expected_status"], "diagnostic status drift")
        if pin["role"] == "legacy_gravity_checksum_skeleton":
            require(
                payload.get("status")
                == "skeleton_passes_forbidden_gravity_path_check_but_missing_component_certificates",
                "legacy checksum skeleton was promoted",
            )
        if pin["role"] == "synthetic_feshbach_fixture":
            require(payload.get("synthetic_fixture") is True, "Feshbach fixture lost synthetic tag")
            require(
                payload.get("physical_source_prediction_ready") is False,
                "synthetic Feshbach fixture was promoted",
            )
            require(
                payload.get("public_clock_gap_promotion_allowed") is False,
                "synthetic clock gap was promoted",
            )


def verify_unit_chart(frontier: Mapping[str, Any]) -> None:
    chart = frontier["exact_unit_chart"]
    expected = {
        "h_J_s": fraction_text(H_SI),
        "nu_Cs_Hz": fraction_text(NU_CS_HZ),
        "GeV_J": fraction_text(GEV_J),
        "c_m_s": fraction_text(C_SI),
    }
    for key, value in expected.items():
        require(chart[key] == value, f"exact SI definition drifted: {key}")
        parse_fraction(chart[key], f"exact_unit_chart/{key}")
    require(chart["source_selection_input"] is False, "SI chart entered source selection")
    require(
        chart["physical_clock_gap_input_present"] is False,
        "physical clock-gap input was forged",
    )


def verify_interval_inversion(frontier: Mapping[str, Any]) -> None:
    receipt = frontier["interval_inversion_receipt"]
    epsilon_lo = parse_fraction(receipt["epsilon_interval"]["lower"], "epsilon/lower")
    epsilon_hi = parse_fraction(receipt["epsilon_interval"]["upper"], "epsilon/upper")
    require(
        epsilon_lo == SYNTHETIC_EPSILON_LO and epsilon_hi == SYNTHETIC_EPSILON_HI,
        "synthetic epsilon theorem fixture drifted",
    )
    require(Fraction(0) < epsilon_lo <= epsilon_hi, "epsilon interval is not positive")
    numerator_j = H_SI * NU_CS_HZ
    expected_j_lo = numerator_j / epsilon_hi
    expected_j_hi = numerator_j / epsilon_lo
    expected_gev_lo = expected_j_lo / GEV_J
    expected_gev_hi = expected_j_hi / GEV_J
    actual_j_lo = parse_fraction(receipt["energy_interval_J"]["lower"], "energy_J/lower")
    actual_j_hi = parse_fraction(receipt["energy_interval_J"]["upper"], "energy_J/upper")
    actual_gev_lo = parse_fraction(
        receipt["energy_interval_GeV"]["lower"],
        "energy_GeV/lower",
    )
    actual_gev_hi = parse_fraction(
        receipt["energy_interval_GeV"]["upper"],
        "energy_GeV/upper",
    )
    require(
        (actual_j_lo, actual_j_hi) == (expected_j_lo, expected_j_hi),
        "joule interval inversion mismatch",
    )
    require(
        (actual_gev_lo, actual_gev_hi) == (expected_gev_lo, expected_gev_hi),
        "GeV interval conversion mismatch",
    )
    require(receipt["synthetic_fixture_only"] is True, "interval fixture was promoted")
    require(
        receipt["physical_application_input_present"] is False,
        "physical application input was forged",
    )
    require(all(receipt["checks"].values()), "stored exact arithmetic check failed")
    require(receipt["checks_pass"] is True, "interval receipt does not pass")


def _acyclic(nodes: list[dict[str, Any]], edges: list[dict[str, str]]) -> bool:
    ids = [node["id"] for node in nodes]
    if len(ids) != len(set(ids)):
        return False
    identifiers = set(ids)
    if any(edge["from"] not in identifiers or edge["to"] not in identifiers for edge in edges):
        return False
    outgoing = {identifier: [] for identifier in identifiers}
    indegree = {identifier: 0 for identifier in identifiers}
    for edge in edges:
        outgoing[edge["from"]].append(edge["to"])
        indegree[edge["to"]] += 1
    queue = sorted(identifier for identifier, degree in indegree.items() if degree == 0)
    seen = 0
    while queue:
        current = queue.pop(0)
        seen += 1
        for target in sorted(outgoing[current]):
            indegree[target] -= 1
            if indegree[target] == 0:
                queue.append(target)
                queue.sort()
    return seen == len(nodes)


def _path_exists(edges: list[dict[str, str]], source: str, target: str) -> bool:
    outgoing: dict[str, list[str]] = {}
    for edge in edges:
        outgoing.setdefault(edge["from"], []).append(edge["to"])
    stack = [source]
    seen: set[str] = set()
    while stack:
        current = stack.pop()
        if current == target:
            return True
        if current in seen:
            continue
        seen.add(current)
        stack.extend(outgoing.get(current, []))
    return False


def verify_provenance(frontier: Mapping[str, Any]) -> None:
    dag = frontier["provenance_dag"]
    nodes = dag["nodes"]
    edges = dag["edges"]
    node_by_id = {node["id"]: node for node in nodes}
    require(_acyclic(nodes, edges), "provenance DAG is cyclic or malformed")
    require(dag["acyclic"] is True, "stored acyclic flag is false")
    for component_id in EXPECTED_COMPONENT_IDS:
        require(
            node_by_id.get(component_id, {}).get("class")
            == "optional_candidate_component",
            f"component provenance class drifted: {component_id}",
        )
        require(
            _path_exists(edges, component_id, "cesium_candidate_clock_gap"),
            f"candidate component does not feed the cesium candidate: {component_id}",
        )
        require(
            not _path_exists(edges, component_id, "dimensionless_clock_gap"),
            f"optional candidate component feeds the generic clock gap: {component_id}",
        )
        require(
            not _path_exists(edges, component_id, "source_energy_interval"),
            f"optional candidate component feeds the source energy interval: {component_id}",
        )
    require(
        node_by_id.get("cesium_candidate_clock_gap", {}).get("class")
        == "blocked_candidate_output",
        "cesium candidate output provenance class drifted",
    )
    require(
        node_by_id.get("exact_si_unit_chart", {}).get("class") == "exact_unit_definition",
        "SI unit chart provenance class drifted",
    )
    require(
        node_by_id.get("route_neutral_dimensionless_clock_observable", {}).get("class")
        == "open_gate",
        "route-neutral dimensionless-clock gate provenance class drifted",
    )
    require(
        node_by_id.get("newton_g_composition", {}).get("class")
        == "downstream_open_composition",
        "Newton-G composition is not typed as downstream",
    )
    require(
        _path_exists(edges, "source_energy_interval", "newton_g_composition"),
        "Newton-G composition does not consume the source energy interval",
    )
    require(
        not _path_exists(edges, "newton_g_composition", "source_energy_interval"),
        "downstream Newton-G composition became a source-energy ancestor",
    )
    protected = dag["protected_outputs"]
    require(
        protected
        == ["dimensionless_clock_gap", "source_energy_interval", "source_g_si_interval"],
        "protected output set or order drifted",
    )
    discovered: list[list[str]] = []
    for diagnostic in EXPECTED_DIAGNOSTIC_ROLES:
        require(
            node_by_id.get(diagnostic, {}).get("class") == "diagnostic_non_source",
            f"diagnostic provenance class drifted: {diagnostic}",
        )
        for target in protected:
            if _path_exists(edges, diagnostic, target):
                discovered.append([diagnostic, target])
    require(not discovered, "diagnostic artifact has a path to a protected output")
    require(dag["diagnostic_to_protected_paths"] == [], "stored diagnostic path list is nonempty")
    discovered_candidate_paths: list[list[str]] = []
    for component_id in EXPECTED_COMPONENT_IDS:
        for target in protected:
            if _path_exists(edges, component_id, target):
                discovered_candidate_paths.append([component_id, target])
    require(
        not discovered_candidate_paths,
        "optional cesium candidate has a path to a generic protected output",
    )
    require(
        dag["optional_candidate_to_generic_paths"] == [],
        "stored optional-candidate path list is nonempty",
    )


def verify_firewall(frontier: Mapping[str, Any], policy: Mapping[str, Any]) -> None:
    firewall = frontier["target_firewall"]
    require(firewall["source_input_paths"] == [POLICY_REL], "source input closure drifted")
    require(firewall["source_payload_count"] == 0, "source payload count is nonzero")
    require(firewall["forbidden_source_token_hits"] == [], "target token hit is present")
    require(
        firewall["diagnostic_paths"]
        == [spec["path"] for spec in policy["diagnostic_references"]],
        "diagnostic path closure drifted",
    )
    require(
        firewall["diagnostics_excluded_from_source_ancestry"] is True,
        "diagnostics entered source ancestry",
    )
    require(
        firewall["exact_si_definitions_source_selection_allowed"] is False,
        "SI definitions entered source selection",
    )
    require(
        firewall["exact_si_definitions_final_unit_chart_allowed"] is True,
        "exact SI unit chart was disabled",
    )
    require(firewall["pass"] is True, "target firewall does not pass")


def verify_frontier(
    frontier: Mapping[str, Any],
    *,
    repo_root: Path = DEFAULT_REPO_ROOT,
) -> dict[str, Any]:
    require(isinstance(frontier, Mapping), "frontier must be an object")
    validate_schema(frontier)
    require(frontier["schema"] == EXPECTED_SCHEMA, "frontier schema mismatch")
    require(frontier["issue"] == 633, "frontier issue mismatch")
    require(frontier["status"] == EXPECTED_STATUS, "frontier status mismatch")
    policy = verify_policy(frontier, repo_root)
    verify_dependency_semantics(frontier)
    verify_candidate_routes(frontier)
    verify_diagnostics(frontier, policy, repo_root)
    verify_unit_chart(frontier)
    verify_interval_inversion(frontier)
    verify_provenance(frontier)
    verify_firewall(frontier, policy)
    require(
        [gate["gate_id"] for gate in frontier["open_gates"]] == EXPECTED_GATE_IDS,
        "open gate set or order drifted",
    )
    require(
        [gate["owner_issues"] for gate in frontier["open_gates"]]
        == [[634], [633], [633]],
        "route-neutral gate ownership drifted",
    )
    require(
        [row["acceptance_index"] for row in frontier["acceptance_progress"]]
        == [1, 2, 3, 4, 5, 6, 7],
        "acceptance row set or order drifted",
    )
    require(
        frontier["downstream_consumers"]
        == [
            {
                "issue": 334,
                "role": "newton_g_composition",
                "required_input": "source_energy_interval",
                "possible_output": "source_g_si_interval",
                "status": "downstream_open_not_issue_633_acceptance_gate",
            }
        ],
        "#334 is not confined to the downstream gravity composition",
    )
    require(
        all(row["status"] != "complete" for row in frontier["acceptance_progress"]),
        "issue acceptance was forged complete",
    )
    require(frontier["dimensionless_clock_gap_emitted"] is False, "clock gap was emitted")
    require(frontier["source_energy_interval_emitted"] is False, "source energy was emitted")
    require(frontier["source_g_si_interval_emitted"] is False, "G interval was emitted")
    require(frontier["physical_promotion_allowed"] is False, "physical promotion was enabled")
    require(all(frontier["checks"].values()), "stored frontier check failed")
    require(frontier["checks_pass"] is True, "frontier checks do not pass")
    body = {key: value for key, value in frontier.items() if key != "frontier_digest"}
    require(
        frontier["frontier_digest"] == canonical_sha256(body),
        "frontier digest mismatch",
    )
    return {
        "status": "PASS",
        "frontier_digest": frontier["frontier_digest"],
        "hard_issue_dependencies": frontier["dependency_semantics"][
            "hard_issue_dependencies"
        ],
        "alternative_routes_allowed": frontier["dependency_semantics"][
            "alternative_routes_allowed"
        ],
        "optional_candidate_route_count": len(frontier["candidate_routes"]),
        "optional_candidate_component_count": len(
            frontier["candidate_routes"][0]["component_contracts"]
        ),
        "physical_source_payload_count": len(frontier["source_payloads"]),
        "physical_promotion_allowed": frontier["physical_promotion_allowed"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=DEFAULT_REPO_ROOT)
    parser.add_argument("--frontier", type=Path, default=DEFAULT_FRONTIER)
    args = parser.parse_args()
    frontier = load_json(args.frontier)
    result = verify_frontier(frontier, repo_root=args.repo_root)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
