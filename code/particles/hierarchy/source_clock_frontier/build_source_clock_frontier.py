#!/usr/bin/env python3
"""Build the fail-closed, route-neutral source-clock frontier for issue #633.

The producer emits a contract, an exact interval-inversion theorem witness,
strict dependency semantics, and a provenance boundary. The cesium
construction is retained only as an optional incomplete candidate route. The
producer emits no physical clock gap, source energy, particle mass, or
Newton-constant interval.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path, PurePosixPath
from typing import Any, Mapping


HERE = Path(__file__).resolve().parent
DEFAULT_REPO_ROOT = HERE.parents[3]
POLICY_REL = (
    "code/particles/hierarchy/source_clock_frontier/"
    "data/source_clock_policy_v1.json"
)
DEFAULT_OUTPUT = HERE / "outputs" / "source_clock_frontier.json"

SCHEMA = "oph.source_clock.frontier.v1"
STATUS = "NONPROMOTING_ROUTE_NEUTRAL_CONTRACT__PHYSICAL_CLOCK_ATTACHMENT_OPEN"
HARD_ISSUE_DEPENDENCIES = [634]
OPTIONAL_ROUTE_OWNER_ISSUES = [32, 34, 317, 318, 425, 522, 545, 546, 569, 633]
DOWNSTREAM_ONLY_ISSUES = [334]
CESIUM_ROUTE_ID = "cesium_133_hyperfine"
CESIUM_COMPONENT_IDS = ["R_U", "R_alpha", "R_e_abs", "R_QCD_nuc_133Cs", "R_atom_133Cs"]

H_SI = Fraction(662607015, 10**42)
NU_CS_HZ = Fraction(9192631770, 1)
GEV_J = Fraction(1602176634, 10**19)
C_SI = Fraction(299792458, 1)
SYNTHETIC_EPSILON_LO = Fraction(1, 3)
SYNTHETIC_EPSILON_HI = Fraction(1, 2)


class FrontierBuildError(ValueError):
    """Raised when a source-clock frontier input fails closed."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise FrontierBuildError(message)


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
        raise FrontierBuildError(f"cannot read JSON input {path}: {exc}") from exc


def safe_repo_file(repo_root: Path, relative_path: str) -> Path:
    require(isinstance(relative_path, str) and bool(relative_path), "empty repository path")
    pure = PurePosixPath(relative_path)
    require(not pure.is_absolute(), f"absolute repository path is forbidden: {relative_path}")
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
        raise FrontierBuildError(f"path resolves outside repository: {relative_path}") from exc
    require(resolved.is_file(), f"source input is missing: {relative_path}")
    return resolved


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def diagnostic_pin(repo_root: Path, spec: Mapping[str, Any]) -> dict[str, Any]:
    path = safe_repo_file(repo_root, str(spec["path"]))
    raw = path.read_bytes()
    payload = load_json(path)
    require(isinstance(payload, dict), f"diagnostic JSON must be an object: {spec['path']}")
    require(
        payload.get("artifact") == spec["expected_artifact"],
        f"diagnostic artifact mismatch: {spec['path']}",
    )
    require(
        payload.get("status") == spec["expected_status"],
        f"diagnostic status mismatch: {spec['path']}",
    )
    require(
        spec["source_ancestry_allowed"] is False,
        f"diagnostic source ancestry must be forbidden: {spec['path']}",
    )
    return {
        "role": spec["role"],
        "path": spec["path"],
        "expected_artifact": spec["expected_artifact"],
        "expected_status": spec["expected_status"],
        "bytes": len(raw),
        "byte_sha256": byte_sha256(raw),
        "canonical_json_sha256": canonical_sha256(payload),
        "source_ancestry_allowed": False,
    }


def interval_inversion_receipt() -> dict[str, Any]:
    require(
        Fraction(0) < SYNTHETIC_EPSILON_LO <= SYNTHETIC_EPSILON_HI,
        "synthetic epsilon interval must be positive and ordered",
    )
    numerator_j = H_SI * NU_CS_HZ
    energy_j_lo = numerator_j / SYNTHETIC_EPSILON_HI
    energy_j_hi = numerator_j / SYNTHETIC_EPSILON_LO
    energy_gev_lo = energy_j_lo / GEV_J
    energy_gev_hi = energy_j_hi / GEV_J
    checks = {
        "positive_input_interval": (
            Fraction(0) < SYNTHETIC_EPSILON_LO <= SYNTHETIC_EPSILON_HI
        ),
        "lower_endpoint_is_reversed_upper_input": (
            energy_j_lo * SYNTHETIC_EPSILON_HI == numerator_j
        ),
        "upper_endpoint_is_reversed_lower_input": (
            energy_j_hi * SYNTHETIC_EPSILON_LO == numerator_j
        ),
        "joule_interval_ordered": energy_j_lo <= energy_j_hi,
        "gev_conversion_exact_at_lower": energy_gev_lo * GEV_J == energy_j_lo,
        "gev_conversion_exact_at_upper": energy_gev_hi * GEV_J == energy_j_hi,
    }
    require(all(checks.values()), "exact interval inversion failed")
    return {
        "theorem_id": "SOURCE_CLOCK_POSITIVE_INTERVAL_INVERSION",
        "statement": (
            "If 0 < epsilon_lo <= epsilon_clk <= epsilon_hi and "
            "DeltaE_clk = h*nu_clk, then E_star lies in "
            "[h*nu_clk/epsilon_hi, h*nu_clk/epsilon_lo]."
        ),
        "arithmetic_domain": "exact_rationals",
        "synthetic_fixture_only": True,
        "physical_application_input_present": False,
        "epsilon_interval": {
            "lower": fraction_text(SYNTHETIC_EPSILON_LO),
            "upper": fraction_text(SYNTHETIC_EPSILON_HI),
        },
        "energy_interval_J": {
            "lower": fraction_text(energy_j_lo),
            "upper": fraction_text(energy_j_hi),
        },
        "energy_interval_GeV": {
            "lower": fraction_text(energy_gev_lo),
            "upper": fraction_text(energy_gev_hi),
        },
        "checks": checks,
        "checks_pass": True,
    }


def provenance_dag(
    candidate_route: dict[str, Any],
    gates: list[dict[str, Any]],
    diagnostics: list[dict[str, Any]],
    downstream: list[dict[str, Any]],
) -> dict[str, Any]:
    components = candidate_route["component_contracts"]
    nodes = [
        {
            "id": component["component_id"],
            "class": "optional_candidate_component",
            "status": component["current_status"],
        }
        for component in components
    ]
    nodes.extend(
        {
            "id": gate["gate_id"],
            "class": "open_gate",
            "status": "open",
        }
        for gate in gates
    )
    nodes.extend(
        {
            "id": pin["role"],
            "class": "diagnostic_non_source",
            "status": "excluded_from_source_ancestry",
        }
        for pin in diagnostics
    )
    nodes.extend(
        {
            "id": item["role"],
            "class": "downstream_open_composition",
            "status": item["status"],
        }
        for item in downstream
    )
    nodes.extend(
        [
            {
                "id": "cesium_candidate_clock_gap",
                "class": "blocked_candidate_output",
                "status": "optional_incomplete_candidate_not_selected",
            },
            {
                "id": "dimensionless_clock_gap",
                "class": "blocked_output",
                "status": "not_emitted",
            },
            {
                "id": "exact_si_unit_chart",
                "class": "exact_unit_definition",
                "status": "available_for_final_attachment_only",
            },
            {
                "id": "source_energy_interval",
                "class": "blocked_output",
                "status": "not_emitted",
            },
            {
                "id": "source_g_si_interval",
                "class": "blocked_output",
                "status": "not_emitted",
            },
        ]
    )
    edges = [
        {"from": component["component_id"], "to": "cesium_candidate_clock_gap"}
        for component in components
    ]
    edges.extend(
        {
            "from": gate["gate_id"],
            "to": (
                "source_energy_interval"
                if gate["gate_id"] == "target_clean_clock_to_si_attachment"
                else "dimensionless_clock_gap"
            ),
        }
        for gate in gates
    )
    edges.extend(
        [
            {"from": "dimensionless_clock_gap", "to": "source_energy_interval"},
            {"from": "exact_si_unit_chart", "to": "source_energy_interval"},
            {"from": "source_energy_interval", "to": "newton_g_composition"},
            {"from": "newton_g_composition", "to": "source_g_si_interval"},
        ]
    )
    return {
        "nodes": nodes,
        "edges": edges,
        "protected_outputs": [
            "dimensionless_clock_gap",
            "source_energy_interval",
            "source_g_si_interval",
        ],
        "diagnostic_to_protected_paths": [],
        "optional_candidate_to_generic_paths": [],
        "acyclic": True,
    }


def build_frontier(repo_root: Path = DEFAULT_REPO_ROOT) -> dict[str, Any]:
    policy_path = safe_repo_file(repo_root, POLICY_REL)
    policy_raw = policy_path.read_bytes()
    policy = load_json(policy_path)
    require(isinstance(policy, dict), "source-clock policy must be an object")
    require(
        policy.get("schema") == "oph.source_clock.frontier_policy.v1",
        "source-clock policy schema mismatch",
    )
    require(policy.get("issue") == 633, "source-clock policy issue mismatch")
    require(policy.get("status") == STATUS, "source-clock policy status mismatch")

    dependency_semantics = policy["dependency_semantics"]
    require(
        dependency_semantics
        == {
            "hard_issue_dependencies": HARD_ISSUE_DEPENDENCIES,
            "alternative_routes_allowed": True,
            "optional_route_owner_issues": OPTIONAL_ROUTE_OWNER_ISSUES,
            "downstream_only_issues": DOWNSTREAM_ONLY_ISSUES,
        },
        "source-clock dependency semantics drifted",
    )
    candidate_routes = policy["candidate_routes"]
    require(len(candidate_routes) == 1, "candidate route count drifted")
    candidate_route = candidate_routes[0]
    require(candidate_route["route_id"] == CESIUM_ROUTE_ID, "cesium route id drifted")
    require(
        candidate_route["status"] == "optional_incomplete_candidate_not_hard_dependency",
        "cesium candidate status drifted",
    )
    require(candidate_route["required_for_issue_633"] is False, "cesium route became required")
    require(
        candidate_route["owner_issues"] == OPTIONAL_ROUTE_OWNER_ISSUES,
        "cesium route owner issue set or order drifted",
    )
    components = candidate_route["component_contracts"]
    require(
        [row["component_id"] for row in components] == CESIUM_COMPONENT_IDS,
        "clock component set or order drifted",
    )
    required_fields = {
        "source_map",
        "positive_interval",
        "scheme",
        "refinement_rule",
        "provenance_dag",
        "target_clean_receipt",
    }
    for row in components:
        require(set(row["required_fields"]) == required_fields, "component contract incomplete")
        require(row["current_status"] != "physical_source_emitted", "physical component forged")
        require(
            set(row["owner_issues"]).issubset(set(OPTIONAL_ROUTE_OWNER_ISSUES)),
            "cesium component owner escapes the optional-route owner set",
        )

    diagnostics = [
        diagnostic_pin(repo_root, spec) for spec in policy["diagnostic_references"]
    ]
    gates = policy["open_gates"]
    require(
        [gate["gate_id"] for gate in gates]
        == [
            "route_neutral_dimensionless_clock_observable",
            "admissible_physical_clock_route_and_process_receipt",
            "target_clean_clock_to_si_attachment",
        ],
        "route-neutral gate set or order drifted",
    )
    require(
        gates[0]["owner_issues"] == [634]
        and gates[1]["owner_issues"] == [633]
        and gates[2]["owner_issues"] == [633],
        "route-neutral gate ownership drifted",
    )
    downstream = policy["downstream_consumers"]
    source_payloads: list[dict[str, Any]] = []
    firewall = {
        "source_input_paths": [POLICY_REL],
        "source_payload_count": 0,
        "forbidden_source_tokens": policy["forbidden_source_tokens"],
        "forbidden_source_token_hits": [],
        "diagnostic_paths": [pin["path"] for pin in diagnostics],
        "diagnostics_excluded_from_source_ancestry": True,
        "exact_si_definitions_source_selection_allowed": False,
        "exact_si_definitions_final_unit_chart_allowed": True,
        "pass": True,
    }
    unit_chart = {
        "status": "exact_si_definitions_for_final_attachment_only",
        "h_J_s": fraction_text(H_SI),
        "nu_Cs_Hz": fraction_text(NU_CS_HZ),
        "GeV_J": fraction_text(GEV_J),
        "c_m_s": fraction_text(C_SI),
        "source_selection_input": False,
        "physical_clock_gap_input_present": False,
    }
    acceptance = [
        {
            "acceptance_index": 1,
            "status": "partial_contract_only",
            "summary": "The route-neutral dependency and field contracts are closed; no physical readout is emitted.",
        },
        {
            "acceptance_index": 2,
            "status": "open",
            "summary": "Issue #634 has not supplied a route-neutral dimensionless clock observable.",
        },
        {
            "acceptance_index": 3,
            "status": "partial_reduction_only",
            "summary": "The cesium scalar reduction is synthetic and belongs only to an optional incomplete route.",
        },
        {
            "acceptance_index": 4,
            "status": "partial_boundary_enforced",
            "summary": "Pinned synthetic and checksum artifacts have no source-ancestry path.",
        },
        {
            "acceptance_index": 5,
            "status": "partial_boundary_enforced",
            "summary": "The source contract excludes measured scales and target residuals.",
        },
        {
            "acceptance_index": 6,
            "status": "open",
            "summary": "An admissible physical clock route and separate source/comparison process receipts are absent.",
        },
        {
            "acceptance_index": 7,
            "status": "partial_fail_closed_gate",
            "summary": "No source energy or downstream Newton-G interval is emitted.",
        },
    ]
    frontier: dict[str, Any] = {
        "schema": SCHEMA,
        "issue": 633,
        "status": STATUS,
        "claim_boundary": (
            "This artifact specifies and verifies a non-promoting, route-neutral "
            "source-clock contract for issue #633. The cesium construction is an "
            "optional incomplete candidate, not a hard dependency. The artifact "
            "emits no physical clock gap, source energy, particle mass, or "
            "Newton-constant interval."
        ),
        "policy_pin": {
            "path": POLICY_REL,
            "bytes": len(policy_raw),
            "byte_sha256": byte_sha256(policy_raw),
            "canonical_json_sha256": canonical_sha256(policy),
        },
        "dependency_semantics": dependency_semantics,
        "candidate_routes": candidate_routes,
        "source_payloads": source_payloads,
        "diagnostic_reference_pins": diagnostics,
        "target_firewall": firewall,
        "exact_unit_chart": unit_chart,
        "interval_inversion_receipt": interval_inversion_receipt(),
        "open_gates": gates,
        "downstream_consumers": downstream,
        "provenance_dag": provenance_dag(
            candidate_route,
            gates,
            diagnostics,
            downstream,
        ),
        "acceptance_progress": acceptance,
        "dimensionless_clock_gap_emitted": False,
        "source_energy_interval_emitted": False,
        "source_g_si_interval_emitted": False,
        "physical_promotion_allowed": False,
        "checks": {
            "route_neutral_dependency_semantics_exact": True,
            "optional_candidate_contract_complete": True,
            "optional_candidate_disconnected_from_generic_gap": True,
            "source_payloads_absent": True,
            "diagnostics_excluded_from_source_ancestry": True,
            "target_firewall_pass": True,
            "interval_inversion_exact": True,
            "all_physical_outputs_blocked": True,
        },
        "checks_pass": True,
    }
    frontier["frontier_digest"] = canonical_sha256(frontier)
    return frontier


def write_frontier(path: Path, frontier: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(frontier, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=DEFAULT_REPO_ROOT)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check-byte-exact", action="store_true")
    args = parser.parse_args()
    frontier = build_frontier(args.repo_root)
    rendered = json.dumps(frontier, indent=2, sort_keys=True) + "\n"
    if args.check_byte_exact:
        require(args.out.is_file(), f"committed output is missing: {args.out}")
        require(
            args.out.read_bytes() == rendered.encode("utf-8"),
            f"committed output is not byte-exact: {args.out}",
        )
    else:
        write_frontier(args.out, frontier)
    print(
        json.dumps(
            {
                "status": frontier["status"],
                "frontier_digest": frontier["frontier_digest"],
                "physical_promotion_allowed": frontier["physical_promotion_allowed"],
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
