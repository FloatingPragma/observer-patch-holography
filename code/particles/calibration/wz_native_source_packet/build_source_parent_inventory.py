#!/usr/bin/env python3
"""Build the fail-closed producer-side frontier for GitHub issue #594.

The output binds finite OPH source parents and records the exact interfaces
that remain open. It emits no action coefficients, pole coordinates, target
comparison, or physical-unit value. The external-SM calculation in issue #593
is treated only as a downstream consumer dependency.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path, PurePosixPath
from typing import Any, Mapping


HERE = Path(__file__).resolve().parent
DEFAULT_REPO_ROOT = HERE.parents[3]
POLICY_REL = (
    "code/particles/calibration/wz_native_source_packet/"
    "data/source_parent_policy_v1.json"
)
DEFAULT_OUTPUT = HERE / "outputs" / "source_parent_inventory.json"

SCHEMA = "oph.wz.source_parent_inventory.v1"
STATUS = "FINITE_SOURCE_PARENTS_BOUND__NATIVE_ACTION_AND_PHYSICAL_ATTACHMENT_OPEN"


class FrontierBuildError(ValueError):
    """Raised when an input to the source frontier fails closed."""


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
        require(not cursor.is_symlink(), f"symlinked source input is forbidden: {relative_path}")
    resolved = cursor.resolve()
    try:
        resolved.relative_to(root)
    except ValueError as exc:
        raise FrontierBuildError(f"path resolves outside repository: {relative_path}") from exc
    require(resolved.is_file(), f"source input is missing: {relative_path}")
    return resolved


def file_pin(repo_root: Path, spec: Mapping[str, Any]) -> tuple[dict[str, Any], Any]:
    path = safe_repo_file(repo_root, str(spec["path"]))
    raw = path.read_bytes()
    payload = load_json(path)
    require(isinstance(payload, dict), f"JSON source input must be an object: {spec['path']}")
    require(
        payload.get("schema") == spec["schema"],
        (
            f"schema mismatch for {spec['path']}: "
            f"expected {spec['schema']!r}, got {payload.get('schema')!r}"
        ),
    )
    return (
        {
            "role": spec["role"],
            "path": spec["path"],
            "schema": spec["schema"],
            "bytes": len(raw),
            "byte_sha256": byte_sha256(raw),
            "canonical_json_sha256": canonical_sha256(payload),
        },
        payload,
    )


def schema_pin(repo_root: Path, spec: Mapping[str, Any]) -> dict[str, Any]:
    path = safe_repo_file(repo_root, str(spec["path"]))
    raw = path.read_bytes()
    payload = load_json(path)
    require(isinstance(payload, dict), f"consumer schema must be an object: {spec['path']}")
    require(
        payload.get("$schema") == "https://json-schema.org/draft/2020-12/schema",
        f"consumer schema is not Draft 2020-12: {spec['path']}",
    )
    return {
        "slot": spec["slot"],
        "path": spec["path"],
        "status": spec["status"],
        "bytes": len(raw),
        "byte_sha256": byte_sha256(raw),
    }


def under_allowed_root(path: str, allowed_roots: list[str]) -> bool:
    pure = PurePosixPath(path)
    return any(
        pure == PurePosixPath(root) or PurePosixPath(root) in pure.parents
        for root in allowed_roots
    )


def _structured_token_hits(value: Any, tokens: list[str], at: str = "$") -> list[str]:
    hits: list[str] = []
    if isinstance(value, dict):
        for key, child in value.items():
            key_text = str(key).lower()
            for token in tokens:
                pattern = rf"(?<![a-z0-9_]){re.escape(token)}(?![a-z0-9_])"
                if re.search(pattern, key_text):
                    hits.append(f"{at}.{key}:{token}")
            hits.extend(_structured_token_hits(child, tokens, f"{at}.{key}"))
    elif isinstance(value, list):
        for index, child in enumerate(value):
            hits.extend(_structured_token_hits(child, tokens, f"{at}[{index}]"))
    elif isinstance(value, str):
        text = value.lower()
        for token in tokens:
            pattern = rf"(?<![a-z0-9_]){re.escape(token)}(?![a-z0-9_])"
            if re.search(pattern, text):
                hits.append(f"{at}:{token}")
    return hits


def expand_binding(
    repo_root: Path,
    spec: Mapping[str, Any],
    *,
    status: str,
) -> tuple[dict[str, Any], list[Any]]:
    pins: list[dict[str, Any]] = []
    payloads: list[Any] = []
    for file_spec in spec["files"]:
        pin, payload = file_pin(repo_root, file_spec)
        pins.append(pin)
        payloads.append(payload)
    return (
        {
            "role": spec["role"],
            "issue": spec["issue"],
            "verifier_id": spec["verifier_id"],
            "status": status,
            "files": pins,
            "usable_exports": spec["usable_exports"],
            "excluded_promotions": spec["excluded_promotions"],
        },
        payloads,
    )


def _acyclic(nodes: list[dict[str, Any]], edges: list[dict[str, str]]) -> bool:
    identifiers = {node["id"] for node in nodes}
    require(len(identifiers) == len(nodes), "source DAG contains duplicate node identifiers")
    require(
        all(edge["from"] in identifiers and edge["to"] in identifiers for edge in edges),
        "source DAG contains a dangling edge",
    )
    outgoing: dict[str, list[str]] = {identifier: [] for identifier in identifiers}
    indegree = {identifier: 0 for identifier in identifiers}
    for edge in edges:
        outgoing[edge["from"]].append(edge["to"])
        indegree[edge["to"]] += 1
    queue = sorted(identifier for identifier, degree in indegree.items() if degree == 0)
    visited: list[str] = []
    while queue:
        current = queue.pop(0)
        visited.append(current)
        for target in sorted(outgoing[current]):
            indegree[target] -= 1
            if indegree[target] == 0:
                queue.append(target)
                queue.sort()
    return len(visited) == len(nodes)


def _forbidden_ancestry_paths(
    nodes: list[dict[str, Any]],
    edges: list[dict[str, str]],
    protected_outputs: list[str],
    forbidden_classes: list[str],
) -> list[list[str]]:
    classes = {node["id"]: node["class"] for node in nodes}
    incoming: dict[str, list[str]] = {node["id"]: [] for node in nodes}
    for edge in edges:
        incoming[edge["to"]].append(edge["from"])
    paths: list[list[str]] = []

    def visit(current: str, trail: list[str]) -> None:
        if current in trail:
            paths.append(list(reversed(trail + [current])))
            return
        if classes[current] in forbidden_classes:
            paths.append(list(reversed(trail + [current])))
            return
        for parent in incoming[current]:
            visit(parent, trail + [current])

    for output in protected_outputs:
        visit(output, [])
    return sorted(paths)


def source_dag(
    positive: list[dict[str, Any]],
    conditional: list[dict[str, Any]],
    interfaces: list[dict[str, Any]],
) -> dict[str, Any]:
    nodes = [
        {
            "id": binding["role"],
            "class": "verified_finite_parent",
            "issue": binding["issue"],
            "status": "verified_finite_parent",
        }
        for binding in positive
    ]
    nodes.extend(
        {
            "id": binding["role"],
            "class": "conditional_context",
            "issue": binding["issue"],
            "status": "conditional_open_interface",
        }
        for binding in conditional
    )
    nodes.extend(
        {
            "id": item["gate_id"],
            "class": "open_interface",
            "issue": item["owner_issues"][0],
            "status": "open",
        }
        for item in interfaces
    )
    nodes.extend(
        [
            {
                "id": "oph_native_dimensionless_packet",
                "class": "candidate_output",
                "issue": 594,
                "status": "blocked",
            },
            {
                "id": "external_qft_pole_consumer",
                "class": "external_validation_consumer",
                "issue": 593,
                "status": "open_dependency",
            },
        ]
    )

    edges = [
        {"from": "screen_carrier", "to": "finite_port_current_algebra"},
        {"from": "finite_port_current_algebra", "to": "finite_matter_module"},
        {"from": "finite_matter_module", "to": "finite_global_form"},
        {"from": "screen_carrier", "to": "family_band_candidate"},
        {"from": "finite_matter_module", "to": "family_band_candidate"},
    ]
    for binding in positive:
        edges.append({"from": binding["role"], "to": "oph_native_dimensionless_packet"})
    for binding in conditional:
        edges.append({"from": binding["role"], "to": "oph_native_dimensionless_packet"})
    for interface in interfaces:
        edges.append({"from": interface["gate_id"], "to": "oph_native_dimensionless_packet"})
    edges.append(
        {"from": "oph_native_dimensionless_packet", "to": "external_qft_pole_consumer"}
    )

    protected = ["oph_native_dimensionless_packet"]
    forbidden_classes = [
        "experimental_target",
        "post_exposure_comparison",
        "calibrated_proxy",
        "inverse_target_adapter",
    ]
    is_acyclic = _acyclic(nodes, edges)
    require(is_acyclic, "constructed source DAG is cyclic")
    forbidden_paths = _forbidden_ancestry_paths(
        nodes,
        edges,
        protected,
        forbidden_classes,
    )
    require(not forbidden_paths, "constructed source DAG contains forbidden target ancestry")
    return {
        "nodes": nodes,
        "edges": edges,
        "protected_outputs": protected,
        "forbidden_node_classes": forbidden_classes,
        "acyclic": True,
        "forbidden_ancestry_paths": [],
    }


def acceptance_map() -> list[dict[str, Any]]:
    return [
        {
            "acceptance_index": 1,
            "status": "open",
            "summary": "finite group and matter types are bound, but no unique OPH action, complete census, or full Yukawa packet is emitted",
            "blocking_gates": [
                "event_and_spacetime_action_parent",
                "common_screen_electroweak_carrier",
                "source_complete_field_census",
                "scalar_higgs_and_fj_coordinate",
                "full_yukawa_operator_and_coefficients",
            ],
        },
        {
            "acceptance_index": 2,
            "status": "open",
            "summary": "v_chart and v_F remain distinct typed coordinates with no equality receipt",
            "blocking_gates": ["scalar_higgs_and_fj_coordinate"],
        },
        {
            "acceptance_index": 3,
            "status": "open",
            "summary": "source-derived running, thresholds, finite maps, Jacobians, masks, and remainders are absent",
            "blocking_gates": ["target_clean_rg_threshold_matching"],
        },
        {
            "acceptance_index": 4,
            "status": "open",
            "summary": "no unique deterministic source point or target-independent joint law and covariance is emitted",
            "blocking_gates": ["unique_source_root_and_joint_law"],
        },
        {
            "acceptance_index": 5,
            "status": "open",
            "summary": "substitution into unchanged validated QFT algorithms waits for the production consumer",
            "blocking_gates": ["validated_qft_consumer"],
        },
        {
            "acceptance_index": 6,
            "status": "open",
            "summary": "the inventory emits no pole coordinates or physical units; source clock attachment is absent",
            "blocking_gates": ["source_operational_clock"],
        },
        {
            "acceptance_index": 7,
            "status": "partial",
            "summary": "source paths and structured ancestry are allowlisted and target-free; hermetic runtime and human-selection receipts remain open",
            "blocking_gates": ["runtime_and_human_target_firewall"],
        },
        {
            "acceptance_index": 8,
            "status": "open",
            "summary": "the full conjunction cannot be replayed before the source packet and production consumer exist",
            "blocking_gates": [
                "unique_source_root_and_joint_law",
                "validated_qft_consumer",
                "runtime_and_human_target_firewall",
            ],
        },
        {
            "acceptance_index": 9,
            "status": "partial",
            "summary": "four finite source parents are hash-bound; the common carrier and physical family attachment remain open",
            "blocking_gates": [
                "common_screen_electroweak_carrier",
                "physical_family_and_matter_pole_attachment",
                "source_complete_field_census",
            ],
        },
    ]


def build_inventory(repo_root: Path) -> dict[str, Any]:
    repo_root = repo_root.resolve()
    policy_path = safe_repo_file(repo_root, POLICY_REL)
    policy_raw = policy_path.read_bytes()
    policy = load_json(policy_path)
    require(policy.get("schema") == "oph.wz.source_parent_policy.v1", "wrong policy schema")
    require(policy.get("issue") == 594, "policy is not for issue #594")
    require(
        policy.get("claim_lane") == "OPH_NATIVE_DIMENSIONLESS_SOURCE_FRONTIER",
        "policy claim lane drifted",
    )

    positive: list[dict[str, Any]] = []
    conditional: list[dict[str, Any]] = []
    scientific_payloads: list[Any] = []
    for spec in policy["positive_parents"]:
        binding, payloads = expand_binding(
            repo_root,
            spec,
            status="verified_finite_parent",
        )
        positive.append(binding)
        scientific_payloads.extend(payloads)
    for spec in policy["conditional_context"]:
        binding, payloads = expand_binding(
            repo_root,
            spec,
            status="conditional_context_only",
        )
        conditional.append(binding)
        scientific_payloads.extend(payloads)

    consumer_schemas = [schema_pin(repo_root, spec) for spec in policy["consumer_schemas"]]
    interfaces = [
        {
            "gate_id": item["gate_id"],
            "owner_issues": item["owner_issues"],
            "required_output": item["required_output"],
            "status": "open",
            "evidence": [],
        }
        for item in policy["open_interfaces"]
    ]

    firewall_policy = policy["firewall_policy"]
    allowed_roots = firewall_policy["allowed_source_roots"]
    resolved_paths = [
        file["path"]
        for binding in positive + conditional
        for file in binding["files"]
    ] + [item["path"] for item in consumer_schemas]
    require(len(resolved_paths) == len(set(resolved_paths)), "duplicate source input path")
    forbidden_paths = set(firewall_policy["forbidden_source_paths"])
    require(not forbidden_paths.intersection(resolved_paths), "forbidden target path mounted")
    require(
        all(under_allowed_root(path, allowed_roots) for path in resolved_paths),
        "source input lies outside the allowlisted roots",
    )
    tokens = [token.lower() for token in firewall_policy["forbidden_structured_tokens"]]
    content_hits: list[str] = []
    for path, payload in zip(
        [
            file["path"]
            for binding in positive + conditional
            for file in binding["files"]
        ],
        scientific_payloads,
        strict=True,
    ):
        content_hits.extend(
            f"{path}:{hit}" for hit in _structured_token_hits(payload, tokens)
        )
    require(not content_hits, f"forbidden structured target content found: {content_hits}")

    finite_seed = [
        {
            "role": binding["role"],
            "issue": binding["issue"],
            "files": [
                {
                    "path": item["path"],
                    "canonical_json_sha256": item["canonical_json_sha256"],
                }
                for item in binding["files"]
            ],
        }
        for binding in positive
    ]
    conditional_seed = [
        {
            "role": binding["role"],
            "issue": binding["issue"],
            "files": [
                {
                    "path": item["path"],
                    "canonical_json_sha256": item["canonical_json_sha256"],
                }
                for item in binding["files"]
            ],
        }
        for binding in conditional
    ]

    inventory: dict[str, Any] = {
        "schema": SCHEMA,
        "issue": 594,
        "claim_lane": "OPH_NATIVE_DIMENSIONLESS_SOURCE_FRONTIER",
        "campaign_classification": "post_exposure_validation",
        "status": STATUS,
        "promotion_allowed": False,
        "policy": {
            "path": POLICY_REL,
            "bytes": len(policy_raw),
            "byte_sha256": byte_sha256(policy_raw),
            "canonical_json_sha256": canonical_sha256(policy),
        },
        "unit_scope": {
            "native_coordinates": "E_star_normalized_dimensionless",
            "physical_clock_required": True,
            "dimensionful_values_present": False,
            "emitted_observables": [],
        },
        "positive_parent_bindings": positive,
        "conditional_context": conditional,
        "coordinate_bridge": {
            "source_coordinate": "v_chart",
            "renormalized_coordinate": "v_F",
            "equality_receipt": None,
            "status": "open",
            "relabel_allowed": False,
        },
        "consumer_contract": {
            "issue": 593,
            "status": "open_dependency",
            "schemas": consumer_schemas,
            "frozen_algorithm_substitution_ready": False,
            "common_subject_digest_ready": False,
        },
        "open_interfaces": interfaces,
        "source_dag": source_dag(positive, conditional, interfaces),
        "target_firewall": {
            "allowed_source_roots": allowed_roots,
            "resolved_source_paths": sorted(resolved_paths),
            "forbidden_source_paths": firewall_policy["forbidden_source_paths"],
            "forbidden_structured_tokens": firewall_policy["forbidden_structured_tokens"],
            "path_and_content_scan_passed": True,
            "comparison_channel_present": False,
            "network_required": False,
            "sealed_input_replay_supported": True,
            "runtime_execution_receipt": "open",
            "human_formula_selection_ancestry": "open",
            "full_gate_satisfied": False,
        },
        "acceptance_map": acceptance_map(),
        "finite_parent_bundle_digest": canonical_sha256(finite_seed),
        "conditional_context_digest": canonical_sha256(conditional_seed),
    }
    inventory["inventory_digest"] = canonical_sha256(inventory)
    return inventory


def write_inventory(path: Path, inventory: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(inventory, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=DEFAULT_REPO_ROOT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check-byte-exact", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    inventory = build_inventory(args.repo_root)
    rendered = json.dumps(inventory, indent=2, sort_keys=True) + "\n"
    if args.check_byte_exact:
        require(args.output.is_file(), f"stored inventory is missing: {args.output}")
        require(
            args.output.read_bytes() == rendered.encode("utf-8"),
            "stored source-parent inventory is not byte-exact",
        )
    else:
        write_inventory(args.output, inventory)
    print(
        json.dumps(
            {
                "status": inventory["status"],
                "positive_parent_count": len(inventory["positive_parent_bindings"]),
                "conditional_context_count": len(inventory["conditional_context"]),
                "open_interface_count": len(inventory["open_interfaces"]),
                "promotion_allowed": inventory["promotion_allowed"],
                "inventory_digest": inventory["inventory_digest"],
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
