#!/usr/bin/env python3
"""Mutation and replay tests for the bounded #630 source frontier."""

from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import subprocess
import sys

import pytest


PACKAGE = Path(__file__).resolve().parents[1]
REPO_ROOT = PACKAGE.parents[3]
BUILD = PACKAGE / "build_higgs_yukawa_source_frontier.py"
CHECK = PACKAGE / "check_higgs_yukawa_source_frontier.py"
OUTPUT = PACKAGE / "outputs" / "higgs_yukawa_source_frontier.json"
POLICY = PACKAGE / "data" / "higgs_yukawa_source_policy_v1.json"

sys.path.insert(0, str(PACKAGE))
from build_higgs_yukawa_source_frontier import (  # noqa: E402
    FrontierError,
    resolve_sources,
    validate_policy,
)


def load_output() -> dict:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def run_checker(path: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(CHECK), "--input", str(path)],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )


def assert_mutation_rejected(
    tmp_path: Path,
    payload: dict,
    expected_code: str,
) -> None:
    path = tmp_path / f"mutated_{expected_code.lower()}.json"
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    result = run_checker(path)
    assert result.returncode != 0
    assert f"{expected_code}:" in result.stdout + result.stderr


def walk_keys(value):
    if isinstance(value, dict):
        for key, child in value.items():
            yield key
            yield from walk_keys(child)
    elif isinstance(value, list):
        for child in value:
            yield from walk_keys(child)


def test_builder_is_deterministic_and_independent_checker_accepts(tmp_path: Path) -> None:
    first = tmp_path / "first.json"
    second = tmp_path / "second.json"
    subprocess.run([sys.executable, str(BUILD), "--output", str(first)], cwd=REPO_ROOT, check=True)
    subprocess.run([sys.executable, str(BUILD), "--output", str(second)], cwd=REPO_ROOT, check=True)
    assert first.read_bytes() == second.read_bytes()
    assert first.read_bytes() == OUTPUT.read_bytes()
    result = run_checker(first)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "PASS" in result.stdout


def test_artifact_emits_no_coefficient_assignment_or_physical_promotion() -> None:
    payload = load_output()
    forbidden = {
        "coefficient_value",
        "coefficient_values",
        "selected_value",
        "selected_values",
        "matrix_entries",
        "physical_mass",
        "physical_masses",
        "prediction",
    }
    assert not (set(walk_keys(payload)) & forbidden)
    assert payload["promotion_allowed"] is False
    assert payload["physical_source_action_emitted"] is False
    assert payload["coefficient_assignments_emitted"] is False
    assert all(
        row["status"] == "not_emitted"
        for row in payload["positive_source_objects"].values()
    )


def test_issue_503_is_non_gating_partial_receipt_ancestry() -> None:
    topology = load_output()["issue_topology"]
    assert topology["open_blocking_dependencies"] == [636, 637, 638]
    assert 503 not in topology["semantic_dependencies"]
    assert topology["non_gating_partial_receipt_ancestry"] == [
        {
            "issue": 503,
            "path": "code/geometry/runs/realized_event_receipt_report.json",
            "closure_required": False,
            "scope": "completed finite E1/E2/E4 screen-sheet receipts only; E3 bulk depth remains open and is not used",
        }
    ]


def test_issue_503_dependency_promotion_fails_closed(tmp_path: Path) -> None:
    payload = load_output()
    payload["issue_topology"]["open_blocking_dependencies"].insert(0, 503)
    assert_mutation_rejected(tmp_path, payload, "SCHEMA")


def test_source_pin_hash_mutation_fails_closed(tmp_path: Path) -> None:
    payload = load_output()
    payload["source_inputs"]["artifacts"][0]["byte_sha256"] = "0" * 64
    assert_mutation_rejected(tmp_path, payload, "PIN_HASH")


def test_target_path_injection_fails_before_resolution(tmp_path: Path) -> None:
    payload = load_output()
    payload["source_inputs"]["artifacts"][0]["path"] = (
        "code/particles/data/particle_reference_values.json"
    )
    assert_mutation_rejected(tmp_path, payload, "TARGET_PATH")


def test_physical_promotion_mutation_fails_closed(tmp_path: Path) -> None:
    payload = load_output()
    payload["promotion_allowed"] = True
    assert_mutation_rejected(tmp_path, payload, "PROMOTION")


def test_coefficient_assignment_key_mutation_fails_closed(tmp_path: Path) -> None:
    payload = load_output()
    payload["conditional_coefficient_space"]["scalar_kinetic"]["selected_value"] = "injected"
    assert_mutation_rejected(tmp_path, payload, "COEFFICIENT_VALUE")


def test_scalar_completion_collapse_fails_closed(tmp_path: Path) -> None:
    payload = load_output()
    witness = payload["nonidentifiability_witnesses"]["scalar_content"]
    witness["completion_B"] = witness["completion_A"]
    assert_mutation_rejected(tmp_path, payload, "SCALAR_WITNESS")


def test_yukawa_completion_collapse_fails_closed(tmp_path: Path) -> None:
    payload = load_output()
    witness = payload["nonidentifiability_witnesses"]["yukawa_matrices"]
    witness["completion_B"] = witness["completion_A"]
    assert_mutation_rejected(tmp_path, payload, "YUKAWA_WITNESS")


def test_fj_coordinate_completion_collapse_fails_closed(tmp_path: Path) -> None:
    payload = load_output()
    witness = payload["nonidentifiability_witnesses"]["v_chart_to_v_F"]
    witness["completion_B"] = witness["completion_A"]
    assert_mutation_rejected(tmp_path, payload, "FJ_WITNESS")


def test_issue_521_exact_separation_mutation_fails_closed(tmp_path: Path) -> None:
    payload = load_output()
    payload["nonidentifiability_witnesses"]["issue_521_exact_separation"][
        "exact_separation"
    ]["Higgs"] = "collapsed"
    assert_mutation_rejected(tmp_path, payload, "SEPARATION")


def test_conditional_dimension_mutation_is_rejected_by_strict_schema(tmp_path: Path) -> None:
    payload = load_output()
    payload["conditional_coefficient_space"]["yukawa_sector"][
        "total_complex_dimension"
    ] = 26
    assert_mutation_rejected(tmp_path, payload, "SCHEMA")


def test_source_status_mutation_fails_closed(tmp_path: Path) -> None:
    payload = load_output()
    payload["source_inputs"]["artifacts"][0]["status"] = "promoted"
    assert_mutation_rejected(tmp_path, payload, "SOURCE_STATUS")


def test_subject_digest_mutation_fails_closed(tmp_path: Path) -> None:
    payload = load_output()
    payload["subject_digest"] = "f" * 64
    assert_mutation_rejected(tmp_path, payload, "SUBJECT_DIGEST")


def test_policy_target_injection_is_rejected_by_producer() -> None:
    policy = json.loads(POLICY.read_text(encoding="utf-8"))
    mutated = deepcopy(policy)
    mutated["source_inputs"][0]["path"] = (
        "code/particles/data/particle_reference_values.json"
    )
    validate_policy(mutated)
    with pytest.raises(FrontierError) as exc_info:
        resolve_sources(mutated)
    assert exc_info.value.code == "TARGET_PATH"


def test_policy_unknown_field_is_rejected() -> None:
    policy = json.loads(POLICY.read_text(encoding="utf-8"))
    mutated = deepcopy(policy)
    mutated["unreviewed_escape_hatch"] = True
    with pytest.raises(FrontierError) as exc_info:
        validate_policy(mutated)
    assert exc_info.value.code == "POLICY_SCHEMA"
