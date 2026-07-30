"""Mutation and replay gates for the issue-647 pre-generation lock."""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
from jsonschema import Draft202012Validator

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = PACKAGE_ROOT.parents[1]
BUILD_PROJECTION = Path("code/invariant_mining/tools/build_source_projection.py")
BUILD_FREEZE = Path("code/invariant_mining/tools/build_pregeneration_freeze.py")
VERIFY = Path(
    "code/invariant_mining/tools/verify_pregeneration_freeze_independent.py"
)


def run(repo_root: Path, script: Path, *arguments: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(script), *arguments],
        cwd=repo_root,
        check=False,
        capture_output=True,
        text=True,
    )


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def write(path: Path, value: dict) -> None:
    path.write_text(
        json.dumps(value, sort_keys=True, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def make_repo(tmp_path: Path) -> Path:
    repo = tmp_path / "repo"
    target_package = repo / "code" / "invariant_mining"
    shutil.copytree(PACKAGE_ROOT, target_package)
    policy = load(PACKAGE_ROOT / "policy" / "pregeneration_policy.json")
    for relative_text in policy["control_artifacts"]:
        relative = Path(relative_text)
        target = repo / relative
        if target.exists():
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(REPO_ROOT / relative, target)
    source_registry = load(PACKAGE_ROOT / "data" / "source_feature_registry.json")
    for feature in source_registry["features"]:
        for artifact in feature["artifacts"]:
            relative = Path(artifact["path"])
            target = repo / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(REPO_ROOT / relative, target)
    return repo


def assert_failed(result: subprocess.CompletedProcess[str], marker: str) -> None:
    assert result.returncode != 0
    assert marker in result.stdout + result.stderr


def test_committed_projection_and_freeze_are_byte_exact() -> None:
    projection = run(REPO_ROOT, BUILD_PROJECTION, "--check")
    freeze = run(REPO_ROOT, BUILD_FREEZE, "--check")
    independent = run(REPO_ROOT, VERIFY)
    assert projection.returncode == 0, projection.stdout + projection.stderr
    assert freeze.returncode == 0, freeze.stdout + freeze.stderr
    assert independent.returncode == 0, independent.stdout + independent.stderr
    assert "SOURCE_PROJECTION_VALID" in projection.stdout
    assert "PREGENERATION_FREEZE_VALID" in freeze.stdout
    assert "REGISTRY_FINALIZATION_LOCK_VALID" in independent.stdout


def test_all_control_documents_match_their_typed_schemas() -> None:
    document_schema = load(
        PACKAGE_ROOT / "schemas" / "pregeneration_documents.schema.json"
    )
    output_schema = load(
        PACKAGE_ROOT / "schemas" / "pregeneration_outputs.schema.json"
    )
    document_validator = Draft202012Validator(document_schema)
    output_validator = Draft202012Validator(output_schema)
    for path in sorted((PACKAGE_ROOT / "data").glob("*.json")) + [
        PACKAGE_ROOT / "policy" / "pregeneration_policy.json"
    ]:
        document_validator.validate(load(path))
    for path in sorted((PACKAGE_ROOT / "outputs").glob("*.json")):
        output_validator.validate(load(path))


def test_source_byte_mutation_is_rejected(tmp_path: Path) -> None:
    repo = make_repo(tmp_path)
    source_path = repo / "Lean" / "Screen" / "A5PortModule.lean"
    source_path.write_bytes(source_path.read_bytes() + b"\n-- mutation\n")
    assert_failed(run(repo, VERIFY), "SOURCE_PROJECTION_DRIFT")


def test_registered_source_omission_is_rejected(tmp_path: Path) -> None:
    repo = make_repo(tmp_path)
    path = repo / "code" / "invariant_mining" / "data" / "source_feature_registry.json"
    payload = load(path)
    payload["features"] = payload["features"][1:]
    write(path, payload)
    assert_failed(run(repo, BUILD_PROJECTION), "required source-feature set drift")


def test_target_injection_is_rejected_before_projection(tmp_path: Path) -> None:
    repo = make_repo(tmp_path)
    path = repo / "code" / "invariant_mining" / "data" / "producer_slots.json"
    payload = load(path)
    payload["slots"][0]["target_value"] = 137.035999
    write(path, payload)
    assert_failed(run(repo, BUILD_PROJECTION), "forbidden target keys")


def test_omitted_nuisance_direction_is_rejected(tmp_path: Path) -> None:
    repo = make_repo(tmp_path)
    path = repo / "code" / "invariant_mining" / "data" / "nuisance_registry.json"
    payload = load(path)
    payload["unresolved_directions"] = payload["unresolved_directions"][:-1]
    write(path, payload)
    assert_failed(run(repo, BUILD_PROJECTION), "required unresolved-nuisance set drift")


def test_completeness_regression_is_rejected(tmp_path: Path) -> None:
    repo = make_repo(tmp_path)
    path = repo / "code" / "invariant_mining" / "data" / "source_feature_registry.json"
    payload = load(path)
    payload["completeness"]["status"] = "BOUNDED_SEED_REGISTRY__NOT_FINAL"
    write(path, payload)
    assert_failed(run(repo, BUILD_PROJECTION), "completeness status drift")


def test_comparison_budget_mutation_is_rejected(tmp_path: Path) -> None:
    repo = make_repo(tmp_path)
    path = repo / "code" / "invariant_mining" / "policy" / "pregeneration_policy.json"
    payload = load(path)
    payload["campaign_comparison_budget"]["comparisons_consumed"] = 1
    write(path, payload)
    assert_failed(run(repo, BUILD_PROJECTION), "campaign comparison already consumed")


def test_exposure_quarantine_preclaim_is_rejected(tmp_path: Path) -> None:
    repo = make_repo(tmp_path)
    path = repo / "code" / "invariant_mining" / "data" / "exposed_data_registry.json"
    payload = load(path)
    payload["surfaces"][0]["quarantine_evidence"] = "claimed-without-custody"
    write(path, payload)
    assert_failed(run(repo, BUILD_PROJECTION), "quarantine evidence may not be pre-claimed")


def test_uncovered_slot_exposure_is_rejected(tmp_path: Path) -> None:
    repo = make_repo(tmp_path)
    path = repo / "code" / "invariant_mining" / "data" / "exposed_data_registry.json"
    payload = load(path)
    for row in payload["surfaces"]:
        if row["applies_to_slot_ids"] == ["a5_angular_rules"]:
            row["applies_to_slot_ids"] = ["observer_overlap_cross_spectra"]
    write(path, payload)
    assert_failed(run(repo, BUILD_PROJECTION), "producer slots without an exposure surface")


def test_direct_n_fallback_reentry_is_rejected(tmp_path: Path) -> None:
    repo = make_repo(tmp_path)
    path = repo / "code" / "invariant_mining" / "data" / "producer_slots.json"
    payload = load(path)
    row = next(
        item
        for item in payload["slots"]
        if item["slot_id"] == "direct_n_capacity_closure"
    )
    row["fallback_eligibility"] = "PENDING_PHYSICALIZATION"
    write(path, payload)
    assert_failed(run(repo, BUILD_PROJECTION), "direct-N fallback re-entry is forbidden")


def test_branch_producer_enablement_is_rejected(tmp_path: Path) -> None:
    repo = make_repo(tmp_path)
    path = repo / "code" / "invariant_mining" / "data" / "producer_slots.json"
    payload = load(path)
    payload["slots"][1]["execution_status"] = "READY_TO_GENERATE"
    write(path, payload)
    assert_failed(run(repo, BUILD_PROJECTION), "producer slot became executable")


def test_expression_grammar_expansion_is_rejected(tmp_path: Path) -> None:
    repo = make_repo(tmp_path)
    path = repo / "code" / "invariant_mining" / "data" / "observable_grammar.json"
    payload = load(path)
    payload["class_ids"].append("measured_target_polynomial")
    write(path, payload)
    assert_failed(run(repo, BUILD_PROJECTION), "required observable grammar class set drift")


def test_comparison_surface_leakage_is_rejected(tmp_path: Path) -> None:
    repo = make_repo(tmp_path)
    leaked = repo / "public_data" / "comparison.json"
    leaked.parent.mkdir(parents=True)
    leaked.write_text('{"measurement": 1}\n', encoding="utf-8")
    path = repo / "code" / "invariant_mining" / "data" / "source_feature_registry.json"
    payload = load(path)
    payload["features"][0]["artifacts"][0]["path"] = "public_data/comparison.json"
    write(path, payload)
    assert_failed(
        run(repo, BUILD_PROJECTION),
        "registered source enters a forbidden comparison surface",
    )


@pytest.mark.parametrize(
    "name",
    [
        "candidate_inventory.json",
        "candidates.json",
        "comparison.json",
        "scores.json",
        "verdict.json",
    ],
)
def test_candidate_and_scoring_outputs_are_rejected(
    tmp_path: Path,
    name: str,
) -> None:
    repo = make_repo(tmp_path)
    rogue = repo / "code" / "invariant_mining" / "outputs" / name
    rogue.write_text("{}\n", encoding="utf-8")
    assert_failed(run(repo, VERIFY), "UNREGISTERED_OUTPUT")


def test_freeze_identifier_mutation_is_rejected(tmp_path: Path) -> None:
    repo = make_repo(tmp_path)
    path = repo / "code" / "invariant_mining" / "outputs" / "pregeneration_freeze.json"
    payload = load(path)
    payload["freeze_id"] = "sha256:" + "0" * 64
    write(path, payload)
    assert_failed(run(repo, VERIFY), "PREGENERATION_FREEZE_DRIFT")
