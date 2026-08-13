"""Adversarial gates for predictive-observation freeze lineage custody."""

from __future__ import annotations

import copy
import json
import shutil
import subprocess
from pathlib import Path

import pytest

import build_architecture_versions as architecture
import build_observation_ledger as observation
import prediction_lineage_custody as custody


def _git(root: Path, *args: str) -> str:
    return subprocess.run(
        ["git", *args],
        cwd=root,
        check=True,
        stdout=subprocess.PIPE,
        text=True,
    ).stdout.strip()


def _inputs() -> tuple[dict, dict, dict, list[dict]]:
    lineages = custody.load_json(custody.LINEAGE_PATH)
    frozen = custody.load_json(custody.FROZEN_PATH)
    architectures = architecture.load_json(architecture.REGISTER_PATH)
    rows = observation.load_ledger(observation.LEDGER_PATH)["rows"]
    return lineages, frozen, architectures, rows


def _event(architectures: dict, predecessor: str, row: dict) -> dict:
    current = architectures["current_version"]
    version = next(row for row in architectures["versions"] if row["id"] == current)
    event = {
        "id": "FZE-001",
        "source_architecture_version": current,
        "source_architecture_snapshot_sha256": version["snapshot_sha256"],
        "observation_row_id": row["id"],
        "observation_contract_sha256": custody.canonical_sha256(
            custody.observation_contract(row)
        ),
        "supersedes_lineage_ids": ["FZ-06"],
        "frozen_utc": "2026-08-13T12:00:00Z",
        "target": "A distinct source-derived ringdown target on a declared estimator.",
        "comparison_protocol": "Use the preregistered event list, estimator, covariance, and blinded comparison procedure.",
        "kill_band": "The declared interval excludes the measured estimator under the registered covariance rule.",
        "precomparison_statement": "This event was committed before any registered comparison data were inspected.",
        "target_payload_sha256": "",
        "predecessor_register_revision": predecessor,
    }
    event["target_payload_sha256"] = custody.canonical_sha256(
        custody.event_payload(event)
    )
    return event


def _copy_current_surface(source: Path, destination: Path) -> None:
    for relative in (
        "claims/frozen_prediction_architecture_lineages.json",
        "claims/frozen_prediction_register.json",
        "tracking/architecture_versions.json",
        "tracking/observation_ledger.json",
    ):
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source / relative, target)


def _two_commit_event_repo(
    tmp_path: Path,
    *,
    origin_status: str | None = None,
    origin_pointer: bool = False,
    origin_audit_pointer: bool = False,
    historical_attained_then_downgraded: bool = False,
) -> tuple[Path, dict, dict, list[dict]]:
    repo = tmp_path / "repo"
    subprocess.run(
        ["git", "clone", "-q", "--shared", str(custody.ROOT), str(repo)],
        check=True,
    )
    _copy_current_surface(custody.ROOT, repo)
    _git(repo, "config", "user.name", "Lineage Test")
    _git(repo, "config", "user.email", "lineage-test@example.invalid")
    _git(repo, "add", "claims/frozen_prediction_architecture_lineages.json")
    _git(repo, "add", "claims/frozen_prediction_register.json")
    _git(repo, "add", "tracking/architecture_versions.json")
    _git(repo, "add", "tracking/observation_ledger.json")
    _git(repo, "commit", "-q", "-m", "install lineage baseline")
    baseline = _git(repo, "rev-parse", "HEAD")

    path = repo / "claims/frozen_prediction_architecture_lineages.json"
    lineages = json.loads(path.read_text(encoding="utf-8"))
    architectures = json.loads(
        (repo / "tracking/architecture_versions.json").read_text(encoding="utf-8")
    )
    ledger_path = repo / "tracking/observation_ledger.json"
    ledger = json.loads(ledger_path.read_text(encoding="utf-8"))
    row = next(item for item in ledger["rows"] if item["id"] == "OL-B4")
    if historical_attained_then_downgraded:
        row["status"] = "attained"
        ledger_path.write_text(json.dumps(ledger, indent=2) + "\n", encoding="utf-8")
        _git(repo, "add", "tracking/observation_ledger.json")
        _git(repo, "commit", "-q", "-m", "record attained row before freeze")
        row["status"] = "owed"
        ledger_path.write_text(json.dumps(ledger, indent=2) + "\n", encoding="utf-8")
        _git(repo, "add", "tracking/observation_ledger.json")
        _git(repo, "commit", "-q", "-m", "downgrade row before freeze")
        baseline = _git(repo, "rev-parse", "HEAD")
    event = _event(architectures, baseline, row)
    lineages["events"].append(event)
    origin_ledger_changed = False
    if origin_status is not None:
        row["status"] = origin_status
        origin_ledger_changed = True
    if origin_pointer:
        row["prediction_event"] = {
            "id": event["id"],
            "target_payload_sha256": event["target_payload_sha256"],
        }
        origin_ledger_changed = True
    if origin_audit_pointer:
        ledger["audit_pointers"][row["id"]] = ["AUD-LAUNDERED"]
        origin_ledger_changed = True
    path.write_text(json.dumps(lineages, indent=2) + "\n", encoding="utf-8")
    _git(repo, "add", "claims/frozen_prediction_architecture_lineages.json")
    if origin_ledger_changed:
        ledger_path.write_text(json.dumps(ledger, indent=2) + "\n", encoding="utf-8")
        _git(repo, "add", "tracking/observation_ledger.json")
    _git(repo, "commit", "-q", "-m", "freeze predictive event")
    origin = _git(repo, "rev-parse", "HEAD")
    origin_blob = _git(
        repo,
        "rev-parse",
        f"{origin}:claims/frozen_prediction_architecture_lineages.json",
    )

    lineages["event_anchors"].append(
        {
            "id": event["id"],
            "origin_revision": origin,
            "event_sha256": custody.event_sha256(event),
            "register_git_blob_sha1": origin_blob,
        }
    )
    row["prediction_event"] = {
        "id": event["id"],
        "target_payload_sha256": event["target_payload_sha256"],
    }
    path.write_text(json.dumps(lineages, indent=2) + "\n", encoding="utf-8")
    ledger_path.write_text(json.dumps(ledger, indent=2) + "\n", encoding="utf-8")
    _git(repo, "add", "claims/frozen_prediction_architecture_lineages.json")
    _git(repo, "add", "tracking/observation_ledger.json")
    _git(repo, "commit", "-q", "-m", "anchor predictive event")
    frozen = json.loads(
        (repo / "claims/frozen_prediction_register.json").read_text(encoding="utf-8")
    )
    rows = ledger["rows"]
    return repo, frozen, architectures, rows


def test_current_lineage_custody_is_valid_and_every_predictive_row_is_mapped() -> None:
    lineages, frozen, architectures, rows = _inputs()
    state = custody.validate(lineages, frozen, architectures, rows)
    predictive_ids = [row["id"] for row in rows if row["rung"] == "predictive"]
    assert [row["observation_row_id"] for row in state["bindings"]] == predictive_ids
    assert state["events"] == []


def test_superseded_void_baseline_cannot_be_revived_as_pending() -> None:
    lineages, frozen, architectures, rows = _inputs()
    fz06 = next(row for row in lineages["rows"] if row["id"] == "FZ-06")
    fz06["lineage_status"] = "unbound_pending_freeze"
    fz06["current_version_eligibility"] = "ineligible_until_version_bound_freeze"
    fz06["transition_policy"] = custody.PENDING_POLICY
    with pytest.raises(SystemExit, match="baseline lineage rows drifted"):
        custody.validate(lineages, frozen, architectures, rows)


def test_predictive_binding_removal_fails_closed() -> None:
    lineages, frozen, architectures, rows = _inputs()
    lineages["predictive_observation_bindings"].pop()
    with pytest.raises(SystemExit, match="exactly match predictive ledger rows"):
        custody.validate(lineages, frozen, architectures, rows)


def test_coherent_candidate_binding_swap_fails_before_first_v2_commit() -> None:
    lineages, frozen, architectures, rows = _inputs()
    h7 = next(
        row
        for row in lineages["predictive_observation_bindings"]
        if row["observation_row_id"] == "OL-H7"
    )
    i3 = next(
        row
        for row in lineages["predictive_observation_bindings"]
        if row["observation_row_id"] == "OL-I3"
    )
    h7["candidate_baseline_lineages"] = ["FZ-07"]
    i3["candidate_baseline_lineages"] = ["FZ-03", "FZ-08"]
    with pytest.raises(SystemExit, match="complete static v2 lineage surface"):
        custody.validate(lineages, frozen, architectures, rows)


def test_legacy_row_cannot_be_presented_as_a_pending_candidate() -> None:
    lineages, frozen, architectures, rows = _inputs()
    binding = next(
        row
        for row in lineages["predictive_observation_bindings"]
        if row["observation_row_id"] == "OL-F4"
    )
    binding["candidate_baseline_lineages"] = ["FZ-11"]
    binding["historical_only_lineages"] = ["FZ-12"]
    binding["no_candidate_reason"] = None
    with pytest.raises(SystemExit, match="candidate FZ-11 is not pending"):
        custody.validate(lineages, frozen, architectures, rows)


def test_unanchored_event_cannot_promote_a_predictive_row() -> None:
    lineages, frozen, architectures, rows = _inputs()
    row = next(item for item in rows if item["id"] == "OL-B4")
    lineages["events"].append(
        _event(architectures, _git(custody.ROOT, "rev-parse", "HEAD"), row)
    )
    state = custody.validate(lineages, frozen, architectures, rows)
    mutated = copy.deepcopy(rows)
    row = next(item for item in mutated if item["id"] == "OL-B4")
    row["status"] = "attained"
    with pytest.raises(SystemExit, match="has no exact event pointer"):
        custody.require_predictive_promotions(mutated, state)


def test_worktree_only_anchor_is_not_operational() -> None:
    lineages, frozen, architectures, rows = _inputs()
    row = next(item for item in rows if item["id"] == "OL-B4")
    event = _event(
        architectures, _git(custody.ROOT, "rev-parse", "HEAD"), row
    )
    lineages["events"].append(event)
    lineages["event_anchors"].append(
        {
            "id": event["id"],
            "origin_revision": _git(custody.ROOT, "rev-parse", "HEAD"),
            "event_sha256": custody.event_sha256(event),
            "register_git_blob_sha1": "0" * 40,
        }
    )
    with pytest.raises(SystemExit, match="must already be committed at HEAD"):
        custody.validate(lineages, frozen, architectures, rows)


def test_event_origin_must_be_on_head_first_parent(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    repo, frozen, architectures, rows = _two_commit_event_repo(tmp_path)
    lineages = json.loads(
        (repo / "claims/frozen_prediction_architecture_lineages.json").read_text(
            encoding="utf-8"
        )
    )
    origin = lineages["event_anchors"][0]["origin_revision"]
    real = custody._is_first_parent_ancestor

    def reject_origin(root: Path, revision: str) -> bool:
        return False if revision == origin else real(root, revision)

    monkeypatch.setattr(custody, "_is_first_parent_ancestor", reject_origin)
    with pytest.raises(SystemExit, match="not on the HEAD first-parent chain"):
        custody.validate(lineages, frozen, architectures, rows, root=repo)


def test_two_commit_event_qualifies_exactly_its_current_predictive_row(
    tmp_path: Path,
) -> None:
    repo, frozen, architectures, rows = _two_commit_event_repo(tmp_path)
    lineages = json.loads(
        (repo / "claims/frozen_prediction_architecture_lineages.json").read_text(
            encoding="utf-8"
        )
    )
    state = custody.validate(
        lineages, frozen, architectures, rows, root=repo
    )
    assert [event["id"] for event in state["qualifying_events_by_row"]["OL-B4"]] == [
        "FZE-001"
    ]
    mutated = copy.deepcopy(rows)
    next(row for row in mutated if row["id"] == "OL-B4")["status"] = "attained"
    custody.require_predictive_promotions(mutated, state)


def test_retroactive_freeze_cannot_launder_an_already_attained_row(
    tmp_path: Path,
) -> None:
    repo, frozen, architectures, rows = _two_commit_event_repo(
        tmp_path, origin_status="attained"
    )
    lineages = json.loads(
        (repo / "claims/frozen_prediction_architecture_lineages.json").read_text(
            encoding="utf-8"
        )
    )
    with pytest.raises(SystemExit, match="origin observation row was already attained"):
        custody.validate(lineages, frozen, architectures, rows, root=repo)


def test_two_step_downgrade_cannot_launder_prior_attainment(tmp_path: Path) -> None:
    repo, frozen, architectures, rows = _two_commit_event_repo(
        tmp_path, historical_attained_then_downgraded=True
    )
    lineages = json.loads(
        (repo / "claims/frozen_prediction_architecture_lineages.json").read_text(
            encoding="utf-8"
        )
    )
    with pytest.raises(SystemExit, match="same observation contract was already attained"):
        custody.validate(lineages, frozen, architectures, rows, root=repo)


def test_event_origin_must_be_pointer_free(tmp_path: Path) -> None:
    repo, frozen, architectures, rows = _two_commit_event_repo(
        tmp_path, origin_pointer=True
    )
    lineages = json.loads(
        (repo / "claims/frozen_prediction_architecture_lineages.json").read_text(
            encoding="utf-8"
        )
    )
    with pytest.raises(SystemExit, match="origin observation row was not pointer-free"):
        custody.validate(lineages, frozen, architectures, rows, root=repo)


def test_event_origin_cannot_have_a_promotion_audit_pointer(tmp_path: Path) -> None:
    repo, frozen, architectures, rows = _two_commit_event_repo(
        tmp_path, origin_audit_pointer=True
    )
    lineages = json.loads(
        (repo / "claims/frozen_prediction_architecture_lineages.json").read_text(
            encoding="utf-8"
        )
    )
    with pytest.raises(SystemExit, match="already had a promotion audit pointer"):
        custody.validate(lineages, frozen, architectures, rows, root=repo)


def test_event_cannot_supersede_lineage_from_another_observation_row() -> None:
    lineages, frozen, architectures, rows = _inputs()
    row = next(item for item in rows if item["id"] == "OL-H7")
    event = _event(
        architectures, _git(custody.ROOT, "rev-parse", "HEAD"), row
    )
    event["supersedes_lineage_ids"] = ["FZ-11"]
    event["target_payload_sha256"] = custody.canonical_sha256(
        custody.event_payload(event)
    )
    lineages["events"].append(event)
    with pytest.raises(SystemExit, match="declared candidate or historical context"):
        custody.validate(lineages, frozen, architectures, rows)


def test_anchored_event_rewrite_fails_against_origin(tmp_path: Path) -> None:
    repo, frozen, architectures, rows = _two_commit_event_repo(tmp_path)
    path = repo / "claims/frozen_prediction_architecture_lineages.json"
    lineages = json.loads(path.read_text(encoding="utf-8"))
    lineages["events"][0]["kill_band"] = "A post-hoc weakened rule."
    lineages["events"][0]["target_payload_sha256"] = custody.canonical_sha256(
        custody.event_payload(lineages["events"][0])
    )
    lineages["event_anchors"][0]["event_sha256"] = custody.event_sha256(
        lineages["events"][0]
    )
    path.write_text(json.dumps(lineages, indent=2) + "\n", encoding="utf-8")
    _git(repo, "add", "claims/frozen_prediction_architecture_lineages.json")
    _git(repo, "commit", "-q", "-m", "rewrite event and matching anchor")
    with pytest.raises(SystemExit, match="origin must be the exact first append"):
        custody.validate(lineages, frozen, architectures, rows, root=repo)


@pytest.mark.parametrize("surface", ["policy", "binding"])
def test_committed_static_lineage_surface_cannot_be_rewritten(
    tmp_path: Path, surface: str
) -> None:
    repo, frozen, architectures, rows = _two_commit_event_repo(tmp_path)
    path = repo / "claims/frozen_prediction_architecture_lineages.json"
    lineages = json.loads(path.read_text(encoding="utf-8"))
    if surface == "policy":
        lineages["policy"] += " Post-hoc exception."
    else:
        lineages["predictive_observation_bindings"][0][
            "no_candidate_reason"
        ] += " Post-hoc exception."
    path.write_text(json.dumps(lineages, indent=2) + "\n", encoding="utf-8")
    _git(repo, "add", "claims/frozen_prediction_architecture_lineages.json")
    _git(repo, "commit", "-q", "-m", f"rewrite static {surface}")
    with pytest.raises(SystemExit, match="complete static v2 lineage surface"):
        custody.validate(lineages, frozen, architectures, rows, root=repo)


def test_prediction_pointer_hash_must_name_the_exact_event(tmp_path: Path) -> None:
    repo, frozen, architectures, rows = _two_commit_event_repo(tmp_path)
    lineages = json.loads(
        (repo / "claims/frozen_prediction_architecture_lineages.json").read_text(
            encoding="utf-8"
        )
    )
    row = next(item for item in rows if item["id"] == "OL-B4")
    row["prediction_event"]["target_payload_sha256"] = "0" * 64
    with pytest.raises(SystemExit, match="prediction_event target hash drifted"):
        custody.validate(lineages, frozen, architectures, rows, root=repo)


def test_later_event_cannot_ambiguously_qualify_an_earlier_pointer(
    tmp_path: Path,
) -> None:
    repo, frozen, architectures, rows = _two_commit_event_repo(tmp_path)
    path = repo / "claims/frozen_prediction_architecture_lineages.json"
    lineages = json.loads(path.read_text(encoding="utf-8"))
    row = next(item for item in rows if item["id"] == "OL-B4")
    later = copy.deepcopy(lineages["events"][0])
    later["id"] = "FZE-002"
    later["target"] = "A different later ringdown target."
    later["predecessor_register_revision"] = _git(repo, "rev-parse", "HEAD")
    later["target_payload_sha256"] = custody.canonical_sha256(
        custody.event_payload(later)
    )
    lineages["events"].append(later)
    row["prediction_event"] = {
        "id": later["id"],
        "target_payload_sha256": later["target_payload_sha256"],
    }
    state = custody.validate(lineages, frozen, architectures, rows, root=repo)
    row["status"] = "attained"
    with pytest.raises(SystemExit, match="lacks exactly one anchored"):
        custody.require_predictive_promotions(rows, state)


def test_event_is_bound_to_immutable_observation_target_contract(
    tmp_path: Path,
) -> None:
    repo, frozen, architectures, rows = _two_commit_event_repo(tmp_path)
    lineages = json.loads(
        (repo / "claims/frozen_prediction_architecture_lineages.json").read_text(
            encoding="utf-8"
        )
    )
    row = next(item for item in rows if item["id"] == "OL-B4")
    row["target"] = "A broader post-hoc target"
    with pytest.raises(SystemExit, match="prediction_event contract differs"):
        custody.validate(lineages, frozen, architectures, rows, root=repo)


def test_event_cannot_acquire_a_posthoc_premise_bridge(tmp_path: Path) -> None:
    repo, frozen, architectures, rows = _two_commit_event_repo(tmp_path)
    lineages = json.loads(
        (repo / "claims/frozen_prediction_architecture_lineages.json").read_text(
            encoding="utf-8"
        )
    )
    row = next(item for item in rows if item["id"] == "OL-B4")
    row["open_premises"] = ["PR-01"]
    with pytest.raises(SystemExit, match="prediction_event contract differs"):
        custody.validate(lineages, frozen, architectures, rows, root=repo)


def test_old_event_remains_historical_after_pointer_is_cleared_and_contract_changes(
    tmp_path: Path,
) -> None:
    repo, frozen, architectures, rows = _two_commit_event_repo(tmp_path)
    lineages = json.loads(
        (repo / "claims/frozen_prediction_architecture_lineages.json").read_text(
            encoding="utf-8"
        )
    )
    row = next(item for item in rows if item["id"] == "OL-B4")
    row["prediction_event"] = None
    row["target"] = "A successor observation contract"
    state = custody.validate(lineages, frozen, architectures, rows, root=repo)
    assert state["qualifying_events_by_row"]["OL-B4"] == []
    assert state["events"][0]["id"] == "FZE-001"


def test_deleted_then_reintroduced_lineage_register_fails_history_guard(
    tmp_path: Path,
) -> None:
    repo, frozen, architectures, rows = _two_commit_event_repo(tmp_path)
    path = repo / "claims/frozen_prediction_architecture_lineages.json"
    payload = path.read_bytes()
    path.unlink()
    _git(repo, "add", "claims/frozen_prediction_architecture_lineages.json")
    _git(repo, "commit", "-q", "-m", "delete lineage custody")
    path.write_bytes(payload)
    _git(repo, "add", "claims/frozen_prediction_architecture_lineages.json")
    _git(repo, "commit", "-q", "-m", "reintroduce lineage custody")
    lineages = json.loads(path.read_text(encoding="utf-8"))
    with pytest.raises(SystemExit, match="disappeared and later reappeared"):
        custody.validate(lineages, frozen, architectures, rows, root=repo)


def test_history_scan_uses_only_path_change_commits(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[tuple[str, ...]] = []

    def fake_git_text(root: Path, *args: str, check: bool = True) -> str:
        calls.append(args)
        return ""

    monkeypatch.setattr(custody, "_git_text", fake_git_text)
    custody._validate_committed_history(
        {"events": [], "event_anchors": []},
        root=Path("/unused"),
        lineage_relative_path="claims/frozen_prediction_architecture_lineages.json",
    )
    assert calls == [
        (
            "log",
            "--first-parent",
            "--full-history",
            "--reverse",
            "--format=%H",
            "HEAD",
            "--",
            "claims/frozen_prediction_architecture_lineages.json",
        )
    ]
