"""Drift, provenance, and reverse-edge gates for the V3 discharge queue."""

from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))

import build_premise_discharge_queue as discharge_queue


def _sources() -> dict:
    return discharge_queue.load_sources()


def _queue() -> tuple[dict, dict]:
    sources = _sources()
    return discharge_queue.build_queue(sources), sources


def _items_by_id(queue: dict) -> dict[str, dict]:
    return {item["premise_id"]: item for item in queue["items"]}


def test_schema_and_exact_per_premise_coverage() -> None:
    queue, sources = _queue()
    assert queue["schema"] == "oph.premise_discharge_queue.v2"
    expected = {
        row["id"]
        for row in sources["rows"]
        if row["disposition"] in discharge_queue.QUEUED_DISPOSITIONS
    }
    actual = [item["premise_id"] for item in queue["items"]]
    assert len(actual) == 53
    assert len(actual) == len(set(actual))
    assert set(actual) == expected
    assert {row["premise_id"] for row in queue["excluded_imports"]} == {
        row["id"] for row in sources["rows"] if row["disposition"] == "import"
    }


def test_items_preserve_exact_av0_rows_and_separate_evaluation_version() -> None:
    queue, sources = _queue()
    source = {row["id"]: row for row in sources["rows"]}
    for item in queue["items"]:
        row = source[item["premise_id"]]
        assert item["premise_statement"] == row["statement"]
        assert item["next_action_or_decision"] == row["notes"]
        assert item["origin_premise_record_sha256"] == discharge_queue._canonical_sha256(row)
        assert item["origin_architecture_version"] == "AV-0"
        assert item["evaluation_architecture_version"] == sources["current_version"]

    successor_sources = copy.deepcopy(sources)
    successor_sources["current_version"] = "AV-1"
    successor = discharge_queue.build_queue(successor_sources)
    assert {item["origin_architecture_version"] for item in successor["items"]} == {"AV-0"}
    assert {item["evaluation_architecture_version"] for item in successor["items"]} == {"AV-1"}


def test_reverse_observation_edges_are_bidirectionally_exact() -> None:
    queue, sources = _queue()
    actual: set[tuple[str, str, str]] = set()
    for item in queue["items"]:
        for edge in item["consumer_edges"]["observation_rows"]:
            actual.add((item["premise_id"], edge["id"], edge["role"]))

    queued_ids = {item["premise_id"] for item in queue["items"]}
    expected = {
        (premise_id, row["id"], role)
        for row in sources["observation_rows"]
        for field, role in (
            ("premises", "consumed_premise"),
            ("open_premises", "open_premise"),
        )
        for premise_id in row[field]
        if premise_id in queued_ids
    }
    assert actual == expected


def test_reverse_claim_edges_are_explicit_and_exact() -> None:
    queue, sources = _queue()
    items = _items_by_id(queue)
    actual = {
        (item["premise_id"], claim_id)
        for item in queue["items"]
        for claim_id in item["consumer_edges"]["claim_ids"]
    }
    expected = {
        (premise_id, claim["claim_id"])
        for claim in sources["claims"]
        for premise_id in set(
            discharge_queue.PREMISE_TOKEN_RE.findall(json.dumps(claim, sort_keys=True))
        )
        if premise_id in items
    }
    assert actual == expected
    assert (
        "OPH-GEOMETRY-INVERSE-SQUARE-SHELL"
        in items["PR-18"]["consumer_edges"]["claim_ids"]
    )


def test_reverse_architecture_and_non_observation_surface_edges_are_exact() -> None:
    queue, sources = _queue()
    items = _items_by_id(queue)
    actual_decisions = {
        (
            item["premise_id"],
            edge["architecture_version"],
            edge["decision_id"],
        )
        for item in queue["items"]
        for edge in item["consumer_edges"]["architecture_decisions"]
    }
    expected_decisions = {
        (premise_id, version["id"], decision["id"])
        for version in sources["architecture"]["versions"]
        for decision in version["protocol_decisions"]
        for premise_id in decision["premises"]
        if premise_id in items
    }
    assert actual_decisions == expected_decisions

    actual_surfaces = {
        (item["premise_id"], edge["lane_issue"], edge["surface"])
        for item in queue["items"]
        for edge in item["consumer_edges"]["non_observation_surfaces"]
    }
    expected_surfaces = {
        (premise_id, lane, surface)
        for (lane, premise_id), surface in discharge_queue.observation_ledger.NON_OBSERVATION_SURFACE_CONSUMERS.items()
        if premise_id in items
    }
    assert actual_surfaces == expected_surfaces


def test_evidence_reference_edges_are_exact_and_not_semantic_inferences() -> None:
    queue, sources = _queue()
    for item in queue["items"]:
        source_row = next(row for row in sources["rows"] if row["id"] == item["premise_id"])
        expected = []
        expected_unreferenced = []
        for path in source_row["evidence"]:
            observation_ids = sorted(
                row["id"]
                for row in sources["observation_rows"]
                if path in row["evidence"]
            )
            claim_ids = sorted(
                claim["claim_id"]
                for claim in sources["claims"]
                if path in claim["evidence"]
            )
            if observation_ids or claim_ids:
                expected.append(
                    {
                        "path": path,
                        "observation_row_ids": observation_ids,
                        "claim_ids": claim_ids,
                        "role": "shared_evidence_reference_not_semantic_premise_inference",
                    }
                )
            else:
                expected_unreferenced.append(path)
        assert item["consumer_edges"]["evidence_references"] == expected
        assert (
            item["consumer_edges"]["unreferenced_evidence_paths"]
            == expected_unreferenced
        )


def test_reverse_edge_coverage_reports_without_inventing_consumers() -> None:
    queue, sources = _queue()
    coverage = queue["reverse_edge_coverage"]
    expected_unannotated = sorted(
        claim["claim_id"]
        for claim in sources["claims"]
        if not discharge_queue.PREMISE_TOKEN_RE.findall(
            json.dumps(claim, sort_keys=True)
        )
    )
    expected_unreferenced = [
        {
            "premise_id": item["premise_id"],
            "paths": item["consumer_edges"]["unreferenced_evidence_paths"],
        }
        for item in queue["items"]
        if item["consumer_edges"]["unreferenced_evidence_paths"]
    ]
    assert coverage["claim_records_without_explicit_premise_ids"] == expected_unannotated
    assert coverage["evidence_paths_without_registered_reverse_reference"] == expected_unreferenced
    assert "No queue item may claim claim-exhaustive replay" in coverage["closure_guard"]
    assert all(
        reference["role"]
        == "shared_evidence_reference_not_semantic_premise_inference"
        for item in queue["items"]
        for reference in item["consumer_edges"]["evidence_references"]
    )


def test_origin_evidence_is_path_resolved_and_commit_pinned() -> None:
    queue, sources = _queue()
    for item in queue["items"]:
        for pin in item["origin_evidence"]:
            assert (discharge_queue.ROOT / pin["path"]).exists()
            assert pin["origin_revision"] == sources["origin_revision"]
            assert pin == discharge_queue._evidence_pin(
                sources["origin_revision"], pin["path"]
            )
            assert pin["git_object_type"] in {"blob", "tree"}
            assert pin["descendant_count"] >= 1


def test_input_custody_pins_all_reverse_edge_sources_and_architecture_snapshot() -> None:
    queue, sources = _queue()
    custody = queue["input_custody"]
    assert custody == sources["input_custody"]
    assert custody["premise_register"]["origin_revision"] == sources["origin_revision"]
    assert custody["origin_architecture_snapshot"] == {
        "id": "AV-0",
        "origin_revision": sources["origin_revision"],
        "snapshot_sha256": sources["origin_version"]["snapshot_sha256"],
        "record_sha256": sources["origin_anchor"]["record_sha256"],
    }
    assert set(custody["reverse_edge_inputs"]) == {
        "observation_ledger",
        "claim_registry",
        "non_observation_surface_map",
    }
    for key in ("premise_register", "architecture_register", "open_issue_snapshot"):
        assert custody[key]["byte_count"] > 0
        assert discharge_queue.DIGEST_RE.fullmatch(custody[key]["sha256"])
        assert discharge_queue.OBJECT_RE.fullmatch(custody[key]["git_blob_sha1"])


def test_exit_paths_follow_disposition_and_use_durable_decision_custody() -> None:
    queue, _ = _queue()
    assert queue["audit_custody_issue"] == 738
    assert queue["issue"] == 739
    assert queue["decision_custody"] == discharge_queue.DECISION_CUSTODY
    assert queue["decision_custody"]["architecture_register"] == (
        "tracking/architecture_versions.json"
    )
    assert queue["decision_custody"]["audit_register"] == (
        "tracking/audit_custody.json"
    )
    assert queue["decision_custody"]["bootstrap_issue"] == 741
    assert "live state is not decision semantics" in queue["policy"]
    for item in queue["items"]:
        expected_paths = (
            ["derive", "axiomatize"]
            if item["disposition"] == "axiomatize"
            else ["derive"]
        )
        assert item["available_exit_paths"] == expected_paths
        assert item["decision_custody"] == discharge_queue.DECISION_CUSTODY
        assert item["state"] == "deferred_open"


def test_triage_order_is_lane_fanout_first_then_premise_id() -> None:
    queue, _ = _queue()
    observed = [
        (
            -len(item["consumer_edges"]["lane_issues"]),
            int(item["premise_id"].split("-", 1)[1]),
        )
        for item in queue["items"]
    ]
    assert observed == sorted(observed)
    assert [item["triage_rank"] for item in queue["items"]] == list(
        range(1, len(queue["items"]) + 1)
    )


@pytest.mark.parametrize(
    ("field", "mutate"),
    [
        ("origin_evidence", lambda item: item["origin_evidence"].pop()),
        (
            "origin_evidence",
            lambda item: item["origin_evidence"][0].update(content_sha256="0" * 64),
        ),
        ("consumer_edges", lambda item: item["consumer_edges"]["lane_issues"].pop()),
        (
            "origin_architecture_version",
            lambda item: item.update(origin_architecture_version="AV-1"),
        ),
        (
            "evaluation_architecture_version",
            lambda item: item.update(evaluation_architecture_version="AV-1"),
        ),
        (
            "next_action_or_decision",
            lambda item: item.update(next_action_or_decision="done"),
        ),
        ("available_exit_paths", lambda item: item.update(available_exit_paths=[])),
        ("state", lambda item: item.update(state="discharged_derived")),
    ],
)
def test_per_premise_custody_mutations_fail_closed(field, mutate) -> None:
    queue, sources = _queue()
    mutant = copy.deepcopy(queue)
    mutate(mutant["items"][0])
    with pytest.raises(SystemExit, match=rf"{field} has drifted"):
        discharge_queue.validate_queue(mutant, sources)


@pytest.mark.parametrize(
    ("field", "mutate"),
    [
        ("policy", lambda queue: queue.update(policy="generation is discharge")),
        (
            "path_contracts",
            lambda queue: queue["path_contracts"]["derive"]["closure_requirements"].clear(),
        ),
        (
            "inherited_v2_residuals",
            lambda queue: queue["inherited_v2_residuals"].pop(),
        ),
        (
            "input_custody",
            lambda queue: queue["input_custody"]["architecture_register"].update(
                sha256="0" * 64
            ),
        ),
        (
            "input_custody",
            lambda queue: queue["input_custody"]["open_issue_snapshot"].update(
                byte_count=0
            ),
        ),
    ],
)
def test_program_level_contract_mutations_fail_closed(field, mutate) -> None:
    queue, sources = _queue()
    mutant = copy.deepcopy(queue)
    mutate(mutant)
    with pytest.raises(SystemExit, match=rf"{field} has drifted"):
        discharge_queue.validate_queue(mutant, sources)


@pytest.mark.parametrize("owner", [738, 739])
def test_missing_live_queue_or_standing_owner_fails_closed(owner: int) -> None:
    sources = copy.deepcopy(_sources())
    sources["open_issues"].remove(owner)
    with pytest.raises(SystemExit, match="owners must both be live"):
        discharge_queue.build_queue(sources)


def test_closed_architecture_bootstrap_issue_does_not_break_queue() -> None:
    sources = copy.deepcopy(_sources())
    sources["open_issues"].discard(741)
    queue = discharge_queue.build_queue(sources)
    assert queue["decision_custody"]["bootstrap_issue"] == 741
    assert {
        item["decision_custody"]["bootstrap_issue"] for item in queue["items"]
    } == {741}


def test_live_origin_register_rewrite_is_rejected(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    current = discharge_queue.premise_register.load_json(
        discharge_queue.PREMISE_REGISTER_PATH
    )
    mutant = copy.deepcopy(current)
    mutant["rows"][0]["statement"] += " Mutated."
    monkeypatch.setattr(
        discharge_queue.premise_register,
        "load_json",
        lambda _: copy.deepcopy(mutant),
    )
    with pytest.raises(SystemExit, match="differs from the anchored AV-0 origin"):
        discharge_queue.load_sources()


def test_stale_open_issue_count_fails_closed() -> None:
    with pytest.raises(SystemExit, match="count does not match"):
        discharge_queue._issue_numbers(
            {
                "open_issue_count": 4,
                "rows": [
                    {"number": 738},
                    {"number": 739},
                    {"number": 741},
                ],
            }
        )


def test_committed_queue_and_markdown_are_current() -> None:
    queue, sources = _queue()
    committed = discharge_queue.load_json(discharge_queue.QUEUE_PATH)
    discharge_queue.validate_queue(committed, sources)
    assert committed == queue
    assert (
        discharge_queue.render(queue).encode("utf-8")
        == discharge_queue.SURFACE_PATH.read_bytes()
    )


def test_check_mode_rejects_stale_surface(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    stale_queue = tmp_path / "premise_discharge_queue.json"
    stale_surface = tmp_path / "PREMISE_DISCHARGE_QUEUE_V3.md"
    stale_queue.write_text("{}\n", encoding="utf-8")
    stale_surface.write_text("stale\n", encoding="utf-8")
    monkeypatch.setattr(discharge_queue, "QUEUE_PATH", stale_queue)
    monkeypatch.setattr(discharge_queue, "SURFACE_PATH", stale_surface)
    monkeypatch.setattr(sys, "argv", ["build_premise_discharge_queue.py", "--check"])
    assert discharge_queue.main() == 1


def test_duplicate_json_key_rejected(tmp_path: Path) -> None:
    duplicate = tmp_path / "duplicate.json"
    duplicate.write_text('{"schema":"one","schema":"two"}', encoding="utf-8")
    with pytest.raises(SystemExit, match="duplicate JSON key"):
        discharge_queue.load_json(duplicate)
