from __future__ import annotations

import copy
import json

import build_open_problem_ledger as ledger_tool
import pytest


def _ledger() -> dict:
    return {
        "artifact": "oph_open_problem_ledger",
        "generated_utc": "2026-07-25T00:00:00Z",
        "repo": ledger_tool.REPO,
        "open_issue_count": 1,
        "closed_out_of_scope_count": 0,
        "worker_policy": {
            "chrome_pro_workers_default": "local_first",
            "max_parallel_workers": 6,
            "launch_condition": "fixture",
            "obstruction_only_result_allowed": False,
        },
        "closed_out_of_scope_records": [],
        "rows": [
            {
                "number": 7,
                "title": "fixture",
                "url": "https://example.invalid/7",
                "labels": [],
                "updated_at": "2026-07-25T00:00:00Z",
                "phase": "fixture",
                "claim_level": "open",
                "blocker": "fixture",
                "closure": "fixture",
                "falsification": "fixture",
                "chrome_policy": "fixture",
            }
        ],
    }


def _write(tmp_path, payload: dict, markdown: str | None = None):
    json_path = tmp_path / "ledger.json"
    md_path = tmp_path / "OPEN_PROBLEMS.md"
    json_path.write_text(json.dumps(payload), encoding="utf-8")
    md_path.write_text(
        markdown
        if markdown is not None
        else ledger_tool.render_markdown(payload) + "\n",
        encoding="utf-8",
    )
    return json_path, md_path


def test_offline_ledger_check_accepts_synchronized_artifacts(tmp_path) -> None:
    paths = _write(tmp_path, _ledger())
    assert ledger_tool.validate_committed_ledger(*paths) == []


def test_offline_ledger_check_rejects_count_and_markdown_drift(tmp_path) -> None:
    payload = copy.deepcopy(_ledger())
    payload["open_issue_count"] = 2
    paths = _write(tmp_path, payload, markdown="stale\n")
    problems = ledger_tool.validate_committed_ledger(*paths)
    assert "open_issue_count does not equal the number of rows" in problems
    assert "OPEN_PROBLEMS.md is out of sync with the ledger JSON" in problems


def test_offline_ledger_check_rejects_duplicate_rows(tmp_path) -> None:
    payload = copy.deepcopy(_ledger())
    payload["rows"].append(copy.deepcopy(payload["rows"][0]))
    payload["open_issue_count"] = 2
    paths = _write(tmp_path, payload)
    problems = ledger_tool.validate_committed_ledger(*paths)
    assert "ledger has duplicate issue numbers" in problems


def test_live_check_ignores_generation_clock_but_not_issue_state(
    tmp_path,
    monkeypatch,
) -> None:
    live_issues = [
        {
            "number": 7,
            "title": "fixture",
            "url": "https://example.invalid/7",
            "labels": [],
            "updatedAt": "2026-07-25T00:00:00Z",
        }
    ]
    monkeypatch.setattr(
        ledger_tool,
        "ISSUE_POLICY",
        {
            7: {
                "phase": "fixture",
                "claim_level": "open",
                "blocker": "fixture",
                "closure": "fixture",
                "falsification": "fixture",
                "chrome_policy": "fixture",
            }
        },
    )
    monkeypatch.setattr(ledger_tool, "CLOSED_OUT_OF_SCOPE_ISSUES", {})
    committed = ledger_tool.build_ledger(live_issues)
    committed["generated_utc"] = "1999-01-01T00:00:00Z"
    json_path, _ = _write(tmp_path, committed)

    assert ledger_tool.compare_committed_to_live(json_path, live_issues) == []

    changed_live = [copy.deepcopy(live_issues[0])]
    changed_live[0]["number"] = 8
    problems = ledger_tool.compare_committed_to_live(json_path, changed_live)
    assert any("open issue membership differs" in problem for problem in problems)


def _v2_issue(
    *,
    number: int = 700,
    title: str = "[E3] Construct the limit",
    labels: tuple[str, ...] = ("track:covariant-net", "size:L", "parked"),
    body: str | None = None,
) -> dict:
    return {
        "number": number,
        "title": title,
        "url": f"https://example.invalid/{number}",
        "labels": [{"name": label} for label in labels],
        "updatedAt": "2026-08-04T00:00:00Z",
        "body": body
        or (
            "Construct the declared limit object without promoting a finite "
            "example to a continuum theorem.\n\n"
            "Deliverables: QFT/Limit.lean and one typed countermodel matrix. "
            "Boundary: finite examples alone do not close the limit. "
            "Depends on: E2 (#693) and A4 (#699). Wave: V2-W3."
        ),
        "milestone": {"title": "V2-W3 Covariant net and gravity"},
    }


def test_v2_policy_comes_from_live_contract_without_scientific_promotion() -> None:
    row = ledger_tool.build_ledger([_v2_issue()])["rows"][0]

    assert row["phase"] == "v2-w3-covariant-net"
    assert row["claim_level"] == "parked stage-gated construction"
    assert "E2 (#693) and A4 (#699)" in row["blocker"]
    assert "Construct the declared limit object" in row["blocker"]
    assert "QFT/Limit.lean" in row["closure"]
    assert "finite examples alone do not close the limit" in row["falsification"]
    assert "Scientific falsifiers" not in row["falsification"]
    assert "Do not launch workers while parked" in row["chrome_policy"]


def test_v2_standing_policy_preserves_custody_semantics() -> None:
    issue = _v2_issue(
        number=695,
        title="[G1] Custody and comparisons",
        labels=("track:custody", "standing"),
        body=(
            "Maintain frozen-prediction custody and sealed comparisons; immutable "
            "registers stay immutable.\n\n"
            "Depends on: none. Standing lane across all waves. Wave: standing."
        ),
    )
    row = ledger_tool.build_ledger([issue])["rows"][0]

    assert row["phase"] == "v2-standing-custody"
    assert row["claim_level"] == "standing custody/comparison control"
    assert "issue state alone is not a scientific verdict" in row["closure"]
    assert "never use a worker to manufacture comparison evidence" in row["chrome_policy"]


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ({"body": "Construct it. Depends on: none."}, "no Wave: contract"),
        (
            {"labels": [{"name": "track:gravity"}, {"name": "size:L"}]},
            "task code E3 conflicts",
        ),
        (
            {"milestone": {"title": "V2-W2 Geometry and time"}},
            "Wave: contract conflicts with milestone",
        ),
        (
            {"body": "Construct it. Depends on: none. Wave: V2-W3 surprise."},
            "malformed Wave: contract",
        ),
    ],
)
def test_v2_contract_validation_fails_closed(mutation, message) -> None:
    issue = _v2_issue()
    issue.update(mutation)
    with pytest.raises(ValueError, match=message):
        ledger_tool.build_ledger([issue])


def test_offline_check_rejects_placeholder_v2_metadata(tmp_path) -> None:
    payload = _ledger()
    row = payload["rows"][0]
    row["labels"] = ["track:foundations", "size:S"]
    row["phase"] = "unclassified"
    row["blocker"] = "Classify blocker from the live issue body."
    row["closure"] = "Add exact closure criterion to this ledger."
    row["falsification"] = "Add exact falsification criterion to this ledger."
    paths = _write(tmp_path, payload)

    problems = ledger_tool.validate_committed_ledger(*paths)
    assert "row 0 has an unclassified V2 phase" in problems
    assert "row 0 retains placeholder V2 blocker" in problems
    assert "row 0 retains placeholder V2 closure" in problems
    assert "row 0 retains placeholder V2 falsification" in problems
