from __future__ import annotations

import copy
import json

import build_open_problem_ledger as ledger_tool


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
