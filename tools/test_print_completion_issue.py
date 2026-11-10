"""Tests for the canonical completion-issue printer."""

import build_completion_program_current as program
import print_completion_issue


def packet():
    return program.load_program()[0]


def test_dependency_checkbox_matches_predecessor_status():
    contract = packet()
    by_id = {stage["id"]: stage for stage in contract["stages"]}
    for stage in contract["stages"]:
        body = print_completion_issue.issue_body(stage, contract)
        for dependency in stage["requires"]:
            predecessor = by_id[dependency]
            mark = "x" if predecessor["status"] == "ACCEPTED" else " "
            line = (f"- [{mark}] {program.issue_link(predecessor['issue'])}: "
                    f"{predecessor['title']}")
            assert line in body
            opposite = f"- [{' ' if mark == 'x' else 'x'}] {program.issue_link(predecessor['issue'])}:"
            assert opposite not in body


def test_printer_only_changes_dependency_checkboxes():
    contract = packet()
    for stage in contract["stages"]:
        body = print_completion_issue.issue_body(stage, contract)
        assert body.replace("- [x] ", "- [ ] ") == program.render_issue(stage, contract)
