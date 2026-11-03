"""Fail-closed fixtures for tools/check_docs_policy.py."""

from pathlib import Path

from tools import check_docs_policy

ROOT = Path(__file__).resolve().parents[1]


def _write(root: Path, relative: str, content: str) -> None:
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def _rules(root: Path) -> set[str]:
    failures, _ = check_docs_policy.check(root)
    return {violation.rule for violation in failures}


def test_dated_file_name_fails(tmp_path: Path) -> None:
    _write(tmp_path, "docs/research/OPH_AUDIT_2026-11-02.md", "text\n")
    assert "dated-file-name" in _rules(tmp_path)


def test_undated_file_name_passes(tmp_path: Path) -> None:
    _write(tmp_path, "docs/research/OPH_AUDIT.md", "text\n")
    assert _rules(tmp_path) == set()


def test_dated_heading_fails(tmp_path: Path) -> None:
    _write(tmp_path, "docs/research/note.md", "# Note\n\n## Space and clocks - 2 November 2026\n")
    assert "dated-heading" in _rules(tmp_path)
    _write(tmp_path, "docs/research/note.md", "# Note\n\n## 2026-11-02 dispatch\n")
    assert "dated-heading" in _rules(tmp_path)


def test_undated_heading_passes(tmp_path: Path) -> None:
    _write(tmp_path, "docs/research/note.md", "# Note\n\n## Space and clocks\n")
    assert _rules(tmp_path) == set()


def test_undeclared_root_document_fails(tmp_path: Path) -> None:
    _write(tmp_path, "docs/PROPOSAL.md", "# Proposal\n")
    assert "root-document" in _rules(tmp_path)


def test_declared_root_document_passes(tmp_path: Path) -> None:
    _write(tmp_path, "docs/COMPLETION_PROGRAM.md", "# Completion program\n")
    assert _rules(tmp_path) == set()


def test_pull_request_reference_fails(tmp_path: Path) -> None:
    _write(tmp_path, "docs/research/note.md", "# Note\n\nMerged in PR #507.\n")
    assert "pull-request-reference" in _rules(tmp_path)
    _write(tmp_path, "docs/research/note.md", "# Note\n\nSee /pull/507.\n")
    assert "pull-request-reference" in _rules(tmp_path)


def test_status_prose_fails(tmp_path: Path) -> None:
    _write(tmp_path, "docs/research/note.md", "# Note\n\n**Audit date:** 2026-11-02\n")
    assert "status-prose" in _rules(tmp_path)
    _write(tmp_path, "docs/research/note.md", "# Note\n\nThe baseline commit is pinned.\n")
    assert "status-prose" in _rules(tmp_path)
    _write(tmp_path, "docs/research/note.md", "# Note\n\ngenerated_utc = 2026-11-02\n")
    assert "status-prose" in _rules(tmp_path)


def test_chronology_heading_fails(tmp_path: Path) -> None:
    _write(tmp_path, "docs/research/note.md", "# Note\n\n## 2. Chronology\n")
    assert "chronology-heading" in _rules(tmp_path)


def test_sequence_heading_passes(tmp_path: Path) -> None:
    _write(tmp_path, "docs/research/note.md", "# Note\n\n## Sequence of checks\n")
    assert _rules(tmp_path) == set()


def test_dated_directory_name_fails(tmp_path: Path) -> None:
    _write(tmp_path, "docs/research/host-clock-2026-11-02/README.md", "# Host clock\n")
    assert "dated-directory-name" in _rules(tmp_path)


def test_public_data_audit_is_exempt(tmp_path: Path) -> None:
    _write(
        tmp_path,
        "docs/research/PUBLIC_DATA_PREDICTION_AUDIT.md",
        "# Audit\n\n## 2026-11-02 comparison\n\nMerged in PR #507.\n",
    )
    assert _rules(tmp_path) == set()


def test_current_tree_passes() -> None:
    failures, exempted = check_docs_policy.check(ROOT)
    assert failures == []
    assert all(violation.path.startswith("docs/") for violation in exempted)


def test_scoped_check_ignores_untouched_violation(tmp_path: Path) -> None:
    _write(tmp_path, "docs/research/dated-2026-11-02.md", "text\n")
    _write(tmp_path, "docs/research/clean.md", "# Clean\n")
    failures, _ = check_docs_policy.check(
        tmp_path, [tmp_path / "docs/research/clean.md"]
    )
    assert failures == []
    failures, _ = check_docs_policy.check(
        tmp_path, [tmp_path / "docs/research/dated-2026-11-02.md"]
    )
    assert {violation.rule for violation in failures} == {"dated-file-name"}


def test_changed_documents_filters_to_the_docs_tree(
    tmp_path: Path, monkeypatch
) -> None:
    _write(tmp_path, "docs/research/note.md", "# Note\n")
    _write(tmp_path, "README.md", "# Root\n")

    class Result:
        returncode = 0
        stdout = "docs/research/note.md\nREADME.md\n"
        stderr = ""

    monkeypatch.setattr(
        check_docs_policy.subprocess, "run", lambda *arguments, **options: Result()
    )
    assert check_docs_policy.changed_documents(tmp_path, "base") == [
        (tmp_path / "docs/research/note.md").resolve()
    ]
