import json
from pathlib import Path
import subprocess

import pytest

from report_grade_changes import CANDIDATES, LEDGER, grade_changes, grades


def write(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data))


def law(path, criteria):
    return {"candidate": path.parent.name,
            "requirements": [{"id": "section", "criteria": [{"id": k, "status": s} for k, s in criteria]}],
            "premises": [{"id": "PR-01", "status": "fulfilled"}],
            "joint": {"one-law": {"status": "open"}}}


def test_revision_to_worktree_handles_changes_additions_removals_and_reordering(tmp_path):
    def git(*args):
        return subprocess.run(["git", *args], cwd=tmp_path, check=True, capture_output=True)
    git("init", "-q")
    ledger = tmp_path / LEDGER
    a, b, c = (tmp_path / CANDIDATES / name / "EVALUATION.json" for name in ("a", "b", "c"))
    write(ledger, {"rows": [{"id": "r1", "status": "Open"}, {"id": "r2", "status": "Partial"}]})
    write(a, law(a, [("A", "partial"), ("B", "open"), ("C", "demonstrated")]))
    write(b, law(b, [("X", "open")]))
    git("add", ".")
    git("-c", "user.name=test", "-c", "user.email=test@example.test", "commit", "-qm", "baseline")
    assert grade_changes(tmp_path, "HEAD") == []
    write(ledger, {"rows": [{"id": "r3", "status": "Partial"}, {"id": "r1", "status": "Open", "what": "edited"}]})
    current = law(a, [("C", "open"), ("B", "open"), ("A", "demonstrated")])
    current["premises"][0]["status"] = "open"
    current["joint"]["one-law"]["status"] = "fulfilled"
    write(a, current)
    b.unlink()
    write(c, law(c, [("X", "partial")]))  # New local deposits may still be untracked.
    assert grade_changes(tmp_path, "HEAD") == [
        ("register", "r2", "Partial", None), ("register", "r3", None, "Partial"),
        ("a/criterion", "A", "partial", "demonstrated"), ("a/criterion", "C", "demonstrated", "open"),
        ("a/joint", "one-law", "open", "fulfilled"), ("a/premise", "PR-01", "fulfilled", "open"),
        ("b/criterion", "X", "open", None), ("b/joint", "one-law", "open", None),
        ("b/premise", "PR-01", "fulfilled", None),
        ("c/criterion", "X", None, "partial"), ("c/joint", "one-law", None, "open"),
        ("c/premise", "PR-01", None, "fulfilled"),
    ]
    with pytest.raises(subprocess.CalledProcessError):
        grade_changes(tmp_path, "nonexistent-base")


def test_duplicate_criterion_or_misnamed_law_is_an_error():
    path = CANDIDATES / "test" / "EVALUATION.json"
    data = law(path, [("OW-01", "open"), ("OW-01", "partial")])
    with pytest.raises(ValueError, match="duplicate status"):
        grades(path, json.dumps(data))
    with pytest.raises(ValueError, match="differs from directory"):
        grades(Path("model/candidate_laws/other/EVALUATION.json"), json.dumps(data))


def test_executable_checks_are_included_in_the_grade_diff():
    path = CANDIDATES / "test" / "EVALUATION.json"
    data = law(path, [("OW-01", "partial")])
    data["checks"] = {"hermitian": {"status": "pass", "detail": "zero residual"}}
    assert grades(path, json.dumps(data))["test/check", "hermitian"] == "pass"
