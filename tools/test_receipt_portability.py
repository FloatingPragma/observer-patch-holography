from __future__ import annotations

import json
from pathlib import Path

import pytest

from check_receipt_portability import check, find_violations


@pytest.mark.parametrize(
    "leaked",
    [
        "/Users/alice/work/oph/receipt.json",
        "/home/alice/work/oph/receipt.json",
        r"C:\Users\Alice\work\oph\receipt.json",
    ],
)
def test_developer_home_paths_fail_closed(tmp_path: Path, leaked: str) -> None:
    receipt = tmp_path / "receipt.json"
    receipt.write_text(
        json.dumps({"artifact": "mutation", "source": {"path": leaked}}),
        encoding="utf-8",
    )
    violations = find_violations([receipt])
    assert len(violations) == 1
    assert violations[0].json_pointer == "/source/path"
    with pytest.raises(ValueError, match="developer-home path leaked"):
        check([receipt])


def test_clone_stable_references_pass(tmp_path: Path) -> None:
    receipt = tmp_path / "receipt.json"
    receipt.write_text(
        json.dumps(
            {
                "repo_relative": "code/particles/runs/flavor/receipt.json",
                "sibling": "oph-workspace://oph-physics-sim/runs/e1/receipt.json",
                "external": "external-file://measurement-pack.json",
            }
        ),
        encoding="utf-8",
    )
    check([receipt])
