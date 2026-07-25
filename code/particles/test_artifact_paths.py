from __future__ import annotations

import sys
from pathlib import Path

CODE_ROOT = Path(__file__).resolve().parents[1]
if str(CODE_ROOT) not in sys.path:
    sys.path.insert(0, str(CODE_ROOT))

from particles.artifact_paths import (
    REPO_ROOT,
    canonical_artifact_ref,
    canonicalize_artifact_paths,
)


def test_repo_path_becomes_repo_relative() -> None:
    path = REPO_ROOT / "code/particles/runs/flavor/receipt.json"
    assert canonical_artifact_ref(path) == "code/particles/runs/flavor/receipt.json"


def test_sibling_sim_path_becomes_stable_workspace_uri() -> None:
    path = REPO_ROOT.parent / "oph-physics-sim/runs/e1/receipt.json"
    assert (
        canonical_artifact_ref(path)
        == "oph-workspace://oph-physics-sim/runs/e1/receipt.json"
    )


def test_recursive_boundary_removes_cross_platform_home_paths() -> None:
    payload = {
        "mac": "/Users/alice/work/input.json",
        "linux": ["/home/alice/work/input.json"],
        "windows": r"C:\Users\Alice\work\input.json",
    }
    normalized = canonicalize_artifact_paths(payload)
    assert normalized == {
        "mac": "external-file://input.json",
        "linux": ["external-file://input.json"],
        "windows": "external-file://input.json",
    }


def test_relative_logical_path_is_unchanged() -> None:
    assert canonical_artifact_ref(Path("code/particles/runs/a.json")) == (
        "code/particles/runs/a.json"
    )
