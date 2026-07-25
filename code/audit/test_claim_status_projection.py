"""Fail-closed tests for the generated #512 cross-surface status projection."""

import importlib.util
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
BUILDER = REPO_ROOT / "tools" / "build_scoreboard.py"

spec = importlib.util.spec_from_file_location("build_scoreboard", BUILDER)
scoreboard = importlib.util.module_from_spec(spec)
sys.modules["build_scoreboard"] = scoreboard
spec.loader.exec_module(scoreboard)


def fixture_sources():
    registry = {
        "release_id": "r-test",
        "claims": [
            {
                "claim_id": "FIX-1",
                "claim_class": "conditional_implication",
                "status": "proved_on_declared_branch",
                "gates": [42],
            }
        ],
    }
    snapshot = {
        "repo": "example/project",
        "rows": [{"number": 42, "title": "open gate"}],
    }
    return registry, snapshot


@pytest.mark.parametrize(
    ("surface", "mutation"),
    [
        ("claim class", lambda r, s: r["claims"][0].update(claim_class="branch_entry")),
        ("descriptive status", lambda r, s: r["claims"][0].update(status="open")),
        ("claim gate", lambda r, s: r["claims"][0].update(gates=[])),
        ("live issue state", lambda r, s: s.update(rows=[])),
        ("release", lambda r, s: r.update(release_id="r-next")),
    ],
)
def test_projection_digest_binds_every_status_input(surface, mutation):
    registry, snapshot = fixture_sources()
    before = scoreboard.projection_digest(
        scoreboard.canonical_projection_payload(registry, snapshot)
    )
    mutation(registry, snapshot)
    after = scoreboard.projection_digest(
        scoreboard.canonical_projection_payload(registry, snapshot)
    )
    assert after != before, surface


def test_projection_replacement_is_exact_and_idempotent():
    registry, snapshot = fixture_sources()
    projection = scoreboard.render_projection(registry, snapshot)
    seed = (
        "# Surface\n\n"
        f"{scoreboard.PROJECTION_START}\nold\n{scoreboard.PROJECTION_END}\n"
    )
    updated = scoreboard.replace_projection(seed, projection, Path("surface.md"))
    assert projection in updated
    assert scoreboard.replace_projection(
        updated, projection, Path("surface.md")
    ) == updated


@pytest.mark.parametrize(
    "bad",
    [
        "# no markers\n",
        (
            f"{scoreboard.PROJECTION_START}\n{scoreboard.PROJECTION_END}\n"
            f"{scoreboard.PROJECTION_START}\n{scoreboard.PROJECTION_END}\n"
        ),
    ],
)
def test_missing_or_duplicate_projection_markers_fail_closed(bad):
    registry, snapshot = fixture_sources()
    projection = scoreboard.render_projection(registry, snapshot)
    with pytest.raises(SystemExit, match="exactly one generated claim-status block"):
        scoreboard.replace_projection(bad, projection, Path("surface.md"))


def test_proof_spine_and_compression_scorecard_carry_identical_live_projection():
    registry, snapshot = scoreboard.source_documents()
    projection = scoreboard.render_projection(registry, snapshot)
    for relative in ["docs/PROOF_SPINE.md", "docs/COMPRESSION_SCORECARD.md"]:
        path = REPO_ROOT / relative
        current = path.read_text(encoding="utf-8")
        assert scoreboard.replace_projection(
            current, projection, Path(relative)
        ) == current
