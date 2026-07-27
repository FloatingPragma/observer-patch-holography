"""Regression tests for the read-only GitHub release-channel checker."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

import check_github_release_channel as checker


def _write(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)


def _fixture(tmp_path: Path) -> tuple[Path, dict, dict]:
    paper = tmp_path / "paper" / "paper.pdf"
    book = tmp_path / checker.BOOK_RELATIVE
    _write(paper, b"%PDF-paper\n")
    _write(book, b"%PDF-book\n")
    manifest_path = tmp_path / checker.DEFAULT_MANIFEST_RELATIVE
    manifest = {
        "release_id": "r-test",
        "book": {
            "built_for_release_id": "r-test",
            "pdf_path": checker.BOOK_RELATIVE.as_posix(),
            "sha256": hashlib.sha256(book.read_bytes()).hexdigest(),
            "size_bytes": book.stat().st_size,
        },
        "papers": {
            "paper": {
                "pdf_path": "paper/paper.pdf",
                "sha256": hashlib.sha256(paper.read_bytes()).hexdigest(),
                "size_bytes": paper.stat().st_size,
            }
        },
        "supplemental_papers": {},
        "extra_papers": {},
    }
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    release_id, contract = checker.expected_assets(
        repo_root=tmp_path,
        manifest_path=manifest_path,
    )
    assets = [
        {
            "name": name,
            "digest": f"sha256:{record['sha256']}",
            "size": record["size"],
        }
        for name, record in sorted(contract.items())
    ]
    release = {
        "tag_name": release_id,
        "draft": False,
        "prerelease": False,
        "assets": assets,
    }
    latest = {"tag_name": release_id}
    return manifest_path, release, latest


def test_matching_release_is_accepted(tmp_path: Path) -> None:
    manifest_path, release, latest = _fixture(tmp_path)
    release_id, contract = checker.expected_assets(
        repo_root=tmp_path,
        manifest_path=manifest_path,
    )
    assert checker.validate_release_payloads(
        release_payload=release,
        latest_payload=latest,
        release_id=release_id,
        expected=contract,
    ) == []


def test_tampered_public_asset_is_rejected(tmp_path: Path) -> None:
    manifest_path, release, latest = _fixture(tmp_path)
    release["assets"][0]["digest"] = "sha256:" + ("0" * 64)
    release_id, contract = checker.expected_assets(
        repo_root=tmp_path,
        manifest_path=manifest_path,
    )
    problems = checker.validate_release_payloads(
        release_payload=release,
        latest_payload=latest,
        release_id=release_id,
        expected=contract,
    )
    assert any("public digest" in problem for problem in problems)


def test_missing_extra_and_nonlatest_release_are_rejected(tmp_path: Path) -> None:
    manifest_path, release, latest = _fixture(tmp_path)
    release["assets"].pop()
    release["assets"].append(
        {"name": "stray.pdf", "digest": "sha256:" + ("1" * 64), "size": 1}
    )
    latest["tag_name"] = "r-other"
    release_id, contract = checker.expected_assets(
        repo_root=tmp_path,
        manifest_path=manifest_path,
    )
    problems = checker.validate_release_payloads(
        release_payload=release,
        latest_payload=latest,
        release_id=release_id,
        expected=contract,
    )
    assert any("latest public GitHub release" in problem for problem in problems)
    assert any("missing assets" in problem for problem in problems)
    assert any("unexpected assets" in problem for problem in problems)


def test_duplicate_asset_names_fail_closed() -> None:
    release = {
        "tag_name": "r-test",
        "draft": False,
        "prerelease": False,
        "assets": [
            {"name": "paper.pdf", "digest": "sha256:" + ("0" * 64), "size": 1},
            {"name": "paper.pdf", "digest": "sha256:" + ("0" * 64), "size": 1},
        ],
    }
    problems = checker.validate_release_payloads(
        release_payload=release,
        latest_payload={"tag_name": "r-test"},
        release_id="r-test",
        expected={"paper.pdf": {"sha256": "0" * 64, "size": 1}},
    )
    assert problems == ["duplicate GitHub release asset name: paper.pdf"]


def test_tampered_local_manifest_digest_is_rejected(tmp_path: Path) -> None:
    manifest_path, _, _ = _fixture(tmp_path)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["papers"]["paper"]["sha256"] = "0" * 64
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    with pytest.raises(checker.ReleaseChannelError, match="does not match local"):
        checker.expected_assets(
            repo_root=tmp_path,
            manifest_path=manifest_path,
        )


def test_stale_book_receipt_is_rejected(tmp_path: Path) -> None:
    manifest_path, _, _ = _fixture(tmp_path)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["book"]["built_for_release_id"] = "r-old"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    with pytest.raises(checker.ReleaseChannelError, match="must match"):
        checker.expected_assets(
            repo_root=tmp_path,
            manifest_path=manifest_path,
        )


@pytest.mark.parametrize(
    "bad_path",
    [
        "../outside.pdf",
        "/tmp/outside.pdf",
        "paper/../outside.pdf",
        r"paper\paper.pdf",
    ],
)
def test_manifest_path_escape_and_nonportable_separators_are_rejected(
    tmp_path: Path,
    bad_path: str,
) -> None:
    manifest_path, _, _ = _fixture(tmp_path)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["papers"]["paper"]["pdf_path"] = bad_path
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    with pytest.raises(checker.ReleaseChannelError):
        checker.expected_assets(
            repo_root=tmp_path,
            manifest_path=manifest_path,
        )


def test_duplicate_portable_asset_basename_is_rejected(tmp_path: Path) -> None:
    manifest_path, _, _ = _fixture(tmp_path)
    duplicate = tmp_path / "extra" / "paper.pdf"
    _write(duplicate, b"%PDF-other\n")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["supplemental_papers"]["other"] = {
        "pdf_path": "extra/paper.pdf",
        "sha256": hashlib.sha256(duplicate.read_bytes()).hexdigest(),
        "size_bytes": duplicate.stat().st_size,
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    with pytest.raises(checker.ReleaseChannelError, match="unique basenames"):
        checker.expected_assets(
            repo_root=tmp_path,
            manifest_path=manifest_path,
        )


def test_lightweight_tag_resolves_to_immutable_commit() -> None:
    commit = "a" * 40

    def unexpected_fetch(_: str) -> dict:
        raise AssertionError("lightweight tag must not fetch an annotated object")

    assert (
        checker.resolve_tag_commit(
            ref_payload={
                "ref": "refs/tags/r-test",
                "object": {"type": "commit", "sha": commit},
            },
            release_id="r-test",
            fetch_tag_object=unexpected_fetch,
        )
        == commit
    )


def test_annotated_tag_is_peeled_and_name_bound() -> None:
    tag_object = "b" * 40
    commit = "c" * 40

    def fetch(sha: str) -> dict:
        assert sha == tag_object
        return {
            "tag": "r-test",
            "object": {"type": "commit", "sha": commit},
        }

    assert (
        checker.resolve_tag_commit(
            ref_payload={
                "ref": "refs/tags/r-test",
                "object": {"type": "tag", "sha": tag_object},
            },
            release_id="r-test",
            fetch_tag_object=fetch,
        )
        == commit
    )


def test_release_tag_commit_mismatch_is_rejected(tmp_path: Path) -> None:
    manifest_path, release, latest = _fixture(tmp_path)
    release_id, contract = checker.expected_assets(
        repo_root=tmp_path,
        manifest_path=manifest_path,
    )
    problems = checker.validate_release_payloads(
        release_payload=release,
        latest_payload=latest,
        release_id=release_id,
        expected=contract,
        tag_commit="a" * 40,
        expected_commit="b" * 40,
    )
    assert any("tag commit differs" in problem for problem in problems)


def test_captured_tag_name_and_commit_are_validated() -> None:
    with pytest.raises(checker.ReleaseChannelError, match="tag name differs"):
        checker.captured_tag_commit(
            {"tag_name": "r-other", "commit_sha": "a" * 40},
            release_id="r-test",
        )
    with pytest.raises(checker.ReleaseChannelError, match="valid commit_sha"):
        checker.captured_tag_commit(
            {"tag_name": "r-test", "commit_sha": "not-a-sha"},
            release_id="r-test",
        )
