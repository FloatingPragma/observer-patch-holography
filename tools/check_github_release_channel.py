#!/usr/bin/env python3
"""Compare the public GitHub Release with the local release bundle (#514).

The paper manifest validates the local source-to-PDF bundle.  This checker
closes the separate publication-channel boundary: the latest public GitHub
Release must use the same release ID and expose byte-identical copies of every
manifest PDF, the book PDF, and the manifest itself.

The default mode reads the public GitHub API.  Tests and offline audits can
provide captured API payloads with ``--release-json`` and ``--latest-json``.
No command in this module creates, edits, or uploads a release.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path, PurePosixPath
from typing import Any, Mapping


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REPO = "FloatingPragma/observer-patch-holography"
DEFAULT_MANIFEST_RELATIVE = Path("paper/paper_release_manifest.json")
BOOK_RELATIVE = Path("book/reverse-engineering-reality-book.pdf")
# Public GitHub Releases intentionally expose only the flagship plus the
# top-level paper/ stack. Adjunct manifest sections support other surfaces.
MANIFEST_SECTIONS = ("papers",)
SHA256_RE = re.compile(r"[0-9a-f]{64}")
GIT_OBJECT_RE = re.compile(r"[0-9a-fA-F]{40,64}")
REPO_SLUG_RE = re.compile(
    r"[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+"
)


class ReleaseChannelError(RuntimeError):
    """Raised when the public release cannot be validated."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ReleaseChannelError(f"cannot read JSON payload {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise ReleaseChannelError(f"JSON payload root must be an object: {path}")
    return payload


def confined_file(
    *,
    repo_root: Path,
    path: Path,
    label: str,
) -> Path:
    """Return one regular file whose resolved path stays below ``repo_root``."""

    try:
        root = repo_root.expanduser().resolve(strict=True)
    except (OSError, RuntimeError) as exc:
        raise ReleaseChannelError(
            f"cannot resolve repository root {repo_root}: {exc}"
        ) from exc
    if not root.is_dir():
        raise ReleaseChannelError(f"repository root is not a directory: {root}")

    candidate = path.expanduser()
    if not candidate.is_absolute():
        candidate = root / candidate
    try:
        resolved = candidate.resolve(strict=True)
        resolved.relative_to(root)
    except (OSError, RuntimeError, ValueError) as exc:
        raise ReleaseChannelError(
            f"{label} must resolve to a file inside the repository: {path}"
        ) from exc
    if not resolved.is_file() or candidate.is_symlink():
        raise ReleaseChannelError(
            f"{label} must be a regular non-symlink file inside the repository: "
            f"{resolved.relative_to(root).as_posix()}"
        )
    return resolved


def manifest_asset_path(
    *,
    repo_root: Path,
    pdf_path: str,
    where: str,
) -> Path:
    """Resolve a canonical POSIX manifest path without platform ambiguity."""

    if "\\" in pdf_path:
        raise ReleaseChannelError(
            f"{where}.pdf_path must use normalized POSIX separators"
        )
    relative = PurePosixPath(pdf_path)
    if (
        relative.is_absolute()
        or not relative.parts
        or any(part in {"", ".", ".."} for part in relative.parts)
        or relative.as_posix() != pdf_path
    ):
        raise ReleaseChannelError(
            f"{where}.pdf_path must be a normalized repository-relative POSIX path"
        )
    return confined_file(
        repo_root=repo_root,
        path=Path(*relative.parts),
        label=f"{where}.pdf_path",
    )


def validated_manifest_asset(
    *,
    repo_root: Path,
    record: Mapping[str, Any],
    where: str,
) -> Path:
    """Validate one manifest receipt and return its confined local artifact."""

    pdf_path = record.get("pdf_path")
    if not isinstance(pdf_path, str) or not pdf_path:
        raise ReleaseChannelError(f"{where} has no pdf_path")
    path = manifest_asset_path(
        repo_root=repo_root,
        pdf_path=pdf_path,
        where=where,
    )
    declared_sha = record.get("sha256")
    if not isinstance(declared_sha, str) or SHA256_RE.fullmatch(declared_sha) is None:
        raise ReleaseChannelError(
            f"{where}.sha256 must be a lowercase 64-digit SHA-256"
        )
    declared_size = record.get("size_bytes")
    if (
        isinstance(declared_size, bool)
        or not isinstance(declared_size, int)
        or declared_size < 0
    ):
        raise ReleaseChannelError(
            f"{where}.size_bytes must be a nonnegative integer"
        )
    try:
        actual_sha = sha256_file(path)
        actual_size = path.stat().st_size
    except OSError as exc:
        raise ReleaseChannelError(
            f"cannot inspect release asset {pdf_path}: {exc}"
        ) from exc
    if declared_sha != actual_sha:
        raise ReleaseChannelError(
            f"{where}: manifest sha256 {declared_sha!r} does not match "
            f"local {actual_sha!r}"
        )
    if declared_size != actual_size:
        raise ReleaseChannelError(
            f"{where}: manifest size_bytes {declared_size!r} does not "
            f"match local {actual_size!r}"
        )
    return path


def expected_assets(
    *,
    repo_root: Path,
    manifest_path: Path,
) -> tuple[str, dict[str, dict[str, int | str]]]:
    """Return the release ID and exact public asset contract."""
    try:
        normalized_root = repo_root.expanduser().resolve(strict=True)
    except (OSError, RuntimeError) as exc:
        raise ReleaseChannelError(
            f"cannot resolve repository root {repo_root}: {exc}"
        ) from exc
    if not normalized_root.is_dir():
        raise ReleaseChannelError(
            f"repository root is not a directory: {normalized_root}"
        )
    normalized_manifest = confined_file(
        repo_root=normalized_root,
        path=manifest_path,
        label="paper release manifest",
    )
    manifest = load_json(normalized_manifest)
    release_id = manifest.get("release_id")
    if not isinstance(release_id, str) or not release_id.strip():
        raise ReleaseChannelError("paper release manifest has no release_id")

    paths: list[Path] = []
    for section_name in MANIFEST_SECTIONS:
        section = manifest.get(section_name)
        if not isinstance(section, dict):
            raise ReleaseChannelError(
                f"paper release manifest section {section_name!r} must be an object"
            )
        for paper_id, record in section.items():
            if not isinstance(record, dict):
                raise ReleaseChannelError(
                    f"{section_name}.{paper_id} must be an object"
                )
            where = f"{section_name}.{paper_id}"
            path = validated_manifest_asset(
                repo_root=normalized_root,
                record=record,
                where=where,
            )
            paths.append(path)

    book = manifest.get("book")
    if not isinstance(book, dict):
        raise ReleaseChannelError(
            "paper release manifest has no canonical book receipt"
        )
    if book.get("built_for_release_id") != release_id:
        raise ReleaseChannelError(
            "book.built_for_release_id must match the manifest release_id"
        )
    if book.get("pdf_path") != BOOK_RELATIVE.as_posix():
        raise ReleaseChannelError(
            f"book.pdf_path must be {BOOK_RELATIVE.as_posix()!r}"
        )
    paths.extend(
        (
            validated_manifest_asset(
                repo_root=normalized_root,
                record=book,
                where="book",
            ),
            normalized_manifest,
        )
    )
    names = [path.name for path in paths]
    if len(names) != len({name.casefold() for name in names}):
        raise ReleaseChannelError(
            "release assets must have portable case-insensitive unique basenames "
            "for GitHub upload"
        )

    contract: dict[str, dict[str, int | str]] = {}
    for path in paths:
        try:
            contract[path.name] = {
                "sha256": sha256_file(path),
                "size": path.stat().st_size,
            }
        except OSError as exc:
            raise ReleaseChannelError(
                f"cannot inspect release asset {path.name}: {exc}"
            ) from exc
    return release_id, contract


def _asset_map(payload: Mapping[str, Any]) -> dict[str, Mapping[str, Any]]:
    assets = payload.get("assets")
    if not isinstance(assets, list):
        raise ReleaseChannelError("GitHub release payload has no assets array")
    result: dict[str, Mapping[str, Any]] = {}
    for index, asset in enumerate(assets):
        if not isinstance(asset, dict):
            raise ReleaseChannelError(f"GitHub release asset {index} is not an object")
        name = asset.get("name")
        if not isinstance(name, str) or not name:
            raise ReleaseChannelError(f"GitHub release asset {index} has no name")
        if name in result:
            raise ReleaseChannelError(f"duplicate GitHub release asset name: {name}")
        result[name] = asset
    return result


def _git_object(payload: Mapping[str, Any], *, where: str) -> tuple[str, str]:
    value = payload.get("object")
    if not isinstance(value, dict):
        raise ReleaseChannelError(f"{where} has no Git object")
    object_type = value.get("type")
    sha = value.get("sha")
    if object_type not in {"commit", "tag"}:
        raise ReleaseChannelError(
            f"{where} has unsupported Git object type {object_type!r}"
        )
    if not isinstance(sha, str) or GIT_OBJECT_RE.fullmatch(sha) is None:
        raise ReleaseChannelError(f"{where} has an invalid Git object SHA")
    return object_type, sha.lower()


def resolve_tag_commit(
    *,
    ref_payload: Mapping[str, Any],
    release_id: str,
    fetch_tag_object,
) -> str:
    """Peel a lightweight or annotated Git tag to its immutable commit."""

    expected_ref = f"refs/tags/{release_id}"
    if ref_payload.get("ref") != expected_ref:
        raise ReleaseChannelError(
            f"Git tag ref differs from the manifest: "
            f"{ref_payload.get('ref')!r} != {expected_ref!r}"
        )
    object_type, sha = _git_object(ref_payload, where=expected_ref)
    if object_type == "commit":
        return sha

    tag_payload = fetch_tag_object(sha)
    if not isinstance(tag_payload, dict):
        raise ReleaseChannelError("annotated Git tag payload root must be an object")
    if tag_payload.get("tag") != release_id:
        raise ReleaseChannelError(
            f"annotated Git tag name differs from the manifest: "
            f"{tag_payload.get('tag')!r} != {release_id!r}"
        )
    peeled_type, peeled_sha = _git_object(
        tag_payload,
        where=f"annotated tag {release_id}",
    )
    if peeled_type != "commit":
        raise ReleaseChannelError(
            f"annotated tag {release_id} does not point directly to a commit"
        )
    return peeled_sha


def captured_tag_commit(
    payload: Mapping[str, Any],
    *,
    release_id: str,
) -> str:
    """Read the normalized tag capture used by offline audits."""

    if payload.get("tag_name") != release_id:
        raise ReleaseChannelError(
            f"captured tag name differs from the manifest: "
            f"{payload.get('tag_name')!r} != {release_id!r}"
        )
    commit = payload.get("commit_sha")
    if not isinstance(commit, str) or GIT_OBJECT_RE.fullmatch(commit) is None:
        raise ReleaseChannelError("captured tag payload has no valid commit_sha")
    return commit.lower()


def validate_release_payloads(
    *,
    release_payload: Mapping[str, Any],
    latest_payload: Mapping[str, Any],
    release_id: str,
    expected: Mapping[str, Mapping[str, int | str]],
    tag_commit: str | None = None,
    expected_commit: str | None = None,
) -> list[str]:
    """Return every release-channel mismatch without mutating remote state."""
    problems: list[str] = []
    if release_payload.get("tag_name") != release_id:
        problems.append(
            "requested release tag differs from the manifest: "
            f"{release_payload.get('tag_name')!r} != {release_id!r}"
        )
    if latest_payload.get("tag_name") != release_id:
        problems.append(
            "latest public GitHub release differs from the manifest: "
            f"{latest_payload.get('tag_name')!r} != {release_id!r}"
        )
    if release_payload.get("draft") is not False:
        problems.append("manifest release is absent or remains a draft")
    if release_payload.get("prerelease") is not False:
        problems.append("manifest release is marked as a prerelease")
    if expected_commit is not None:
        normalized_expected = expected_commit.lower()
        if GIT_OBJECT_RE.fullmatch(normalized_expected) is None:
            problems.append("expected commit is not a valid Git object SHA")
        elif tag_commit is None:
            problems.append("public release tag commit was not resolved")
        elif tag_commit.lower() != normalized_expected:
            problems.append(
                "public release tag commit differs from the requested commit: "
                f"{tag_commit.lower()} != {normalized_expected}"
            )

    try:
        actual = _asset_map(release_payload)
    except ReleaseChannelError as exc:
        return problems + [str(exc)]

    expected_names = set(expected)
    actual_names = set(actual)
    missing = sorted(expected_names - actual_names)
    extra = sorted(actual_names - expected_names)
    if missing:
        problems.append(f"public release is missing assets: {missing}")
    if extra:
        problems.append(f"public release has unexpected assets: {extra}")

    for name in sorted(expected_names & actual_names):
        contract = expected[name]
        asset = actual[name]
        expected_digest = f"sha256:{contract['sha256']}"
        if asset.get("digest") != expected_digest:
            problems.append(
                f"{name}: public digest {asset.get('digest')!r} "
                f"does not match local {expected_digest!r}"
            )
        if asset.get("size") != contract["size"]:
            problems.append(
                f"{name}: public size {asset.get('size')!r} "
                f"does not match local {contract['size']!r}"
            )
    return problems


def github_json(url: str, *, token: str | None = None) -> dict[str, Any]:
    headers = {
        "Accept": "application/vnd.github+json",
        "User-Agent": "oph-release-channel-integrity-check",
        "X-GitHub-Api-Version": "2022-11-28",
    }
    if token:
        headers["Authorization"] = f"Bearer {token}"
    request = urllib.request.Request(url, headers=headers)
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            payload = json.load(response)
    except (
        OSError,
        ValueError,
        UnicodeError,
        urllib.error.URLError,
        json.JSONDecodeError,
    ) as exc:
        raise ReleaseChannelError(f"GitHub API request failed for {url}: {exc}") from exc
    if not isinstance(payload, dict):
        raise ReleaseChannelError(f"GitHub API returned a non-object payload for {url}")
    return payload


def public_tag_commit(
    *,
    api_base: str,
    release_id: str,
    token: str | None,
) -> str:
    """Resolve the public Git tag without trusting release.target_commitish.

    GitHub can report a branch name such as ``main`` in ``target_commitish``
    even when an existing annotated release tag points to an older immutable
    commit. The Git Data ref and annotated-tag objects are the binding source.
    """

    encoded_release = urllib.parse.quote(release_id, safe="")
    ref_payload = github_json(
        f"{api_base}/git/ref/tags/{encoded_release}",
        token=token,
    )

    def fetch_tag_object(sha: str) -> Mapping[str, Any]:
        return github_json(f"{api_base}/git/tags/{sha}", token=token)

    return resolve_tag_commit(
        ref_payload=ref_payload,
        release_id=release_id,
        fetch_tag_object=fetch_tag_object,
    )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=ROOT)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--repo", default=DEFAULT_REPO)
    parser.add_argument("--release-json", type=Path)
    parser.add_argument("--latest-json", type=Path)
    parser.add_argument(
        "--tag-json",
        type=Path,
        help=(
            "offline normalized tag capture with tag_name and commit_sha; "
            "requires --expected-commit"
        ),
    )
    parser.add_argument(
        "--expected-commit",
        help=(
            "require the peeled public Git tag to point to this explicit commit; "
            "release.target_commitish is never used as a commit binding"
        ),
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        try:
            repo_root = args.repo_root.expanduser().resolve(strict=True)
        except (OSError, RuntimeError) as exc:
            raise ReleaseChannelError(
                f"cannot resolve repository root {args.repo_root}: {exc}"
            ) from exc
        if not repo_root.is_dir():
            raise ReleaseChannelError(
                f"repository root is not a directory: {repo_root}"
            )
        manifest_path = confined_file(
            repo_root=repo_root,
            path=args.manifest or DEFAULT_MANIFEST_RELATIVE,
            label="paper release manifest",
        )
        if REPO_SLUG_RE.fullmatch(args.repo) is None:
            raise ReleaseChannelError(
                "repository must use the exact owner/name form"
            )
        release_id, contract = expected_assets(
            repo_root=repo_root,
            manifest_path=manifest_path,
        )
        fixture_mode = args.release_json is not None or args.latest_json is not None
        if fixture_mode:
            if args.release_json is None or args.latest_json is None:
                raise ReleaseChannelError(
                    "--release-json and --latest-json must be supplied together"
                )
            release_payload = load_json(args.release_json)
            latest_payload = load_json(args.latest_json)
        else:
            api_base = f"https://api.github.com/repos/{args.repo}"
            token = os.environ.get("GITHUB_TOKEN") or os.environ.get("GH_TOKEN")
            encoded_release = urllib.parse.quote(release_id, safe="")
            release_payload = github_json(
                f"{api_base}/releases/tags/{encoded_release}",
                token=token,
            )
            latest_payload = github_json(
                f"{api_base}/releases/latest",
                token=token,
            )

        tag_commit: str | None = None
        if args.expected_commit is not None:
            if GIT_OBJECT_RE.fullmatch(args.expected_commit) is None:
                raise ReleaseChannelError(
                    "--expected-commit must be a 40- or 64-digit Git object SHA"
                )
            if fixture_mode:
                if args.tag_json is None:
                    raise ReleaseChannelError(
                        "--tag-json is required with fixture payloads and "
                        "--expected-commit"
                    )
                tag_commit = captured_tag_commit(
                    load_json(args.tag_json),
                    release_id=release_id,
                )
            else:
                tag_commit = public_tag_commit(
                    api_base=api_base,
                    release_id=release_id,
                    token=token,
                )
        elif args.tag_json is not None:
            raise ReleaseChannelError(
                "--tag-json is meaningful only with --expected-commit"
            )
        problems = validate_release_payloads(
            release_payload=release_payload,
            latest_payload=latest_payload,
            release_id=release_id,
            expected=contract,
            tag_commit=tag_commit,
            expected_commit=args.expected_commit,
        )
    except ReleaseChannelError as exc:
        problems = [str(exc)]

    if problems:
        print("GitHub release channel FAILED:")
        for problem in problems:
            print(f"  - {problem}")
        return 1

    print(
        f"GitHub release channel OK: {release_id}, "
        f"{len(contract)} byte-exact assets"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
