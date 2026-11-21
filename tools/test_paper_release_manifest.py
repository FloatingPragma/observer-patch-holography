#!/usr/bin/env python3
"""Regression tests for tools/validate_paper_release_manifest.py (issue #514).

Run with: python -m pytest tools/test_paper_release_manifest.py
"""
from __future__ import annotations

import copy
import json
import os
import re
import shutil
import subprocess
import textwrap
from pathlib import Path, PureWindowsPath

import pytest

import bump_paper_release as bumper
import generate_paper_release_manifest as generator
import refresh_paper_release as refresher
import validate_paper_release_manifest as validator

REPO_ROOT = Path(__file__).resolve().parent.parent
MANIFEST = REPO_ROOT / "paper" / "paper_release_manifest.json"


def _base() -> dict:
    return json.loads(MANIFEST.read_text(encoding="utf-8"))


def _write(tmp_path: Path, manifest: dict) -> Path:
    path = tmp_path / "manifest.json"
    path.write_text(json.dumps(manifest), encoding="utf-8")
    return path


def test_committed_manifest_matches_source_set(tmp_path: Path) -> None:
    # The manifest membership is derived from build_tex_papers.py, so the
    # committed manifest must validate with no problems.
    assert validator.validate(MANIFEST) == []


def test_expected_sets_are_derived_not_fixed() -> None:
    sections = validator.expected_sections()
    assert set(sections["papers"]) == set(validator.source.RELEASE_TRACKED)
    assert set(sections["extra_papers"]) == set(
        validator.source.RELEASED_ADJUNCT_PAPERS
    )
    assert set(sections["supplemental_papers"]) == set(validator.source.PAPERS) - set(validator.source.RELEASE_TRACKED)


def test_manifest_generation_timestamp_is_release_derived() -> None:
    assert (
        generator.deterministic_generated_at("July 25, 2026")
        == "2026-07-25T00:00:00Z"
    )


def test_manifest_paths_are_posix_on_windows() -> None:
    assert (
        generator.manifest_path(PureWindowsPath(r"paper\example.pdf"))
        == "paper/example.pdf"
    )


@pytest.mark.parametrize(
    "message",
    [
        "response code 429 Too Many Requests",
        'failed to retrieve "ec-lmtt8.tfm" from the network',
        "connection reset by peer",
        "operation timed out",
    ],
)
def test_tectonic_transient_fetch_failures_are_retryable(message: str) -> None:
    result = subprocess.CompletedProcess(
        ["tectonic"],
        1,
        stdout="",
        stderr=message,
    )
    assert refresher.paper_sources.transient_fetch_failure(result)


def test_tectonic_tex_errors_are_not_retryable() -> None:
    result = subprocess.CompletedProcess(
        ["tectonic"],
        1,
        stdout="",
        stderr="undefined control sequence",
    )
    assert not refresher.paper_sources.transient_fetch_failure(result)


def test_manifest_binds_canonical_book_bytes_and_release() -> None:
    manifest = _base()
    book = manifest["book"]
    book_path = REPO_ROOT / book["pdf_path"]

    assert book["pdf_path"] == "book/reverse-engineering-reality-book.pdf"
    assert book["built_for_release_id"] == manifest["release_id"]
    assert book["sha256"] == generator.sha256(book_path)
    assert book["size_bytes"] == book_path.stat().st_size


def test_generator_carries_book_builder_receipt_across_release_bump() -> None:
    previous = _base()
    carried = generator.carried_book_manifest_entry(
        REPO_ROOT,
        previous,
        "r-next",
    )

    assert carried == previous["book"]
    assert carried["built_for_release_id"] == previous["release_id"]


def test_book_builder_updates_manifest_receipt(tmp_path: Path) -> None:
    repo = tmp_path / "repo"
    (repo / "paper").mkdir(parents=True)
    (repo / "book").mkdir()
    (repo / "paper/release_info.tex").write_text(
        r"\newcommand{\OPHPaperReleaseID}{r-next}" + "\n",
        encoding="utf-8",
    )
    (repo / "paper/paper_release_manifest.json").write_text(
        json.dumps({"release_id": "r-next", "book": {}}),
        encoding="utf-8",
    )
    book_path = repo / generator.BOOK_PDF_RELATIVE
    book_path.write_bytes(b"%PDF-rebuilt-book")

    generator.update_book_manifest_entry(repo)

    manifest = json.loads(
        (repo / "paper/paper_release_manifest.json").read_text(encoding="utf-8")
    )
    assert manifest["book"] == {
        "built_for_release_id": "r-next",
        "pdf_path": "book/reverse-engineering-reality-book.pdf",
        "sha256": generator.sha256(book_path),
        "size_bytes": book_path.stat().st_size,
    }


def test_generator_rejects_changed_pdf_behind_existing_tag() -> None:
    tagged = _base()
    current = copy.deepcopy(tagged)
    paper_id = next(iter(current["papers"]))
    current["papers"][paper_id]["sha256"] = "0" * 64

    with pytest.raises(SystemExit, match="immutable Git tag.*--preview"):
        generator.enforce_tag_immutability(tagged, current)


def test_generator_preview_accepts_changed_pdf_behind_existing_tag() -> None:
    tagged = _base()
    current = copy.deepcopy(tagged)
    paper_id = next(iter(current["papers"]))
    current["papers"][paper_id]["sha256"] = "0" * 64

    generator.enforce_generation_policy(
        tagged,
        tagged,
        current,
        preview=True,
    )


def test_generator_publication_candidate_rejects_changed_tagged_pdf() -> None:
    tagged = _base()
    current = copy.deepcopy(tagged)
    paper_id = next(iter(current["papers"]))
    current["papers"][paper_id]["sha256"] = "0" * 64

    with pytest.raises(SystemExit, match="immutable Git tag"):
        generator.enforce_generation_policy(
            tagged,
            tagged,
            current,
            preview=False,
        )


def test_generator_rejects_removed_pdf_behind_existing_tag() -> None:
    tagged = _base()
    current = copy.deepcopy(tagged)
    removed = next(iter(current["papers"]))
    current["papers"].pop(removed)

    with pytest.raises(SystemExit, match="immutable Git tag") as exc:
        generator.enforce_tag_immutability(tagged, current)

    assert f"papers.{removed}" in str(exc.value)


def test_generator_requires_release_bump_for_removed_pdf() -> None:
    previous = _base()
    current = copy.deepcopy(previous)
    removed = next(iter(current["papers"]))
    current["papers"].pop(removed)

    with pytest.raises(SystemExit, match="manifest entries changed") as exc:
        generator.enforce_release_bump(previous, current)

    assert f"papers.{removed}" in str(exc.value)


def test_generator_accepts_unchanged_pdf_behind_existing_tag() -> None:
    tagged = _base()
    generator.enforce_tag_immutability(tagged, copy.deepcopy(tagged))


def test_generator_accepts_new_untagged_release() -> None:
    current = _base()
    current["release_id"] = "r-next"
    generator.enforce_tag_immutability(_base(), current)


def test_publication_validation_rejects_existing_local_release_id(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[list[str]] = []

    def fake_run(argv, **kwargs):
        calls.append(argv)

        class Result:
            returncode = 0
            stdout = "tag-object\n"
            stderr = ""

        return Result()

    monkeypatch.setattr(validator.subprocess, "run", fake_run)
    problems: list[str] = []
    validator.check_publication_release_id(
        REPO_ROOT,
        "r-test",
        problems,
        remote="origin",
    )

    assert any("requires a new release ID" in problem for problem in problems)
    assert calls == [["git", "rev-parse", "--verify", "--quiet", "refs/tags/r-test"]]


def test_publication_validation_accepts_unused_local_and_remote_release_id(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    results = iter((1, 2))

    def fake_run(argv, **kwargs):
        class Result:
            returncode = next(results)
            stdout = ""
            stderr = ""

        return Result()

    monkeypatch.setattr(validator.subprocess, "run", fake_run)
    problems: list[str] = []
    validator.check_publication_release_id(
        REPO_ROOT,
        "r-next",
        problems,
        remote="origin",
    )

    assert problems == []


def test_publication_validation_rejects_existing_remote_release_id(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[list[str]] = []

    def fake_run(argv, **kwargs):
        calls.append(argv)

        class Result:
            returncode = 1 if argv[1] == "rev-parse" else 0
            stdout = "" if returncode else "tag-object\n"
            stderr = ""

        return Result()

    monkeypatch.setattr(validator.subprocess, "run", fake_run)
    problems: list[str] = []
    validator.check_publication_release_id(
        REPO_ROOT,
        "r-test",
        problems,
        remote="origin",
    )

    assert any("exists on 'origin'" in problem for problem in problems)
    assert [call[1] for call in calls] == ["rev-parse", "ls-remote"]


def test_publication_validation_fails_closed_on_remote_error(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fake_run(argv, **kwargs):
        class Result:
            returncode = 1 if argv[1] == "rev-parse" else 128
            stdout = ""
            stderr = "" if returncode == 1 else "network unavailable"

        return Result()

    monkeypatch.setattr(validator.subprocess, "run", fake_run)
    problems: list[str] = []
    validator.check_publication_release_id(
        REPO_ROOT,
        "r-test",
        problems,
        remote="origin",
    )

    assert problems == [
        "could not verify that release ID 'r-test' is unused on "
        "'origin': network unavailable"
    ]


def test_manifest_generator_has_no_legacy_same_release_alias() -> None:
    with pytest.raises(SystemExit):
        generator.parse_args(["--allow-same-release"])


def test_pdf_release_check_can_be_skipped_only_for_preview() -> None:
    with pytest.raises(SystemExit):
        generator.parse_args(["--skip-pdf-release-check"])

    args = generator.parse_args(["--preview", "--skip-pdf-release-check"])
    assert args.preview is True
    assert args.skip_pdf_release_check is True


def test_rejects_missing_paper(tmp_path: Path) -> None:
    manifest = _base()
    removed = next(iter(manifest["papers"]))
    manifest["papers"].pop(removed)
    problems = validator.validate(_write(tmp_path, manifest))
    assert any("missing" in p and removed in p for p in problems)


def test_rejects_missing_or_stale_book_receipt(tmp_path: Path) -> None:
    manifest = _base()
    manifest.pop("book")
    problems = validator.validate(_write(tmp_path, manifest))
    assert any("book receipt is missing" in problem for problem in problems)

    manifest = _base()
    manifest["book"]["sha256"] = "0" * 64
    manifest["book"]["built_for_release_id"] = "r-old"
    problems = validator.validate(_write(tmp_path, manifest))
    assert any("book: sha256 mismatch" in problem for problem in problems)
    assert any("built_for_release_id" in problem for problem in problems)


def test_rejects_unexpected_paper(tmp_path: Path) -> None:
    manifest = _base()
    manifest["papers"]["not_a_source_paper"] = {"pdf_path": "paper/x.pdf", "sha256": "x", "size_bytes": 1}
    problems = validator.validate(_write(tmp_path, manifest))
    assert any("unexpected" in p and "not_a_source_paper" in p for p in problems)


def test_rejects_absent_artifact(tmp_path: Path) -> None:
    manifest = _base()
    paper_id = next(iter(manifest["papers"]))
    manifest["papers"][paper_id] = {"pdf_path": "paper/DOES_NOT_EXIST.pdf", "sha256": "x", "size_bytes": 1}
    problems = validator.validate(_write(tmp_path, manifest))
    assert any("missing on disk" in p for p in problems)


def test_rejects_cross_paper_artifact_mapping(tmp_path: Path) -> None:
    """A valid digest for the wrong paper cannot satisfy source-derived membership."""
    manifest = _base()
    left, right = list(manifest["papers"])[:2]
    manifest["papers"][left] = copy.deepcopy(manifest["papers"][right])
    problems = validator.validate(_write(tmp_path, manifest))
    assert any(left in p and "source derives" in p for p in problems)


def test_rejects_sha256_mismatch(tmp_path: Path) -> None:
    # A listed PDF whose content differs from its declared digest (silent rebuild /
    # swap / tamper) must be rejected even when the path exists.
    manifest = _base()
    paper_id = next(iter(manifest["papers"]))
    manifest["papers"][paper_id]["sha256"] = "0" * 64
    problems = validator.validate(_write(tmp_path, manifest))
    assert any("sha256 mismatch" in p and paper_id in p for p in problems)


def test_rejects_size_mismatch(tmp_path: Path) -> None:
    # A truncated / regrown artifact whose byte count differs from the manifest is rejected.
    manifest = _base()
    paper_id = next(iter(manifest["papers"]))
    manifest["papers"][paper_id]["size_bytes"] = "1"
    problems = validator.validate(_write(tmp_path, manifest))
    assert any("size_bytes mismatch" in p and paper_id in p for p in problems)


def test_generator_allows_repository_only_extra_pdf(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Focused extra PDFs remain in the repository outside the release set."""
    paper_dir = tmp_path / "paper"
    extra_dir = tmp_path / "extra"
    paper_dir.mkdir()
    extra_dir.mkdir()
    expected_pdf = paper_dir / "expected.pdf"
    expected_pdf.write_bytes(b"%PDF expected")
    manifest = {
        "papers": {
            "expected": {
                "pdf_path": "paper/expected.pdf",
                "sha256": "unused-by-this-check",
                "size_bytes": expected_pdf.stat().st_size,
            }
        },
        "supplemental_papers": {},
        "extra_papers": {},
    }
    monkeypatch.setattr(generator, "NON_TEX_SOURCE_PDFS", {})

    generator.verify_no_stray_pdfs(tmp_path, manifest)
    (extra_dir / "not_implied_by_any_source.pdf").write_bytes(b"%PDF stray")

    generator.verify_no_stray_pdfs(tmp_path, manifest)


def test_generator_rejects_missing_registered_non_tex_source(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    (tmp_path / "paper").mkdir()
    (tmp_path / "extra").mkdir()
    monkeypatch.setattr(
        generator,
        "NON_TEX_SOURCE_PDFS",
        {
            Path("extra/derived.pdf"): Path("extra/derived/build.sh"),
        },
    )

    with pytest.raises(SystemExit, match="registered source.*missing"):
        generator.verify_no_stray_pdfs(
            tmp_path,
            {
                "papers": {},
                "supplemental_papers": {},
                "extra_papers": {},
            },
        )


def test_generator_rejects_missing_registered_non_tex_output(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    (tmp_path / "paper").mkdir()
    source = tmp_path / "extra" / "derived" / "build.sh"
    source.parent.mkdir(parents=True)
    source.write_text("#!/bin/sh\n", encoding="utf-8")
    monkeypatch.setattr(
        generator,
        "NON_TEX_SOURCE_PDFS",
        {
            Path("extra/derived.pdf"): Path("extra/derived/build.sh"),
        },
    )

    with pytest.raises(SystemExit, match="registered non-TeX output is missing"):
        generator.verify_no_stray_pdfs(
            tmp_path,
            {
                "papers": {},
                "supplemental_papers": {},
                "extra_papers": {},
            },
        )


def test_validator_rejects_missing_registered_non_tex_output(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    (tmp_path / "paper").mkdir()
    source_path = tmp_path / "extra" / "derived" / "build.sh"
    source_path.parent.mkdir(parents=True)
    source_path.write_text("#!/bin/sh\n", encoding="utf-8")
    monkeypatch.setattr(validator, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(
        validator.source,
        "NON_TEX_SOURCE_PDFS",
        {
            Path("extra/derived.pdf"): Path("extra/derived/build.sh"),
        },
    )
    problems: list[str] = []
    validator.check_release_surface(
        {
            "papers": {},
            "supplemental_papers": {},
            "extra_papers": {},
        },
        problems,
    )
    assert any("registered non-TeX output is missing" in p for p in problems)


@pytest.mark.parametrize(
    ("case", "diagnostic"),
    [
        ("matching", None),
        ("untracked", "Generated PDF is not committed"),
        ("ignored", "Generated PDF is not committed"),
        ("untracked_whitespace", "Generated PDF is not committed"),
        ("changed_tracked", "The committed PDFs or manifest do not match"),
        ("stale_book", "The committed book PDF does not match"),
    ],
)
def test_actual_preview_artifact_step_rejects_missing_or_stale_commits(
    tmp_path: Path, case: str, diagnostic: str | None
) -> None:
    """Execute the workflow's own Bash block, without rebuilding any real PDF."""
    workflow = (
        REPO_ROOT / ".github" / "workflows" / "ci.yml"
    ).read_text(encoding="utf-8")
    marker = (
        "      - name: Require committed PDFs to match the source rebuild\n"
        "        run: |\n"
    )
    assert workflow.count(marker) == 1
    script = textwrap.dedent(
        re.split(r"\n      - name:|\n\n  \S", workflow.split(marker, 1)[1], maxsplit=1)[0]
    )
    env = os.environ.copy()
    bash = shutil.which("bash")
    if os.name == "nt":
        # Execute this Ubuntu workflow regression with Git Bash, never WSL Bash.
        git = shutil.which("git")
        assert git is not None, "Git is required for the preview fixture"
        git_root = Path(git).resolve().parent.parent
        candidates = [git_root / "bin/bash.exe", git_root / "usr/bin/bash.exe"]
        bash_path = next((p for p in candidates if p.is_file()), None)
        if bash_path is None:
            pytest.skip("Git Bash is required to replay the Ubuntu workflow on Windows")
        bash = str(bash_path)
        env["PATH"] = os.pathsep.join(
            [str(git_root / "usr/bin"), str(git_root / "bin"), env["PATH"]]
        )
    assert bash is not None, "Bash is required to replay the preview workflow"

    repo = tmp_path / "repo"
    repo.mkdir()
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    env["RUNNER_TEMP"] = scratch.as_posix()

    def git_run(*args: str) -> None:
        subprocess.run(
            ["git", *args], cwd=repo, env=env, check=True,
            capture_output=True, text=True, encoding="utf-8",
        )

    git_run("init", "--quiet")
    git_run("config", "user.name", "Preview fixture")
    git_run("config", "user.email", "preview@example.invalid")
    git_run("config", "commit.gpgSign", "false")
    git_run("config", "core.autocrlf", "false")
    hooks = tmp_path / "empty-hooks"
    hooks.mkdir()
    git_run("config", "core.hooksPath", hooks.as_posix())
    for directory in ("paper", "flagship", "extra", "cosmology", "book"):
        (repo / directory).mkdir()
    for directory in ("paper", "flagship", "extra", "cosmology"):
        (repo / directory / "existing paper.pdf").write_bytes(b"%PDF-existing\n")
    book = repo / "book/reverse-engineering-reality-book.pdf"
    book.write_bytes(b"%PDF-book\n")
    (scratch / "book.pdf").write_bytes(book.read_bytes())
    if case == "ignored":
        (repo / ".gitignore").write_text("cosmology/generated.pdf\n", encoding="utf-8")
    if case.startswith("untracked") or case == "ignored":
        (repo / "cosmology/generated.tex").write_text("new source\n", encoding="utf-8")
    git_run("add", ".")
    git_run("commit", "--quiet", "-m", "committed preview fixture")

    if case in ("untracked", "ignored"):
        (repo / "cosmology/generated.pdf").write_bytes(b"%PDF-generated\n")
    elif case == "untracked_whitespace":
        (repo / "cosmology/generated paper.pdf").write_bytes(b"%PDF-generated\n")
    elif case == "changed_tracked":
        (repo / "cosmology/existing paper.pdf").write_bytes(b"%PDF-rebuilt\n")
    elif case == "stale_book":
        (scratch / "book.pdf").write_bytes(b"%PDF-current-book\n")

    result = subprocess.run(
        [bash, "--noprofile", "--norc", "-e", "-c", script],
        cwd=repo, env=env, capture_output=True, text=True, encoding="utf-8",
    )
    output = result.stdout + result.stderr
    assert result.returncode == (0 if diagnostic is None else 1), output
    if diagnostic is not None:
        assert diagnostic in output
    assert book.read_bytes() == b"%PDF-book\n", "the gate must not refresh the book"


def test_refresh_wrapper_dispatches_preview_by_default() -> None:
    args = refresher.parse_args([])
    generator_command, validator_command = refresher.manifest_commands(
        "python",
        publication=args.publication,
    )

    assert generator_command == [
        "python",
        "tools/generate_paper_release_manifest.py",
        "--preview",
    ]
    assert validator_command == [
        "python",
        "tools/validate_paper_release_manifest.py",
    ]


def test_refresh_wrapper_dispatches_explicit_preview() -> None:
    args = refresher.parse_args(["--preview"])
    generator_command, validator_command = refresher.manifest_commands(
        "python",
        publication=args.publication,
    )

    assert args.preview is True
    assert args.publication is False
    assert generator_command[-1] == "--preview"
    assert "--publication" not in validator_command


def test_refresh_wrapper_dispatches_strict_publication() -> None:
    args = refresher.parse_args(["--publication"])
    generator_command, validator_command = refresher.manifest_commands(
        "python",
        publication=args.publication,
    )

    assert args.preview is False
    assert args.publication is True
    assert "--preview" not in generator_command
    assert validator_command[-1] == "--publication"


def test_refresh_wrapper_builds_book_inside_publication_mode(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[tuple[str, list[str]]] = []
    monkeypatch.setattr(
        "sys.argv",
        ["refresh_paper_release.py", "--publication", "--no-gate"],
    )
    monkeypatch.setattr(
        refresher,
        "run",
        lambda description, argv: calls.append((description, argv)),
    )

    assert refresher.main() == 0
    descriptions = [description for description, _argv in calls]
    assert descriptions == [
        "build all registered papers",
        "regenerate the release manifest",
        "build the canonical book and stamp its release receipt",
        "validate the release manifest",
    ]
    assert calls[2][1][-1] == "tools/build_book_pdf.py"


def test_refresh_wrapper_modes_are_mutually_exclusive() -> None:
    with pytest.raises(SystemExit):
        refresher.parse_args(["--preview", "--publication"])


def test_bump_helper_points_to_strict_publication_wrapper() -> None:
    assert bumper.NEXT_STEP == (
        "Next: run python3 tools/refresh_paper_release.py --publication"
    )


def test_reproduce_marks_release_channel_integrity_check_as_publication_only() -> None:
    guide = (REPO_ROOT / "REPRODUCE.md").read_text(encoding="utf-8")
    normalized = " ".join(guide.split())
    assert (
        "The manually dispatched `Release Channel Integrity` workflow is a "
        "post-publication integrity check."
    ) in normalized
    assert "Do not use it to validate a same-release preview." in normalized


def test_content_inventory_names_only_tracked_sources() -> None:
    """Every path and paper the content inventory names has to resolve."""
    import yaml

    inventory = yaml.safe_load(
        (REPO_ROOT / "paper" / "content_inventory.yml").read_text(encoding="utf-8")
    )
    roots = ("paper", "extra", "cosmology", "flagship")
    for fragment, entry in inventory.get("shared_fragments", {}).items():
        assert (REPO_ROOT / fragment).is_file(), fragment
        relative = fragment.split("/", 1)[1]
        for paper in entry.get("included_by", []):
            source = next(
                (
                    REPO_ROOT / root / f"{paper}.tex"
                    for root in roots
                    if (REPO_ROOT / root / f"{paper}.tex").is_file()
                ),
                None,
            )
            assert source is not None, f"{fragment} included_by: {paper}"
            assert relative in source.read_text(encoding="utf-8"), (
                f"{paper} does not input {relative}"
            )
    for block in ("papers", "focused_papers"):
        for name in inventory.get(block, {}) or {}:
            assert any(
                (REPO_ROOT / root / f"{name}.tex").is_file() for root in roots
            ), f"{block}: {name}"
