#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import build_tex_papers as paper_sources


RELEASE_INFO_RELATIVE = Path("paper/release_info.tex")
OUTPUT_RELATIVE = Path("paper/paper_release_manifest.json")
BOOK_PDF_RELATIVE = Path("book/reverse-engineering-reality-book.pdf")
RELEASED_COSMOLOGY_TEX = ()
RELEASE_TRACKED_PDFS = {
    "from_observer_consensus_to_standard_physics": Path(
        "flagship/from_observer_consensus_to_standard_physics.pdf"
    ),
    "deriving_the_particle_zoo_from_observer_consistency": Path(
        "paper/deriving_the_particle_zoo_from_observer_consistency.pdf"
    ),
    "observers_are_all_you_need": Path("paper/observers_are_all_you_need.pdf"),
    "paradise_as_fixed_point_consensus": Path("paper/paradise_as_fixed_point_consensus.pdf"),
    "reality_as_consensus_protocol": Path("paper/reality_as_consensus_protocol.pdf"),
    "recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency": Path(
        "paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.pdf"
    ),
    "deriving_standard_model_gauge_structure_from_observer_overlap_consistency": Path(
        "paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.pdf"
    ),
    "screen_microphysics_and_observer_synchronization": Path(
        "paper/screen_microphysics_and_observer_synchronization.pdf"
    ),
}
SUPPLEMENTAL_RELEASE_PDFS = {}

# Kept as a module-level alias for the focused regression tests. The canonical
# registry lives beside the source-derived TeX inventory.
NON_TEX_SOURCE_PDFS = paper_sources.NON_TEX_SOURCE_PDFS


def main() -> int:
    args = parse_args()
    repo_root = Path(__file__).resolve().parent.parent
    verify_lean_theorem_count(repo_root)
    release_info = (repo_root / RELEASE_INFO_RELATIVE).read_text(encoding="utf-8")
    release_id = extract_macro(release_info, "OPHPaperReleaseID")
    release_date = extract_macro(release_info, "OPHPaperReleaseDate")
    output_path = repo_root / OUTPUT_RELATIVE
    previous_manifest = load_existing_manifest(output_path)

    manifest = {
        "release_id": release_id,
        "released_at": release_date,
        # Deterministic release metadata: a clean rebuild must not dirty the
        # manifest merely because it ran at a different wall-clock time.
        "generated_at": deterministic_generated_at(release_date),
        "papers": {},
        "supplemental_papers": {},
        "extra_papers": {},
    }
    manifest["book"] = carried_book_manifest_entry(
        repo_root,
        previous_manifest,
        release_id,
    )
    fill_section(repo_root, manifest["papers"], RELEASE_TRACKED_PDFS)
    fill_section(repo_root, manifest["supplemental_papers"], SUPPLEMENTAL_RELEASE_PDFS)
    fill_section(repo_root, manifest["extra_papers"], discover_extra_pdfs(repo_root))

    tagged_manifest = None
    if not args.preview:
        tagged_manifest = load_tagged_manifest(repo_root, release_id)
    enforce_generation_policy(
        tagged_manifest,
        previous_manifest,
        manifest,
        preview=args.preview,
    )
    verify_no_stray_pdfs(repo_root, manifest)
    verify_pdf_release_lines(repo_root, manifest, args.skip_pdf_release_check)
    output_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    mode = "preview" if args.preview else "publication candidate"
    print(f"{output_path} ({mode})")
    return 0


def verify_lean_theorem_count(repo_root: Path) -> None:
    """Refuse to build the manifest while a stated Lean theorem-count floor drifts."""
    check = repo_root / "tools" / "check_lean_theorem_count.py"
    result = subprocess.run([sys.executable, str(check)], cwd=repo_root)
    if result.returncode != 0:
        raise SystemExit(
            "Lean theorem-count floor drifted; run "
            "`python3 tools/check_lean_theorem_count.py --fix` and re-run"
        )


def discover_extra_pdfs(repo_root: Path) -> dict[str, Path]:
    discovered: dict[str, Path] = {}
    for tex_path in paper_sources.EXTRA_PAPERS.values():
        pdf_path = tex_path.with_suffix(".pdf")
        if not pdf_path.is_file():
            raise SystemExit(
                f"missing PDF for extra paper {tex_path.name}: {pdf_path}. "
                "Run python3 tools/build_tex_papers.py --extra-only before regenerating the manifest."
            )
        discovered[tex_path.stem] = pdf_path.relative_to(repo_root)
    return discovered


def verify_no_stray_pdfs(repo_root: Path, manifest: dict) -> None:
    """Reject release-surface PDFs that no source implies.

    Membership is derived from the sources: the curated release set for
    ``paper/`` and the ``extra/*.tex`` glob for ``extra/``. Any other PDF in
    those directories would ship unhashed and unaudited, so its presence
    fails the manifest build (issue #514).
    """
    expected = {
        str(Path(payload["pdf_path"]))
        for payload in manifest_pdf_entries(manifest).values()
    }
    for pdf, source in NON_TEX_SOURCE_PDFS.items():
        if not (repo_root / source).is_file():
            raise SystemExit(
                f"registered source for {pdf} is missing: {source}. "
                "Fix the NON_TEX_SOURCE_PDFS registry."
            )
        if not (repo_root / pdf).is_file():
            raise SystemExit(
                f"registered non-TeX output is missing: {pdf}. "
                f"Rebuild it from {source} before regenerating the manifest."
            )
        expected.add(str(pdf))
    actual = {
        str(pdf.relative_to(repo_root))
        for pdf in (repo_root / "paper").glob("*.pdf")
    }
    actual.update(
        str(pdf.relative_to(repo_root))
        for tex_path in paper_sources.EXTRA_PAPERS.values()
        if (pdf := repo_root / "extra" / tex_path.with_suffix(".pdf").name).is_file()
    )
    strays = sorted(actual - expected)
    if strays:
        raise SystemExit(
            "stray PDFs on the release surface are not implied by any source: "
            f"{', '.join(strays)}. Remove them or register their source before "
            "regenerating the manifest."
        )


def fill_section(repo_root: Path, section: dict, pdfs: dict[str, Path]) -> None:
    for paper_id, relative_path in pdfs.items():
        pdf_path = repo_root / relative_path
        if not pdf_path.is_file():
            raise SystemExit(f"missing release PDF for {paper_id}: {pdf_path}")
        section[paper_id] = {
            "pdf_path": str(relative_path),
            "sha256": sha256(pdf_path),
            "size_bytes": pdf_path.stat().st_size,
        }


def book_manifest_entry(repo_root: Path, release_id: str) -> dict:
    book_path = repo_root / BOOK_PDF_RELATIVE
    if not book_path.is_file():
        raise SystemExit(
            f"missing canonical book PDF: {book_path}. "
            "Run python3 tools/build_book_pdf.py before regenerating the manifest."
        )
    return {
        "built_for_release_id": release_id,
        "pdf_path": BOOK_PDF_RELATIVE.as_posix(),
        "sha256": sha256(book_path),
        "size_bytes": book_path.stat().st_size,
    }


def carried_book_manifest_entry(
    repo_root: Path,
    previous_manifest: dict | None,
    release_id: str,
) -> dict:
    """Carry the last builder-stamped book receipt across a paper rebuild.

    A release bump must not relabel old book bytes. The book builder replaces
    this entry only after producing the canonical book for the selected
    release. The missing-entry branch is a one-time migration for manifests
    created before the book joined the canonical bundle.
    """

    if previous_manifest is None or not isinstance(previous_manifest.get("book"), dict):
        return book_manifest_entry(repo_root, release_id)
    payload = dict(previous_manifest["book"])
    expected = book_manifest_entry(
        repo_root,
        str(payload.get("built_for_release_id", "")).strip(),
    )
    if payload != expected:
        raise SystemExit(
            "the canonical book PDF and its manifest receipt disagree. "
            "Run python3 tools/build_book_pdf.py before regenerating the paper manifest."
        )
    return payload


def update_book_manifest_entry(repo_root: Path) -> None:
    """Stamp the canonical book bytes for the manifest's selected release."""

    output_path = repo_root / OUTPUT_RELATIVE
    manifest = load_existing_manifest(output_path)
    if manifest is None:
        raise SystemExit(
            f"missing {output_path}; generate the paper manifest before building "
            "the canonical book"
        )
    release_info = (repo_root / RELEASE_INFO_RELATIVE).read_text(encoding="utf-8")
    source_release_id = extract_macro(release_info, "OPHPaperReleaseID")
    manifest_release_id = str(manifest.get("release_id", "")).strip()
    if manifest_release_id != source_release_id:
        raise SystemExit(
            "paper/release_info.tex and paper/paper_release_manifest.json "
            f"disagree on release ID: {source_release_id} vs {manifest_release_id}"
        )
    manifest["book"] = book_manifest_entry(repo_root, source_release_id)
    output_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Write the current release manifest for the synced paper PDFs and "
            "canonical book receipt."
        ),
    )
    parser.add_argument(
        "--preview",
        action="store_true",
        help=(
            "Write a local review manifest for the current PDF bytes without "
            "requiring a release bump. Preview generation is allowed even when "
            "the current release ID is already tagged; publication still requires "
            "the strict, no-flag mode."
        ),
    )
    parser.add_argument(
        "--skip-pdf-release-check",
        action="store_true",
        help="Skip checking that each local PDF exposes the visible release line.",
    )
    args = parser.parse_args(argv)
    if args.skip_pdf_release_check and not args.preview:
        parser.error("--skip-pdf-release-check may be used only with --preview")
    return args


def extract_macro(text: str, macro_name: str) -> str:
    pattern = re.compile(r"\\newcommand\{\\%s\}\{([^}]*)\}" % re.escape(macro_name))
    match = pattern.search(text)
    if not match:
        raise SystemExit(f"missing macro {macro_name} in release info")
    return match.group(1).strip()


def deterministic_generated_at(release_date: str) -> str:
    try:
        parsed = datetime.strptime(release_date, "%B %d, %Y")
    except ValueError as exc:
        raise SystemExit(
            f"paper release date must use 'Month D, YYYY': {release_date}"
        ) from exc
    return parsed.strftime("%Y-%m-%dT00:00:00Z")


def load_existing_manifest(path: Path) -> dict | None:
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def load_tagged_manifest(repo_root: Path, release_id: str) -> dict | None:
    """Read the manifest committed at an existing local release tag.

    Same-release refreshes are useful while preparing a release. Once the
    release ID has a Git tag, that tag is the immutable local publication
    boundary and its PDF hashes cannot be replaced under the same ID.
    """

    tag_ref = f"refs/tags/{release_id}"
    exists = subprocess.run(
        ["git", "rev-parse", "--verify", "--quiet", tag_ref],
        cwd=repo_root,
        capture_output=True,
        text=True,
        check=False,
    )
    if exists.returncode != 0:
        return None

    tagged_path = f"{tag_ref}:{OUTPUT_RELATIVE.as_posix()}"
    result = subprocess.run(
        ["git", "show", tagged_path],
        cwd=repo_root,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
    )
    if result.returncode != 0:
        detail = result.stderr.strip() or f"git show exited {result.returncode}"
        raise SystemExit(
            f"could not read the release manifest at {tag_ref}: {detail}"
        )
    try:
        payload = json.loads(result.stdout)
    except json.JSONDecodeError as exc:
        raise SystemExit(
            f"release manifest at {tag_ref} is not valid JSON: {exc}"
        ) from exc
    if not isinstance(payload, dict):
        raise SystemExit(f"release manifest at {tag_ref} is not a JSON object")
    return payload


def enforce_tag_immutability(tagged_manifest: dict | None, manifest: dict) -> None:
    """Reject replacement of PDF bytes behind an existing release tag."""

    if tagged_manifest is None:
        return

    tagged_release_id = str(tagged_manifest.get("release_id", "")).strip()
    current_release_id = str(manifest.get("release_id", "")).strip()
    if tagged_release_id != current_release_id:
        return

    changed_papers = changed_manifest_entries(tagged_manifest, manifest)
    if not changed_papers:
        return

    changed_list = ", ".join(sorted(changed_papers))
    raise SystemExit(
        "PDF manifest entries differ from immutable Git tag "
        f"{current_release_id}: {changed_list}. "
        "Run python3 tools/bump_paper_release.py, rebuild every paper PDF, "
        "and regenerate the publication candidate. Use --preview only for "
        "local review; a preview cannot be published under the tagged ID."
    )


def enforce_generation_policy(
    tagged_manifest: dict | None,
    previous_manifest: dict | None,
    manifest: dict,
    *,
    preview: bool,
) -> None:
    """Apply preview or publication-candidate release policy.

    A preview manifest describes the local bytes under review. It is not a
    publication claim, so it may reuse the visible release line from an
    existing tag. Publication-candidate generation keeps both historical
    protections: tagged bytes are immutable and changed bytes require a new
    release ID.
    """

    if preview:
        return
    enforce_tag_immutability(tagged_manifest, manifest)
    enforce_release_bump(previous_manifest, manifest)


def enforce_release_bump(previous_manifest: dict | None, manifest: dict) -> None:
    if previous_manifest is None:
        return

    previous_release_id = str(previous_manifest.get("release_id", "")).strip()
    current_release_id = str(manifest.get("release_id", "")).strip()
    if previous_release_id != current_release_id:
        return

    changed_papers = changed_manifest_entries(previous_manifest, manifest)
    if not changed_papers:
        return

    changed_list = ", ".join(sorted(changed_papers))
    raise SystemExit(
        "PDF manifest entries changed for the current release ID "
        f"{current_release_id}, but the release was not bumped first: {changed_list}. "
        "Run python3 tools/bump_paper_release.py, rebuild all current paper PDFs, and rerun this command. "
        "Use --preview only for local review under the unchanged release line."
    )


def verify_pdf_release_lines(repo_root: Path, manifest: dict, skip_pdf_release_check: bool) -> None:
    if skip_pdf_release_check:
        return

    release_id = str(manifest["release_id"]).strip()
    missing_release_line: list[str] = []
    tool_failures: list[str] = []
    for paper_id, payload in manifest_pdf_entries(manifest).items():
        pdf_path = repo_root / payload["pdf_path"]
        contains_release = pdf_contains_text(pdf_path, f"Paper release: {release_id}")
        if contains_release is True:
            continue
        if contains_release is False:
            missing_release_line.append(paper_id)
        else:
            tool_failures.append(f"{paper_id}: {contains_release}")

    if missing_release_line:
        raise SystemExit(
            "Local PDFs do not expose the current visible release line "
            f"{release_id}: {', '.join(sorted(missing_release_line))}. "
            "Rebuild all current paper PDFs after bumping paper/release_info.tex, then rerun this command. "
            "Every release bump must propagate to every paper PDF so the release IDs stay in sync."
        )

    if tool_failures:
        failures = "; ".join(tool_failures)
        raise SystemExit(
            "Could not verify the visible release line in the local PDFs: "
            f"{failures}. Install pdftotext or rerun the preview with "
            "--preview --skip-pdf-release-check."
        )


def pdf_contains_text(path: Path, needle: str) -> bool | str:
    if shutil.which("pdftotext") is None:
        return "pdftotext not installed"
    try:
        result = subprocess.run(
            ["pdftotext", str(path), "-"],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            check=False,
        )
    except OSError as exc:
        return f"pdftotext failed: {exc}"

    if result.returncode != 0:
        stderr = result.stderr.strip() or f"exit code {result.returncode}"
        return f"pdftotext failed: {stderr}"
    return needle in result.stdout


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def manifest_pdf_entries(manifest: dict) -> dict:
    entries: dict = {}
    for section_name in ("papers", "supplemental_papers", "extra_papers"):
        entries.update(manifest.get(section_name, {}))
    return entries


def manifest_pdf_entries_by_section(manifest: dict) -> dict[str, dict]:
    entries: dict[str, dict] = {}
    if "book" in manifest:
        entries["book"] = manifest["book"]
    for section_name in ("papers", "supplemental_papers", "extra_papers"):
        section = manifest.get(section_name, {})
        for paper_id, payload in section.items():
            entries[f"{section_name}.{paper_id}"] = payload
    return entries


def changed_manifest_entries(baseline: dict, current: dict) -> list[str]:
    """Return added, removed, moved, or changed PDF manifest entries."""

    baseline_entries = manifest_pdf_entries_by_section(baseline)
    current_entries = manifest_pdf_entries_by_section(current)
    return sorted(
        entry_id
        for entry_id in set(baseline_entries) | set(current_entries)
        if baseline_entries.get(entry_id) != current_entries.get(entry_id)
    )


if __name__ == "__main__":
    raise SystemExit(main())
